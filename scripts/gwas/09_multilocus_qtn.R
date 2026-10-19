###############################################################################
# 09_multilocus_qtn.R -- multi-locus QTN detection (3VmrMLM family)
#
# Title:    Rice-skills multi-locus GWAS: QTN / QEI (QTN x environment) /
#           QQI (QTN x QTN epistasis) via the 3VmrMLM method family
# Description:
#   Wraps the 3VmrMLM-family engines, which estimate ALL genetic effects
#   (additive + dominance, and their interaction counterparts) while fully
#   controlling all polygenic backgrounds via a compressed variance
#   component mixed model, then select QTNs with a mrMLM-style two-stage
#   framework (genome-wide scan -> empirical-Bayes multi-locus ->
#   likelihood-ratio test). Compared with the single-locus engines in
#   03_association.R this detects more small-effect and rare QTNs.
#
#   Engines (installed at RUNTIME, never vendored; see
#   docs/references/3vmrmlm_multilocus_qtn.md):
#     auto       : use Fast3VmrMLM if it loads, else IIIVmrMLM (default)
#     IIIVmrMLM  : CRAN, GPL-3; the 3VmrMLM reference implementation
#                  (method = Single_env / Multi_env / Epistasis). Loads on
#                  any R 4.x.
#     Fast3VmrMLM: GitHub, GPL-3 per package DESCRIPTION; big-data engine
#                  for million-marker datasets. Prebuilt zips target
#                  R 4.3; on newer R versions the DLL may fail to load, in
#                  which case this wrapper reports the issue and you can
#                 still run the method-equivalent IIIVmrMLM.
#
# Usage:
#   Rscript scripts/gwas/09_multilocus_qtn.R \
#       --genotype results/gwas/01_qc/qc_geno \
#       --phenotype data/phenotype.csv \
#       --method single_env --trait 1:3 \
#       --pop-structure results/gwas/01_qc/pca.csv --pop-type PC \
#       --outdir results/gwas/09_multilocus
#
# Inputs:
#   --genotype      PLINK binary prefix (.bed/.bim/.fam), OR a SNP x sample
#                   dosage CSV (0/1/2, missing = empty) when --snp-info is
#                   also given (converted to a temporary PLINK set)
#   --snp-info      snp_id, chr, pos  (only with dosage CSV)
#   --phenotype     CSV: sample_id + numeric trait columns. For
#                   --method multi_env supply trait-major wide columns:
#                   one column per trait x environment (trait1's
#                   environments first, then trait2's ...)
#   --trait         1-based indices into the trait columns, e.g. "1:3" or
#                   "1,2" (default: all)
#   --n-en          multi_env only: comma-separated environment counts per
#                   selected trait, e.g. "2,2" (2 traits x 2 environments)
#   --pop-structure optional CSV: PC scores (--pop-type PC) or a Q matrix
#                   from ADMIXTURE/STRUCTURE (--pop-type Q)
#   --svpal         scan-stage P-value threshold (default 0.01; guidance:
#                   N < 4000 -> 0.01, N >= 4000 -> 1e-5)
#   --svrad         de-collinearity search radius in bp (default 2e4)
#   --svmlod        LOD threshold for suggested QTNs (default 3)
#   --engine        auto | IIIVmrMLM | Fast3VmrMLM
#   --draw-plot     also produce the engine's Manhattan plots (PDF)
#   --outdir        output directory
#   --demo          synthetic end-to-end self-test (no input files needed)
#
# Outputs (in --outdir):
#   qtn_results.csv            normalised QTN/QEI/QQI table per trait
#   association_results.csv    snp_id/p table, feeds 05_plot_gwas.R
#   engine_raw/<trait>_*       engine-native outputs (K, midresult,
#                              result.xlsx) copied out of the hidden
#                              .result-*_detection/ folder
#   run_info.txt / summary.txt
#
# Notes:
#   - Upstream engines are NOT bundled (GPL-3; IIIVmrMLM/FastBiCmrMLM ship
#     no acceptable LICENSE for vendoring). Method distilled from the
#     papers cited in the method note; papers themselves are not copied.
#   - Final significance: Bonferroni P = 0.05 / n_markers is the paper's
#     recommended threshold; --svpal only controls which markers enter the
#     multi-locus stage.
#   - PLINK v1.9 on PATH is required only for the epistasis module of
#     Fast3VmrMLM (it shells out to plink); IIIVmrMLM needs no PLINK.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
source(file.path(dirname(script_file), "..", "utils", "genomics_utils.R"))
source(file.path(dirname(script_file), "..", "utils", "vmrml_utils.R"))
rice_set_seed(42)

parse_args <- function(args) {
  opts <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) stop("Unexpected argument '", a, "'. Use --key value or --flag.", call. = FALSE)
    key <- sub("^--", "", a)
    if (i == length(args) || grepl("^--", args[i + 1L])) { opts[[key]] <- TRUE; i <- i + 1L }
    else { opts[[key]] <- args[i + 1L]; i <- i + 2L }
  }
  opts
}
opt <- function(opts, key, default = NULL) {
  v <- opts[[key]]; if (is.null(v)) return(default)
  if (isTRUE(v)) return(if (is.logical(default) && length(default) == 1L) default else TRUE)
  v
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

method <- tolower(opt(opts, "method", "single_env"))
engine <- tolower(opt(opts, "engine", "auto"))
genotype_path <- opt(opts, "genotype")
snp_info_path <- opt(opts, "snp-info")
phenotype_path <- opt(opts, "phenotype")
trait_spec <- opt(opts, "trait", NULL)
n_en_spec <- opt(opts, "n-en", NULL)
ps_path <- opt(opts, "pop-structure", NULL)
ps_type <- toupper(opt(opts, "pop-type", "PC"))
svpal <- suppressWarnings(as.numeric(opt(opts, "svpal", "0.01")))
svrad <- suppressWarnings(as.numeric(opt(opts, "svrad", "20000")))
svmlod <- suppressWarnings(as.numeric(opt(opts, "svmlod", "3")))
draw_plot <- isTRUE(opts[["draw-plot"]]) || identical(opt(opts, "draw-plot", FALSE), TRUE)
outdir <- opt(opts, "outdir", "results/gwas/09_multilocus")
demo_mode <- isTRUE(opts[["demo"]])

if (!method %in% c("single_env", "multi_env", "epistasis")) {
  stop("--method must be single_env / multi_env / epistasis; got '", method, "'.", call. = FALSE)
}
if (!engine %in% c("auto", "iiivmrml", "fast3vmrml")) {
  stop("--engine must be auto / IIIVmrMLM / Fast3VmrMLM; got '", engine, "'.", call. = FALSE)
}
if (!ps_type %in% c("PC", "Q")) stop("--pop-type must be PC or Q.", call. = FALSE)

parse_trait_spec <- function(spec, n_avail) {
  if (is.null(spec)) return(seq_len(n_avail))
  spec <- gsub("\\s", "", spec)
  if (grepl(":", spec)) {
    p <- as.integer(strsplit(spec, ":")[[1]])
    if (length(p) != 2 || any(is.na(p))) stop("Bad --trait range '", spec, "'.", call. = FALSE)
    return(seq(p[1], p[2]))
  }
  v <- as.integer(strsplit(spec, ",")[[1]])
  if (any(is.na(v))) stop("Bad --trait list '", spec, "'.", call. = FALSE)
  v
}

need_pkgs("jsonlite")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir)                       # absolute from here on
engine_raw <- file.path(outdir, "engine_raw")
if (!dir.exists(engine_raw)) dir.create(engine_raw, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 0) Demo mode: synthetic data + engine self-test
## ---------------------------------------------------------------------------
demo_dir <- NULL
if (demo_mode) {
  message("== DEMO MODE: generating synthetic PLINK dataset (seed = 42) ==")
  demo_dir <- file.path(outdir, "_demo_work")
  dir.create(demo_dir, showWarnings = FALSE)
  set.seed(42)
  n <- 200; m <- 400
  chr <- rep(1:4, each = m / 4)
  pos <- ave(seq_len(m), chr, FUN = function(x) x * 5000)
  snp_id <- paste0("S", seq_len(m))
  snp_info <- data.frame(snp_id = snp_id, chr = chr, pos = pos,
                         stringsAsFactors = FALSE)
  af <- runif(m, 0.2, 0.5)
  geno <- matrix(rbinom(n * m, 2, rep(af, each = n)), nrow = n, ncol = m,
                 dimnames = list(paste0("Ind", seq_len(n)), snp_id))
  qtl_idx <- c(120, 260)                                   # 2 planted QTNs
  y <- 2.8 * geno[, qtl_idx[1]] + 1.8 * geno[, qtl_idx[2]] +
       scale(rnorm(n)) * 2 + rnorm(n, 0, 2)
  pheno <- data.frame(sample_id = rownames(geno), trait1 = y,
                      trait2 = scale(rnorm(n, 100, 10)),
                      stringsAsFactors = FALSE)
  if (method == "multi_env") {
    # widen trait1 into 2 environments (+ a QEI on the planted QTN)
    pheno <- data.frame(
      sample_id = pheno$sample_id,
      trait1__E1 = pheno$trait1 - 3,
      trait1__E2 = pheno$trait1 + 3 + 0.6 * geno[, qtl_idx[1]],
      trait2__E1 = pheno$trait2, trait2__E2 = pheno$trait2 + 1.5,
      stringsAsFactors = FALSE)
  }
  genotype_path <- file.path(demo_dir, "demo_geno")
  write_plink_set(genotype_path, geno, snp_info)
  phenotype_path <- file.path(demo_dir, "demo_pheno.csv")
  env_labels <- if (method == "multi_env") c("E1", "E2", "E1", "E2") else NULL
  write_engine_pheno(pheno, phenotype_path, env_labels = env_labels)
  # single_env: --trait selects trait columns; multi_env: --trait selects the
  # two traits, each contributing n_en columns (env-per-trait layout).
  trait_spec <- if (method == "multi_env") "1:2" else NULL
  n_en_spec <- if (method == "multi_env") "2,2" else NULL
  message(sprintf("  synthetic: %d individuals x %d SNPs, planted QTNs: %s",
                  n, m, paste(snp_id[qtl_idx], collapse = ", ")))
  message(sprintf("  files: %s(.bed/.bim/.fam) + %s",
                  basename(genotype_path), basename(phenotype_path)))
} else {
  if (is.null(genotype_path) || is.null(phenotype_path)) {
    stop("Required: --genotype and --phenotype (or --demo).\n",
         "Hint: see the Usage block at the top of this script.", call. = FALSE)
  }
}

## ---------------------------------------------------------------------------
## 1) Resolve genotype to a PLINK prefix
## ---------------------------------------------------------------------------
plink_prefix <- NULL
if (all(file.exists(paste0(genotype_path, c(".bed", ".bim", ".fam"))))) {
  plink_prefix <- genotype_path
  message(sprintf("PLINK set found: %s(.bed/.bim/.fam)", genotype_path))
} else if (file.exists(genotype_path)) {
  if (is.null(snp_info_path)) {
    stop("--genotype is not a PLINK set; a dosage CSV requires --snp-info ",
         "(snp_id, chr, pos).", call. = FALSE)
  }
  gt <- read_genotype(genotype_path)
  snp_info <- read_snp_info(snp_info_path)
  conv_dir <- file.path(outdir, "_plink_convert")
  if (!dir.exists(conv_dir)) dir.create(conv_dir, recursive = TRUE, showWarnings = FALSE)
  plink_prefix <- file.path(conv_dir, "converted")
  write_plink_set(plink_prefix, gt, snp_info)
  message(sprintf("Converted dosage CSV (%d x %d) to PLINK set: %s",
                  nrow(gt), ncol(gt), plink_prefix))
} else {
  stop("Genotype not found: '", genotype_path, "' (expected a PLINK prefix or a CSV).",
       call. = FALSE)
}

## ---------------------------------------------------------------------------
## 2) Phenotype in engine format
## ---------------------------------------------------------------------------
pheno_raw <- utils::read.csv(phenotype_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
if (ncol(pheno_raw) < 2) stop("Phenotype CSV needs sample_id + >= 1 trait column.", call. = FALSE)
trait_cols <- setdiff(colnames(pheno_raw), colnames(pheno_raw)[1])
n_all <- length(trait_cols)
if (method == "multi_env") {
  # --trait selects TRAITS; every trait's trait-major columns go to the engine
  if (is.null(n_en_spec)) {
    stop("--method multi_env requires --n-en (environment counts per trait), e.g. --n-en 2,2.",
         call. = FALSE)
  }
  n_en <- as.integer(strsplit(gsub("\\s", "", n_en_spec), ",")[[1]])
  n_trait <- length(n_en)
  if (sum(n_en) != n_all) {
    stop("--n-en sums to ", sum(n_en), " but the phenotype has ", n_all,
         " trait columns (one column per trait x environment, trait-major).",
         call. = FALSE)
  }
  if (is.null(trait_spec)) {
    trait_spec <- paste0("1:", n_trait)
  }
  trait_idx <- parse_trait_spec(trait_spec, n_trait)
  if (length(trait_idx) != n_trait) {
    stop("multi_env analyses ALL selected traits jointly; give --n-en with one count ",
         "per trait and --trait covering the same set (demo uses 1:2 / 2,2).",
         call. = FALSE)
  }
  pheno_engine <- pheno_raw                                   # all columns
} else {
  trait_idx <- parse_trait_spec(trait_spec, n_all)
  if (any(trait_idx < 1) || any(trait_idx > n_all)) {
    stop("--trait out of range; available trait columns: ",
         paste(trait_cols, collapse = ", "), call. = FALSE)
  }
  pheno_engine <- pheno_raw[, c(1, trait_idx + 1)]
}
colnames(pheno_engine)[1] <- "<Phenotype>"
phe_engine_path <- file.path(outdir, "phenotype_engine_format.csv")
utils::write.csv(pheno_engine, phe_engine_path, row.names = FALSE, quote = FALSE)

fam_ids <- scan(paste0(plink_prefix, ".fam"), what = character(), quiet = TRUE)
fam_ids <- fam_ids[seq(2, length(fam_ids), by = 6)]          # IID column
n_shared <- sum(fam_ids %in% pheno_raw[[1]])
message(sprintf("Samples: %d in genotype, %d shared with phenotype.",
                length(fam_ids), n_shared))
if (n_shared < 50) {
  stop("Fewer than 50 shared samples. Hint: check that phenotype sample IDs ",
       "match the PLINK .fam IIDs.", call. = FALSE)
}

## ---------------------------------------------------------------------------
## 3) Engine selection
## ---------------------------------------------------------------------------
use_fast <- FALSE
if (engine == "fast3vmrml") {
  chk <- engine_ready("Fast3VmrMLM",
    "Prebuilt zips target R 4.3; on newer R the DLL may fail to load. ",
    "Install: see docs/references/3vmrmlm_multilocus_qtn.md (runtime install, GPL-3).")
  if (!chk$ok) stop(chk$message, call. = FALSE)
  use_fast <- TRUE
} else if (engine == "auto" || engine == "iiivmrml") {
  chk <- engine_ready("IIIVmrMLM",
    "Install from CRAN: install.packages('IIIVmrMLM') (GPL-3, runtime install).")
  if (engine == "auto" && chk$ok) {
    chk2 <- engine_ready("Fast3VmrMLM")
    use_fast <- chk2$ok
    if (!use_fast) {
      message("Fast3VmrMLM not loadable (its prebuilt DLL targets R 4.3); ",
              "using IIIVmrMLM (method-equivalent 3VmrMLM implementation).")
    }
  }
  if (!chk$ok) {
    if (engine == "auto") {
      stop("Neither IIIVmrMLM nor Fast3VmrMLM is installed.\n",
           "  install.packages('IIIVmrMLM')            # CRAN, GPL-3\n",
           "  # Fast3VmrMLM (GitHub, GPL-3) for million-marker datasets:\n",
           "  # see docs/references/3vmrmlm_multilocus_qtn.md", call. = FALSE)
    }
    stop(chk$message, call. = FALSE)
  }
} else {
  stop("unreachable", call. = FALSE)
}
message(sprintf("Engine: %s | method: %s", if (use_fast) "Fast3VmrMLM" else "IIIVmrMLM", method))

ps_arg <- NULL; ps_type_arg <- NULL
if (!is.null(ps_path)) {
  if (!file.exists(ps_path)) stop("Population structure file not found: ", ps_path, call. = FALSE)
  ps_arg <- normalizePath(ps_path); ps_type_arg <- ps_type
}

## ---------------------------------------------------------------------------
## 4) Run
## ---------------------------------------------------------------------------
engine_method <- c(single_env = "Single_env", multi_env = "Multi_env",
                   epistasis = "Epistasis")[[method]]   # IIIVmrMLM spelling
run_dir <- file.path(outdir, "_engine_run")
if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
wd0 <- getwd()
on.exit(setwd(wd0), add = TRUE)
setwd(run_dir)
geno_rel <- normalizePath(plink_prefix, mustWork = FALSE)
genofile_rel <- if (.Platform$OS.type == "windows") gsub("\\\\", "/", geno_rel) else geno_rel
phe_rel <- gsub("\\\\", "/", normalizePath(phe_engine_path, mustWork = FALSE))

t0 <- Sys.time()
status <- "ok"
# Engine warnings (e.g. "No QTN is detected in <trait>") are informative, not
# fatal: demote them to messages so one trait without QTNs cannot kill the run.
quiet_warn <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    message("engine warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}
if (!use_fast) {
  message("Running IIIVmrMLM (", engine_method, ") ...")
  r <- tryCatch({
    suppressMessages(library(IIIVmrMLM))
    if (method == "multi_env") {
      quiet_warn(IIIVmrMLM(fileGen = genofile_rel, filePhe = phe_rel,
                fileKin = NULL, filePS = ps_arg, PopStrType = ps_type_arg,
                fileCov = NULL, method = engine_method, trait = trait_idx,
                n.en = n_en, SearchRadius = svrad, svpal = svpal,
                DrawPlot = draw_plot, dir = ".", Chr_name_com = NULL))
    } else {
      quiet_warn(IIIVmrMLM(fileGen = genofile_rel, filePhe = phe_rel,
                fileKin = NULL, filePS = ps_arg, PopStrType = ps_type_arg,
                fileCov = NULL, method = engine_method, trait = trait_idx,
                SearchRadius = svrad, svpal = svpal, DrawPlot = draw_plot,
                dir = ".", Chr_name_com = NULL))
    }
    "ok"
  }, error = function(e) paste("RUN FAIL:", conditionMessage(e)))
  status <- r
} else {
  message("Running Fast3VmrMLM ...")
  r <- tryCatch({
    suppressMessages(library(Fast3VmrMLM))
    if (method == "single_env") {
      Fast3VmrMLM(fileGen = genofile_rel, filePhe = phe_rel, filePS = ps_arg,
                  PopStrType = if (is.null(ps_arg)) "PC" else ps_type,
                  fileOut = "./", trait = trait_idx,
                  svrad = svrad, svpal = svpal, svmlod = svmlod,
                  DrawPlot = draw_plot, Plotformat = "pdf")
    } else if (method == "multi_env") {
      Fast3VmrMLM_MEJA(fileGen = genofile_rel, filePhe = phe_rel,
                       filePS = ps_arg, PopStrType = if (is.null(ps_arg)) "PC" else ps_type,
                       fileOut = "./", trait = trait_idx, n_en = n_en,
                       svrad = svrad, svpal = svpal, svmlod = svmlod,
                       DrawPlot = draw_plot, Plotformat = "pdf")
    } else {
      # epistasis: Fast3VmrMLM shells out to PLINK v1.9
      plink_exe <- Sys.which("plink")
      if (length(plink_exe) == 0 || plink_exe == "") {
        stop("Fast3VmrMLM epistasis requires PLINK v1.9 on PATH ",
             "(plinkdir parameter). Not found.\n",
             "  Alternative: --engine IIIVmrMLM runs Epistasis mode with no PLINK needed.",
             call. = FALSE)
      }
      m_markers <- length(readLines(paste0(plink_prefix, ".bim"), warn = FALSE))
      svpal_qqi <- if (n_shared >= 4000) 0.05 / (m_markers * (m_markers - 1) / 2) else 1e-5
      message(sprintf("epistasis scan: %d markers, %d individuals -> svpal=%g",
                      m_markers, n_shared, svpal_qqi))
      Fast3VmrMLM_QQI(fileGen = genofile_rel, filePhe = phe_rel, filePS = ps_arg,
                      PopStrType = if (is.null(ps_arg)) "PC" else ps_type,
                      fileOut = "./", genoType = "SNP", trait = trait_idx,
                      thinMarker = TRUE, LDthin_threshold = 0.8,
                      LDthin_max.dist = 200e3,
                      plinkdir = dirname(plink_exe), plinkepi1 = 1e-5,
                      svpal = svpal_qqi, svpal_ad = 1e-3,
                      svrad = 1e4, svmlod = svmlod, svmlod_QQI = svmlod)
    }
    "ok"
  }, error = function(e) paste("RUN FAIL:", conditionMessage(e)))
  status <- r
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (is.character(status) && grepl("^RUN FAIL", status)) {
  stop(status, "\nHint: check the engine's own console output above.", call. = FALSE)
}

## ---------------------------------------------------------------------------
## 5) Collect outputs (engine writes a hidden .result-*_detection/ dir)
## ---------------------------------------------------------------------------
hidden_dirs <- list.files(".", pattern = "^\\.result-", all.files = TRUE,
                          include.dirs = TRUE, full.names = TRUE)
collected <- character(0)
for (hd in hidden_dirs) {
  fs <- list.files(hd, full.names = TRUE)
  for (f in fs) {
    dest <- file.path(engine_raw, basename(hd), basename(f))
    if (!dir.exists(dirname(dest))) dir.create(dirname(dest), recursive = TRUE,
                                               showWarnings = FALSE)
    ok <- file.copy(f, dest, overwrite = TRUE)
    if (ok) collected <- c(collected, dest)
  }
}
setwd(wd0)
if (length(collected) == 0L) {
  stop("The engine produced no output files. ", engine_method,
       " did not write results -- check the console above for engine warnings ",
       "(e.g. wrong method name, empty phenotype).", call. = FALSE)
}

## ---------------------------------------------------------------------------
## 6) Normalised summary + plot-ready association table
## ---------------------------------------------------------------------------
qtn_tbl <- data.frame()
assoc_tbl <- data.frame()
for (f in collected) {
  if (!grepl("_result\\.xlsx$", f)) next
  trait_tag <- sub("_result\\.xlsx$", "", basename(f))
  sheets <- xlsx_sheet_names(f)
  for (sh in sheets) {
    df <- suppressWarnings(try(openxlsx::read.xlsx(f, sheet = sh), silent = TRUE))
    if (inherits(df, "try-error") || is.null(df) || nrow(df) == 0) next
    if (!"Marker" %in% names(df)) next
    q <- normalise_qtns(df, trait_label = paste0(trait_tag, "/", sh))
    qtn_tbl <- rbind(qtn_tbl, q)
    assoc_tbl <- rbind(assoc_tbl, data.frame(
      snp_id = q$marker, chr = q$chr, pos = q$pos, p = q$p,
      engine = if (use_fast) "Fast3VmrMLM" else "IIIVmrMLM",
      stringsAsFactors = FALSE))
  }
}
if (nrow(qtn_tbl)) {
  utils::write.csv(qtn_tbl, file.path(outdir, "qtn_results.csv"), row.names = FALSE)
  utils::write.csv(assoc_tbl, file.path(outdir, "association_results.csv"),
                   row.names = FALSE)
} else {
  utils::write.csv(qtn_tbl, file.path(outdir, "qtn_results.csv"), row.names = FALSE)
}

summary_lines <- c(
  sprintf("method: %s", method),
  sprintf("engine: %s", if (use_fast) "Fast3VmrMLM" else "IIIVmrMLM"),
  sprintf("genotype: %s", plink_prefix),
  sprintf("phenotype: %s", phenotype_path),
  sprintf("traits analysed: %s", paste(trait_cols[trait_idx], collapse = ", ")),
  sprintf("shared samples: %d", n_shared),
  sprintf("svpal: %g | svrad: %g | svmlod: %g", svpal, svrad, svmlod),
  sprintf("QTNs/QEIs/QQIs in results: %d", nrow(qtn_tbl)),
  sprintf("elapsed: %.1f s", elapsed),
  "",
  "Significance guidance (3VmrMLM papers): use the Bonferroni threshold",
  "P = 0.05 / n_markers for final QTNs; --svpal only controls which markers",
  "enter the multi-locus (machine-learning) stage.",
  "",
  "Plot: Rscript scripts/gwas/05_plot_gwas.R --assoc association_results.csv \\",
  "       --snp-info <snp_id,chr,pos csv> --outdir <dir>"
)
writeLines(summary_lines, file.path(outdir, "summary.txt"))

write_run_info(file.path(outdir, "run_info.txt"), list(
  script = "scripts/gwas/09_multilocus_qtn.R",
  method = method, engine = if (use_fast) "Fast3VmrMLM" else "IIIVmrMLM",
  genotype = plink_prefix, phenotype = phenotype_path,
  traits = paste(trait_cols[trait_idx], collapse = ","),
  n_shared = n_shared, svpal = svpal, svrad = svrad, svmlod = svmlod,
  pop_structure = if (is.null(ps_path)) "none" else paste0(ps_path, " (", ps_type, ")"),
  elapsed_sec = round(elapsed, 1),
  demo = demo_mode
))

cat("\n")
cat("=== 09_multilocus_qtn.R finished ===\n")
cat(sprintf("method=%s engine=%s traits=%d shared=%d elapsed=%.1fs\n",
            method, if (use_fast) "Fast3VmrMLM" else "IIIVmrMLM",
            length(trait_idx), n_shared, elapsed))
cat(sprintf("QTNs/QEIs/QQIs: %d -> %s\n", nrow(qtn_tbl),
            file.path(outdir, "qtn_results.csv")))
cat(sprintf("outputs: %s\n", outdir))
