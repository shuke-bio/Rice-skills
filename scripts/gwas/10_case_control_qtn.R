###############################################################################
# 10_case_control_qtn.R -- case-control (binary trait) GWAS via FastBiCmrMLM
#
# Title:    Rice-skills binary/case-control QTN detection (FastBiCmrMLM)
# Description:
#   Wraps FastBiCmrMLM, a fast multi-locus mixed LOGISTIC model for
#   binary (case-control) traits -- e.g. disease resistance / susceptibility
#   scores in crops, or case-control phenotypes in human cohorts. Like the
#   3VmrMLM family it estimates all effects (additive + dominance) and fully
#   controls all polygenic backgrounds, but with a logistic link suitable
#   for 0/1 traits; LD markers can be replaced by bin haplotypes
#   (FastBiCmrMLM-Hap), enabling rare-variant and small-alpha detection.
#
#   Upstream: Wang JT, Chang XY, Zhao Q, Zhang YM (2024) Briefings in
#   Bioinformatics 25(4): bbae290. Software github.com/YuanmingZhang65/
#   FastBiCmrMLM -- NO open-source LICENSE file, so it is NEVER vendored:
#   this wrapper installs it at runtime from the official zip and records
#   the provenance. *** The shipped build is Linux-only. ***
#
#   Algorithm choice by sample size n:
#     FastBiCmrMLM        n < 1,000        (primary algorithm)
#     FastBiCmrMLM-Time   1,000 <= n < 20,000   (~7 min per WTCCC1 disease)
#     FastBiCmrMLM-RAM    n >= 20,000      (memory-efficient; 500k x 1M
#                                         SNPs in ~13.5 h / 10 CPU / 141 GB)
#     FastBiCmrMLM-Hap    bin haplotypes   (LD-marker blocks, c/NumHap)
#
# Usage (Linux / WSL):
#   Rscript scripts/gwas/10_case_control_qtn.R \
#       --genotype results/gwas/01_qc/qc_geno \
#       --phenotype data/binary_trait.csv \
#       --trait blast_resistance \
#       --method auto \
#       --outdir results/gwas/10_case_control
#
# Inputs:
#   --genotype      PLINK binary prefix (.bed/.bim/.fam), OR a SNP x sample
#                   dosage CSV (0/1/2, missing = empty) with --snp-info
#   --snp-info      snp_id, chr, pos  (only with dosage CSV)
#   --phenotype     CSV: sample_id + binary trait column(s) (0/1); support
#                   continuous or categorical phenotypes via --cut-to-binary
#   --trait         name (or comma-separated names) of binary trait columns
#   --cut-to-binary cut a continuous phenotype at this threshold
#                   (<= threshold -> 0 control, > threshold -> 1 case)
#   --method        auto | FastBiCmrMLM | FastBiCmrMLM-Time |
#                   FastBiCmrMLM-RAM | FastBiCmrMLM-Hap
#   --pop-structure CSV of PCs (PCA-type structure), optional
#   --search-radius de-collinearity window in bp (default 1e6)
#   --svpal         scan-stage P threshold (default 1e-5)
#   --lod-threshold LOD for suggested variants (default 3)
#   --outdir        output directory
#   --demo          synthetic self-test (Linux only; other OSes skip with
#                   a clear message)
#
# Outputs (in --outdir):
#   qtn_results.csv          normalised result table
#   association_results.csv  snp_id/p table, feeds 05_plot_gwas.R
#   engine_raw/              engine-native workbook(s)
#   run_info.txt / summary.txt
#
# Notes:
#   - Windows users: FastBiCmrMLM ships no Windows build; run this script
#     under Linux or WSL2. For quantitative traits use
#     scripts/gwas/09_multilocus_qtn.R (IIIVmrMLM works on Windows).
#   - The upstream software has no LICENSE file; using it is at your own
#     responsibility (cite the paper). Nothing is bundled into this repo.
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
method <- opt(opts, "method", "auto")
genotype_path <- opt(opts, "genotype")
snp_info_path <- opt(opts, "snp-info")
phenotype_path <- opt(opts, "phenotype")
trait <- opt(opts, "trait")
cut_binary <- suppressWarnings(as.numeric(opt(opts, "cut-to-binary", "NA")))
ps_path <- opt(opts, "pop-structure", NULL)
search_radius <- suppressWarnings(as.numeric(opt(opts, "search-radius", "1e6")))
svpal <- suppressWarnings(as.numeric(opt(opts, "svpal", "1e-5")))
lod_thr <- suppressWarnings(as.numeric(opt(opts, "lod-threshold", "3")))
outdir <- opt(opts, "outdir", "results/gwas/10_case_control")
demo_mode <- isTRUE(opts[["demo"]])

valid_methods <- c("auto", "FastBiCmrMLM", "FastBiCmrMLM-Time",
                   "FastBiCmrMLM-RAM", "FastBiCmrMLM-Hap")
if (!method %in% valid_methods) {
  stop("--method must be one of ", paste(valid_methods, collapse = " / "),
       "; got '", method, "'.", call. = FALSE)
}

need_pkgs("jsonlite")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir)

## ---------------------------------------------------------------------------
## 0) OS gate (FastBiCmrMLM is Linux-only) + engine check
## ---------------------------------------------------------------------------
on_linux <- .Platform$OS.type == "unix"
if (!on_linux) {
  msg <- c("FastBiCmrMLM ships a Linux-only build (no Windows binary).",
           "This script must run under Linux or WSL2.",
           "Quantitative-trait multi-locus GWAS on Windows: use",
           "  scripts/gwas/09_multilocus_qtn.R  (IIIVmrMLM, CRAN binary)")
  if (demo_mode) {
    message("== DEMO MODE (skipped) ==")
    for (m in msg) message("  ", m)
    writeLines(c("status: SKIPPED (non-Linux host)",
                 "FastBiCmrMLM is Linux-only; run under Linux/WSL2 to exercise this route."),
               file.path(outdir, "run_info.txt"))
    quit(status = 0)
  }
  stop(paste(msg, collapse = "\n"), call. = FALSE)
}

chk <- engine_ready("FastBiCmrMLM",
  paste("Install (runtime, from the official GitHub zip):",
        "download.file('https://github.com/YuanmingZhang65/FastBiCmrMLM/raw/main/FastBiCmrMLM_Linux.zip',",
        "  'FastBiCmrMLM_Linux.zip', mode='wb'); unzip(...);",
        "  install.packages('FastBiCmrMLM', repos=NULL)",
        "(GPL-style research software; no LICENSE file -- use at your own",
        "responsibility, cite Wang et al. 2024, Brief Bioinform 25: bbae290.)",
        sep = "\n  "))
if (!chk$ok) stop(chk$message, call. = FALSE)

## ---------------------------------------------------------------------------
## 1) Demo data (Linux)
## ---------------------------------------------------------------------------
demo_dir <- NULL
if (demo_mode) {
  message("== DEMO MODE: synthetic binary-trait PLINK dataset (seed = 42) ==")
  demo_dir <- file.path(outdir, "_demo_work")
  dir.create(demo_dir, showWarnings = FALSE)
  set.seed(42)
  n <- 600; m <- 300
  chr <- rep(1:3, each = m / 3)
  pos <- ave(seq_len(m), chr, FUN = function(x) x * 4000)
  snp_id <- paste0("S", seq_len(m))
  snp_info <- data.frame(snp_id = snp_id, chr = chr, pos = pos,
                         stringsAsFactors = FALSE)
  af <- runif(m, 0.1, 0.4)
  geno <- matrix(rbinom(n * m, 2, rep(af, each = n)), nrow = n, ncol = m,
                 dimnames = list(paste0("Ind", seq_len(n)), snp_id))
  qtl_idx <- c(88, 190)
  lin <- 1.2 * geno[, qtl_idx[1]] + 0.9 * geno[, qtl_idx[2]] + rnorm(n, 0, 1.5)
  binary <- as.integer(lin > median(lin))              # balanced case/control
  pheno <- data.frame(sample_id = rownames(geno), blast_resistance = binary,
                      stringsAsFactors = FALSE)
  genotype_path <- file.path(demo_dir, "demo_geno")
  write_plink_set(genotype_path, geno, snp_info)
  phenotype_path <- file.path(demo_dir, "demo_pheno.csv")
  utils::write.csv(pheno, phenotype_path, row.names = FALSE, quote = FALSE)
  trait <- "blast_resistance"
  message(sprintf("  synthetic: %d individuals x %d SNPs, planted QTNs: %s; cases=%d controls=%d",
                  n, m, paste(snp_id[qtl_idx], collapse = ", "),
                  sum(binary), n - sum(binary)))
} else {
  if (is.null(genotype_path) || is.null(phenotype_path) || is.null(trait)) {
    stop("Required: --genotype, --phenotype, --trait (or --demo).", call. = FALSE)
  }
}

## ---------------------------------------------------------------------------
## 2) Genotype -> PLINK prefix
## ---------------------------------------------------------------------------
if (all(file.exists(paste0(genotype_path, c(".bed", ".bim", ".fam"))))) {
  plink_prefix <- genotype_path
} else if (file.exists(genotype_path)) {
  if (is.null(snp_info_path)) stop("Dosage CSV needs --snp-info (snp_id, chr, pos).", call. = FALSE)
  gt <- read_genotype(genotype_path); si <- read_snp_info(snp_info_path)
  conv_dir <- file.path(outdir, "_plink_convert")
  dir.create(conv_dir, recursive = TRUE, showWarnings = FALSE)
  plink_prefix <- file.path(conv_dir, "converted")
  write_plink_set(plink_prefix, gt, si)
} else {
  stop("Genotype not found: ", genotype_path, call. = FALSE)
}

## ---------------------------------------------------------------------------
## 3) Phenotype -> engine format
## ---------------------------------------------------------------------------
pheno_raw <- utils::read.csv(phenotype_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
if (!trait %in% colnames(pheno_raw)) {
  stop("Trait '", trait, "' not in phenotype. Available: ",
       paste(setdiff(colnames(pheno_raw), colnames(pheno_raw)[1]), collapse = ", "),
       call. = FALSE)
}
y <- pheno_raw[[trait]]
if (!is.na(cut_binary)) {
  y <- as.integer(y > cut_binary)
  message(sprintf("Cut continuous trait at %g -> case=%d control=%d",
                  cut_binary, sum(y == 1), sum(y == 0)))
}
if (!all(y %in% c(0, 1))) {
  stop("Trait '", trait, "' must be binary (0/1). Use --cut-to-binary for ",
       "continuous phenotypes.", call. = FALSE)
}
n_ind <- length(unique(pheno_raw[[1]]))
message(sprintf("Binary trait '%s': %d cases / %d controls (n = %d).",
                trait, sum(y == 1), sum(y == 0), n_ind))

eng_method <- method
if (method == "auto") {
  eng_method <- if (n_ind >= 20000) "FastBiCmrMLM-RAM" else
                if (n_ind >= 1000) "FastBiCmrMLM-Time" else "FastBiCmrMLM"
}
message(sprintf("Algorithm: %s (auto-selected by n = %d)", eng_method, n_ind))

phe_engine <- data.frame("<Phenotype>" = pheno_raw[[1]], "trait1" = y,
                         stringsAsFactors = FALSE)
phe_engine_path <- file.path(outdir, "phenotype_engine_format.csv")
utils::write.csv(phe_engine, phe_engine_path, row.names = FALSE, quote = FALSE)

## ---------------------------------------------------------------------------
## 4) Run
## ---------------------------------------------------------------------------
run_dir <- file.path(outdir, "_engine_run")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
wd0 <- getwd(); on.exit(setwd(wd0), add = TRUE); setwd(run_dir)
geno_abs <- gsub("\\\\", "/", normalizePath(plink_prefix, mustWork = FALSE))
phe_abs <- gsub("\\\\", "/", normalizePath(phe_engine_path, mustWork = FALSE))
ps_abs <- if (is.null(ps_path)) NULL else
  gsub("\\\\", "/", normalizePath(ps_path, mustWork = FALSE))

t0 <- Sys.time()
r <- tryCatch({
  suppressMessages(library(FastBiCmrMLM))
  FastBiCmrMLM(fileGen = geno_abs, filePhe = phe_abs, fileGRM = NULL,
               filePS = ps_abs, PopStrType = "PCA", fileOut = "./",
               method = eng_method, trait = 1, SearchRadius = search_radius,
               svpal = svpal, LODThreshold = lod_thr, CutToBinary = NULL,
               SubforGRM = FALSE, SubNum = 1e4, ChunkNum = 5e4, fastSPA = FALSE)
  "ok"
}, error = function(e) paste("RUN FAIL:", conditionMessage(e)))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (is.character(r) && grepl("^RUN FAIL", r)) stop(r, call. = FALSE)

## ---------------------------------------------------------------------------
## 5) Collect + normalise
## ---------------------------------------------------------------------------
engine_raw <- file.path(outdir, "engine_raw")
dir.create(engine_raw, recursive = TRUE, showWarnings = FALSE)
new_files <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = TRUE)
new_files <- new_files[!grepl("_engine_run$", new_files)]
collected <- character(0)
for (f in new_files) {
  dest <- file.path(engine_raw, basename(f))
  if (file.exists(dest)) next
  if (file.copy(f, dest, overwrite = TRUE)) collected <- c(collected, dest)
}
setwd(wd0)

qtn_tbl <- data.frame()
for (f in collected) {
  if (!grepl("\\.xlsx$", f)) next
  sh <- try(xlsx_sheet_names(f), silent = TRUE)
  if (inherits(sh, "try-error")) next
  for (s in sh) {
    df <- suppressWarnings(try(openxlsx::read.xlsx(f, sheet = s), silent = TRUE))
    if (inherits(df, "try-error") || is.null(df) || nrow(df) == 0) next
    if (!"Marker" %in% names(df)) next
    qtn_tbl <- rbind(qtn_tbl, normalise_qtns(df, trait_label = s))
  }
}
if (nrow(qtn_tbl)) {
  utils::write.csv(qtn_tbl, file.path(outdir, "qtn_results.csv"), row.names = FALSE)
  assoc <- data.frame(snp_id = qtn_tbl$marker, chr = qtn_tbl$chr,
                      pos = qtn_tbl$pos, p = qtn_tbl$p, lod = qtn_tbl$lod,
                      engine = "FastBiCmrMLM", stringsAsFactors = FALSE)
  utils::write.csv(assoc, file.path(outdir, "association_results.csv"),
                   row.names = FALSE)
}

summary_lines <- c(
  sprintf("engine: FastBiCmrMLM (%s)", eng_method),
  sprintf("trait: %s (binary)", trait),
  sprintf("genotype: %s", plink_prefix),
  sprintf("n = %d (%d cases / %d controls)", n_ind, sum(y == 1), sum(y == 0)),
  sprintf("search_radius: %g | svpal: %g | LODThreshold: %g",
          search_radius, svpal, lod_thr),
  sprintf("QTNs in results: %d", nrow(qtn_tbl)),
  sprintf("elapsed: %.1f s", elapsed))
writeLines(summary_lines, file.path(outdir, "summary.txt"))

write_run_info(file.path(outdir, "run_info.txt"), list(
  script = "scripts/gwas/10_case_control_qtn.R",
  engine = "FastBiCmrMLM", method = eng_method, trait = trait,
  genotype = plink_prefix, phenotype = phenotype_path,
  n_individuals = n_ind, cases = sum(y == 1), controls = sum(y == 0),
  search_radius = search_radius, svpal = svpal, lod_threshold = lod_thr,
  elapsed_sec = round(elapsed, 1), demo = demo_mode,
  note = "upstream has no LICENSE file; runtime install only, not vendored"))

cat("\n=== 10_case_control_qtn.R finished ===\n")
cat(sprintf("engine=FastBiCmrMLM method=%s cases=%d controls=%d elapsed=%.1fs\n",
            eng_method, sum(y == 1), sum(y == 0), elapsed))
cat(sprintf("QTNs: %d -> %s\n", nrow(qtn_tbl), file.path(outdir, "qtn_results.csv")))
