###############################################################################
# 11_multiparent_qtl.R -- multi-parent population QTL mapping (mppQTL)
#
# Title:    Rice-skills multi-parent population QTL mapping (mppQTL)
# Description:
#   Wraps mppQTL, a multi-locus QTL method for multi-parent populations
#   (MAGIC / NAM / ROAM): step 1 scans every marker with a single-locus
#   linear mixed model (founder design probabilities via qtl2 machinery),
#   selects peaks with group-lasso; step 2 fits a multi-locus model on the
#   candidate loci, estimates effects by EM empirical Bayes, and calls
#   significant QTLs by likelihood-ratio test -- a small-effect-friendly
#   complement to biparental R/QTL workflows.
#
#   Upstream: Li G, Zhou YH, Li HF, Zhang YM (2023) Computational and
#   Structural Biotechnology Journal 21: 2241-2252. Software
#   github.com/YuanmingZhang65/IIIVmrMLM (mppQTL_1.0.0.zip) and
#   github.com/YuanmingZhang65/mppQTL -- GPL-3; NEVER vendored, runtime
#   install only.
#
# Usage:
#   Rscript scripts/gwas/11_multiparent_qtl.R \
#       --population magic \
#       --genotype data/4way_magic_genotype.csv \
#       --phenotype data/4way_magic_phenotype.csv \
#       --map data/4way_magic_map.csv \
#       --cross-info data/4way_magic_cross.csv \
#       --crosstype riself4 --map-function haldane \
#       --cores 4 --outdir results/gwas/11_multiparent
#
# Inputs (native mppQTL CSV formats):
#   --population  magic | nam | roam
#   --genotype    magic/roam: single genotype CSV (rows = markers, columns =
#                 founder allele letters then one allelic letter per line).
#                 nam: a FOLDER with one <family>_genotype.csv per family.
#   --phenotype   CSV: rows = Trait x Environment, columns = lines
#                 (numeric); first two columns are Trait and Environment.
#   --map         CSV: Marker, Chromosome, Genetic_Position_cM,
#                 Physical_Position_bp (genetic position optional for NAM/ROAM)
#   --cross-info  magic/roam only: CSV with ID + one founder index column per
#                 founder (e.g. founder_A..founder_D holding 1..4)
#   --crosstype   riself4/8/16, risib4/8/16, magic19, nam5/26, roam14, ...
#   --map-function haldane | kosambi | c-f | morgan (MAGIC/ROAM)
#   --p-value     final LRT P threshold (default: 0.05 / s, s = candidate QTLs)
#   --cores       CPU cores for the inner parallel scan
#   --outdir      output directory
#   --demo        synthetic MAGIC (riself4) end-to-end self-test
#
# Outputs (in --outdir):
#   qtl_results.csv     per trait x environment significant QTL table
#   engine_raw/         mppQTL-native workbooks (Trait_Env*.xlsx) + Recode/
#   run_info.txt / summary.txt
#
# Notes:
#   - mppQTL is GPL-3 and not bundled; install at runtime from the official
#     zip (see docs/references/3vmrmlm_multilocus_qtn.md section on mppQTL).
#   - Requires the qtl2 package (CRAN) and its dependencies.
#   - --demo uses the same 4-founder selfing MAGIC layout as the upstream
#     example, generated synthetically with a planted QTL.
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
population <- tolower(opt(opts, "population", "magic"))
genotype_path <- opt(opts, "genotype")
phenotype_path <- opt(opts, "phenotype")
map_path <- opt(opts, "map")
cross_path <- opt(opts, "cross-info", NULL)
crosstype <- opt(opts, "crosstype", NULL)
map_function <- opt(opts, "map-function", NULL)
p_value <- suppressWarnings(as.numeric(opt(opts, "p-value", "NA")))
cores <- suppressWarnings(as.integer(opt(opts, "cores", "2")))
outdir <- opt(opts, "outdir", "results/gwas/11_multiparent")
demo_mode <- isTRUE(opts[["demo"]])

if (!population %in% c("magic", "nam", "roam")) {
  stop("--population must be magic / nam / roam; got '", population, "'.", call. = FALSE)
}
need_pkgs("jsonlite")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir)

chk <- engine_ready("mppQTL",
  paste("Install (runtime, GPL-3):",
        "download.file('https://github.com/YuanmingZhang65/IIIVmrMLM/raw/main/mppQTL_1.0.0.zip',",
        "  'mppQTL.zip', mode='wb'); install.packages('mppQTL.zip', repos=NULL)",
        "Requires CRAN package qtl2. Cite Li et al. 2023, CSBJ 21: 2241-2252.",
        sep = "\n  "))
if (!chk$ok) stop(chk$message, call. = FALSE)
if (!requireNamespace("qtl2", quietly = TRUE)) {
  stop("mppQTL depends on the CRAN package 'qtl2': install.packages('qtl2').",
       call. = FALSE)
}

## ---------------------------------------------------------------------------
## 1) Demo data: synthetic 4-founder MAGIC (riself4), 3 chr x 75 markers,
##    60 lines, 2 traits x 2 environments, planted QTL on chr 2
## ---------------------------------------------------------------------------
demo_dir <- NULL
if (demo_mode) {
  message("== DEMO MODE: synthetic 4-way MAGIC (selfing) population ==")
  demo_dir <- file.path(outdir, "_demo_work")
  dir.create(demo_dir, showWarnings = FALSE)
  set.seed(42)
  founders <- LETTERS[1:4]
  n_chr <- 4L; n_marker <- 60L; n_line <- 150L
  # NOTE: mppQTL's recode_MAGIC recognises only A/T/C/G nucleotide letters;
  # any other symbol (e.g. "B") is not recoded and breaks the design matrix.
  allele_letters <- c("A", "T", "C", "G")
  markers <- data.frame(
    Marker = paste0("D", rep(1:n_chr, each = n_marker), "M", 1:n_marker),
    Chromosome = rep(1:n_chr, each = n_marker),
    Genetic_Position_cM = 0,
    Physical_Position_bp = 0,
    stringsAsFactors = FALSE)
  # increasing physical positions per chromosome (qtl2 needs real positions
  # to insert pseudo-markers and compute recombination fractions)
  markers$Physical_Position_bp <- ave(markers$Chromosome, markers$Chromosome,
                                      FUN = function(x) seq_along(x) * 50000)
  markers$Genetic_Position_cM <- markers$Physical_Position_bp / 1e4   # 1 cM / Mb
  # each marker carries two allele letters; every founder gets one of them.
  # Guarantee polymorphism (both alleles present among founders): a marker
  # where all four founders share one allele yields a degenerate design
  # matrix and crashes qtl2's genoprob step.
  pair_letters <- t(replicate(nrow(markers), sample(allele_letters, 2, replace = FALSE)))
  f_alleles <- t(apply(pair_letters, 1, function(r) {
    k <- sample(2:3, 1)                       # 2/2 or 3/1 split
    sample(c(rep(r[1], k), rep(r[2], 4 - k)))
  }))                                          # m x 4
  colnames(f_alleles) <- paste0("founder_", founders)
  # each line: founder order per chromosome block (1..4). qtl2's
  # riself4 cross_info expects each row to be a PERMUTATION of 1:n_founder
  # (the crossing order), not an arbitrary assignment.
  line_ids <- paste0("Line", seq_len(n_line))
  cross_perm <- t(replicate(n_line, sample(1:4)))
  cross <- data.frame(ID = line_ids,
                      setNames(as.data.frame(cross_perm),
                               paste0("founder_", founders)),
                      check.names = FALSE)
  # line alleles: markers (rows) x lines (cols). Walk each chromosome with
  # occasional ancestry switches (recombination) so the planted QTL is
  # localised rather than smeared over a whole founder block.
  line_alleles <- sapply(seq_len(n_line), function(i) {
    founder_idx <- integer(nrow(markers))
    for (k in seq_len(n_chr)) {
      rows <- which(markers$Chromosome == k)
      cur <- cross[i, k + 1]
      for (r in rows) {
        if (runif(1) < 0.02) cur <- sample(1:4, 1)
        founder_idx[r] <- cur
      }
    }
    mapply(function(r, f) f_alleles[r, f], seq_len(nrow(markers)), founder_idx)
  })
  colnames(line_alleles) <- line_ids
  geno_df <- cbind(markers[, "Marker", drop = FALSE], f_alleles, line_alleles)
  genotype_path <- file.path(demo_dir, "demo_magic_genotype.csv")
  utils::write.csv(geno_df, genotype_path, row.names = FALSE, quote = FALSE)
  map_path <- file.path(demo_dir, "demo_magic_map.csv")
  utils::write.csv(markers, map_path, row.names = FALSE, quote = FALSE)
  cross_path <- file.path(demo_dir, "demo_magic_cross.csv")
  utils::write.csv(cross, cross_path, row.names = FALSE, quote = FALSE)
  # phenotype: Trait x Environment rows x lines; chr2 marker D2M30 is the QTL.
  # line_alleles is (markers x lines); the QTL ROW gives the allele per line.
  qtl_marker <- "D2M30"; qi <- which(markers$Marker == qtl_marker)
  qtl_effect <- sample(c(-1, 1), 1) * 2.5
  # effect allele = the first of the QTL marker's two founder alleles
  qtl_pair <- as.character(unlist(f_alleles[qi, ]))
  carrier <- line_alleles[qi, ] == qtl_pair[1]
  lin <- scale(as.numeric(carrier)) * qtl_effect + rnorm(n_line, 0, 1)
  phe <- data.frame(
    Trait = rep("Trait1", 2), Environment = c("Env_1", "Env_2"),
    stringsAsFactors = FALSE)
  for (j in seq_len(n_line)) phe[[line_ids[j]]] <- c(lin[j] - 2, lin[j] + 2)
  # Trait2 also carries a weaker signal at the same locus: a purely noisy
  # trait yields an empty QTL table, and mppQTL's internal write.xlsx()
  # fails on NULL results (upstream limitation).
  phe2 <- phe
  phe2$Trait <- "Trait2"
  for (j in seq_len(n_line)) phe2[[line_ids[j]]] <- 50 + 0.6 * lin[j] + rnorm(2, 0, 2)
  phe <- rbind(phe, phe2)
  phenotype_path <- file.path(demo_dir, "demo_magic_phenotype.csv")
  utils::write.csv(phe, phenotype_path, row.names = FALSE, quote = FALSE)
  crosstype <- "riself4"; map_function <- "haldane"
  message(sprintf("  synthetic: %d lines x %d markers, %d chr; planted QTL: %s",
                  n_line, nrow(markers), n_chr, qtl_marker))
} else {
  for (p in list(genotype_path, phenotype_path, map_path)) {
    if (is.null(p) || !file.exists(p)) {
      stop("Required existing inputs: --genotype, --phenotype, --map (or --demo).",
           call. = FALSE)
    }
  }
  if (is.null(crosstype)) {
    stop("--crosstype is required (e.g. riself4 / risib8 / magic19 / nam5).",
         call. = FALSE)
  }
  if (population %in% c("magic", "roam")) {
    if (is.null(cross_path) || !file.exists(cross_path)) {
      stop(population, " populations require --cross-info (founder composition file).",
           call. = FALSE)
    }
    if (is.null(map_function)) map_function <- "haldane"
  }
}

## ---------------------------------------------------------------------------
## 2) Run mppQTL
## ---------------------------------------------------------------------------
run_dir <- file.path(outdir, "_engine_run")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
wd0 <- getwd(); on.exit(setwd(wd0), add = TRUE); setwd(run_dir)
abs <- function(p) gsub("\\\\", "/", normalizePath(as.character(p), mustWork = TRUE))
geno_abs <- abs(genotype_path); phe_abs <- abs(phenotype_path); map_abs <- abs(map_path)
cross_abs <- if (is.null(cross_path)) NULL else abs(cross_path)

t0 <- Sys.time()
r <- tryCatch({
  suppressMessages(library(mppQTL))
  # NOTE: sea() silently fails when P is NA (it expects NULL for the
  # default 0.05/s threshold); never pass NA through.
  mppQTL(genotype = geno_abs, phenotype = phe_abs, map = map_abs,
         crosstype = crosstype, cross_info = cross_abs,
         P = if (is.na(p_value)) NULL else p_value,
         cores = if (is.na(cores)) NULL else cores,
         dire = "./", itera = 200, eps = 1e-8,
         map_function = map_function, step = NULL)
  "ok"
}, error = function(e) paste("RUN FAIL:", conditionMessage(e)))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (is.character(r) && grepl("^RUN FAIL", r)) stop(r, call. = FALSE)

## ---------------------------------------------------------------------------
## 3) Collect + normalise
## ---------------------------------------------------------------------------
engine_raw <- file.path(outdir, "engine_raw")
dir.create(engine_raw, recursive = TRUE, showWarnings = FALSE)
new_files <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = TRUE)
for (f in new_files) {
  dest <- file.path(engine_raw, basename(f))
  if (!file.exists(dest)) file.copy(f, dest, overwrite = TRUE)
}
setwd(wd0)

# engine writes Trait_Env*.xlsx per trait x environment, each with a "GWS"
# (genome-wide scan) and a "QTL" (significant QTL) sheet. Keep them separate:
# qtl_results.csv = significant QTLs; gws_scan.csv = full scan.
qtl_tbl <- data.frame(); gws_tbl <- data.frame()
for (f in list.files(engine_raw, pattern = "\\.xlsx$", full.names = TRUE)) {
  tag <- tools::file_path_sans_ext(basename(f))
  sh <- try(xlsx_sheet_names(f), silent = TRUE)
  if (inherits(sh, "try-error")) next
  for (s in sh) {
    df <- suppressWarnings(try(openxlsx::read.xlsx(f, sheet = s), silent = TRUE))
    if (inherits(df, "try-error") || is.null(df) || nrow(df) == 0) next
    if (!"Marker" %in% names(df)) next
    pick <- function(...) {
      cand <- c(...); hit <- cand[cand %in% names(df)]
      if (length(hit) == 0) return(rep(NA, nrow(df)))
      df[[hit[1]]]
    }
    row <- data.frame(
      trait_env = tag, sheet = s, marker = df$Marker,
      chr = pick("Chromosome", "chr"),
      pos_bp = pick("Physical_Position_bp", "Position.(bp)", "pos"),
      pos_cm = pick("Genetic_Position_cM"),
      lod = pick("LOD", "LOD.score"),
      p = pick("minus_log10_p", "P-value", "P.value"),
      stringsAsFactors = FALSE)
    if (toupper(s) == "QTL") qtl_tbl <- rbind(qtl_tbl, row) else gws_tbl <- rbind(gws_tbl, row)
  }
}
if (nrow(qtl_tbl)) utils::write.csv(qtl_tbl, file.path(outdir, "qtl_results.csv"), row.names = FALSE)
if (nrow(gws_tbl)) utils::write.csv(gws_tbl, file.path(outdir, "gws_scan.csv"), row.names = FALSE)

summary_lines <- c(
  sprintf("population: %s (%s)", population, crosstype),
  sprintf("genotype: %s", genotype_path),
  sprintf("markers: %s", map_path),
  sprintf("engines outputs: %d workbook(s)", length(list.files(engine_raw, pattern = "\\.xlsx$"))),
  sprintf("elapsed: %.1f s", elapsed),
  "Significant QTLs are called by likelihood-ratio test inside mppQTL;",
  "P (final threshold) defaults to 0.05 / number of candidate loci.")
writeLines(summary_lines, file.path(outdir, "summary.txt"))

write_run_info(file.path(outdir, "run_info.txt"), list(
  script = "scripts/gwas/11_multiparent_qtl.R",
  engine = "mppQTL (GPL-3, runtime install)", population = crosstype,
  genotype = genotype_path, phenotype = phenotype_path, map = map_path,
  cross_info = if (is.null(cross_path)) "none" else cross_path,
  map_function = if (is.null(map_function)) "none" else map_function,
  elapsed_sec = round(elapsed, 1), demo = demo_mode))

cat("\n=== 11_multiparent_qtl.R finished ===\n")
cat(sprintf("crosstype=%s elapsed=%.1fs\n", crosstype, elapsed))
cat(sprintf("outputs: %s\n", outdir))
