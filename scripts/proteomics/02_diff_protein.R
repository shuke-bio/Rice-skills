###############################################################################
# 02_diff_protein.R -- differential protein abundance analysis
#
# Title:    Rice-skills differential proteomics step
# Description:
#   Two-group differential abundance analysis on the preprocessed (log2,
#   median-normalised) protein matrix from 01_preprocess.py:
#     engine limma : moderated t-statistics (recommended; install limma)
#     engine welch : Welch t-test per protein + BH FDR (built-in fallback)
#   Optional batch covariate is supported by the limma engine
#   (design ~ 0 + group + batch).
#
# Usage:
#   Rscript scripts/proteomics/02_diff_protein.R \
#       --matrix results/proteomics/01_preprocess/proteins_filtered.csv \
#       --metadata demo_data/proteomics/metadata.csv \
#       --ref control --alt salt \
#       --outdir results/proteomics/02_diff
#
# Inputs:
#   --matrix    protein x sample matrix CSV (log2-scale, first column
#               protein_id; output of 01_preprocess.py)
#   --metadata  sample metadata CSV (first column sample_id; needs group col)
#   --group-col metadata column with the two groups (default "group")
#   --ref / --alt   the two group labels (alt contrasted against ref;
#               log2FC = mean(alt) - mean(ref))
#   --batch-col optional metadata column for batch adjustment (limma engine)
#   --lf2fc     |log2FC| cutoff for the significant list (default 1)
#   --fdr       FDR cutoff (default 0.05)
#   --outdir    output directory (default results/proteomics/02_diff)
# Outputs:
#   protein_diff_all.csv   all proteins: log2FC, stat, p, fdr, group means
#   protein_diff_sig.csv   significant subset (|log2FC| >= cutoff & fdr <= cutoff)
#   diff_summary.txt       counts up/down, engine used, cutoff info
#   run_info.txt
#
# Notes:
#   - limma install: BiocManager::install("limma") — strongly recommended
#     (moderated statistics are far more stable for small n).
#   - The Welch fallback tests each protein independently; results are
#     serviceable but less powerful than limma.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
source(file.path(dirname(script_file), "..", "utils", "genomics_utils.R"))
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
  if (isTRUE(v)) return(if (is.logical(default) && length(default) == 1L) default else TRUE); v
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
matrix_path <- opt(opts, "matrix")
metadata_path <- opt(opts, "metadata")
group_col <- opt(opts, "group-col", "group")
ref_group <- opt(opts, "ref")
alt_group <- opt(opts, "alt")
batch_col <- opt(opts, "batch-col")
lf2fc <- suppressWarnings(as.numeric(opt(opts, "lf2fc", "1")))
fdr_cut <- suppressWarnings(as.numeric(opt(opts, "fdr", "0.05")))
outdir <- opt(opts, "outdir", "results/proteomics/02_diff")

if (any(vapply(list(matrix_path, metadata_path, ref_group, alt_group), is.null, logical(1)))) {
  stop("Required: --matrix, --metadata, --ref, --alt.\n",
       "Hint: see the Usage block at the top of this script.", call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills differential proteomics ==")

# ---------------------------------------------------------------------------
# 1) Load and align
# ---------------------------------------------------------------------------
mat_df <- read_csv_safe(matrix_path, "protein matrix")
protein_ids <- as.character(mat_df[[1]])
X <- as.matrix(mat_df[, -1, drop = FALSE])
suppressWarnings(storage.mode(X) <- "double")
rownames(X) <- protein_ids

meta <- read_pheno(metadata_path)
if (!group_col %in% colnames(meta)) {
  stop("Metadata has no column '", group_col, "'. Available: ",
       paste(colnames(meta), collapse = ", "), call. = FALSE)
}
shared <- intersect(colnames(X), meta$sample_id)
if (length(shared) < 4) stop("Fewer than 4 shared samples between matrix and metadata.", call. = FALSE)
X <- X[, shared, drop = FALSE]
meta <- meta[match(shared, meta$sample_id), , drop = FALSE]
grp <- as.character(meta[[group_col]])

ref_i <- which(grp == ref_group)
alt_i <- which(grp == alt_group)
if (length(ref_i) < 2 || length(alt_i) < 2) {
  stop("Need >= 2 samples per group. Found ref='", ref_group, "': ", length(ref_i),
       ", alt='", alt_group, "': ", length(alt_i), ".", call. = FALSE)
}
message(sprintf("  proteins: %d | samples: %d (%s n=%d vs %s n=%d)",
                nrow(X), length(shared), ref_group, length(ref_i), alt_group, length(alt_i)))

# require no NA for the tested proteins (01_preprocess imputes by default)
na_rows <- rowSums(is.na(X))
usable <- na_rows == 0
if (!all(usable)) message(sprintf("  dropping %d proteins with NA (re-run 01 with --impute)", sum(!usable)))
X <- X[usable, , drop = FALSE]
if (nrow(X) < 10) stop("Fewer than 10 complete proteins remain.", call. = FALSE)

mean_ref <- rowMeans(X[, ref_i, drop = FALSE])
mean_alt <- rowMeans(X[, alt_i, drop = FALSE])
log2fc <- mean_alt - mean_ref

# ---------------------------------------------------------------------------
# 2) Statistics: limma (preferred) or Welch t-test
# ---------------------------------------------------------------------------
use_limma <- requireNamespace("limma", quietly = TRUE)
if (use_limma) {
  message("  engine: limma (moderated t)")
  group_f <- factor(grp, levels = c(ref_group, alt_group))
  if (!is.null(batch_col) && batch_col %in% colnames(meta)) {
    batch_f <- factor(as.character(meta[[batch_col]]))
    design <- stats::model.matrix(~ 0 + group_f + batch_f)
    colnames(design)[seq_len(2)] <- c("ref", "alt")
    fit <- limma::lmFit(X, design)
    cont <- limma::makeContrasts(contrasts = "alt - ref", levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  } else {
    design <- stats::model.matrix(~ 0 + group_f)
    colnames(design) <- c("ref", "alt")
    cont <- limma::makeContrasts(contrasts = "alt - ref", levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(X, design), cont))
  }
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
  stat <- tt$t; pval <- tt$P.Value; fdr <- tt$adj.P.Val
} else {
  message("  engine: Welch t-test (built-in fallback; install limma for moderated stats)")
  stat <- rep(NA_real_, nrow(X)); pval <- stat; fdr <- stat
  for (i in seq_len(nrow(X))) {
    y1 <- X[i, ref_i]; y2 <- X[i, alt_i]
    if (stats::sd(y1) == 0 && stats::sd(y2) == 0) {
      stat[i] <- 0; pval[i] <- 1
    } else {
      tt <- tryCatch(stats::t.test(y2, y1), error = function(e) NULL)
      stat[i] <- if (is.null(tt)) NA_real_ else unname(tt$statistic)
      pval[i] <- if (is.null(tt)) NA_real_ else tt$p.value
    }
  }
  fdr <- stats::p.adjust(pval, method = "BH")
}

res <- data.frame(
  protein_id = rownames(X),
  log2FC = log2fc,
  mean_ref = mean_ref, mean_alt = mean_alt,
  stat = stat, p = pval, fdr = fdr,
  engine = if (use_limma) "limma" else "welch",
  stringsAsFactors = FALSE
)
res <- res[order(res$p), , drop = FALSE]
utils::write.csv(res, file.path(outdir, "protein_diff_all.csv"), row.names = FALSE)

sig <- res[!is.na(res$fdr) & res$fdr <= fdr_cut & abs(res$log2FC) >= lf2fc, , drop = FALSE]
utils::write.csv(sig, file.path(outdir, "protein_diff_sig.csv"), row.names = FALSE)

up <- sum(sig$log2FC > 0); down <- sum(sig$log2FC < 0)
sum_lines <- c(
  "Rice-skills -- differential proteomics summary",
  sprintf("contrast          : %s vs %s (log2FC = alt - ref)", alt_group, ref_group),
  sprintf("engine            : %s", if (use_limma) "limma" else "welch t + BH"),
  sprintf("proteins tested   : %d", nrow(res)),
  sprintf("significant       : %d (up %d / down %d)", nrow(sig), up, down),
  sprintf("cutoffs           : |log2FC| >= %g, FDR <= %g", lf2fc, fdr_cut)
)
writeLines(sum_lines, file.path(outdir, "diff_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

write_run_info_r(outdir, "02_diff_protein.R", params = list(
  matrix = matrix_path, metadata = metadata_path, ref = ref_group, alt = alt_group,
  group_col = group_col, batch_col = batch_col, engine = if (use_limma) "limma" else "welch",
  n_tested = nrow(res), n_sig = nrow(sig), lf2fc = lf2fc, fdr = fdr_cut
))
message("Next step: Rscript scripts/proteomics/03_plot_proteome.R --diff ",
        file.path(outdir, "protein_diff_all.csv"), " --matrix ", matrix_path, " ...")
message("Done. Outputs written to: ", normalizePath(outdir))
