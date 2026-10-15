###############################################################################
# 02_atac_chip_analysis.R -- differential peak analysis for ATAC-seq / ChIP-seq
#
# Title:    Rice-skills chromatin peak analysis step
# Description:
#   Differential accessibility (ATAC-seq, DAR) / differential binding
#   (ChIP-seq, DHR) analysis from a peak count matrix:
#     CPM filtering -> CPM + log2 normalisation -> sample PCA ->
#     differential test (limma if installed, else Wilcoxon + BH) ->
#     volcano + top-peak heatmap figures.
#   Peak counts come from upstream peak calling (MACS2 etc.; see
#   docs/epigenomics_usage.md for the upstream command templates).
#
# Usage:
#   Rscript scripts/epigenomics/02_atac_chip_analysis.R \
#       --matrix demo_data/epigenome/peak_counts.csv \
#       --metadata demo_data/epigenome/metadata.csv \
#       --ref control --alt salt \
#       --outdir results/epigenome/02_dar
#
# Inputs:
#   --matrix    peak x sample count matrix CSV: peak_id, chr, start, end,
#               then samples; integer counts (raw fragment counts).
#   --metadata  sample metadata CSV (sample_id, group, ...)
#   --group-col metadata group column (default "group")
#   --ref/--alt two group labels (log2FC = alt - ref, on CPM)
#   --min-cpm   CPM threshold (default 1)
#   --min-samples  peaks kept if CPM >= min-cpm in >= this many samples
#               (default: half the samples)
#   --lf2fc     |log2FC| cutoff for the significant list (default 1)
#   --fdr       FDR cutoff (default 0.05)
#   --outdir    output directory (default results/epigenome/02_dar)
# Outputs:
#   peaks_filtered.csv     filtered + CPM-log2 matrix used for testing
#   dar_all.csv            all peaks: log2FC (CPM), stat, p, fdr
#   dar_significant.csv    significant peaks (up/down in alt)
#   pca_peaks.*            PCA of samples (group-coloured)
#   dar_volcano.*          volcano plot
#   dar_heatmap.*          top-peak heatmap
#   dar_summary.txt / run_info.txt
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
min_cpm <- suppressWarnings(as.numeric(opt(opts, "min-cpm", "1")))
min_samples <- suppressWarnings(as.integer(opt(opts, "min-samples", NA)))
lf2fc <- suppressWarnings(as.numeric(opt(opts, "lf2fc", "1")))
fdr_cut <- suppressWarnings(as.numeric(opt(opts, "fdr", "0.05")))
outdir <- opt(opts, "outdir", "results/epigenome/02_dar")

if (any(vapply(list(matrix_path, metadata_path, ref_group, alt_group), is.null, logical(1)))) {
  stop("Required: --matrix, --metadata, --ref, --alt. See --help.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills ATAC/ChIP differential peak analysis ==")

# ---------------------------------------------------------------------------
# 1) Load + filter + normalise
# ---------------------------------------------------------------------------
pdf_df <- read_csv_safe(matrix_path, "peak count matrix")
meta_cols <- intersect(c("peak_id", "chr", "start", "end"), colnames(pdf_df))
if (length(meta_cols) < 1) stop("Matrix needs a peak_id column.", call. = FALSE)
peak_ids <- as.character(pdf_df[[meta_cols[1]]])
C <- as.matrix(pdf_df[, setdiff(colnames(pdf_df), meta_cols), drop = FALSE])
suppressWarnings(storage.mode(C) <- "double")
rownames(C) <- peak_ids
C[is.na(C)] <- 0

meta <- read_pheno(metadata_path)
if (!group_col %in% colnames(meta)) stop("Metadata has no '", group_col, "' column.", call. = FALSE)
shared <- intersect(colnames(C), meta$sample_id)
if (length(shared) < 4) stop("Fewer than 4 shared samples.", call. = FALSE)
C <- C[, shared, drop = FALSE]
meta <- meta[match(shared, meta$sample_id), , drop = FALSE]
grp <- as.character(meta[[group_col]])
ref_i <- which(grp == ref_group); alt_i <- which(grp == alt_group)
if (length(ref_i) < 2 || length(alt_i) < 2) {
  stop("Need >= 2 samples per group.", call. = FALSE)
}

# CPM
lib <- colSums(C)
cpm <- sweep(C, 2, pmax(lib, 1e-9), "/") * 1e6
keep_n <- if (is.na(min_samples)) max(2, floor(length(shared) / 2)) else min_samples
keep <- rowSums(cpm >= min_cpm) >= keep_n
message(sprintf("  peaks: %d -> %d after CPM >= %g in >= %d samples",
                nrow(C), sum(keep), min_cpm, keep_n))
if (sum(keep) < 20) stop("Fewer than 20 peaks remain after filtering.", call. = FALSE)
cpm <- cpm[keep, , drop = FALSE]
lg <- log2(cpm + 1)

# ---------------------------------------------------------------------------
# 2) Differential test
# ---------------------------------------------------------------------------
use_limma <- requireNamespace("limma", quietly = TRUE)
log2fc <- rowMeans(lg[, alt_i, drop = FALSE]) - rowMeans(lg[, ref_i, drop = FALSE])
if (use_limma) {
  message("  engine: limma")
  gf <- factor(grp, levels = c(ref_group, alt_group))
  design <- stats::model.matrix(~ 0 + gf)
  colnames(design) <- c("ref", "alt")
  cont <- limma::makeContrasts(contrasts = "alt - ref", levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(lg, design), cont))
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
  stat <- tt$t; pval <- tt$P.Value; fdr <- tt$adj.P.Val
} else {
  message("  engine: Wilcoxon + BH (install limma for moderated stats)")
  stat <- rep(NA_real_, nrow(lg)); pval <- stat
  for (i in seq_len(nrow(lg))) {
    y1 <- lg[i, ref_i]; y2 <- lg[i, alt_i]
    if (all(y1 == y1[1]) && all(y2 == y2[1])) { pval[i] <- 1; next }
    tt <- suppressWarnings(tryCatch(stats::wilcox.test(y2, y1), error = function(e) NULL))
    pval[i] <- if (is.null(tt)) NA_real_ else tt$p.value
  }
  fdr <- stats::p.adjust(pval, method = "BH")
}

res <- data.frame(peak_id = rownames(lg), log2FC = log2fc, stat = stat,
                  p = pval, fdr = fdr, stringsAsFactors = FALSE)
for (cc in c("chr", "start", "end")) if (cc %in% colnames(pdf_df)) {
  res[[cc]] <- pdf_df[[cc]][keep]
}
res <- res[order(res$p), , drop = FALSE]
utils::write.csv(res, file.path(outdir, "dar_all.csv"), row.names = FALSE)

sig <- res[!is.na(res$fdr) & res$fdr <= fdr_cut & abs(res$log2FC) >= lf2fc, , drop = FALSE]
utils::write.csv(sig, file.path(outdir, "dar_significant.csv"), row.names = FALSE)
out_peaks <- data.frame(peak_id = rownames(lg), lg, check.names = FALSE)
utils::write.csv(out_peaks, file.path(outdir, "peaks_filtered.csv"), row.names = FALSE)

n_up <- sum(sig$log2FC > 0); n_dn <- sum(sig$log2FC < 0)
sum_lines <- c(
  "Rice-skills -- differential peak summary",
  sprintf("contrast        : %s vs %s", alt_group, ref_group),
  sprintf("engine          : %s", if (use_limma) "limma" else "wilcox + BH"),
  sprintf("peaks tested    : %d", nrow(res)),
  sprintf("significant     : %d (open/up %d / closed/down %d)", nrow(sig), n_up, n_dn),
  sprintf("cutoffs         : |log2FC| >= %g, FDR <= %g", lf2fc, fdr_cut)
)
writeLines(sum_lines, file.path(outdir, "dar_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

# ---------------------------------------------------------------------------
# 3) Figures
# ---------------------------------------------------------------------------
pal <- rice_palette()
col_up <- nature_semantic("up"); col_dn <- nature_semantic("down"); col_ns <- nature_semantic("not_significant")

# PCA
pc <- stats::prcomp(t(lg), center = TRUE, scale. = FALSE)
ve <- 100 * (pc$sdev^2 / sum(pc$sdev^2))
df_pc <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = grp)
p_pca <- ggplot2::ggplot(df_pc, ggplot2::aes(PC1, PC2, colour = group)) +
  ggplot2::geom_point(size = 1.7, alpha = 0.9) +
  scale_colour_groups() +
  ggplot2::labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
                title = "Peak-count PCA") +
  theme_nature()
save_nature_figure(p_pca, file.path(outdir, "pca_peaks"), width_mm = 89, height_mm = 75)

# volcano
res$call <- ifelse(!is.na(res$fdr) & res$fdr <= fdr_cut & res$log2FC >= lf2fc, "up",
             ifelse(!is.na(res$fdr) & res$fdr <= fdr_cut & res$log2FC <= -lf2fc, "down", "ns"))
p_vol <- ggplot2::ggplot(res, ggplot2::aes(log2FC, -log10(pmax(p, 1e-300)))) +
  ggplot2::geom_point(ggplot2::aes(colour = call), size = 0.7, alpha = 0.7) +
  ggplot2::scale_colour_manual(
    values = c(up = col_up, down = col_dn, ns = col_ns),
    breaks = c("up", "down", "ns"),
    labels = c(sprintf("Up (%d)", n_up), sprintf("Down (%d)", n_dn), "n.s.")) +
  ggplot2::geom_vline(xintercept = c(-lf2fc, lf2fc), linetype = "dashed",
                      colour = "grey55", linewidth = 0.3) +
  ggplot2::labs(x = expression(log[2]~"FC (CPM)"), y = expression(-log[10](italic(p))),
                colour = NULL, title = "Differential peaks") +
  theme_nature()
save_nature_figure(p_vol, file.path(outdir, "dar_volcano"), width_mm = 89, height_mm = 80)

# top-peak heatmap
top_n <- min(50, nrow(sig))
if (top_n >= 2) {
  ord <- c(ref_i, alt_i)
  Z <- lg[sig$peak_id[seq_len(top_n)], ord, drop = FALSE]
  Z <- t(scale(t(Z))); Z[is.na(Z)] <- 0
  df_hm <- data.frame(
    peak = rep(rownames(Z), times = ncol(Z)),
    sample = rep(colnames(Z), each = nrow(Z)),
    z = as.vector(Z))
  df_hm$peak <- factor(df_hm$peak, levels = rev(rownames(Z)))
  df_hm$sample <- factor(df_hm$sample, levels = colnames(Z))
  p_hm <- ggplot2::ggplot(df_hm, ggplot2::aes(sample, peak, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = diverging_palette()[1], mid = "white",
                                  high = diverging_palette()[7], midpoint = 0,
                                  name = "z", oob = scales::squish, limits = c(-2, 2)) +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("Top %d differential peaks", top_n)) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 4),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
  save_nature_figure(p_hm, file.path(outdir, "dar_heatmap"), width_mm = 120, height_mm = 110)
}

write_run_info_r(outdir, "02_atac_chip_analysis.R", params = list(
  matrix = matrix_path, metadata = metadata_path, ref = ref_group, alt = alt_group,
  min_cpm = min_cpm, n_tested = nrow(res), n_sig = nrow(sig),
  engine = if (use_limma) "limma" else "wilcox"
))
message("Done. Outputs written to: ", normalizePath(outdir))
