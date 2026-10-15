###############################################################################
# 01_methylation_analysis.R -- region-level DNA methylation analysis (DMR)
#
# Title:    Rice-skills methylation analysis step
# Description:
#   Differential methylation analysis on a region-level methylation matrix
#   (WGBS-derived methylation levels per region per sample):
#     - global methylation level QC (per-sample distributions);
#     - two-group differential methylation regions (DMR) with t-test /
#       Wilcoxon + BH FDR and a |delta| cutoff;
#     - volcano + DMR heatmap figures (shared Nature figure spec).
#   The region matrix can be produced upstream by any WGBS pipeline
#   (Bismark CX_report -> region aggregation; see docs/epigenomics_usage.md),
#   or from public rice methylome datasets (e.g. RiceENCODE).
#
# Usage:
#   Rscript scripts/epigenomics/01_methylation_analysis.R \
#       --matrix demo_data/epigenome/methylation_regions.csv \
#       --metadata demo_data/epigenome/metadata.csv \
#       --ref control --alt salt \
#       --outdir results/epigenome/01_dmr
#
# Inputs:
#   --matrix    region x sample methylation matrix CSV: region_id, chr, start,
#               end [, context], then samples; values = methylation fraction
#               (0-1) or percentage (0-100; auto-detected).
#   --metadata  sample metadata CSV (sample_id, group, ...)
#   --group-col metadata group column (default "group")
#   --ref/--alt two group labels (delta = alt - ref)
#   --delta     |delta| cutoff for the DMR list (default 0.1, fraction scale)
#   --fdr       FDR cutoff (default 0.05)
#   --test      t | wilcox (default t)
#   --outdir    output directory (default results/epigenome/01_dmr)
# Outputs:
#   dmr_all.csv       all regions: mean_ref, mean_alt, delta, p, fdr
#   dmr_significant.csv  DMR list (|delta| >= cutoff & fdr <= cutoff)
#   methyl_summary.txt   global levels + DMR counts (hyper/hypo)
#   meth_distribution.* methylation level by sample (boxplot, group-coloured)
#   dmr_volcano.*        delta vs -log10(p)
#   dmr_heatmap.*        top DMR heatmap (z-scored delta matrix)
#   run_info.txt
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
delta_cut <- suppressWarnings(as.numeric(opt(opts, "delta", "0.1")))
fdr_cut <- suppressWarnings(as.numeric(opt(opts, "fdr", "0.05")))
test_method <- tolower(opt(opts, "test", "t"))
outdir <- opt(opts, "outdir", "results/epigenome/01_dmr")

if (any(vapply(list(matrix_path, metadata_path, ref_group, alt_group), is.null, logical(1)))) {
  stop("Required: --matrix, --metadata, --ref, --alt. See --help.", call. = FALSE)
}
if (!test_method %in% c("t", "wilcox")) stop("--test must be t | wilcox.", call. = FALSE)
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills methylation (DMR) analysis ==")

# ---------------------------------------------------------------------------
# 1) Load
# ---------------------------------------------------------------------------
mdf <- read_csv_safe(matrix_path, "methylation matrix")
meta_cols <- intersect(c("region_id", "chr", "start", "end", "context"), colnames(mdf))
if (length(meta_cols) < 1) stop("Matrix needs a region_id column.", call. = FALSE)
region_ids <- as.character(mdf[[meta_cols[1]]])
X <- as.matrix(mdf[, setdiff(colnames(mdf), meta_cols), drop = FALSE])
suppressWarnings(storage.mode(X) <- "double")
rownames(X) <- region_ids

# scale: fraction (0-1) or percent (0-100)
scale_pct <- FALSE
if (max(X, na.rm = TRUE) > 1.5) {
  scale_pct <- TRUE
  message("  percentage-scale values detected -> converted to fraction (0-1)")
  X <- X / 100
  delta_cut <- delta_cut  # cutoffs stay on the fraction scale by convention
}

meta <- read_pheno(metadata_path)
if (!group_col %in% colnames(meta)) stop("Metadata has no '", group_col, "' column.", call. = FALSE)
shared <- intersect(colnames(X), meta$sample_id)
if (length(shared) < 4) stop("Fewer than 4 shared samples.", call. = FALSE)
X <- X[, shared, drop = FALSE]
meta <- meta[match(shared, meta$sample_id), , drop = FALSE]
grp <- as.character(meta[[group_col]])
ref_i <- which(grp == ref_group); alt_i <- which(grp == alt_group)
if (length(ref_i) < 2 || length(alt_i) < 2) {
  stop("Need >= 2 samples per group (ref: ", length(ref_i), ", alt: ", length(alt_i), ").", call. = FALSE)
}
message(sprintf("  regions: %d | samples: %d (%s n=%d vs %s n=%d)",
                nrow(X), length(shared), ref_group, length(ref_i), alt_group, length(alt_i)))

# ---------------------------------------------------------------------------
# 2) Per-region statistics
# ---------------------------------------------------------------------------
mean_ref <- rowMeans(X[, ref_i, drop = FALSE], na.rm = TRUE)
mean_alt <- rowMeans(X[, alt_i, drop = FALSE], na.rm = TRUE)
delta <- mean_alt - mean_ref

pval <- rep(NA_real_, nrow(X))
for (i in seq_len(nrow(X))) {
  y1 <- X[i, ref_i]; y2 <- X[i, alt_i]
  y1 <- y1[!is.na(y1)]; y2 <- y2[!is.na(y2)]
  if (length(y1) < 2 || length(y2) < 2) next
  pval[i] <- tryCatch({
    if (test_method == "t") stats::t.test(y2, y1)$p.value
    else suppressWarnings(stats::wilcox.test(y2, y1)$p.value)
  }, error = function(e) NA_real_)
}
fdr <- stats::p.adjust(pval, method = "BH")
coverage <- rowMeans(!is.na(X))

res <- data.frame(
  region_id = region_ids, mean_ref = mean_ref, mean_alt = mean_alt,
  delta = delta, p = pval, fdr = fdr, coverage = round(coverage, 3),
  stringsAsFactors = FALSE
)
# attach coordinates when present
for (cc in c("chr", "start", "end", "context")) {
  if (cc %in% colnames(mdf)) res[[cc]] <- mdf[[cc]]
}
res <- res[order(res$p), , drop = FALSE]
utils::write.csv(res, file.path(outdir, "dmr_all.csv"), row.names = FALSE)

sig <- res[!is.na(res$fdr) & res$fdr <= fdr_cut & abs(res$delta) >= delta_cut, , drop = FALSE]
utils::write.csv(sig, file.path(outdir, "dmr_significant.csv"), row.names = FALSE)
n_hyper <- sum(sig$delta > 0); n_hypo <- sum(sig$delta < 0)

global <- colMeans(X, na.rm = TRUE)
sum_lines <- c(
  "Rice-skills -- methylation summary",
  sprintf("groups            : %s vs %s (delta = alt - ref)", alt_group, ref_group),
  sprintf("regions tested    : %d", nrow(res)),
  sprintf("global methylation: %s mean %.3f", ref_group, mean(global[grp == ref_group])),
  sprintf("                    %s mean %.3f", alt_group, mean(global[grp == alt_group])),
  sprintf("DMRs (%s): %d total | hyper %d / hypo %d",
          sprintf("|delta| >= %g, FDR <= %g", delta_cut, fdr_cut), nrow(sig), n_hyper, n_hypo)
)
writeLines(sum_lines, file.path(outdir, "methyl_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

# ---------------------------------------------------------------------------
# 3) Figures
# ---------------------------------------------------------------------------
pal <- rice_palette()
col_hyper <- nature_semantic("up"); col_hypo <- nature_semantic("down")
col_ns <- nature_semantic("not_significant")

# (a) per-sample methylation distribution
df_dist <- data.frame(
  sample = rep(colnames(X), each = nrow(X)),
  meth = as.vector(X),
  group = rep(grp, each = nrow(X))
)
df_dist$sample <- factor(df_dist$sample, levels = colnames(X))
p_dist <- ggplot2::ggplot(df_dist, ggplot2::aes(x = sample, y = meth, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = 0.15, linewidth = 0.25, colour = "grey30") +
  scale_fill_groups() +
  ggplot2::labs(x = NULL, y = "Methylation level (fraction)", title = "Per-sample methylation") +
  theme_nature() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 4.5))
save_nature_figure(p_dist, file.path(outdir, "meth_distribution"), width_mm = 120, height_mm = 70)

# (b) volcano
res$call <- ifelse(!is.na(res$fdr) & res$fdr <= fdr_cut & res$delta >= delta_cut, "hyper",
             ifelse(!is.na(res$fdr) & res$fdr <= fdr_cut & res$delta <= -delta_cut, "hypo", "ns"))
p_vol <- ggplot2::ggplot(res, ggplot2::aes(x = delta, y = -log10(pmax(p, 1e-300)))) +
  ggplot2::geom_point(ggplot2::aes(colour = call), size = 0.7, alpha = 0.7) +
  ggplot2::scale_colour_manual(
    values = c(hyper = col_hyper, hypo = col_hypo, ns = col_ns),
    labels = c(sprintf("Hyper (%d)", n_hyper), sprintf("Hypo (%d)", n_hypo), "n.s."),
    breaks = c("hyper", "hypo", "ns")) +
  ggplot2::geom_vline(xintercept = c(-delta_cut, delta_cut), linetype = "dashed",
                      colour = "grey55", linewidth = 0.3) +
  ggplot2::labs(x = expression(Delta~"methylation (alt - ref)"),
                y = expression(-log[10](italic(p))), colour = NULL,
                title = "Differential methylation regions") +
  theme_nature()
save_nature_figure(p_vol, file.path(outdir, "dmr_volcano"), width_mm = 89, height_mm = 80)

# (c) top-DMR heatmap (delta, samples ordered by group)
top_n <- min(50, nrow(sig))
if (top_n >= 2) {
  top_res <- head(sig, top_n)
  ord <- c(ref_i, alt_i)
  Z <- X[top_res$region_id, ord, drop = FALSE]
  df_hm <- data.frame(
    region = rep(rownames(Z), times = ncol(Z)),
    sample = rep(colnames(Z), each = nrow(Z)),
    meth = as.vector(Z),
    group = rep(grp[ord], each = nrow(Z))
  )
  df_hm$region <- factor(df_hm$region, levels = rev(sig$region_id[seq_len(top_n)]))
  df_hm$sample <- factor(df_hm$sample, levels = colnames(Z))
  p_hm <- ggplot2::ggplot(df_hm, ggplot2::aes(x = sample, y = region, fill = meth)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = pal$sequential6[length(pal$sequential6)],
                                 name = "meth", limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf("Top %d DMR methylation", top_n)) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 4),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
  save_nature_figure(p_hm, file.path(outdir, "dmr_heatmap"), width_mm = 120, height_mm = 110)
}

write_run_info_r(outdir, "01_methylation_analysis.R", params = list(
  matrix = matrix_path, metadata = metadata_path, ref = ref_group, alt = alt_group,
  test = test_method, delta_cut = delta_cut, fdr_cut = fdr_cut,
  n_dmr = nrow(sig), n_hyper = n_hyper, n_hypo = n_hypo
))
message("Done. Outputs written to: ", normalizePath(outdir))
