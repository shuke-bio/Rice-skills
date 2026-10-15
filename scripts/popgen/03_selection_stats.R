###############################################################################
# 03_selection_stats.R -- selection scan: Hudson Fst and nucleotide diversity
#
# Title:    Rice-skills selection-statistics step
# Description:
#   Computes per-SNP Hudson Fst between two groups (ratio-of-averages
#   estimator), genome-wide Fst, windowed Fst, and windowed nucleotide
#   diversity pi per group. The genome-wide Fst scatter (top windows
#   highlighted) localises candidate regions of divergent selection
#   (e.g. indica/japonica differentiation).
#
# Usage:
#   Rscript scripts/popgen/03_selection_stats.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --phenotype demo_data/genomics/phenotype.csv --group-col group \
#       --outdir results/popgen/03_selection --window 100000
#
# Notes:
#   - Group labels come from --group-col in the phenotype table; the first two
#     levels are used unless --group1 / --group2 are given.
#   - For >2 groups run the script once per pairwise comparison.
# Outputs:
#   fst_per_snp.csv / fst_windows.csv / pi_windows.csv
#   fst_genome_scatter.{pdf,svg,tiff} / pi_plot.{pdf,svg,tiff}
#   fig_selection_overview.{pdf,svg,tiff} (patchwork)
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
genotype_path <- opt(opts, "genotype")
snp_info_path <- opt(opts, "snp-info")
phenotype_path <- opt(opts, "phenotype")
group_col <- opt(opts, "group-col", "group")
group1 <- opt(opts, "group1")
group2 <- opt(opts, "group2")
outdir <- opt(opts, "outdir", "results/popgen/03_selection")
window_bp <- suppressWarnings(as.numeric(opt(opts, "window", "100000")))
maf_cut <- suppressWarnings(as.numeric(opt(opts, "maf", "0.05")))

if (is.null(genotype_path) || is.null(snp_info_path) || is.null(phenotype_path)) {
  stop("Required: --genotype, --snp-info, --phenotype (with the group column).",
       call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
pheno <- read_pheno(phenotype_path)
if (!group_col %in% colnames(pheno)) {
  stop("Column '", group_col, "' not found in the phenotype table. Available: ",
       paste(colnames(pheno), collapse = ", "), call. = FALSE)
}
al <- align_samples_geno(gt$geno, pheno)
ct <- table(al$pheno[[group_col]])
if (length(ct) < 2) stop("The group column has fewer than 2 levels.", call. = FALSE)
levels_g <- names(ct)
g1 <- if (!is.null(group1)) group1 else levels_g[1]
g2 <- if (!is.null(group2)) group2 else levels_g[2]
if (!all(c(g1, g2) %in% levels_g)) {
  stop("Group labels not found. Available levels: ", paste(levels_g, collapse = ", "),
       call. = FALSE)
}

samples1 <- al$pheno$sample_id[al$pheno[[group_col]] == g1]
samples2 <- al$pheno$sample_id[al$pheno[[group_col]] == g2]
message(sprintf("Fst comparison: '%s' (n=%d) vs '%s' (n=%d).",
                g1, length(samples1), g2, length(samples2)))
if (min(length(samples1), length(samples2)) < 3) {
  stop("Each group needs at least 3 samples for a meaningful Fst.", call. = FALSE)
}

fl <- filter_snps(al$geno, snp_info, min_maf = maf_cut, max_missing = 0.1)
message(sprintf("SNPs after filter: %d.", ncol(fl$geno)))

# --- Hudson Fst ----------------------------------------------------------------
fst <- hudson_fst(fl$geno, samples1, samples2)
fst <- merge(fst, fl$snp_info[, c("snp_id", "chr", "pos")], by = "snp_id")
fst <- fst[order(fst$chr, fst$pos), , drop = FALSE]
gw_fst <- sum(fst$fst_num, na.rm = TRUE) / sum(fst$fst_den, na.rm = TRUE)
message(sprintf("Genome-wide Hudson Fst (%s vs %s): %.4f", g1, g2, gw_fst))
utils::write.csv(fst, file.path(outdir, "fst_per_snp.csv"), row.names = FALSE)

# windowed Fst = ratio of summed numerator/denominator per window
fst$window <- floor(fst$pos / window_bp)
fst$one <- 1
win_num <- stats::aggregate(cbind(fst_num, fst_den, one) ~ chr + window, data = fst, FUN = sum)
fst_win <- data.frame(
  chr = as.character(win_num$chr),
  start = win_num$window * window_bp,
  end = win_num$window * window_bp + window_bp - 1,
  mid = win_num$window * window_bp + window_bp / 2,
  fst = win_num$fst_num / win_num$fst_den,
  n_snp = as.integer(win_num$one)
)
fst_win <- fst_win[order(fst_win$chr, fst_win$start), , drop = FALSE]
utils::write.csv(fst_win, file.path(outdir, "fst_windows.csv"), row.names = FALSE)

# --- nucleotide diversity per group -------------------------------------------
pi1 <- site_pi(fl$geno[samples1, , drop = FALSE])
pi2 <- site_pi(fl$geno[samples2, , drop = FALSE])
pi_w1 <- window_summarise(fl$snp_info, pi1, window_bp)
pi_w2 <- window_summarise(fl$snp_info, pi2, window_bp)
pi_all <- rbind(
  data.frame(group = g1, pi_w1, stringsAsFactors = FALSE),
  data.frame(group = g2, pi_w2, stringsAsFactors = FALSE)
)
utils::write.csv(pi_all, file.path(outdir, "pi_windows.csv"), row.names = FALSE)

# --- figures -------------------------------------------------------------------
# Fst scatter: windows above the 99th percentile are candidate divergent-
# selection regions (annotated in the title with their count).
q99 <- stats::quantile(fst_win$fst, 0.99, na.rm = TRUE)
n_out <- sum(fst_win$fst >= q99, na.rm = TRUE)
scatter_cols <- rep(rice_palette()$series12[c(1, 5)], length.out = length(unique(fst_win$chr)))
names(scatter_cols) <- unique(fst_win$chr)
p_fst <- ggplot2::ggplot(fst_win, ggplot2::aes(x = mid / 1e6, y = fst, colour = chr)) +
  ggplot2::geom_point(size = 0.9, alpha = 0.8) +
  ggplot2::scale_colour_manual(values = scatter_cols, guide = "none") +
  ggplot2::geom_point(data = fst_win[fst_win$fst >= q99, , drop = FALSE],
                      colour = nature_semantic("up"), size = 1.4) +
  ggplot2::geom_hline(yintercept = gw_fst, linetype = "dashed",
                      colour = nature_semantic("reference_line"), linewidth = 0.4) +
  ggplot2::labs(x = "Chromosome position (Mb)", y = expression(F[ST]),
                title = sprintf("Fst scan (%s vs %s)", g1, g2),
                subtitle = sprintf("Genome-wide Fst = %.4f; top 1%% windows = %d",
                                   gw_fst, n_out)) +
  theme_nature()
save_nature_figure(p_fst, file.path(outdir, "fst_genome_scatter"),
                   width_mm = 180, height_mm = 80)

# pi plot: lower diversity in one group often marks a sweep in that group.
pi_all$group <- factor(pi_all$group, levels = c(g1, g2))
p_pi <- ggplot2::ggplot(pi_all, ggplot2::aes(x = mid / 1e6, y = value, colour = group)) +
  ggplot2::geom_line(linewidth = 0.5, alpha = 0.9) +
  scale_colour_groups(name = group_col) +
  ggplot2::labs(x = "Chromosome position (Mb)", y = expression(pi),
                title = "Nucleotide diversity") +
  theme_nature()
save_nature_figure(p_pi, file.path(outdir, "pi_plot"), width_mm = 180, height_mm = 70)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_ov <- patchwork::wrap_plots(p_fst, p_pi, ncol = 1, heights = c(1, 1)) +
    patchwork::plot_annotation(title = "Selection overview")
  save_nature_figure(p_ov, file.path(outdir, "fig_selection_overview"),
                     width_mm = 180, height_mm = 150)
} else {
  message("Note: patchwork not installed - overview figure skipped.")
}

write_run_info_r(
  outdir, "03_selection_stats.R",
  params = list(genotype = genotype_path, snp_info = snp_info_path,
                phenotype = phenotype_path, group_col = group_col,
                group1 = g1, group2 = g2, window = window_bp, maf = maf_cut,
                genome_wide_fst = round(gw_fst, 5), seed = RICE_SEED),
  packages = c("ggplot2", "patchwork")
)
message("Selection-statistics step complete. Results in ", outdir)
