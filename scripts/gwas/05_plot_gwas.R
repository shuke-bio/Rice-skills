###############################################################################
# 05_plot_gwas.R -- publication figures for GWAS results (Nature theme)
#
# Title:    Rice-skills GWAS plotting step
# Description:
#   Produces the classic GWAS figures from association_results.csv:
#     - Manhattan plot (Bonferroni / custom threshold line, significant points
#       in the semantic 'up' colour);
#     - QQ plot with the 95% confidence band;
#     - a combined overview figure (patchwork) when patchwork is installed.
#
# Usage:
#   Rscript scripts/gwas/05_plot_gwas.R \
#       --assoc results/gwas/03_association/association_results.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --outdir results/gwas/05_figures
#
# Outputs (each PDF + SVG + TIFF @300 dpi):
#   manhattan.{pdf,svg,tiff} / qq.{pdf,svg,tiff}
#   fig_gwas_overview.{pdf,svg,tiff}   (requires patchwork; else two files)
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
assoc_path <- opt(opts, "assoc")
snp_info_path <- opt(opts, "snp-info")
outdir <- opt(opts, "outdir", "results/gwas/05_figures")
p_threshold <- opt(opts, "p-threshold", NULL)
plot_title <- opt(opts, "title", "GWAS")
trait_label <- opt(opts, "trait", "")

if (is.null(assoc_path) || is.null(snp_info_path)) {
  stop("Required: --assoc (association_results.csv) and --snp-info.\n",
       "Hint: run scripts/gwas/03_association.R first.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

assoc <- read_csv_safe(assoc_path, "Association results")
check_columns(assoc, c("snp_id", "p"), "Association results")
snp_info <- read_snp_info(snp_info_path)

if (is.null(p_threshold)) {
  p_threshold <- 0.05 / nrow(assoc)
} else {
  p_threshold <- as.numeric(p_threshold)
}
n_sig <- sum(assoc$p <= p_threshold, na.rm = TRUE)

title_suffix <- if (nzchar(trait_label)) paste0(" - ", trait_label) else ""

# --- Manhattan ---------------------------------------------------------------
p_man <- manhattan_plot(assoc, snp_info, threshold = p_threshold,
                        title = paste0("Manhattan plot", title_suffix))
save_nature_figure(p_man, file.path(outdir, "manhattan"),
                   width_mm = 180, height_mm = 80)

# --- QQ ----------------------------------------------------------------------
p_qq <- qq_plot(assoc, title = paste0("QQ plot", title_suffix))
save_nature_figure(p_qq, file.path(outdir, "qq"), width_mm = 89, height_mm = 89)

# --- Combined overview --------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_overview <- patchwork::wrap_plots(p_man, p_qq, ncol = 1,
                                      heights = c(2, 1.6)) +
    patchwork::plot_annotation(
      title = paste0("GWAS overview", title_suffix),
      subtitle = sprintf("%d SNPs tested; %d significant (p <= %.3g)",
                         nrow(assoc), n_sig, p_threshold)
    )
  save_nature_figure(p_overview, file.path(outdir, "fig_gwas_overview"),
                     width_mm = 180, height_mm = 170)
  message("Combined overview figure saved (patchwork).")
} else {
  message("Note: patchwork not installed - overview figure skipped ",
          "(manhattan and qq files are still written). ",
          "Install with install.packages(\"patchwork\").")
}

write_run_info_r(
  outdir, "05_plot_gwas.R",
  params = list(assoc = assoc_path, snp_info = snp_info_path,
                p_threshold = p_threshold, n_sig = n_sig,
                title = plot_title, seed = RICE_SEED),
  packages = c("ggplot2", "patchwork")
)
message(sprintf("GWAS figures complete: %d significant SNPs at p <= %.3g. Results in %s",
                n_sig, p_threshold, outdir))
