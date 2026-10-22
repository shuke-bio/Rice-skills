###############################################################################
# 02_ld_decay.R -- linkage-disequilibrium decay curve
#
# Title:    Rice-skills LD-decay step
# Description:
#   Computes r^2 between SNP pairs binned by physical distance and draws the
#   classic LD-decay curve. The distance at which the fitted curve crosses
#   r^2 = 0.2 (configurable) is reported - it sets the plausible candidate-gene
#   search window for GWAS.
#
# Usage:
#   Rscript scripts/popgen/02_ld_decay.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --outdir results/popgen/02_ld --max-dist 500000
#
# Outputs:
#   ld_decay.csv, run_info.txt
#   ld_decay_curve.{pdf,tiff,png}
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
outdir <- opt(opts, "outdir", "results/popgen/02_ld")
max_dist <- suppressWarnings(as.numeric(opt(opts, "max-dist", "500000")))
maf_cut <- suppressWarnings(as.numeric(opt(opts, "maf", "0.05")))
r2_threshold <- suppressWarnings(as.numeric(opt(opts, "r2-threshold", "0.2")))

if (is.null(genotype_path) || is.null(snp_info_path)) {
  stop("Required: --genotype and --snp-info.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
fl <- filter_snps(gt$geno, snp_info, min_maf = maf_cut, max_missing = 0.1)
message(sprintf("LD input: %d samples x %d SNPs.", nrow(fl$geno), ncol(fl$geno)))

ld <- ld_decay(fl$geno, fl$snp_info, max_dist_bp = max_dist)
if (nrow(ld) == 0) {
  stop("No SNP pairs within --max-dist. Hint: increase --max-dist or check ",
       "that snp_info positions are in bp.", call. = FALSE)
}
utils::write.csv(ld, file.path(outdir, "ld_decay.csv"), row.names = FALSE)

# Half-decay distance (r^2 crossing the threshold from a loess fit).
half_decay_bp <- NA_real_
fit <- tryCatch(stats::loess(mean_r2 ~ dist, data = ld, span = 0.9), error = function(e) NULL)
if (!is.null(fit)) {
  grid <- seq(min(ld$dist), max(ld$dist), length.out = 400)
  pred <- stats::predict(fit, newdata = data.frame(dist = grid))
  below <- which(pred <= r2_threshold)
  if (length(below) > 0) half_decay_bp <- grid[below[1]]
}
msg <- if (is.na(half_decay_bp)) {
  sprintf("LD did not decay below r2 = %g within %s bp (fitting range too short).",
          r2_threshold, format(max_dist, big.mark = ""))
} else {
  sprintf("LD decays to r2 = %g at ~%s bp (candidate-search window guide).",
          r2_threshold, format(round(half_decay_bp / 1000) * 1000, big.mark = ""))
}
message(msg)

# Figure: decay curve with the r2 threshold and half-decay marker.
p <- ggplot2::ggplot(ld, ggplot2::aes(x = dist / 1000, y = mean_r2)) +
  ggplot2::geom_point(size = 1.1, colour = rice_palette()$series12[1], alpha = 0.8) +
  {if (!is.null(fit)) ggplot2::geom_smooth(method = "loess", formula = "y ~ x",
                                           span = 0.9, se = FALSE,
                                           colour = rice_palette()$series12[2],
                                           linewidth = 0.7)} +
  ggplot2::geom_hline(yintercept = r2_threshold, linetype = "dashed",
                      colour = nature_semantic("reference_line"), linewidth = 0.4) +
  {if (!is.na(half_decay_bp)) ggplot2::geom_vline(xintercept = half_decay_bp / 1000,
                                                  linetype = "dotted",
                                                  colour = nature_semantic("highlight"),
                                                  linewidth = 0.5)} +
  ggplot2::labs(x = "Physical distance (kb)", y = expression(Mean~r^2),
                title = "LD decay", subtitle = msg) +
  theme_nature()
save_nature_figure(p, file.path(outdir, "ld_decay_curve"), width_mm = 110, height_mm = 85)

write_run_info_r(
  outdir, "02_ld_decay.R",
  params = list(genotype = genotype_path, snp_info = snp_info_path,
                max_dist = max_dist, maf = maf_cut, r2_threshold = r2_threshold,
                half_decay_bp = half_decay_bp, seed = RICE_SEED),
  packages = c("ggplot2")
)
message("LD-decay step complete. Results in ", outdir)
