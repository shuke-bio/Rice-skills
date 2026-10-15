###############################################################################
# 02_colocalization.R -- GWAS x QTL colocalization assessment (region-level)
#
# Title:    Rice-skills colocalization step
# Description:
#   Practical region-level colocalization between GWAS association signals
#   and molecular-trait QTLs (eQTL/mQTL; input = 01_qtl_mapping.R output):
#     1. per-feature overlay plot of GWAS vs QTL -log10(p) tracks across the
#        candidate region (the classic colocalization figure);
#     2. lead-SNP distance between the two signals;
#     3. |Z|-profile correlation within the region as a colocalization
#        statistic, with a three-tier call:
#          colocalized (r >= 0.8 and lead distance <= --lead-window)
#          partial      (r >= 0.5)
#          distinct     (otherwise)
#   For publication-grade Bayesian colocalization (coloc/PP4) run the R
#   package `coloc` on the same region — this script's calls are the fast,
#   dependency-free screen that tells you which loci deserve it.
#
# Usage:
#   Rscript scripts/qtl/02_colocalization.R \
#       --gwas results/gwas/03_association/association_results.csv \
#       --eqtl results/qtl/01_map/qtl_results.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --region Chr05:26700000-26750000 \
#       --outdir results/qtl/02_coloc
#
#   # restrict to one feature of interest
#   Rscript scripts/qtl/02_colocalization.R ... --feature GENE0032
#
# Inputs:
#   --gwas        GWAS results CSV (snp_id, p; beta/se optional)
#   --eqtl        QTL results CSV (feature_id, snp_id, p; 01_qtl_mapping output)
#   --snp-info    snp_id, chr, pos
#   --region      "chr:start-end" (required: defines the candidate region)
#   --feature     optional: single feature_id (default: every feature with
#                 enough QTL SNPs inside the region)
#   --min-snps    minimum SNPs with QTL results in the region (default 5)
#   --lead-window lead-vs-lead distance for the colocalized call (default 1e5)
#   --r-strong / --r-partial  correlation cutoffs (default 0.8 / 0.5)
#   --outdir      output directory (default results/qtl/02_coloc)
# Outputs:
#   coloc_summary.csv     per feature: leads, distance, r, call, n_snps
#   region_tracks.csv     merged GWAS/QTL -log10p tracks (for custom plots)
#   coloc_tracks.*        overlay figure per feature (top 6 by r)
#   coloc_summary.txt
#   run_info.txt
#
# Notes:
#   - QTL results from 01_qtl_mapping.R only retain pairs passing the output
#     cutoff; SNPs absent there are treated as p = 1 (flattened track), which
#     is conservative for the correlation.
#   - The correlation uses |Z| = -log10(p)-based profiles; sign information
#     (effect direction) is not required for the screen.
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
gwas_path <- opt(opts, "gwas")
eqtl_path <- opt(opts, "eqtl")
snp_info_path <- opt(opts, "snp-info")
region_str <- opt(opts, "region")
feature_arg <- opt(opts, "feature")
min_snps <- suppressWarnings(as.integer(opt(opts, "min-snps", "5")))
lead_window <- suppressWarnings(as.numeric(opt(opts, "lead-window", "1e5")))
r_strong <- suppressWarnings(as.numeric(opt(opts, "r-strong", "0.8")))
r_partial <- suppressWarnings(as.numeric(opt(opts, "r-partial", "0.5")))
outdir <- opt(opts, "outdir", "results/qtl/02_coloc")

if (any(vapply(list(gwas_path, eqtl_path, snp_info_path, region_str), is.null, logical(1)))) {
  stop("Required: --gwas, --eqtl, --snp-info, --region. See --help.", call. = FALSE)
}
m <- regmatches(region_str, regexec("^(.+):(\\d+)-(\\d+)$", region_str))[[1]]
if (length(m) != 4) stop("Cannot parse --region. Expected chr:start-end.", call. = FALSE)
chr_r <- m[2]; start_r <- as.numeric(m[3]); end_r <- as.numeric(m[4])
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message(sprintf("== Rice-skills colocalization (%s:%d-%d) ==", chr_r, start_r, end_r))

# ---------------------------------------------------------------------------
# 1) Region tracks: GWAS + QTL merged on the SNP grid
# ---------------------------------------------------------------------------
info <- read_snp_info(snp_info_path)
gwas <- read_csv_safe(gwas_path, "GWAS results")
check_columns(gwas, c("snp_id", "p"), "GWAS results")
eqtl <- read_csv_safe(eqtl_path, "QTL results")
check_columns(eqtl, c("feature_id", "snp_id", "p"), "QTL results")

grid <- info[info$chr == chr_r & info$pos >= start_r & info$pos <= end_r, , drop = FALSE]
if (nrow(grid) < 10) stop("Fewer than 10 SNPs in the region.", call. = FALSE)

g_region <- gwas[match(grid$snp_id, gwas$snp_id), , drop = FALSE]
g_p <- ifelse(is.na(g_region$p), 1, g_region$p)
g_neglog10 <- -log10(pmax(g_p, 1e-300))
g_lead_i <- which.max(g_neglog10)

tracks <- data.frame(snp_id = grid$snp_id, chr = grid$chr, pos = grid$pos,
                     gwas_p = g_p, gwas_neglog10p = g_neglog10, stringsAsFactors = FALSE)

feats <- if (!is.null(feature_arg)) feature_arg else unique(eqtl$feature_id)
message(sprintf("  region SNPs: %d | GWAS lead: %s (-log10p %.2f) | features: %d",
                nrow(grid), grid$snp_id[g_lead_i], g_neglog10[g_lead_i], length(feats)))

# ---------------------------------------------------------------------------
# 2) Per-feature assessment
# ---------------------------------------------------------------------------
rows <- list(); track_list <- list(); plot_data <- list()
for (fid in feats) {
  q_f <- eqtl[eqtl$feature_id == fid, , drop = FALSE]
  if (!nrow(q_f)) next
  qi <- match(grid$snp_id, q_f$snp_id)
  q_p <- ifelse(is.na(qi), 1, q_f$p[qi])
  if (sum(q_p < 1) < min_snps) next
  q_neglog10 <- -log10(pmax(q_p, 1e-300))
  q_lead_i <- which.max(q_neglog10)

  # |Z|-style profile correlation (use -log10p tracks; flat QTL tails are p=1)
  use <- g_p < 1 | q_p < 1
  r <- if (sum(use) >= 5) suppressWarnings(stats::cor(g_neglog10[use], q_neglog10[use])) else NA_real_
  if (is.na(r)) next

  dist <- abs(grid$pos[q_lead_i] - grid$pos[g_lead_i])
  call <- if (!is.na(r) && r >= r_strong && dist <= lead_window) "colocalized"
          else if (!is.na(r) && r >= r_partial) "partial"
          else "distinct"

  rows[[length(rows) + 1L]] <- data.frame(
    feature_id = fid,
    gwas_lead = grid$snp_id[g_lead_i], gwas_lead_p = g_p[g_lead_i],
    qtl_lead = grid$snp_id[q_lead_i], qtl_lead_p = q_p[q_lead_i],
    lead_distance_bp = dist, r_profile = round(r, 3),
    n_qtl_snps = sum(q_p < 1), call = call, stringsAsFactors = FALSE)
  track_list[[length(track_list) + 1L]] <- data.frame(
    snp_id = grid$snp_id, pos = grid$pos, feature = fid, qtl_neglog10p = q_neglog10)
  plot_data[[length(plot_data) + 1L]] <- list(
    feature = fid, r = r, df = data.frame(pos = grid$pos,
                                          gwas = g_neglog10, qtl = q_neglog10),
    glead = grid$pos[g_lead_i], qlead = grid$pos[q_lead_i])
  tracks[[fid]] <- q_neglog10
}

if (!length(rows)) {
  message("No feature had >= min-snps QTL results inside the region.")
  write_run_info_r(outdir, "02_colocalization.R", params = list(
    gwas = gwas_path, eqtl = eqtl_path, region = region_str, n_features = 0))
  quit(save = "no", status = 0)
}
res <- do.call(rbind, rows)
res <- res[order(match(res$call, c("colocalized", "partial", "distinct")), -res$r_profile), , drop = FALSE]
utils::write.csv(res, file.path(outdir, "coloc_summary.csv"), row.names = FALSE)

track_df <- if (length(track_list)) do.call(rbind, track_list) else NULL
if (!is.null(track_df)) {
  tracks_long <- merge(track_df,
                       data.frame(snp_id = tracks$snp_id, gwas_neglog10p = tracks$gwas_neglog10p),
                       by = "snp_id")
  utils::write.csv(tracks_long, file.path(outdir, "region_tracks.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 3) Overlay figures (top 6 by r)
# ---------------------------------------------------------------------------
pal <- rice_palette()
plot_data <- plot_data[order(-vapply(plot_data, function(x) x$r, numeric(1)))]
for (pd in head(plot_data, 6)) {
  q_name <- paste("QTL:", pd$feature)
  df_long <- rbind(
    data.frame(pos = pd$df$pos, neglog10p = pd$df$gwas, track = "GWAS"),
    data.frame(pos = pd$df$pos, neglog10p = pd$df$qtl, track = q_name))
  track_vals <- c("GWAS" = pal$groups5[1])
  track_vals[q_name] <- pal$groups5[2]
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = pos / 1e6, y = neglog10p)) +
    ggplot2::geom_line(ggplot2::aes(colour = track), linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(colour = track), size = 0.6) +
    ggplot2::geom_vline(xintercept = pd$glead / 1e6, linetype = "dashed",
                        colour = pal$groups5[1], linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = pd$qlead / 1e6, linetype = "dashed",
                        colour = pal$groups5[2], linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = track_vals) +
    ggplot2::labs(x = sprintf("%s position (Mb)", chr_r), y = expression(-log[10](italic(p))),
                  colour = NULL,
                  title = sprintf("%s — r = %.2f", pd$feature, pd$r)) +
    theme_nature()
  save_nature_figure(p, file.path(outdir, paste0("coloc_tracks_", pd$feature)),
                     width_mm = 130, height_mm = 70)
}

n_coloc <- sum(res$call == "colocalized"); n_part <- sum(res$call == "partial")
sum_lines <- c(
  "Rice-skills -- colocalization summary",
  sprintf("region            : %s", region_str),
  sprintf("features assessed : %d", nrow(res)),
  sprintf("colocalized       : %d (r >= %g, lead distance <= %g bp)", n_coloc, r_strong, lead_window),
  sprintf("partial           : %d", n_part),
  sprintf("distinct          : %d", nrow(res) - n_coloc - n_part),
  "note              : fast screen; run Bayesian coloc for PP4 on colocalized loci"
)
writeLines(sum_lines, file.path(outdir, "coloc_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

write_run_info_r(outdir, "02_colocalization.R", params = list(
  gwas = gwas_path, eqtl = eqtl_path, region = region_str,
  n_features = nrow(res), n_colocalized = n_coloc, n_partial = n_part
))
message("Done. Outputs written to: ", normalizePath(outdir))
