#!/usr/bin/env Rscript
###############################################################################
# 02_spectral_indices.R -- vegetation/water spectral index computation
#
# Title:    Rice-skills spectral phenotyping step
# Description:
#   Computes standard vegetation and water indices from reflectance spectra
#   stored as a wide table (sample_id + wavelength columns). Wavelength
#   columns are auto-detected from name patterns like "R670", "X800",
#   "800nm" or plain "800"; specific bands can be overridden.
#   Also produces index correlation and distribution figures.
#
# Usage:
#   Rscript scripts/phenomics/02_spectral_indices.R \
#       --spectra demo_data/phenomics/spectra.csv \
#       --outdir results/phenomics/02_indices
#
# Inputs:
#   --spectra   CSV: first column sample_id; remaining columns = reflectance
#               at wavelengths (column names must contain the wavelength).
#   --nir/--red/--blue/--re   override band column names (default: auto by
#               wavelength 800 / 670 / 470 / 730 nm; nearest available +/- 15 nm)
#   --indices   comma-separated subset of: NDVI,EVI,SAVI,RVI,NDRE,WDRVI,GNDVI
#               (default: all)
#   --outdir    output directory (default results/phenomics/02_indices)
# Outputs:
#   spectral_indices.csv   sample_id + one column per computed index
#   index_correlation.*    correlation heatmap of indices
#   index_distribution.*   distribution of each index
#   run_info.txt
#
# Notes:
#   - Formulas (R = red, NIR, B = blue, RE = red-edge):
#       NDVI  = (NIR - R) / (NIR + R)
#       EVI   = 2.5 * (NIR - R) / (NIR + 6R - 7.5B + 1)
#       SAVI  = 1.5 * (NIR - R) / (NIR + R + 0.5)
#       RVI   = NIR / R
#       NDRE  = (NIR - RE) / (NIR + RE)
#       WDRVI = (0.1*NIR - R) / (0.1*NIR + R)
#       GNDVI = (NIR - G) / (NIR + G)   (G: 550 nm)
#   - Negative denominators are guarded; results are unitless indices.
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
spectra_path <- opt(opts, "spectra")
outdir <- opt(opts, "outdir", "results/phenomics/02_indices")
idx_arg <- opt(opts, "indices")

if (is.null(spectra_path)) stop("Required: --spectra <table.csv>. See --help.", call. = FALSE)
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills spectral indices ==")

d <- read_csv_safe(spectra_path, "spectral table")
sample_ids <- as.character(d[[1]])
wl_cols <- colnames(d)[-1]

# parse wavelengths from column names
wl <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", wl_cols)))
if (all(is.na(wl))) {
  stop("Could not parse wavelengths from column names. Use names containing the ",
       "wavelength, e.g. R670 / X800 / 730nm.", call. = FALSE)
}
find_band <- function(target, tol = 15) {
  cand <- which(abs(wl - target) <= tol)
  if (!length(cand)) return(NA_integer_)
  cand[which.min(abs(wl[cand] - target))]
}

band_override <- function(arg, target) {
  v <- opt(opts, arg)
  if (!is.null(v)) {
    if (!v %in% wl_cols) stop("--", arg, " '", v, "' not found among spectral columns.", call. = FALSE)
    return(match(v, wl_cols))
  }
  find_band(target)
}
i_nir <- band_override("nir", 800)
i_red <- band_override("red", 670)
i_blue <- band_override("blue", 470)
i_re <- band_override("re", 730)
i_green <- band_override("green", 550)
if (is.na(i_nir) || is.na(i_red)) {
  stop("NIR (800 nm) and Red (670 nm) bands are required (use --nir / --red to set columns).", call. = FALSE)
}
message(sprintf("  bands: NIR=%s Red=%s Blue=%s RE=%s Green=%s",
                wl_cols[i_nir], wl_cols[i_red],
                if (is.na(i_blue)) "-" else wl_cols[i_blue],
                if (is.na(i_re)) "-" else wl_cols[i_re],
                if (is.na(i_green)) "-" else wl_cols[i_green]))

nir <- d[[wl_cols[i_nir]]]; red <- d[[wl_cols[i_red]]]
blue <- if (is.na(i_blue)) NULL else d[[wl_cols[i_blue]]]
re <- if (is.na(i_re)) NULL else d[[wl_cols[i_re]]]
green <- if (is.na(i_green)) NULL else d[[wl_cols[i_green]]]

safe_div <- function(num, den) {
  den[abs(den) < 1e-9] <- NA_real_
  num / den
}

out <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
all_idx <- c("NDVI", "EVI", "SAVI", "RVI", "NDRE", "WDRVI", "GNDVI")
want <- if (is.null(idx_arg)) all_idx else trimws(strsplit(idx_arg, ",")[[1]])
bad <- setdiff(want, all_idx); if (length(bad)) stop("Unknown index(es): ", paste(bad, collapse = ", "), call. = FALSE)

if ("NDVI" %in% want)  out$NDVI  <- safe_div(nir - red, nir + red)
if ("EVI" %in% want) {
  if (is.null(blue)) message("  EVI skipped (blue band not found)")
  else out$EVI <- safe_div(2.5 * (nir - red), nir + 6 * red - 7.5 * blue + 1)
}
if ("SAVI" %in% want)  out$SAVI  <- safe_div(1.5 * (nir - red), nir + red + 0.5)
if ("RVI" %in% want)   out$RVI   <- safe_div(nir, red)
if ("NDRE" %in% want) {
  if (is.null(re)) message("  NDRE skipped (red-edge band not found)")
  else out$NDRE <- safe_div(nir - re, nir + re)
}
if ("WDRVI" %in% want) out$WDRVI <- safe_div(0.1 * nir - red, 0.1 * nir + red)
if ("GNDVI" %in% want) {
  if (is.null(green)) message("  GNDVI skipped (green band not found)")
  else out$GNDVI <- safe_div(nir - green, nir + green)
}

utils::write.csv(out, file.path(outdir, "spectral_indices.csv"), row.names = FALSE)
idx_cols <- setdiff(colnames(out), "sample_id")
message(sprintf("  indices computed: %s", paste(idx_cols, collapse = ", ")))

# figures
pal <- rice_palette()
if (length(idx_cols)) {
  # distribution
  long <- do.call(rbind, lapply(idx_cols, function(cc) {
    data.frame(index = cc, value = out[[cc]])
  }))
  p_dist <- ggplot2::ggplot(long, ggplot2::aes(x = index, y = value)) +
    ggplot2::geom_boxplot(fill = pal$groups5[1], alpha = 0.85, linewidth = 0.3,
                          outlier.size = 0.4, colour = "grey25") +
    ggplot2::facet_wrap(~ index, scales = "free", nrow = 1) +
    ggplot2::labs(x = NULL, y = "Index value", title = "Spectral index distributions") +
    theme_nature() + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  save_nature_figure(p_dist, file.path(outdir, "index_distribution"),
                     width_mm = 160, height_mm = 65)

  if (length(idx_cols) >= 2) {
    cm <- stats::cor(out[, idx_cols, drop = FALSE], use = "pairwise.complete.obs")
    df_c <- data.frame(
      a = rep(rownames(cm), times = ncol(cm)),
      b = rep(colnames(cm), each = nrow(cm)),
      r = as.vector(cm))
    p_cor <- ggplot2::ggplot(df_c, ggplot2::aes(a, b, fill = r)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 2) +
      ggplot2::scale_fill_gradient2(low = diverging_palette()[1], mid = "white",
                                    high = diverging_palette()[7], midpoint = 0,
                                    limits = c(-1, 1), name = "r") +
      ggplot2::labs(x = NULL, y = NULL, title = "Index correlations") +
      theme_nature()
    save_nature_figure(p_cor, file.path(outdir, "index_correlation"),
                       width_mm = 100, height_mm = 90)
  }
}

write_run_info_r(outdir, "02_spectral_indices.R", params = list(
  spectra = spectra_path, indices = paste(idx_cols, collapse = ", "),
  n_samples = length(sample_ids)
))
message("Done. Outputs written to: ", normalizePath(outdir))
