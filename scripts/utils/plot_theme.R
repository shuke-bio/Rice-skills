###############################################################################
# plot_theme.R -- Rice-skills shared Nature-style plotting theme (R backend)
#
# Title:    Rice-skills Nature-style ggplot2 theme and figure helpers
# Description:
#   Provides a low-saturation, colorblind-friendly palette and a ggplot2 theme
#   that follow Nature-series illustration conventions (Arial 8 pt base,
#   bold 9 pt titles, 0.5 pt lines, PDF/TIFF/PNG >=300 dpi vector output).
#   The palette is loaded from config/nature_palette.json so that the R
#   (ggplot2/ComplexHeatmap) and Python (matplotlib) backends stay visually
#   identical. All code is implemented from scratch; no third-party theme
#   source is copied.
#
# Usage:
#   source("scripts/utils/plot_theme.R")   # from the Rice-skills root
#   p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point() +
#        theme_nature()
#   save_nature_figure(p, "results/fig1_volcano", width_mm = 89, height_mm = 70)
#
# Inputs:  config/nature_palette.json (auto-located; override with the
#          RICE_PALETTE environment variable or the palette_file argument)
# Outputs: figures in PDF + TIFF + PNG (300 dpi)
#
# Notes:
#   - Works with ggplot2 >= 3.4 (linewidth aesthetic API).
#   - Arial is used when available; otherwise the device default is used.
###############################################################################

RICE_SEED <- 42L

# ---------------------------------------------------------------------------
# Palette location and loading
# ---------------------------------------------------------------------------

#' Locate the shared palette JSON file.
#'
#' Search order: explicit path -> RICE_PALETTE env var -> config/ subfolder of
#' the current working directory and its parents (up to 4 levels).
#'
#' @param path Optional explicit path to nature_palette.json.
#' @return Absolute path to the palette file.
rice_find_palette <- function(path = NULL) {
  if (!is.null(path) && nzchar(path)) {
    if (!file.exists(path)) {
      stop("Palette file not found: ", path,
           "\nHint: pass a valid path or leave NULL to auto-locate config/nature_palette.json.",
           call. = FALSE)
    }
    return(normalizePath(path, mustWork = TRUE))
  }
  env <- Sys.getenv("RICE_PALETTE", unset = "")
  if (nzchar(env) && file.exists(env)) {
    return(normalizePath(env, mustWork = TRUE))
  }
  d <- getwd()
  for (i in 0:4) {
    candidate <- file.path(d, "config", "nature_palette.json")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate config/nature_palette.json.",
       "\nHint: run scripts from the Rice-skills repository root, or set the",
       "\nRICE_PALETTE environment variable to the absolute JSON path.", call. = FALSE)
}

#' Load the shared palette (cached).
#'
#' @param palette_file Optional explicit path; NULL uses the cached/default file.
#' @return A list with semantic colours, categorical series, diverging and
#'   sequential ramps, typography and figure geometry.
rice_palette <- local({
  cache <- NULL
  function(palette_file = NULL) {
    if (is.null(palette_file) && !is.null(cache)) return(cache)
    pal <- jsonlite::fromJSON(rice_find_palette(palette_file), simplifyVector = TRUE)
    if (is.null(palette_file)) cache <<- pal
    pal
  }
})

# ---------------------------------------------------------------------------
# Font handling
# ---------------------------------------------------------------------------

#' Resolve the font family used for figure text.
#'
#' Arial is used when it is registered as a Type 1 font for the PDF/TIFF/PNG
#' devices; otherwise Helvetica is used. Helvetica is a PDF base-14 font that
#' is always available and metrically compatible with Arial, which guarantees
#' that text is embedded (and remains editable) in vector output on every
#' platform, including Windows.
#'
#' @param pal Palette list from rice_palette().
#' @return Font family string usable in ggplot2 base_family.
.rice_resolve_family <- function(pal) {
  fam <- pal$typography$font_family
  ps <- tryCatch(grDevices::postscriptFonts(), error = function(e) NULL)
  if (!is.null(ps) && fam %in% names(ps)) return(fam)
  "Helvetica"
}

# ---------------------------------------------------------------------------
# ggplot2 theme and scales
# ---------------------------------------------------------------------------

#' Nature-style ggplot2 theme.
#'
#' @param base_size Base font size in points (default from palette: 8).
#' @param ... Additional arguments passed to ggplot2::theme().
#' @return A ggplot2 theme object.
theme_nature <- function(base_size = NULL, ...) {
  pal <- rice_palette()
  bs <- if (is.null(base_size)) pal$typography$base_size_pt else base_size
  fam <- .rice_resolve_family(pal)
  # Note: theme elements are merged with '+' (no need for ggplot2's %+replace%
  # infix, which would require the package to be attached).
  th <- ggplot2::theme_classic(base_size = bs, base_family = fam) +
    ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = pal$typography$axis_line_width_pt,
                                      colour = "grey30"),
    axis.ticks = ggplot2::element_line(linewidth = pal$typography$line_width_pt,
                                       colour = "grey30"),
    axis.ticks.length = ggplot2::unit(1.5, "pt"),
    plot.title = ggplot2::element_text(size = pal$typography$title_size_pt, face = "bold",
                                       hjust = 0, margin = ggplot2::margin(b = 4)),
    plot.subtitle = ggplot2::element_text(size = bs, hjust = 0, colour = "grey40",
                                          margin = ggplot2::margin(b = 6)),
    axis.title = ggplot2::element_text(size = bs, colour = "grey20"),
    axis.text = ggplot2::element_text(size = pal$typography$tick_size_pt, colour = "grey20"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 1.5)),
    axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 1.5)),
    legend.text = ggplot2::element_text(size = pal$typography$legend_size_pt),
    legend.title = ggplot2::element_text(size = pal$typography$legend_size_pt),
    legend.key.size = ggplot2::unit(4, "mm"),
    legend.key.spacing = ggplot2::unit(0.8, "mm"),
    legend.spacing = ggplot2::unit(0.8, "mm"),
    panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = ggplot2::unit(2, "mm"),
    strip.text = ggplot2::element_text(size = bs, face = "bold", colour = "grey20"),
    strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
    plot.margin = ggplot2::margin(3, 4, 3, 1, "mm"),
    ...
  )
  th
}

#' Discrete colour scale using the 12-colour categorical series.
scale_colour_nature <- function(..., aesthetics = "colour") {
  ggplot2::scale_colour_manual(values = rice_palette()$series12, aesthetics = aesthetics, ...)
}

#' Discrete fill scale using the 12-colour categorical series.
scale_fill_nature <- function(..., aesthetics = "fill") {
  ggplot2::scale_fill_manual(values = rice_palette()$series12, aesthetics = aesthetics, ...)
}

#' Discrete colour scale using the 5-colour group series (sample groups).
scale_colour_groups <- function(..., aesthetics = "colour") {
  ggplot2::scale_colour_manual(values = rice_palette()$groups5, aesthetics = aesthetics, ...)
}

#' Discrete fill scale using the 5-colour group series (sample groups).
scale_fill_groups <- function(..., aesthetics = "fill") {
  ggplot2::scale_fill_manual(values = rice_palette()$groups5, aesthetics = aesthetics, ...)
}

#' Diverging colour ramp (blue-white-red, muted) for correlation heatmaps.
diverging_palette <- function() rice_palette()$diverging7

#' Sequential colour ramp (single-hue blue) for intensity heatmaps.
sequential_palette <- function() rice_palette()$sequential6

#' Semantic colours: "up", "down", "not_significant", "highlight", "reference_line".
nature_semantic <- function(name) {
  pal <- rice_palette()
  if (!name %in% names(pal$semantic)) {
    stop("Unknown semantic colour '", name, "'. Available: ",
         paste(names(pal$semantic), collapse = ", "), call. = FALSE)
  }
  unname(pal$semantic[[name]])
}

# ---------------------------------------------------------------------------
# Figure geometry and saving
# ---------------------------------------------------------------------------

#' Convert millimetres to inches.
mm_to_inch <- function(mm) mm / 25.4

#' Save a figure in Nature-compatible formats (PDF + TIFF + PNG 300 dpi).
#'
#' Accepts ggplot2 objects and ComplexHeatmap objects (Heatmap / HeatmapList).
#'
#' @param plot A ggplot object or ComplexHeatmap object.
#' @param file_stem Output path without extension.
#' @param width_mm Figure width in mm (default 89 mm = single column).
#' @param height_mm Figure height in mm (default 89 mm).
#' @param formats Character vector of formats (default from palette).
#' @param ... Extra arguments passed to the graphics device.
#' @return Character vector of saved file paths (invisible).
save_nature_figure <- function(plot, file_stem, width_mm = NULL, height_mm = NULL,
                               formats = NULL, ...) {
  pal <- rice_palette()
  w_mm <- if (is.null(width_mm)) pal$figure$single_column_mm else width_mm
  h_mm <- if (is.null(height_mm)) pal$figure$single_column_mm else height_mm
  fmts <- if (is.null(formats)) pal$figure$formats else formats
  dpi <- pal$figure$dpi_raster
  stem <- normalizePath(file_stem, mustWork = FALSE)
  out_dir <- dirname(stem)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  w_in <- mm_to_inch(w_mm)
  h_in <- mm_to_inch(h_mm)
  is_gg <- inherits(plot, "ggplot")
  is_ch <- any(c("Heatmap", "HeatmapList") %in% class(plot))
  if (!is_gg && !is_ch) {
    stop("save_nature_figure() expects a ggplot object or a ComplexHeatmap object; got: ",
         paste(class(plot), collapse = ", "), call. = FALSE)
  }
  saved <- character(0)
  for (f in fmts) {
    f <- tolower(f)
    out <- paste0(stem, ".", if (f %in% c("tif", "tiff")) "tiff" else f)
    if (is_gg) {
      if (f %in% c("tif", "tiff")) {
        if (requireNamespace("ragg", quietly = TRUE)) {
          ggplot2::ggsave(out, plot, device = ragg::agg_tiff, width = w_in, height = h_in,
                          units = "in", res = dpi, ...)
        } else {
          grDevices::tiff(out, width = w_in, height = h_in, units = "in", res = dpi,
                          compression = "lzw")
          print(plot)
          grDevices::dev.off()
        }
      } else {
        ggplot2::ggsave(out, plot, width = w_in, height = h_in, units = "in", device = f, ...)
      }
    } else {
      if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop("Saving a ComplexHeatmap object requires the ComplexHeatmap package.",
             "\nHint: BiocManager::install('ComplexHeatmap').", call. = FALSE)
      }
      dev <- switch(f,
                    pdf = "pdf", svg = "svg", png = "png",
                    tiff = "tiff", tif = "tiff",
                    stop("Unsupported format for ComplexHeatmap objects: ", f, call. = FALSE))
      devfun <- get(dev, envir = asNamespace("ComplexHeatmap"))
      args <- list(file = out, width = w_in, height = h_in)
      if (f %in% c("tif", "tiff", "png")) args$res <- dpi
      args <- c(args, list(...))
      do.call(devfun, args)
      ComplexHeatmap::draw(plot)
      grDevices::dev.off()
    }
    saved <- c(saved, out)
  }
  invisible(saved)
}

# ---------------------------------------------------------------------------
# Package availability check (shared by all R scripts)
# ---------------------------------------------------------------------------

#' Check that required packages are installed; stop with install hints.
#'
#' @param pkgs Character vector of package names.
#' @return TRUE (invisible) when all packages are available.
need_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    bioc <- c("DESeq2", "edgeR", "clusterProfiler", "org.Osativa.eg.db", "ReactomePA",
              "fgsea", "WGCNA", "MOFA2", "mixOmics", "ComplexHeatmap", "MetaboAnalystR",
              "imputeLCMD")
    hint <- if (any(missing %in% bioc)) {
      paste0("\nHint: install Bioconductor packages with:\n  Rscript -e 'BiocManager::install(c(\"",
             paste(missing, collapse = "\", \""), "\"))'")
    } else {
      paste0("\nHint: install CRAN packages with:\n  install.packages(c(\"",
             paste(missing, collapse = "\", \""), "\"))")
    }
    stop("Missing required R package(s): ", paste(missing, collapse = ", "), hint,
         "\nSee docs/installation.md for the full dependency list.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Reproducibility and I/O helpers shared by all R scripts
# ---------------------------------------------------------------------------

#' Set the analysis seed (default 42) and record it for run_info.txt.
rice_set_seed <- function(seed = RICE_SEED) {
  assign("RICE_SEED", seed, envir = .GlobalEnv)
  set.seed(seed)
  invisible(seed)
}

#' Friendly CSV reader with clear error messages.
#'
#' @param path File path.
#' @param what Human-readable description used in the error message.
#' @return A data.frame.
read_csv_safe <- function(path, what = "input file") {
  if (!file.exists(path)) {
    stop(what, " not found: ", path,
         "\nHint: check the path and run the previous pipeline step first.", call. = FALSE)
  }
  df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) {
                   stop("Failed to parse ", what, ": ", path, "\n", conditionMessage(e),
                        "\nHint: expected a comma-separated file with a header row.",
                        call. = FALSE)
                 })
  if (ncol(df) < 2) {
    stop(what, " must contain at least an ID column and one sample column: ", path, call. = FALSE)
  }
  df
}

#' Check that required columns are present; stop with a friendly message.
check_columns <- function(df, required, df_name = "table") {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing required column(s): ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Write run_info.txt documenting parameters and package versions.
#'
#' @param outdir Output directory.
#' @param script_name Script file name (for provenance).
#' @param params Named list of parameters.
#' @param packages Character vector of key package names to version-check.
#' @param extra Extra note lines.
#' @return Path to run_info.txt (invisible).
write_run_info_r <- function(outdir, script_name, params = list(),
                             packages = character(0), extra = NULL) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  seed <- get0("RICE_SEED", envir = .GlobalEnv, ifnotfound = RICE_SEED)
  lines <- c(
    sprintf("script: %s", script_name),
    sprintf("run_time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    sprintf("R_version: %s", R.version.string),
    sprintf("seed: %s", seed),
    "parameters:"
  )
  if (length(params) > 0) {
    lines <- c(lines, vapply(names(params), function(n) {
      sprintf("  %s = %s", n, paste(unlist(params[[n]]), collapse = ","))
    }, character(1)))
  }
  if (length(packages) > 0) {
    ip <- utils::installed.packages()
    lines <- c(lines, "packages:", vapply(packages, function(p) {
      ver <- if (p %in% rownames(ip)) ip[p, "Version"] else "NOT INSTALLED"
      sprintf("  %s = %s", p, ver)
    }, character(1)))
  }
  if (!is.null(extra)) lines <- c(lines, "notes:", paste0("  ", extra))
  out <- file.path(outdir, "run_info.txt")
  writeLines(lines, out)
  invisible(out)
}
