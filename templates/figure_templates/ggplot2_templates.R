###############################################################################
# ggplot2_templates.R -- Rice-skills publication figure templates (R backend)
#
# Title:    Rice-skills ggplot2 figure template library
# Description:
#   Copy-and-adapt starting points for every publication figure of the
#   Rice-skills pipeline. Each template_*() function returns a ready-to-use
#   ggplot2 object (or a pheatmap object for template_heatmap); the caller is
#   responsible for saving it with save_nature_figure() / save_pheatmap_nature()
#   from the shared theme (PDF + SVG + TIFF at 300 dpi).
#
#   All colours are taken from config/nature_palette.json through the shared
#   theme API (rice_palette(), nature_semantic(), diverging_palette(),
#   sequential_palette(), scale_*_nature(), scale_*_groups()); no colour is
#   hard-coded here, so the R and Python (matplotlib) backends stay visually
#   identical.
#
#   The finished scripts under scripts/ are generated from these templates:
#   copy the template function into the analysis script, wire it to the
#   pipeline CSV outputs, and keep the theme calls untouched.
#
# IMPORTANT: source the shared theme BEFORE this file. From the repository
# root:
#   source("scripts/utils/plot_theme.R")
#   source("templates/figure_templates/ggplot2_templates.R")
# or, when working from this folder:
#   source("../../scripts/utils/plot_theme.R")
#   source("ggplot2_templates.R")
#
# Dependencies: ggplot2 (required); ggrepel (optional, falls back to
#   geom_text); pheatmap (for template_heatmap). Everything else is base R.
###############################################################################

# ---------------------------------------------------------------------------
# Internal helpers (not part of the public template API)
# ---------------------------------------------------------------------------

#' Normalise a `significant` column to the factor levels Up / DOWN / NS.
#'
#' Accepts the common spellings produced by DESeq2 wrappers and by the
#' pipeline scripts: "Up", "UP", "up", "DOWN", "Down", "down", "NS", "ns",
#' "NOT_SIGNIFICANT", ... Anything that is not up/down is treated as NS.
#'
#' @param sig Character vector of significance labels.
#' @return A factor with levels c("Up", "DOWN", "NS").
.rice_sig_factor <- function(sig) {
  s <- toupper(trimws(as.character(sig)))
  out <- rep("NS", length(s))
  out[grepl("^UP", s)] <- "Up"
  out[grepl("^DOWN", s)] <- "DOWN"
  factor(out, levels = c("Up", "DOWN", "NS"))
}

#' Discrete colour scale for the semantic significance colours.
#' @param ... Passed to ggplot2::scale_colour_manual().
#' @return A ggplot2 discrete colour scale.
scale_colour_significant <- function(...) {
  ggplot2::scale_colour_manual(
    values = c(Up = nature_semantic("up"),
               DOWN = nature_semantic("down"),
               NS = nature_semantic("not_significant")),
    ...)
}

#' Discrete fill scale for the semantic significance colours.
#' @param ... Passed to ggplot2::scale_fill_manual().
#' @return A ggplot2 discrete fill scale.
scale_fill_significant <- function(...) {
  ggplot2::scale_fill_manual(
    values = c(Up = nature_semantic("up"),
               DOWN = nature_semantic("down"),
               NS = nature_semantic("not_significant")),
    ...)
}

#' Derive the `significant` labels from the significance thresholds.
#'
#' Used when the input table carries no pre-computed `significant` column.
#'
#' @param df Data frame with `log2FoldChange` and `padj` columns.
#' @param log2fc_cutoff |log2FC| threshold.
#' @param padj_cutoff Adjusted p-value threshold.
#' @return Character vector of "Up" / "DOWN" / "NS".
.rice_derive_sig <- function(df, log2fc_cutoff, padj_cutoff) {
  if (!"padj" %in% names(df)) {
    stop("Cannot derive 'significant': the input has neither a 'significant' ",
         "column nor a 'padj' column.", call. = FALSE)
  }
  ifelse(df$padj < padj_cutoff,
         ifelse(df$log2FoldChange >= log2fc_cutoff, "Up",
                ifelse(df$log2FoldChange <= -log2fc_cutoff, "DOWN", "NS")),
         "NS")
}

#' Prepare a differential-expression table for volcano/MA plots.
#'
#' Guarantees the presence of a normalised `sig` factor column.
#'
#' @param df Data frame with at least gene_id and log2FoldChange.
#' @param log2fc_cutoff |log2FC| threshold (only used when deriving labels).
#' @param padj_cutoff Adjusted p-value threshold (only used when deriving).
#' @return The input data frame with an added `sig` factor column.
.rice_prepare_de <- function(df, log2fc_cutoff = 1, padj_cutoff = 0.05) {
  check_columns(df, c("gene_id", "log2FoldChange"), "DE table")
  d <- df
  if (!"significant" %in% names(d)) {
    d$significant <- .rice_derive_sig(d, log2fc_cutoff, padj_cutoff)
  }
  d$sig <- .rice_sig_factor(d$significant)
  d
}

#' Text layer for gene labels: ggrepel when available, geom_text otherwise.
#' @param labelled_df Data frame restricted to the genes to label.
#' @param size Text size in points.
#' @return A ggplot2 layer.
.rice_label_layer <- function(labelled_df, size = 2.2) {
  ink <- "grey20"  # theme ink colour (matches theme_nature() axis text)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = labelled_df, ggplot2::aes(label = gene_id),
      size = size, colour = ink, min.segment.length = 0, seed = RICE_SEED,
      max.overlaps = 20, box.padding = 0.3, point.padding = 0.2,
      show.legend = FALSE)
  } else {
    ggplot2::geom_text(
      data = labelled_df, ggplot2::aes(label = gene_id),
      size = size, colour = ink, hjust = 0, vjust = -0.5, show.legend = FALSE)
  }
}

#' Keep the top `top_n` rows per group (or overall) by a score column.
#'
#' Base-R implementation (no dplyr dependency) used by the enrichment
#' templates to keep multi-category plots readable.
#'
#' @param d Data frame.
#' @param score_col Numeric column used for ranking (higher = kept).
#' @param top_n Number of rows to keep per group (NULL or <= 0 keeps all).
#' @param group_col Optional grouping column (e.g. "Category").
#' @return The filtered data frame.
.rice_top_n <- function(d, score_col, top_n, group_col = "Category") {
  if (is.null(top_n) || top_n <= 0 || nrow(d) <= top_n) return(d)
  ord <- order(-d[[score_col]])
  if (!is.null(group_col) && group_col %in% names(d)) {
    keep <- unlist(lapply(split(ord, d[[group_col]][ord]),
                          utils::head, n = top_n), use.names = FALSE)
  } else {
    keep <- utils::head(ord, top_n)
  }
  d[sort(keep), , drop = FALSE]
}

#' Return the first of `candidates` present in `df` (NULL when none).
#' @param df Data frame.
#' @param candidates Character vector of candidate column names.
#' @return The matched column name or NULL.
.rice_pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

#' Annotation colour list for pheatmap from the shared palette.
#'
#' Uses the 5 group colours first (groups5) and extends with series12 when a
#' column has more than five levels.
#'
#' @param annotation_df Annotation data frame (or NULL).
#' @return Named list of named colour vectors (level -> colour), or NULL.
.rice_annotation_colors <- function(annotation_df) {
  if (is.null(annotation_df) || ncol(annotation_df) == 0) return(NULL)
  pal <- rice_palette()
  cols <- pal$series12  # the first five entries are exactly groups5
  lapply(annotation_df, function(x) {
    levs <- sort(unique(as.character(x)))
    if (length(levs) > length(cols)) {
      warning("More annotation levels (", length(levs), ") than palette ",
              "colours (", length(cols), "); colours are recycled.", call. = FALSE)
    }
    setNames(rep(cols, length.out = length(levs)), levs)
  })
}

#' Save a pheatmap object in Nature-compatible formats (PDF + SVG + TIFF 300 dpi).
#'
#' pheatmap objects cannot be passed to save_nature_figure() (which accepts
#' ggplot2 and ComplexHeatmap objects only), so this helper opens the devices
#' directly and draws the pheatmap gtable.
#'
#' @param ph A pheatmap object returned by template_heatmap().
#' @param file_stem Output path without extension.
#' @param width_mm Figure width in mm (default: single column, 89 mm).
#' @param height_mm Figure height in mm (default: single column, 89 mm).
#' @param ... Extra arguments passed to the graphics device.
#' @return Character vector of saved file paths (invisible).
save_pheatmap_nature <- function(ph, file_stem, width_mm = NULL, height_mm = NULL,
                                 ...) {
  pal <- rice_palette()
  w_mm <- if (is.null(width_mm)) pal$figure$single_column_mm else width_mm
  h_mm <- if (is.null(height_mm)) pal$figure$single_column_mm else height_mm
  w_in <- mm_to_inch(w_mm)
  h_in <- mm_to_inch(h_mm)
  stem <- normalizePath(file_stem, mustWork = FALSE)
  out_dir <- dirname(stem)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saved <- character(0)
  for (f in pal$figure$formats) {
    out <- paste0(stem, ".", if (f %in% c("tif", "tiff")) "tiff" else tolower(f))
    if (f %in% c("tif", "tiff")) {
      grDevices::tiff(out, width = w_in, height = h_in, units = "in",
                      res = pal$figure$dpi_raster, compression = "lzw")
    } else {
      switch(tolower(f),
             pdf = grDevices::pdf(out, width = w_in, height = h_in),
             svg = grDevices::svg(out, width = w_in, height = h_in),
             stop("Unsupported format for pheatmap objects: ", f, call. = FALSE))
    }
    grid::grid.draw(ph$gtable)
    grDevices::dev.off()
    saved <- c(saved, out)
  }
  invisible(saved)
}

# ---------------------------------------------------------------------------
# Public templates -- each returns a ggplot2 object (save with
# save_nature_figure()) except template_heatmap() (save with
# save_pheatmap_nature()).
# ---------------------------------------------------------------------------

#' Volcano plot of differential expression.
#'
#' Points are coloured by the semantic significance colours (up = red,
#' down = blue, not significant = grey). Dashed grey reference lines mark the
#' log2FC and adjusted-p thresholds. The top `label_top` significant genes are
#' ranked by |log2FC| (ties broken by padj) and labelled with ggrepel when it
#' is installed (geom_text fallback otherwise).
#'
#' @param df Data frame with columns gene_id, log2FoldChange, padj and
#'   (optionally) significant ("Up"/"DOWN"/"NS"); when `significant` is absent
#'   it is derived from the cutoffs.
#' @param log2fc_cutoff |log2FC| threshold for the vertical lines (default 1).
#' @param padj_cutoff Adjusted p-value threshold for the horizontal line
#'   (default 0.05).
#' @param label_top Number of genes to label (0 disables labelling).
#' @param title Plot title.
#' @return A ggplot2 object.
#'
#' Usage:
#'   de <- read.csv("results/transcriptome/02_de/de_all_salt_vs_control.csv")
#'   p <- template_volcano(de)
#'   save_nature_figure(p, "results/transcriptome/07_figures/fig_volcano",
#'                      width_mm = 89, height_mm = 80)
template_volcano <- function(df, log2fc_cutoff = 1, padj_cutoff = 0.05,
                             label_top = 5, title = "Volcano plot") {
  d <- .rice_prepare_de(df, log2fc_cutoff, padj_cutoff)
  d$y <- -log10(pmax(d$padj, 1e-300))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = log2FoldChange, y = y, colour = sig)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
                        linetype = 2, linewidth = 0.4,
                        colour = nature_semantic("reference_line")) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                        linetype = 2, linewidth = 0.4,
                        colour = nature_semantic("reference_line")) +
    scale_colour_significant(name = NULL) +
    ggplot2::labs(title = title,
                  x = expression(log[2]~fold~change),
                  y = expression(-log[10]~adjusted~italic(P))) +
    theme_nature()
  if (label_top > 0) {
    sig_rows <- d[d$sig != "NS", , drop = FALSE]
    if (nrow(sig_rows) > 0) {
      ord <- order(-abs(sig_rows$log2FoldChange), sig_rows$padj)
      top <- sig_rows[utils::head(ord, label_top), , drop = FALSE]
      p <- p + .rice_label_layer(top)
    }
  }
  p
}

#' MA plot (log2 fold change vs mean expression).
#'
#' x-axis is log10(baseMean); genes with baseMean == 0 (removed by the QC
#' filter in practice) are not plotted. Colours follow the semantic palette;
#' a dashed grey line marks y = 0 and the |log2FC| cutoffs.
#'
#' @param df Data frame with columns gene_id, baseMean, log2FoldChange and
#'   (optionally) significant or padj.
#' @param log2fc_cutoff |log2FC| threshold for the vertical lines (default 1).
#' @param padj_cutoff Adjusted p-value threshold (only used when deriving
#'   labels; default 0.05).
#' @param title Plot title.
#' @return A ggplot2 object.
#'
#' Usage:
#'   de <- read.csv("results/transcriptome/02_de/de_all_salt_vs_control.csv")
#'   p <- template_ma(de)
#'   save_nature_figure(p, "results/transcriptome/07_figures/fig_ma",
#'                      width_mm = 89, height_mm = 80)
template_ma <- function(df, log2fc_cutoff = 1, padj_cutoff = 0.05,
                        title = "MA plot") {
  check_columns(df, c("gene_id", "baseMean", "log2FoldChange"), "DE table (MA)")
  d <- .rice_prepare_de(df, log2fc_cutoff, padj_cutoff)
  d$x <- ifelse(d$baseMean > 0, log10(d$baseMean), NA_real_)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = log2FoldChange, colour = sig)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.85, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4,
                        colour = nature_semantic("reference_line")) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
                        linetype = 2, linewidth = 0.4,
                        colour = nature_semantic("reference_line")) +
    scale_colour_significant(name = NULL) +
    ggplot2::labs(title = title,
                  x = expression(log[10]~base~mean),
                  y = expression(log[2]~fold~change)) +
    theme_nature()
}

#' Normal-theory confidence ellipse coordinates per group.
#'
#' ggplot2's stat_ellipse() (>= 4.0) refuses groups with three or fewer points
#' ("Too few points to calculate an ellipse"), which is the common case for
#' three-replicate designs; this helper computes the covariance ellipse
#' directly and works for every group with at least three samples. Groups with
#' fewer points (or a singular covariance) are skipped with a message.
#'
#' @param d Data frame with the x/y columns and a grouping column.
#' @param xvar Name of the x column.
#' @param yvar Name of the y column.
#' @param group_col Name of the grouping column (a factor).
#' @param level Confidence level (default 0.95).
#' @return Data frame with columns grp, x, y, colour; NULL when no group has
#'   enough points.
.rice_ellipse_df <- function(d, xvar, yvar, group_col, level = 0.95) {
  chi2 <- stats::qchisq(level, df = 2)
  cols <- rice_palette()$series12  # same order as scale_colour_nature()
  levs <- levels(d[[group_col]])
  pieces <- list()
  skipped <- 0L
  for (g in levs) {
    sub <- d[!is.na(d[[group_col]]) & d[[group_col]] == g, , drop = FALSE]
    ok <- nrow(sub) >= 3 &&
      all(stats::complete.cases(sub[c(xvar, yvar)]))
    if (ok) {
      xy <- cbind(as.numeric(sub[[xvar]]), as.numeric(sub[[yvar]]))
      covm <- stats::cov(xy)
      ok <- all(is.finite(covm)) &&
        all(eigen(covm, symmetric = TRUE, only.values = TRUE)$values > 0)
    }
    if (!ok) {
      skipped <- skipped + 1L
      next
    }
    chol <- chol(covm)
    ang <- seq(0, 2 * pi, length.out = 61)
    circle <- cbind(cos(ang), sin(ang)) * sqrt(chi2)
    pts <- sweep(circle %*% chol, 2, colMeans(xy), "+")
    idx <- match(g, levs)
    pieces[[length(pieces) + 1L]] <- data.frame(
      grp = g, x = pts[, 1], y = pts[, 2],
      colour = cols[(idx - 1L) %% length(cols) + 1L],
      stringsAsFactors = FALSE)
  }
  if (length(pieces) == 0) return(NULL)
  if (skipped > 0L) {
    message("Skipped confidence ellipse(s) for ", skipped,
            " group(s) with fewer than 3 samples.")
  }
  do.call(rbind, pieces)
}

#' PCA score plot with per-group 95% confidence ellipses.
#'
#' `meta_df` is indexed by sample_id (rownames) or carries a sample_id column.
#' Colours use the categorical series (whose first five entries are the group
#' colours); `shape_col` optionally maps a second metadata column to shapes.
#' Numeric metadata columns (e.g. batch = 1/2) are converted to factors
#' automatically. The dashed ellipse is a normal-theory 95% confidence ellipse
#' computed from the sample covariance (see .rice_ellipse_df()); set
#' `ellipse = FALSE` to draw points only.
#'
#' @param scores_df Data frame with columns sample_id, PC1, PC2.
#' @param meta_df Metadata data frame indexed by sample_id.
#' @param color_col Name of the metadata column used for point colours.
#' @param shape_col Optional name of the metadata column used for shapes.
#' @param ellipse Logical; draw 95% confidence ellipses (default TRUE).
#' @param title Plot title.
#' @return A ggplot2 object.
#'
#' Usage:
#'   scores <- read.csv("results/transcriptome/01_qc/pca_scores.csv")
#'   meta <- read.csv("metadata.csv", row.names = 1)
#'   p <- template_pca(scores, meta, color_col = "group", shape_col = "batch")
#'   save_nature_figure(p, "results/transcriptome/07_figures/fig_pca_qc",
#'                      width_mm = 89, height_mm = 80)
template_pca <- function(scores_df, meta_df, color_col, shape_col = NULL,
                         ellipse = TRUE, title = "PCA score plot") {
  check_columns(scores_df, c("sample_id", "PC1", "PC2"), "scores_df")
  if (!color_col %in% names(meta_df)) {
    stop("color_col '", color_col, "' not found in meta_df. Available columns: ",
         paste(names(meta_df), collapse = ", "), call. = FALSE)
  }
  if (!is.null(shape_col) && !shape_col %in% names(meta_df)) {
    stop("shape_col '", shape_col, "' not found in meta_df. Available columns: ",
         paste(names(meta_df), collapse = ", "), call. = FALSE)
  }
  meta2 <- meta_df
  if (!"sample_id" %in% names(meta2)) meta2$sample_id <- rownames(meta_df)
  d <- merge(scores_df, meta2, by = "sample_id", all.x = TRUE, sort = FALSE)
  # Discrete scales require discrete variables: read.csv() turns e.g. batch
  # into an integer column, which ggplot2 (>= 4.0) refuses to map to shape.
  d[[color_col]] <- factor(d[[color_col]])
  if (!is.null(shape_col)) d[[shape_col]] <- factor(d[[shape_col]])
  if (anyNA(d[[color_col]])) {
    warning("Some samples in scores_df have no metadata row; their colour ",
            "value is missing.", call. = FALSE)
  }
  aes_point <- ggplot2::aes(x = PC1, y = PC2, colour = .data[[color_col]])
  if (!is.null(shape_col)) {
    aes_point <- ggplot2::aes(x = PC1, y = PC2, colour = .data[[color_col]],
                              shape = .data[[shape_col]])
  }
  p <- ggplot2::ggplot(d, aes_point) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(title = title, x = "PC1", y = "PC2",
                  colour = color_col, shape = shape_col) +
    scale_colour_nature(name = color_col) +
    theme_nature()
  if (!is.null(shape_col)) {
    p <- p + ggplot2::scale_shape_discrete(name = shape_col)
  }
  if (isTRUE(ellipse)) {
    ell <- .rice_ellipse_df(d, xvar = "PC1", yvar = "PC2", group_col = color_col)
    if (!is.null(ell)) {
      p <- p + ggplot2::geom_path(
        data = ell,
        ggplot2::aes(x = x, y = y, group = grp, colour = I(colour)),
        linetype = 2, linewidth = 0.4, inherit.aes = FALSE)
    }
  }
  p
}

#' Enrichment dot plot (GO / KEGG / Reactome style).
#'
#' x = GeneRatio (or Count when GeneRatio is absent), y = reorder(Description),
#' size = Count, colour = -log10(p.adjust) on the sequential ramp. When the
#' table has a `Category` column the plot is faceted per category and the
#' `top_n` most significant terms of each category are shown.
#'
#' @param enrich_df Enrichment result table (clusterProfiler-style columns:
#'   Description, p.adjust, plus GeneRatio/Count and optional Category).
#' @param top_n Maximum number of terms per category (NULL keeps all).
#' @param title Plot title.
#' @return A ggplot2 object.
#'
#' Usage:
#'   ego <- read.csv("results/transcriptome/03_enrichment/go_bp.csv")
#'   p <- template_dotplot(ego, top_n = 15)
#'   save_nature_figure(p, "results/transcriptome/07_figures/fig_go_bp",
#'                      width_mm = 89, height_mm = 120)
template_dotplot <- function(enrich_df, top_n = 20, title = "Enrichment dot plot") {
  check_columns(enrich_df, c("Description", "p.adjust"), "enrich_df")
  d <- enrich_df
  xvar <- .rice_pick_col(d, c("GeneRatio", "Count"))
  if (is.null(xvar)) {
    stop("enrich_df needs a 'GeneRatio' or 'Count' column for the x-axis.",
         call. = FALSE)
  }
  svar <- .rice_pick_col(d, c("Count", "GeneRatio"))
  d$neglog10p <- -log10(pmax(as.numeric(d$p.adjust), 1e-300))
  d <- .rice_top_n(d, "neglog10p", top_n, group_col = "Category")
  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data[[xvar]], y = reorder(Description, .data[[xvar]]),
                    size = .data[[svar]], colour = neglog10p)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradientn(colours = sequential_palette(),
                                    name = expression(-log[10]~adjusted~italic(P))) +
    ggplot2::scale_size_continuous(name = svar) +
    ggplot2::labs(title = title, x = xvar, y = NULL) +
    theme_nature()
  if ("Category" %in% names(d) && length(unique(d$Category)) > 1) {
    p <- p + ggplot2::facet_wrap(~ Category, scales = "free_y", ncol = 1)
  }
  p
}

#' Enrichment bubble plot (GSEA style).
#'
#' x = NES (or a fold-enrichment / ratio column), y = Description,
#' size = gene-set size (Count / set_size when available), colour =
#' -log10(p.adjust) on the sequential ramp.
#'
#' @param enrich_df Enrichment result table; the x column is auto-detected
#'   from c("NES", "fold_enrichment", "FoldEnrichment", "foldEnrichment",
#'   "GeneRatio", "Count"), the size column from c("set_size", "setSize",
#'   "set size", "Count", "GeneRatio").
#' @param top_n Maximum number of terms per category (NULL keeps all).
#' @param title Plot title.
#' @return A ggplot2 object.
#'
#' Usage:
#'   gsea <- read.csv("results/transcriptome/04_gsea/gsea_results.csv")
#'   p <- template_bubble(gsea, top_n = 15)
#'   save_nature_figure(p, "results/transcriptome/07_figures/fig_gsea",
#'                      width_mm = 89, height_mm = 120)
template_bubble <- function(enrich_df, top_n = 20, title = "Enrichment bubble plot") {
  check_columns(enrich_df, c("Description"), "enrich_df")
  pcol <- .rice_pick_col(enrich_df, c("p.adjust", "padj", "pvalue", "p_value",
                                      "Pvalue", "P.adjust"))
  if (is.null(pcol)) {
    stop("enrich_df needs a p-value column (e.g. 'p.adjust').", call. = FALSE)
  }
  xvar <- .rice_pick_col(enrich_df, c("NES", "fold_enrichment", "FoldEnrichment",
                                      "foldEnrichment", "GeneRatio", "Count"))
  if (is.null(xvar)) {
    stop("enrich_df needs an x-axis column (e.g. 'NES').", call. = FALSE)
  }
  svar <- .rice_pick_col(enrich_df, c("set_size", "setSize", "set size",
                                      "Count", "GeneRatio"))
  xlab <- switch(xvar,
                 NES = "Normalised enrichment score",
                 fold_enrichment = , FoldEnrichment = , foldEnrichment = "Fold enrichment",
                 xvar)
  d <- enrich_df
  d$neglog10p <- -log10(pmax(as.numeric(d[[pcol]]), 1e-300))
  d <- .rice_top_n(d, "neglog10p", top_n, group_col = "Category")
  if (is.null(svar)) {
    aes <- ggplot2::aes(x = .data[[xvar]],
                        y = reorder(Description, .data[[xvar]]),
                        colour = neglog10p)
  } else {
    aes <- ggplot2::aes(x = .data[[xvar]],
                        y = reorder(Description, .data[[xvar]]),
                        colour = neglog10p, size = .data[[svar]])
  }
  p <- ggplot2::ggplot(d, aes) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradientn(colours = sequential_palette(),
                                    name = expression(-log[10]~adjusted~italic(P))) +
    ggplot2::labs(title = title, x = xlab, y = NULL, size = svar) +
    theme_nature()
  if ("Category" %in% names(d) && length(unique(d$Category)) > 1) {
    p <- p + ggplot2::facet_wrap(~ Category, scales = "free_y", ncol = 1)
  }
  p
}

#' Heatmap template (pheatmap wrapper).
#'
#' `type = "diverging"` uses the blue-white-red ramp (correlation matrices);
#' `type = "sequential"` uses the single-hue blue ramp (expression / abundance
#' matrices). Annotation colours are taken from the group palette (groups5,
#' extended by series12). If `mat` is a symmetric correlation matrix and
#' `annotation_df` is given, the annotations are drawn on rows and columns.
#'
#' NOTE: pheatmap objects must be saved with save_pheatmap_nature() because
#' save_nature_figure() only accepts ggplot2 / ComplexHeatmap objects. When
#' ComplexHeatmap is installed you can switch this template to
#' ComplexHeatmap::Heatmap() and then save with save_nature_figure() directly.
#'
#' @param mat Numeric matrix with row/column names (genes x samples or a
#'   symmetric correlation matrix).
#' @param annotation_df Optional annotation data frame whose rownames match
#'   colnames(mat) (and rownames(mat) for symmetric matrices).
#' @param type "diverging" (default) or "sequential".
#' @param fontsize Base font size for pheatmap labels (default 8).
#' @param ... Extra arguments passed to pheatmap::pheatmap().
#' @return A pheatmap object.
#'
#' Usage:
#'   mat <- as.matrix(read.csv("results/transcriptome/02_de/vst_matrix.csv",
#'                             row.names = 1))
#'   meta <- read.csv("metadata.csv", row.names = 1)
#'   ph <- template_heatmap(mat, annotation_df = meta["group"])
#'   save_pheatmap_nature(ph, "results/transcriptome/07_figures/fig_heatmap_de",
#'                        width_mm = 120, height_mm = 180)
template_heatmap <- function(mat, annotation_df = NULL,
                             type = c("diverging", "sequential"),
                             fontsize = 8, ...) {
  type <- match.arg(type)
  need_pkgs("pheatmap")
  ramp <- grDevices::colorRampPalette(
    if (identical(type, "diverging")) diverging_palette() else sequential_palette()
  )(99)
  args <- list(mat = mat, color = ramp, fontsize = fontsize,
               annotation_col = annotation_df,
               annotation_colors = .rice_annotation_colors(annotation_df))
  # Symmetric (correlation) matrices: annotate rows and columns alike.
  if (isSymmetric(mat) && !is.null(annotation_df)) {
    args$annotation_row <- annotation_df
  }
  args <- c(args, list(...))
  do.call(pheatmap::pheatmap, args)
}

#' Generic horizontal bar chart.
#'
#' y is reordered by x so the largest bar ends up on top; `fill` optionally
#' maps a grouping column to the categorical palette.
#'
#' @param df Data frame with the x and y columns.
#' @param x Name of the numeric value column.
#' @param y Name of the category column.
#' @param fill Optional name of the grouping column used for bar colours.
#' @param title Plot title (NULL for none).
#' @param xlab x-axis label (defaults to the column name).
#' @return A ggplot2 object.
#'
#' Usage:
#'   cls <- read.csv("results/metabolome/04_classification/class_counts.csv")
#'   p <- template_bar(cls, x = "n_metabolites", y = "class")
#'   save_nature_figure(p, "results/metabolome/07_figures/fig_class_bar",
#'                      width_mm = 89, height_mm = 80)
template_bar <- function(df, x, y, fill = NULL, title = NULL, xlab = NULL) {
  check_columns(df, c(x, y), "df (bar input)")
  p <- ggplot2::ggplot(
    df, ggplot2::aes(x = .data[[x]], y = reorder(.data[[y]], .data[[x]])))
  if (is.null(fill)) {
    p <- p + ggplot2::geom_col(fill = rice_palette()$series12[[1]], width = 0.7)
  } else {
    check_columns(df, fill, "df (bar input)")
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data[[fill]]), width = 0.7) +
      scale_fill_nature(name = fill)
  }
  p + ggplot2::labs(title = title, x = if (is.null(xlab)) x else xlab, y = NULL) +
    theme_nature()
}
