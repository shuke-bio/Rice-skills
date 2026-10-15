###############################################################################
# 07_plot_metabolome.R -- publication figure assembly for the metabolome module
#
# Title:    Metabolome overview figures rebuilt from steps 01-06 (step 07)
# Description:
#   Re-assembles the metabolome results into four publication-ready overview
#   figures, re-drawing everything from the saved CSV tables so that the
#   figures stay reproducible even when the individual steps are re-run with
#   different parameters:
#     fig1_preprocess_overview  PCA (normalised) + per-sample boxplot +
#                               missing-value overview
#     fig2_oplsda_overview      OPLS-DA scores + S-plot + top-30 VIP
#     fig3_biomarker_overview   ROC curves (top 5) + random-forest importance
#     fig4_network_overview     metabolite correlation network
#   Missing upstream results are skipped with a message; when patchwork is not
#   installed, the individual panels are saved separately instead of combined.
#
# Usage:
#   Rscript scripts/metabolome/07_plot_metabolome.R \
#       --results-dir results/metabolome --outdir results/metabolome/07_figures
#
# Inputs:
#   --results-dir  Directory containing the per-step output folders
#                  (01_preprocess, 02_oplsda, 05_biomarker, 06_network).
#
# Outputs (written to --outdir):
#   fig1_preprocess_overview.*   (PDF/SVG/TIFF)
#   fig2_oplsda_overview.*
#   fig3_biomarker_overview.*
#   fig4_network_overview.*
#   run_info.txt
#
# Notes:
#   - Every panel carries a subtitle stating its scientific conclusion.
#   - The network figure is rebuilt from metabolite_edges.csv /
#     metabolite_nodes.csv; ggraph is used when available, otherwise the igraph
#     base device (PDF/SVG/TIFF via grDevices).
#   - patchwork combines the multi-panel figures; without it, panels are saved
#     as figN_..._panelA/B/C.*.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
rice_set_seed(42)

# ---------------------------------------------------------------------------
# Command-line parsing (standard Rice-skills implementation)
# ---------------------------------------------------------------------------

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
  if (isTRUE(v)) return(if (is.logical(default)) default else TRUE); v
}

need_pkgs(c("ggplot2", "ggrepel"))

# ---------------------------------------------------------------------------
# Result discovery helpers
# ---------------------------------------------------------------------------

ANNOTATION_COLS <- c("metabolite_id", "metabolite_name", "class", "kegg_id")

#' Read a CSV if it exists; return NULL (with a message) otherwise.
read_if_exists <- function(path, what) {
  if (!file.exists(path)) {
    message("Skipping '", what, "': file not found (", path, ").")
    return(NULL)
  }
  read_csv_safe(path, what)
}

#' Identify the numeric sample columns of a metabolite matrix.
sample_columns <- function(df) {
  names(df)[vapply(df, is.numeric, logical(1))]
}

#' Add 95% confidence ellipses when every group has >= 4 points.
#'
#' ggplot2::stat_ellipse needs more than 3 points per group; with the typical
#' 3-replicate metabolome design the ellipses are omitted with a message.
maybe_add_ellipses <- function(p) {
  grp <- p$data$group
  if (is.null(grp)) return(p)
  sizes <- table(grp)
  if (min(sizes) >= 4L) {
    p + ggplot2::stat_ellipse(ggplot2::aes(group = group), level = 0.95,
                              type = "norm", linewidth = 0.5)
  } else {
    message("Confidence ellipses omitted (stat_ellipse needs >= 4 samples per ",
            "group; smallest group: ", min(sizes), ").")
    p
  }
}

#' Coerce a 'significant' column to logical, tolerating "TRUE"/"True"/"T"/1.
as_significant <- function(x) {
  if (is.logical(x)) return(x)
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(FALSE, length(x_chr))
  out[x_chr %in% c("true", "t", "yes", "1")] <- TRUE
  out[is.na(x_chr) | !nzchar(x_chr)] <- FALSE
  out
}

#' Save a list of panels: combined via patchwork when available, else separately.
save_panels <- function(panels, stem, width_mm, height_mm, ncol, widths) {
  keep <- !vapply(panels, is.null, logical(1))
  panels <- panels[keep]
  if (!length(panels)) {
    message("No panels available for ", stem, "; nothing saved.")
    return(invisible(NULL))
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- Reduce(`+`, panels) +
      patchwork::plot_layout(ncol = ncol, widths = widths)
    save_nature_figure(combined, stem, width_mm = width_mm, height_mm = height_mm)
  } else {
    letters <- LETTERS[seq_along(panels)]
    for (i in seq_along(panels)) {
      save_nature_figure(panels[[i]], paste0(stem, "_panel", letters[i]),
                         width_mm = 89, height_mm = height_mm)
    }
    message("patchwork is not installed: saved ", length(panels), " panels separately (",
            paste(paste0(stem, "_panel", letters), collapse = ", "), ").")
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Fig 1 panels (01_preprocess)
# ---------------------------------------------------------------------------

panel_pca_after <- function(res_dir) {
  scores <- read_if_exists(file.path(res_dir, "01_preprocess", "pca_scores.csv"), "PCA scores")
  if (is.null(scores) || !all(c("PC1", "PC2", "group") %in% colnames(scores))) {
    if (!is.null(scores)) message("pca_scores.csv lacks PC1/PC2/group columns; skipping the PCA panel.")
    return(NULL)
  }
  p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, colour = group)) +
    ggplot2::geom_point(size = 2)
  p <- maybe_add_ellipses(p)
  p + scale_colour_groups() +
    ggplot2::labs(title = "PCA after normalisation",
                  subtitle = "Sum + log2 + pareto scaling separates the groups along PC1",
                  x = "PC1", y = "PC2") +
    theme_nature()
}

panel_boxplot <- function(res_dir, group_of) {
  mat <- read_if_exists(file.path(res_dir, "01_preprocess", "metabolome_preprocessed.csv"),
                        "preprocessed matrix")
  if (is.null(mat)) return(NULL)
  sample_cols <- setdiff(sample_columns(mat), ANNOTATION_COLS)
  if (!length(sample_cols)) {
    message("No numeric sample columns in metabolome_preprocessed.csv; skipping the boxplot panel.")
    return(NULL)
  }
  vals <- as.matrix(mat[, sample_cols, drop = FALSE])
  long <- data.frame(sample = rep(sample_cols, each = nrow(vals)),
                     value = as.numeric(vals), stringsAsFactors = FALSE)
  long$group <- if (!is.null(group_of)) unname(group_of[long$sample]) else "all"
  ggplot2::ggplot(long, ggplot2::aes(x = sample, y = value, fill = group)) +
    ggplot2::geom_boxplot(width = 0.6, outlier.size = 0.3, outlier.alpha = 0.4) +
    scale_fill_groups() +
    ggplot2::labs(title = "Normalised abundance",
                  subtitle = "Pareto scaling equalises metabolite variance across samples",
                  x = NULL, y = "Normalised value") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

panel_missing <- function(res_dir, group_of) {
  summary <- read_if_exists(file.path(res_dir, "01_preprocess", "missing_value_summary.csv"),
                            "missing-value summary")
  if (is.null(summary) || !all(c("level", "id", "missing_ratio") %in% colnames(summary))) {
    if (!is.null(summary)) message("missing_value_summary.csv has unexpected columns; skipping.")
    return(NULL)
  }
  srows <- summary[summary$level == "sample", , drop = FALSE]
  if (!nrow(srows)) return(NULL)
  srows$group <- if (!is.null(group_of)) unname(group_of[srows$id]) else "all"
  ggplot2::ggplot(srows, ggplot2::aes(x = id, y = missing_ratio, fill = group)) +
    ggplot2::geom_col(width = 0.7) +
    scale_fill_groups() +
    ggplot2::labs(title = "Missing values per sample",
                  subtitle = "Missingness is metabolite-driven (MNAR); half-min imputation applied",
                  x = NULL, y = "Missing fraction") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# ---------------------------------------------------------------------------
# Fig 2 panels (02_oplsda)
# ---------------------------------------------------------------------------

panel_oplsda_scores <- function(res_dir) {
  scores <- read_if_exists(file.path(res_dir, "02_oplsda", "oplsda_scores.csv"), "OPLS-DA scores")
  if (is.null(scores) || !all(c("t1", "to1", "group") %in% colnames(scores))) {
    if (!is.null(scores)) message("oplsda_scores.csv lacks t1/to1/group; skipping the score panel.")
    return(NULL)
  }
  p <- ggplot2::ggplot(scores, ggplot2::aes(x = t1, y = to1, colour = group)) +
    ggplot2::geom_point(size = 2)
  p <- maybe_add_ellipses(p)
  p + scale_colour_groups() +
    ggplot2::labs(title = "OPLS-DA scores",
                  subtitle = "Groups separate on the predictive component t1; to1 = within-group spread",
                  x = "t1 (predictive)", y = "to1 (orthogonal)") +
    theme_nature()
}

panel_splot <- function(res_dir) {
  loadings <- read_if_exists(file.path(res_dir, "02_oplsda", "oplsda_loadings.csv"),
                             "OPLS-DA loadings")
  vip <- read_if_exists(file.path(res_dir, "02_oplsda", "vip_scores.csv"), "VIP scores")
  if (is.null(loadings) || is.null(vip) ||
      !all(c("metabolite_id", "p1", "p_corr") %in% colnames(loadings))) {
    if (!is.null(loadings) && !is.null(vip)) {
      message("oplsda_loadings.csv/vip_scores.csv columns incomplete; skipping the S-plot panel.")
    }
    return(NULL)
  }
  df <- merge(loadings, vip[, c("metabolite_id", "VIP", "log2FC", "significant"), drop = FALSE],
              by = "metabolite_id", all.x = TRUE)
  df$significant <- as_significant(df$significant)
  df$log2FC[is.na(df$log2FC)] <- 0
  df$cls <- ifelse(!df$significant, "ns", ifelse(df$log2FC >= 0, "up", "down"))
  df$cls <- factor(df$cls, levels = c("up", "down", "ns"))
  col_map <- c(up = nature_semantic("up"), down = nature_semantic("down"),
               ns = nature_semantic("not_significant"))
  top_idx <- head(order(-df$VIP), 8)
  ggplot2::ggplot(df, ggplot2::aes(x = p1, y = p_corr, colour = cls)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, colour = nature_semantic("reference_line"),
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, colour = nature_semantic("reference_line"),
                        linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = col_map, aesthetics = "colour", name = "Class",
                                 labels = c(up = "up", down = "down", ns = "ns")) +
    ggrepel::geom_label_repel(data = df[top_idx, ], ggplot2::aes(label = metabolite_id),
                              size = 1.8, colour = "grey20", max.overlaps = 20,
                              min.segment.length = 0) +
    ggplot2::labs(title = "S-plot",
                  subtitle = "Outer arms = high covariance + reliable correlation (biomarker candidates)",
                  x = "p(1) (loading)", y = "p(corr)") +
    theme_nature()
}

panel_vip <- function(res_dir, vip_cut = 1) {
  vip <- read_if_exists(file.path(res_dir, "02_oplsda", "vip_scores.csv"), "VIP scores")
  if (is.null(vip) || !all(c("metabolite_id", "VIP") %in% colnames(vip))) {
    if (!is.null(vip)) message("vip_scores.csv lacks metabolite_id/VIP; skipping the VIP panel.")
    return(NULL)
  }
  df <- utils::head(vip[order(-vip$VIP), ], 30)
  df$metabolite_id <- factor(df$metabolite_id, levels = rev(df$metabolite_id))
  df$above <- df$VIP >= vip_cut
  ggplot2::ggplot(df, ggplot2::aes(x = metabolite_id, y = VIP, fill = above)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = vip_cut, colour = nature_semantic("reference_line"),
                        linewidth = 0.5, linetype = "dashed") +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = nature_semantic("up"), `FALSE` = nature_semantic("not_significant")),
      name = NULL, labels = c(`TRUE` = paste0("VIP >= ", vip_cut),
                              `FALSE` = paste0("VIP < ", vip_cut))) +
    ggplot2::labs(title = "Top 30 VIP",
                  subtitle = "VIP >= 1 metabolites drive the group separation",
                  x = NULL, y = "VIP") +
    theme_nature()
}

# ---------------------------------------------------------------------------
# Fig 3 panels (05_biomarker)
# ---------------------------------------------------------------------------

panel_roc <- function(res_dir, n_show = 5) {
  roc <- read_if_exists(file.path(res_dir, "05_biomarker", "biomarker_roc.csv"), "biomarker ROC table")
  coords <- read_if_exists(file.path(res_dir, "05_biomarker", "roc_coordinates.csv"),
                           "ROC coordinates")
  if (is.null(roc) || is.null(coords) ||
      !all(c("metabolite_id", "AUC") %in% colnames(roc)) ||
      !all(c("metabolite_id", "specificity", "sensitivity") %in% colnames(coords))) {
    if (!is.null(roc) && !is.null(coords)) {
      message("biomarker_roc.csv/roc_coordinates.csv columns incomplete; skipping the ROC panel.")
    }
    return(NULL)
  }
  top <- utils::head(roc[order(-roc$AUC), ], min(n_show, nrow(roc)))
  sub <- coords[coords$metabolite_id %in% top$metabolite_id, , drop = FALSE]
  auc_map <- setNames(round(top$AUC, 3), top$metabolite_id)
  sub$label <- paste0(sub$metabolite_id, " (AUC ", auc_map[sub$metabolite_id], ")")
  ggplot2::ggplot(sub, ggplot2::aes(x = 1 - specificity, y = sensitivity, colour = label)) +
    ggplot2::geom_path(linewidth = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = nature_semantic("reference_line"),
                         linewidth = 0.4, linetype = "dashed") +
    scale_colour_nature(name = "Metabolite") +
    ggplot2::labs(title = "Biomarker ROC curves",
                  subtitle = "AUC near 1 = near-perfect single-metabolite discrimination",
                  x = "1 - specificity", y = "Sensitivity") +
    ggplot2::coord_equal() +
    theme_nature()
}

panel_rf <- function(res_dir, top_n = 20) {
  imp <- read_if_exists(file.path(res_dir, "05_biomarker", "rf_importance.csv"),
                        "random-forest importance")
  if (is.null(imp) || !all(c("metabolite_id", "MeanDecreaseGini") %in% colnames(imp))) {
    if (!is.null(imp)) message("rf_importance.csv columns incomplete; skipping the RF panel.")
    return(NULL)
  }
  df <- utils::head(imp[order(-imp$MeanDecreaseGini), ], top_n)
  df$metabolite_id <- factor(df$metabolite_id, levels = rev(df$metabolite_id))
  ggplot2::ggplot(df, ggplot2::aes(x = metabolite_id, y = MeanDecreaseGini)) +
    ggplot2::geom_col(fill = rice_palette()$series12[1], width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Random forest importance",
                  subtitle = "Multivariate ranking (MeanDecreaseGini) of all metabolites",
                  x = NULL, y = "MeanDecreaseGini") +
    theme_nature()
}

# ---------------------------------------------------------------------------
# Fig 4 (06_network): rebuilt from the edge/node tables
# ---------------------------------------------------------------------------

build_network_figure <- function(res_dir, outdir) {
  edges <- read_if_exists(file.path(res_dir, "06_network", "metabolite_edges.csv"), "network edges")
  nodes <- read_if_exists(file.path(res_dir, "06_network", "metabolite_nodes.csv"), "network nodes")
  if (is.null(edges) || is.null(nodes)) return(invisible(NULL))
  if (!nrow(edges) || !all(c("from", "to", "r") %in% colnames(edges))) {
    message("Empty or malformed edge table; skipping fig4_network_overview.")
    return(invisible(NULL))
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("igraph is not installed; skipping fig4_network_overview.")
    return(invisible(NULL))
  }
  g <- igraph::graph_from_data_frame(edges[, c("from", "to"), drop = FALSE],
                                     directed = FALSE,
                                     vertices = data.frame(name = nodes$metabolite_id,
                                                           stringsAsFactors = FALSE))
  v_match <- match(igraph::V(g)$name, nodes$metabolite_id)
  g <- igraph::set_vertex_attr(g, "degree", value = nodes$degree[v_match])
  g <- igraph::set_vertex_attr(g, "is_hub", value = nodes$is_hub[v_match])
  igraph::E(g)$r <- edges$r
  lay <- igraph::layout_with_fr(g)
  stem <- file.path(outdir, "fig4_network_overview")
  if (requireNamespace("ggraph", quietly = TRUE)) {
    lay_df <- ggraph::create_layout(g, layout = "fr")
    p <- ggraph::ggraph(lay_df) +
      ggraph::geom_edge_link(ggplot2::aes(colour = r), alpha = 0.55, linewidth = 0.4) +
      ggraph::geom_node_point(ggplot2::aes(size = degree), colour = rice_palette()$series12[1]) +
      ggraph::geom_node_text(ggplot2::aes(label = ifelse(is_hub, name, "")),
                             size = 2.0, colour = "grey20", repel = TRUE) +
      ggplot2::scale_colour_gradient2(low = nature_semantic("down"),
                                      mid = diverging_palette()[4],
                                      high = nature_semantic("up"),
                                      name = "Correlation (r)") +
      ggplot2::scale_size_continuous(range = c(1.5, 6), name = "Degree") +
      ggplot2::labs(title = "Metabolite correlation network",
                    subtitle = paste0("|r| >= threshold and padj-significant; ",
                                      igraph::vcount(g), " metabolites, ",
                                      igraph::ecount(g), " edges; hubs labelled"),
                    x = NULL, y = NULL) +
      theme_nature() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    save_nature_figure(p, stem, width_mm = 89, height_mm = 89)
  } else {
    deg <- nodes$degree
    rng <- max(deg) - min(deg)
    v_size <- if (rng > 0) 2 + 8 * (deg - min(deg)) / rng else rep(4, length(deg))
    hub <- nodes$is_hub
    for (fmt in c("pdf", "svg", "tiff")) {
      out <- paste0(stem, ".", fmt)
      if (fmt == "tiff") {
        grDevices::tiff(out, width = mm_to_inch(89), height = mm_to_inch(89),
                        units = "in", res = 300, compression = "lzw")
      } else if (fmt == "pdf") {
        grDevices::pdf(out, width = mm_to_inch(89), height = mm_to_inch(89))
      } else {
        grDevices::svg(out, width = mm_to_inch(89), height = mm_to_inch(89))
      }
      igraph::plot.igraph(g, layout = lay,
                   vertex.size = v_size, vertex.color = rice_palette()$series12[1],
                   vertex.frame.color = NA,
                   vertex.label = ifelse(hub, nodes$metabolite_id, NA),
                   vertex.label.cex = 0.55, vertex.label.dist = 0.6,
                   vertex.label.color = "grey20",
                   edge.color = ifelse(igraph::E(g)$r > 0,
                                       nature_semantic("up"), nature_semantic("down")),
                   edge.width = 0.5, main = "Metabolite correlation network")
      grDevices::dev.off()
    }
    message("ggraph is not installed: fig4_network_overview saved with the igraph base device.")
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  res_dir <- opt(opts, "results-dir", "results/metabolome")
  outdir <- opt(opts, "outdir", file.path(res_dir, "07_figures"))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Results directory: ", normalizePath(res_dir))
  message("Output directory: ", normalizePath(outdir))

  # Group lookup (for colouring the fig1 panels), taken from the PCA scores.
  pca <- read_if_exists(file.path(res_dir, "01_preprocess", "pca_scores.csv"), "PCA scores")
  group_of <- NULL
  if (!is.null(pca) && all(c("sample_id", "group") %in% colnames(pca))) {
    group_of <- setNames(as.character(pca$group), pca$sample_id)
  }

  # ---- Fig 1: preprocessing overview -------------------------------------------
  fig1_panels <- list(panel_pca_after(res_dir), panel_boxplot(res_dir, group_of),
                      panel_missing(res_dir, group_of))
  if (any(!vapply(fig1_panels, is.null, logical(1)))) {
    save_panels(fig1_panels, file.path(outdir, "fig1_preprocess_overview"),
                width_mm = 180, height_mm = 78, ncol = 3, widths = c(1, 1.15, 1))
  } else {
    message("01_preprocess results not found; fig1_preprocess_overview skipped.")
  }

  # ---- Fig 2: OPLS-DA overview ---------------------------------------------------
  fig2_panels <- list(panel_oplsda_scores(res_dir), panel_splot(res_dir),
                      panel_vip(res_dir))
  if (any(!vapply(fig2_panels, is.null, logical(1)))) {
    save_panels(fig2_panels, file.path(outdir, "fig2_oplsda_overview"),
                width_mm = 180, height_mm = 88, ncol = 3, widths = c(1, 1, 1.1))
  } else {
    message("02_oplsda results not found; fig2_oplsda_overview skipped.")
  }

  # ---- Fig 3: biomarker overview ---------------------------------------------------
  fig3_panels <- list(panel_roc(res_dir), panel_rf(res_dir))
  if (any(!vapply(fig3_panels, is.null, logical(1)))) {
    save_panels(fig3_panels, file.path(outdir, "fig3_biomarker_overview"),
                width_mm = 180, height_mm = 80, ncol = 2, widths = c(1, 1))
  } else {
    message("05_biomarker results not found; fig3_biomarker_overview skipped.")
  }

  # ---- Fig 4: network overview ------------------------------------------------------
  build_network_figure(res_dir, outdir)

  write_run_info_r(
    outdir,
    script_name = "07_plot_metabolome.R",
    params = list(results_dir = res_dir, outdir = outdir),
    packages = c("ggplot2", "ggrepel", "patchwork", "igraph", "ggraph"),
    extra = c(
      "Figures are re-drawn from the CSV outputs of steps 01-06; missing steps are skipped.",
      "Multi-panel figures use patchwork when installed (otherwise panels are saved separately)."
    )
  )
  message("Done. Results in ", normalizePath(outdir))
  invisible(NULL)
}

main()
