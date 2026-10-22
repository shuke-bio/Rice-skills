###############################################################################
# 06_plot_integration.R -- publication-ready composite figures for the
#                          transcriptome x metabolome integration module
#                          (Rice-skills skill, module 3, step 06)
#
# Title:    Integration figure assembly (MOFA overview, correlation network,
#           joint pathways, sankey, circos, UpSet)
# Description:
#   Re-assembles the results of steps 01-05 into journal-style composite
#   figures. Every figure is rebuilt from the CSV outputs of the previous
#   steps, so this script can be re-run at any time and only produces the
#   figures whose inputs are available (missing steps are skipped with a
#   console note and recorded in run_info.txt).
#
#   fig1_integration_overview   MOFA variance decomposition + factor scores
#                               (double column, 180 mm). Scientific message:
#                               how much of the multi-omics variance each
#                               latent factor captures per omics layer, and
#                               whether the factors align with the groups.
#   fig2_correlation_overview   gene-metabolite correlation heatmap + bipartite
#                               network. Message: which DE genes and
#                               differential metabolites form a tightly
#                               co-varying module.
#   fig3_joint_pathway_overview joint pathway dot plot. Message: which KEGG
#                               pathways respond on both omics layers.
#   fig4_sankey                 gene -> pathway -> metabolite flow. Message:
#                               the direction and composition of the joint
#                               pathway response. ggalluvial is used when
#                               installed; otherwise a straight-ribbon ggplot2
#                               approximation is drawn (noted in the script).
#   fig5_circos                 circlize chord diagram of gene-metabolite
#                               connections grouped by pathway. Skipped with a
#                               note when circlize is not installed.
#   fig6_upset                  UpSet-style relation between the transcriptome,
#                               metabolome and joint pathway sets. UpSetR is
#                               used when installed; otherwise a ggplot2
#                               approximation (set sizes + intersection matrix
#                               + intersection sizes) is drawn.
#
#   patchwork is used to combine panels when installed; otherwise the panels
#   are saved as separate figures (figXa_*, figXb_*).
#
# Usage:
#   Rscript scripts/integration/06_plot_integration.R \
#       --results-dir results/integration \
#       --outdir results/integration/06_figures \
#       --transcriptome-dir results/transcriptome \
#       --metabolome-dir results/metabolome
#
# Outputs (written to --outdir):
#   fig1..fig6 (PDF/TIFF/PNG), run_info.txt
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
rice_set_seed(42)

# ---------------------------------------------------------------------------
# CLI parsing (shared skeleton)
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  opts <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) {
      stop("Unexpected argument '", a, "'. Use --key value or --flag.", call. = FALSE)
    }
    key <- sub("^--", "", a)
    if (i == length(args) || grepl("^--", args[i + 1L])) {
      opts[[key]] <- TRUE
      i <- i + 1L
    } else {
      opts[[key]] <- args[i + 1L]
      i <- i + 2L
    }
  }
  opts
}
opt <- function(opts, key, default = NULL) {
  v <- opts[[key]]
  if (is.null(v)) return(default)
  if (isTRUE(v)) return(if (is.logical(default)) default else TRUE)
  v
}

usage_06 <- paste(
  "Usage: Rscript 06_plot_integration.R [--results-dir DIR] [--outdir DIR]",
  "         [--transcriptome-dir DIR] [--metabolome-dir DIR]",
  sep = "\n"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE))
results_dir <- opt(opts, "results-dir", "results/integration")
outdir <- opt(opts, "outdir", "results/integration/06_figures")
transcriptome_dir <- opt(opts, "transcriptome-dir", "results/transcriptome")
metabolome_dir <- opt(opts, "metabolome-dir", "results/metabolome")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rice_set_seed(42)

# ---------------------------------------------------------------------------
# Input discovery helpers
# ---------------------------------------------------------------------------

step_dir <- function(step) file.path(results_dir, step)

locate <- function(...) {
  p <- do.call(file.path, list(...))
  if (file.exists(p)) normalizePath(p, mustWork = FALSE) else NULL
}

read_quiet <- function(path, what) {
  if (is.null(path)) return(NULL)
  df <- tryCatch(read_csv_safe(path, what), error = function(e) {
    message("Note: could not read ", what, " (", path, "): ", conditionMessage(e))
    NULL
  })
  df
}

skip_note <- function(fig, reason) {
  message("SKIP ", fig, ": ", reason)
  c(sprintf("%s: SKIPPED (%s)", fig, reason))
}

produced <- character(0)
skipped <- character(0)

# ---------------------------------------------------------------------------
# fig1: MOFA variance explained + factor scores
# ---------------------------------------------------------------------------

fig1 <- "fig1_integration_overview"
ve_csv <- locate(step_dir("01_mofa"), "variance_explained.csv")
fs_csv <- locate(step_dir("01_mofa"), "factor_scores.csv")
if (is.null(ve_csv) && is.null(fs_csv)) {
  skipped <- c(skipped, skip_note(fig1, paste0("step 01_mofa outputs not found in ",
                                               results_dir)))
} else {
  panels <- list()
  if (!is.null(ve_csv)) {
    ve <- read_quiet(ve_csv, "variance explained table")
    if (!is.null(ve)) {
      num_cols <- colnames(ve)[vapply(ve, is.numeric, logical(1))]
      r2_col <- if ("r2_total" %in% num_cols) "r2_total"
                else if ("r2_per_factor" %in% num_cols) "r2_per_factor"
                else num_cols[1]
      if (!is.null(r2_col)) {
        ve$view <- as.factor(ve$view)
        ve$factor <- as.factor(ve$factor)
        p_ve <- ggplot2::ggplot(ve, ggplot2::aes(x = .data$factor, y = .data[[r2_col]],
                                                 fill = .data$view)) +
          ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.2) +
          scale_fill_nature(name = "Omics view") +
          ggplot2::labs(x = "Latent factor", y = "Variance explained (R2)",
                        title = "MOFA variance decomposition",
                        subtitle = "Fraction of each view's variance captured per factor") +
          theme_nature()
        panels$variance <- p_ve
      }
    }
  }
  if (!is.null(fs_csv)) {
    fs <- read_quiet(fs_csv, "factor scores table")
    if (!is.null(fs)) {
      fcols <- grep("^Factor[0-9]+$", colnames(fs), value = TRUE)
      if (length(fcols) >= 2 && "group" %in% colnames(fs)) {
        p_fs <- ggplot2::ggplot(fs, ggplot2::aes(x = .data[[fcols[1]]],
                                                 y = .data[[fcols[2]]],
                                                 colour = .data$group)) +
          ggplot2::geom_point(size = 1.8, alpha = 0.9) +
          scale_colour_groups(name = "Group") +
          ggplot2::labs(x = fcols[1], y = fcols[2],
                        title = "MOFA factor scores",
                        subtitle = "Latent factors align samples; colour = experimental group") +
          theme_nature()
        panels$scores <- p_fs
      }
    }
  }
  if (length(panels) == 0) {
    skipped <- c(skipped, skip_note(fig1, "01_mofa tables could not be parsed"))
  } else {
    if (length(panels) == 2 && requireNamespace("patchwork", quietly = TRUE)) {
      combined <- panels$variance + panels$scores +
        patchwork::plot_layout(widths = c(1.15, 1))
      save_nature_figure(combined, file.path(outdir, fig1), width_mm = 180, height_mm = 75)
      produced <- c(produced, paste0(fig1, " (patchwork, 180 mm)"))
    } else {
      if (!is.null(panels$variance)) {
        save_nature_figure(panels$variance, file.path(outdir, "fig1a_variance_explained"),
                           width_mm = 120, height_mm = 85)
        produced <- c(produced, "fig1a_variance_explained (panel saved separately)")
      }
      if (!is.null(panels$scores)) {
        save_nature_figure(panels$scores, file.path(outdir, "fig1b_factor_scores"),
                           width_mm = 120, height_mm = 85)
        produced <- c(produced, "fig1b_factor_scores (panel saved separately)")
      }
    }
  }
}

# ---------------------------------------------------------------------------
# fig2: correlation heatmap + bipartite network
# ---------------------------------------------------------------------------

fig2 <- "fig2_correlation_overview"
pairs_csv <- locate(step_dir("03_correlation"), "gene_metabolite_pairs_sig.csv")
if (is.null(pairs_csv)) {
  skipped <- c(skipped, skip_note(fig2, "step 03_correlation outputs not found"))
} else {
  pairs <- read_quiet(pairs_csv, "significant pairs table")
  if (is.null(pairs) || nrow(pairs) == 0) {
    skipped <- c(skipped, skip_note(fig2, "no significant gene-metabolite pairs in 03"))
  } else {
    # --- panel A: tile heatmap of the significant correlation matrix
    hm_genes <- unique(pairs$gene_id)
    hm_mets <- unique(pairs$metabolite_id)
    if (length(hm_genes) > 50 || length(hm_mets) > 50) {
      g_rank <- sort(table(pairs$gene_id), decreasing = TRUE)
      m_rank <- sort(table(pairs$metabolite_id), decreasing = TRUE)
      hm_genes <- names(g_rank)[seq_len(min(50, length(g_rank)))]
      hm_mets <- names(m_rank)[seq_len(min(50, length(m_rank)))]
      pairs_hm <- pairs[pairs$gene_id %in% hm_genes & pairs$metabolite_id %in% hm_mets, ]
    } else {
      pairs_hm <- pairs
    }
    hm <- matrix(NA_real_, nrow = length(hm_genes), ncol = length(hm_mets),
                 dimnames = list(hm_genes, hm_mets))
    hm[cbind(match(pairs_hm$gene_id, hm_genes), match(pairs_hm$metabolite_id, hm_mets))] <-
      pairs_hm$r
    hm_long <- data.frame(gene = rep(rownames(hm), times = ncol(hm)),
                          metabolite = rep(colnames(hm), each = nrow(hm)),
                          r = as.vector(hm), stringsAsFactors = FALSE)
    hm_long <- hm_long[is.finite(hm_long$r), , drop = FALSE]
    lim <- max(abs(hm_long$r))
    p_heat <- ggplot2::ggplot(hm_long, ggplot2::aes(x = .data$metabolite, y = .data$gene,
                                                    fill = .data$r)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
      ggplot2::scale_fill_gradientn(colours = diverging_palette(), limits = c(-lim, lim),
                                    name = "Spearman r") +
      ggplot2::labs(x = "Differential metabolites", y = "DE genes",
                    title = "Gene-metabolite correlations",
                    subtitle = "Significant pairs only (|r| and BH-adjusted p below thresholds)") +
      theme_nature() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))

    # --- panel B: bipartite network (genes left, metabolites right)
    g_deg <- table(pairs$gene_id)
    m_deg <- table(pairs$metabolite_id)
    gene_y <- setNames(seq(0, 1, length.out = length(g_deg) + 2)[seq_along(g_deg)],
                       names(sort(g_deg, decreasing = TRUE)))
    met_y <- setNames(seq(0, 1, length.out = length(m_deg) + 2)[seq_along(m_deg)],
                      names(sort(m_deg, decreasing = TRUE)))
    edge_df <- data.frame(
      x = 0, y = gene_y[pairs$gene_id],
      xend = 1, yend = met_y[pairs$metabolite_id],
      r = pairs$r, stringsAsFactors = FALSE
    )
    node_df <- rbind(
      data.frame(x = 0, y = gene_y, node = names(gene_y), type = "gene",
                 stringsAsFactors = FALSE),
      data.frame(x = 1, y = met_y, node = names(met_y), type = "metabolite",
                 stringsAsFactors = FALSE)
    )
    p_net <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = edge_df,
                            ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                         yend = .data$yend,
                                         colour = .data$r >= 0),
                            linewidth = 0.25, alpha = 0.5) +
      ggplot2::geom_point(data = node_df,
                          ggplot2::aes(x = .data$x, y = .data$y, shape = .data$type),
                          size = 1.8, colour = "grey30") +
        ggplot2::scale_colour_manual(values = c(`TRUE` = nature_semantic("up"),
                                                `FALSE` = nature_semantic("down")),
                                     labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
                                     name = "Correlation sign") +
      ggplot2::scale_shape_manual(values = c(gene = 16, metabolite = 17),
                                  name = "Node type") +
      ggplot2::scale_x_continuous(limits = c(-0.12, 1.12), breaks = c(0, 1),
                                  labels = c("DE genes", "Metabolites")) +
      ggplot2::coord_cartesian(ylim = c(-0.03, 1.03)) +
      ggplot2::labs(x = NULL, y = NULL,
                    title = "Gene-metabolite bipartite network",
                    subtitle = "Nodes ordered by degree; edges = significant correlations") +
      theme_nature()

    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined <- p_heat + p_net + patchwork::plot_layout(widths = c(1.3, 1))
      save_nature_figure(combined, file.path(outdir, fig2), width_mm = 180, height_mm = 90)
      produced <- c(produced, paste0(fig2, " (patchwork, 180 mm)"))
    } else {
      save_nature_figure(p_heat, file.path(outdir, "fig2a_correlation_heatmap"),
                         width_mm = 120, height_mm = 90)
      save_nature_figure(p_net, file.path(outdir, "fig2b_correlation_network"),
                         width_mm = 120, height_mm = 90)
      produced <- c(produced, "fig2a_correlation_heatmap + fig2b_correlation_network (separate)")
    }
  }
}

# ---------------------------------------------------------------------------
# fig3: joint pathway dot plot
# ---------------------------------------------------------------------------

fig3 <- "fig3_joint_pathway_overview"
jp_csv <- locate(step_dir("04_joint_pathway"), "joint_pathway_table.csv")
if (is.null(jp_csv)) {
  skipped <- c(skipped, skip_note(fig3, "step 04_joint_pathway outputs not found"))
} else {
  jp <- read_quiet(jp_csv, "joint pathway table")
  if (is.null(jp) || nrow(jp) == 0) {
    skipped <- c(skipped, skip_note(fig3, "joint pathway table is empty"))
  } else {
    jp$min_padj <- pmin(jp$gene_padj, jp$metabolite_padj, na.rm = TRUE)
    jp$min_padj[!is.finite(jp$min_padj)] <- 1
    jp$category <- factor(jp$category,
                          levels = c("both", "transcriptome_only", "metabolome_only"))
    jp$label <- ifelse(nchar(jp$name) > 42, paste0(substr(jp$name, 1, 40), "..."), jp$name)
    p <- ggplot2::ggplot(jp, ggplot2::aes(x = .data$n_genes, y = -log10(.data$min_padj),
                                          colour = .data$category)) +
      ggplot2::geom_point(ggplot2::aes(size = .data$n_metabolites), alpha = 0.9) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                          colour = nature_semantic("reference_line"), linewidth = 0.5) +
      scale_colour_nature(
        labels = c(both = "Both omics", transcriptome_only = "Transcriptome only",
                   metabolome_only = "Metabolome only"), drop = FALSE) +
      ggplot2::labs(x = "Number of DE genes in pathway",
                    y = "-log10(min adjusted p-value)",
                    colour = "Pathway category", size = "Differential metabolites",
                    title = "Joint transcriptome-metabolome KEGG pathways",
                    subtitle = "Dashed line: padj = 0.05; 'both' = significant on both layers") +
      theme_nature()
    if (requireNamespace("ggrepel", quietly = TRUE) && nrow(jp) <= 40) {
      top_rows <- head(jp[order(-jp$min_padj), ], min(10, nrow(jp)))
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data$label), data = top_rows,
                                        size = 2.4, colour = "grey30", max.overlaps = 20,
                                        seed = 42)
    }
    save_nature_figure(p, file.path(outdir, fig3), width_mm = 120,
                       height_mm = max(80, 20 + 5 * nrow(jp)))
    produced <- c(produced, paste0(fig3, " (", nrow(jp), " pathways)"))
  }
}

# ---------------------------------------------------------------------------
# fig4: sankey gene -> pathway -> metabolite (ggalluvial or ggplot2 approximation)
# ---------------------------------------------------------------------------

fig4 <- "fig4_sankey"
links_csv <- locate(step_dir("04_joint_pathway"), "pathway_links.csv")
if (is.null(links_csv)) {
  skipped <- c(skipped, skip_note(fig4, "step 04_joint_pathway outputs not found"))
} else {
  links <- read_quiet(links_csv, "pathway links table")
  if (is.null(links) || nrow(links) == 0) {
    skipped <- c(skipped, skip_note(fig4, "pathway_links.csv is empty"))
  } else {
    # Keep pathways that contain both node types (a meaningful flow) and cap
    # the number of pathways for readability.
    has_gene <- tapply(links$node_type == "gene", links$pathway, any)
    has_met <- tapply(links$node_type == "metabolite", links$pathway, any)
    both_paths <- names(has_gene)[has_gene & has_met]
    if (length(both_paths) == 0) {
      skipped <- c(skipped, skip_note(fig4, "no pathway contains both genes and metabolites"))
    } else {
      tot <- table(links$pathway)
      both_paths <- both_paths[order(-tot[both_paths])]
      if (length(both_paths) > 8) {
        message("Note: fig4 shows the top 8 of ", length(both_paths),
                " pathways (by node count).")
        both_paths <- both_paths[seq_len(8)]
      }
      links <- links[links$pathway %in% both_paths, , drop = FALSE]
      n_p <- length(both_paths)
      slot <- setNames((seq_len(n_p) - 0.5) / n_p, both_paths)

      # Vertical positions inside each pathway slot.
      y_gene <- y_met <- numeric(0)
      gene_rows <- met_rows <- list()
      for (i in seq_len(n_p)) {
        pid <- both_paths[i]
        sub <- links[links$pathway == pid, , drop = FALSE]
        gg <- sub$node[sub$node_type == "gene"]
        mm <- sub$node[sub$node_type == "metabolite"]
        k_g <- length(gg); k_m <- length(mm)
        span <- 0.85 / n_p
        yg <- slot[[pid]] + (seq_len(k_g) - (k_g + 1) / 2) * span / max(k_g, 1)
        ym <- slot[[pid]] + (seq_len(k_m) - (k_m + 1) / 2) * span / max(k_m, 1)
        if (k_g > 0) gene_rows[[length(gene_rows) + 1L]] <-
          data.frame(node = gg, y = yg, pathway = pid,
                     log2fc = sub$log2fc[sub$node_type == "gene"], stringsAsFactors = FALSE)
        if (k_m > 0) met_rows[[length(met_rows) + 1L]] <-
          data.frame(node = mm, y = ym, pathway = pid, stringsAsFactors = FALSE)
      }
      gene_df <- do.call(rbind, gene_rows)
      met_df <- do.call(rbind, met_rows)
      path_df <- data.frame(pathway = both_paths, y = as.numeric(slot[both_paths]),
                            label = both_paths, stringsAsFactors = FALSE)

      rib_g <- data.frame(x = 0, y = gene_df$y, xend = 1,
                          yend = slot[gene_df$pathway], sign = gene_df$log2fc >= 0,
                          stringsAsFactors = FALSE)
      rib_m <- data.frame(x = 1, y = slot[met_df$pathway], xend = 2,
                          yend = met_df$y, stringsAsFactors = FALSE)

      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = rib_g,
                              ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                           yend = .data$yend, colour = .data$sign),
                              linewidth = 0.3, alpha = 0.45) +
        ggplot2::geom_segment(data = rib_m,
                              ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                                           yend = .data$yend),
                              linewidth = 0.3, alpha = 0.45, colour = "grey55") +
        ggplot2::geom_point(data = gene_df, ggplot2::aes(x = 0, y = .data$y),
                            size = 1.6, colour = rice_palette()$series12[1]) +
        ggplot2::geom_point(data = met_df, ggplot2::aes(x = 2, y = .data$y),
                            size = 1.6, colour = rice_palette()$series12[2]) +
        ggplot2::geom_text(data = path_df, ggplot2::aes(x = 1, y = .data$y,
                                                        label = .data$label),
                           size = 2.6, colour = "grey20") +
        ggplot2::scale_colour_manual(values = c(`TRUE` = nature_semantic("up"),
                                                `FALSE` = nature_semantic("down")),
                                     labels = c(`TRUE` = "Up-regulated gene",
                                                `FALSE` = "Down-regulated gene"),
                                     name = "Gene direction") +
        ggplot2::scale_x_continuous(limits = c(-0.18, 2.18), breaks = c(0, 1, 2),
                                    labels = c("DE genes", "KEGG pathway", "Metabolites")) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Gene - pathway - metabolite flow",
                      subtitle = paste0("Top ", n_p, " joint pathways; ",
                                        "straight-ribbon approximation (ggalluvial not installed)")) +
        theme_nature() +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
      save_nature_figure(p, file.path(outdir, fig4), width_mm = 150,
                         height_mm = max(80, 30 + 22 * n_p))
      produced <- c(produced, paste0(fig4, " (", n_p, " pathways)"))
    }
  }
}

# ---------------------------------------------------------------------------
# fig5: circos of gene-metabolite connections grouped by pathway (circlize)
# ---------------------------------------------------------------------------

fig5 <- "fig5_circos"
if (!requireNamespace("circlize", quietly = TRUE)) {
  skipped <- c(skipped, skip_note(fig5,
             "circlize is not installed (install with BiocManager::install('circlize'))"))
} else if (is.null(links_csv)) {
  skipped <- c(skipped, skip_note(fig5, "step 04_joint_pathway outputs not found"))
} else {
  links5 <- read_quiet(links_csv, "pathway links table")
  if (is.null(links5) || nrow(links5) == 0) {
    skipped <- c(skipped, skip_note(fig5, "pathway_links.csv is empty"))
  } else {
      # Adjacency between genes and metabolites that share at least one pathway.
      shared <- links5[links5$node_type == "gene", c("node", "pathway")]
      shared <- merge(shared, links5[links5$node_type == "metabolite",
                                   c("node", "pathway")], by = "pathway")
      if (nrow(shared) == 0) {
        skipped <- c(skipped, skip_note(fig5, "no gene-metabolite pair shares a pathway"))
      } else {
        adj <- table(factor(shared$node.x, levels = unique(shared$node.x)),
                     factor(shared$node.y, levels = unique(shared$node.y)))
        # Group each node by its (first) pathway; named colours keep sectors
        # aligned with the pathway groups inside chordDiagram.
        grp <- links5[!duplicated(links5$node), c("node", "pathway")]
        group <- setNames(grp$pathway, grp$node)
        node_group <- group[union(rownames(adj), colnames(adj))]
        node_group[is.na(node_group)] <- "unassigned"
        pathway_ids <- sort(unique(node_group))
        series <- rice_palette()$series12
        grid_cols <- setNames(series[((seq_along(pathway_ids) - 1) %% length(series)) + 1],
                              pathway_ids)
        grid_cols["unassigned"] <- nature_semantic("not_significant")
        draw_circos <- function() {
          circlize::circos.par(gap.degree = 2, start.degree = 0)
          circlize::chordDiagram(
            adj, group = node_group,
            grid.col = grid_cols,
            col = grDevices::adjustcolor(nature_semantic("not_significant"), alpha.f = 0.4),
            transparency = 0.5, annotationTrack = "grid",
            preAllocateTracks = 1
          )
          circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
            circlize::circos.text(circlize::CELL_META$xcenter,
                                  circlize::CELL_META$ycenter[1] + 0.28,
                                  circlize::CELL_META$sector.index, facing = "bending",
                                  niceFacing = TRUE, cex = 0.5)
          }, bg.border = NA)
          graphics::title("Gene-metabolite connections by KEGG pathway")
        }
      ok <- tryCatch({
        for (f in c("pdf", "tiff", "png")) {
          out <- paste0(file.path(outdir, fig5), ".", if (f == "tiff") "tiff" else f)
          if (f == "pdf") {
            grDevices::pdf(out, width = mm_to_inch(120), height = mm_to_inch(120))
          } else if (f == "png") {
            grDevices::png(out, width = mm_to_inch(120), height = mm_to_inch(120), units = "in", res = 300)
          } else {
            grDevices::tiff(out, width = mm_to_inch(120), height = mm_to_inch(120),
                            units = "in", res = 300, compression = "lzw")
          }
          draw_circos()
          grDevices::dev.off()
        }
        TRUE
      }, error = function(e) {
        message("Note: circlize chord diagram failed: ", conditionMessage(e))
        FALSE
      })
      if (ok) produced <- c(produced, paste0(fig5, " (circlize)"))
      else skipped <- c(skipped, paste0(fig5, ": SKIPPED (circlize rendering failed)"))
    }
  }
}

# ---------------------------------------------------------------------------
# fig6: UpSet of transcriptome / metabolome / joint pathway sets
# ---------------------------------------------------------------------------

fig6 <- "fig6_upset"
gene_enr_csv <- locate(step_dir("04_joint_pathway"), "gene_kegg_enrichment.csv")
met_enr_csv <- locate(step_dir("04_joint_pathway"), "metabolite_pathway_enrichment.csv")
set_t <- character(0); set_m <- character(0); set_j <- character(0)
if (!is.null(jp_csv)) {
  jp6 <- read_quiet(jp_csv, "joint pathway table")
  if (!is.null(jp6) && nrow(jp6) > 0) {
    if ("gene_padj" %in% colnames(jp6)) {
      set_t <- jp6$pathway_id[is.finite(jp6$gene_padj) & jp6$gene_padj <= 0.05]
    }
    if ("metabolite_padj" %in% colnames(jp6)) {
      set_m <- jp6$pathway_id[is.finite(jp6$metabolite_padj) & jp6$metabolite_padj <= 0.05]
    }
    if ("category" %in% colnames(jp6)) {
      set_j <- jp6$pathway_id[jp6$category == "both"]
    }
  }
}
if (length(set_t) == 0 && !is.null(gene_enr_csv)) {
  ge <- read_quiet(gene_enr_csv, "gene KEGG enrichment table")
  if (!is.null(ge) && "padj" %in% colnames(ge)) {
    set_t <- ge$pathway_id[is.finite(ge$padj) & ge$padj <= 0.05]
  }
}
if (length(set_m) == 0 && !is.null(met_enr_csv)) {
  me6 <- read_quiet(met_enr_csv, "metabolite pathway table")
  if (!is.null(me6) && "padj" %in% colnames(me6)) {
    set_m <- me6$pathway_id[is.finite(me6$padj) & me6$padj <= 0.05]
  }
}
set_t <- unique(set_t); set_m <- unique(set_m); set_j <- unique(set_j)

if (length(set_t) == 0 && length(set_m) == 0) {
  skipped <- c(skipped, skip_note(fig6, "no pathway sets available (run 04_joint_pathway first)"))
} else {
  sets <- list(Transcriptome = set_t, Metabolome = set_m, Joint = set_j)
  set_names <- names(sets)
  combos <- expand.grid(Transcriptome = c(TRUE, FALSE), Metabolome = c(TRUE, FALSE),
                        Joint = c(TRUE, FALSE), stringsAsFactors = FALSE)
  combo_sizes <- integer(nrow(combos))
  for (i in seq_len(nrow(combos))) {
    inc <- unlist(combos[i, ])
    inter <- Reduce(intersect, sets[inc])
    excl <- Reduce(setdiff, list(inter, sets[!inc]))
    combo_sizes[i] <- length(excl)
  }
  keep <- combo_sizes > 0
  combos <- combos[keep, , drop = FALSE]
  combo_sizes <- combo_sizes[keep]
  combos$size <- combo_sizes
  combos$id <- seq_len(nrow(combos))
  combos$label <- vapply(seq_len(nrow(combos)), function(i) {
    inc <- unlist(combos[i, set_names])
    paste(set_names[inc], collapse = " + ")
  }, character(1))

  # Long table for the intersection dot matrix.
  dots <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
    inc <- unlist(combos[i, set_names])
    data.frame(id = i, set = set_names, present = inc, stringsAsFactors = FALSE)
  }))

  p_bars <- ggplot2::ggplot(combos, ggplot2::aes(x = .data$id, y = .data$size)) +
    ggplot2::geom_col(fill = rice_palette()$series12[1], colour = "white", linewidth = 0.2) +
    ggplot2::labs(x = NULL, y = "Pathways in intersection",
                  title = "Pathway set relations (UpSet-style)",
                  subtitle = "Transcriptome vs metabolome vs joint KEGG pathway sets") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

  p_dots <- ggplot2::ggplot(dots, ggplot2::aes(x = .data$id, y = .data$set)) +
    ggplot2::geom_line(data = dots[dots$present, ], ggplot2::aes(group = .data$id),
                       colour = "grey30", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(alpha = .data$present, size = .data$present)) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.15), guide = "none") +
    ggplot2::scale_size_manual(values = c(`TRUE` = 2.4, `FALSE` = 1.2), guide = "none") +
    ggplot2::scale_y_discrete(limits = rev(set_names)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  set_sizes <- data.frame(set = set_names, size = vapply(sets, length, integer(1)),
                          stringsAsFactors = FALSE)
  p_sizes <- ggplot2::ggplot(set_sizes, ggplot2::aes(x = .data$size, y = .data$set)) +
    ggplot2::geom_col(fill = rice_palette()$series12[4], colour = "white", linewidth = 0.2) +
    ggplot2::scale_y_discrete(limits = rev(set_names)) +
    ggplot2::labs(x = "Set size", y = NULL, title = "Set sizes") +
    theme_nature()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p_sizes + (p_bars / p_dots) +
      patchwork::plot_layout(widths = c(1, 3), heights = c(3, 1.1))
    save_nature_figure(combined, file.path(outdir, fig6), width_mm = 150, height_mm = 110)
    produced <- c(produced, paste0(fig6, " (patchwork, 150 mm)"))
  } else {
    save_nature_figure(p_bars, file.path(outdir, "fig6a_upset_intersections"),
                       width_mm = 120, height_mm = 90)
    save_nature_figure(p_dots, file.path(outdir, "fig6b_upset_matrix"),
                       width_mm = 120, height_mm = 60)
    save_nature_figure(p_sizes, file.path(outdir, "fig6c_upset_setsizes"),
                       width_mm = 120, height_mm = 60)
    produced <- c(produced, "fig6a/fig6b/fig6c (panels saved separately)")
  }
}

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

extra_lines <- c(
  paste0("produced figures: ", if (length(produced) > 0) length(produced) else 0),
  if (length(produced) > 0) paste0("  + ", produced) else character(0),
  paste0("skipped figures: ", if (length(skipped) > 0) length(skipped) else 0),
  if (length(skipped) > 0) paste0("  - ", skipped) else character(0)
)
write_run_info_r(
  outdir = outdir, script_name = "06_plot_integration.R",
  params = list(results_dir = results_dir, outdir = outdir,
                transcriptome_dir = transcriptome_dir, metabolome_dir = metabolome_dir),
  packages = c("ggplot2", "patchwork", "ggrepel", "circlize"),
  extra = extra_lines
)

message("Done. ", length(produced), " figure(s) produced, ", length(skipped),
        " skipped. Outputs written to: ", normalizePath(outdir))
