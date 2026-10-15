###############################################################################
# 07_plot_transcriptome.R -- publication figure assembly for the transcriptome
#
# Title:    Publication-grade transcriptome overview figures
# Description:
#   Re-draws combined overview figures from the results of the transcriptome
#   pipeline (steps 01/02/03/05/06). Steps whose outputs are missing are
#   skipped with a message, so the script can run on a partial results tree.
#   All panels are rebuilt from the saved tables (no intermediate R objects
#   are required), which keeps the figures reproducible from results/ alone.
#
# Usage:
#   Rscript scripts/transcriptomics/07_plot_transcriptome.R \
#       --results-dir results/transcriptome \
#       --outdir results/transcriptome/07_figures
#   # restrict to specific contrasts (default: auto-discovered de_sig_*.csv):
#   Rscript scripts/transcriptomics/07_plot_transcriptome.R --contrasts salt_vs_control
#
# Inputs (read from --results-dir):
#   01_qc/qc_metrics.csv, 01_qc/counts_filtered.csv
#   02_de/de_all_<contrast>.csv, 02_de/de_sig_<contrast>.csv,
#   02_de/de_summary.csv, 02_de/normalized_counts_vst.csv
#   03_enrichment/go_BP_<UP|DOWN>.csv, 03_enrichment/kegg_<UP|DOWN>.csv
#   05_wgcna/module_trait_correlation.csv, 05_wgcna/gene_modules.csv,
#   05_wgcna/hub_genes.csv
#   06_tf/tf_family_stats.csv
#
# Outputs (written to --outdir):
#   fig1_qc_overview.*          QC overview (library size, detected genes, PCA,
#                               sample correlation) - or separate panels
#   fig2_de_overview.*          volcano + MA + DE count panels per contrast
#   fig3_enrichment_overview.*  GO BP + KEGG dot-plot overview
#   fig4_wgcna_overview.*       module-trait heatmap + module sizes + hub genes
#   fig5_tf_overview.*          top-20 TF family bar chart
#   run_info.txt                provenance (parameters, versions, seed)
#
# Notes:
#   - patchwork is used to combine panels into single figures; when it is not
#     installed, each panel is saved as its own figure (fig1a_*, fig2a_*, ...)
#     and a message is printed.
#   - Volcano/MA thresholds are fixed to the pipeline defaults (padj 0.05,
#     |log2FC| 1); adjust the individual analysis scripts to change them.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
rice_set_seed(42)

# ---------------------------------------------------------------------------
# Command-line parsing (base R)
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

opts <- parse_args(commandArgs(trailingOnly = TRUE))

results_dir <- opt(opts, "results-dir", "results/transcriptome")
outdir <- opt(opts, "outdir", "results/transcriptome/07_figures")
contrast_arg <- opt(opts, "contrasts")

need_pkgs(c("ggplot2", "ggrepel", "jsonlite"))
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!has_patchwork) {
  message("NOTE: patchwork is not installed; combined figures will be saved ",
          "as separate panels (fig<N><panel>_*).")
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

qc_dir <- file.path(results_dir, "01_qc")
de_dir <- file.path(results_dir, "02_de")
enr_dir <- file.path(results_dir, "03_enrichment")
wgcna_dir <- file.path(results_dir, "05_wgcna")
tf_dir <- file.path(results_dir, "06_tf")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Read a count matrix CSV (first column = gene_id) into a numeric matrix.
read_counts_matrix <- function(path) {
  df <- read_csv_safe(path, "count matrix")
  gene_ids <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_ids
  mat
}

#' Group colours from the shared palette (groups5, extended with series12
#' when more than five groups are present).
group_palette <- function(n) {
  pal <- rice_palette()
  groups <- pal$groups5
  if (n <= length(groups)) return(groups[seq_len(n)])
  extra <- setdiff(pal$series12, groups)
  c(groups, extra)[seq_len(n)]
}

#' Combine plots with patchwork, or save them as separate panel figures.
combine_or_separate <- function(plots, stems, width_mm, height_mm, ncol = 2) {
  if (has_patchwork) {
    combined <- do.call(patchwork::wrap_plots, c(plots, list(ncol = ncol)))
    save_nature_figure(combined, file.path(outdir, stems[1]),
                       width_mm = width_mm, height_mm = height_mm)
  } else {
    for (i in seq_along(plots)) {
      save_nature_figure(plots[[i]], file.path(outdir, stems[i + 1L]),
                         width_mm = width_mm, height_mm = height_mm)
    }
  }
}

#' Parse a GeneRatio string ("12/187") into a numeric value.
parse_gene_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) == 2L) {
      num <- suppressWarnings(as.numeric(p[1]))
      den <- suppressWarnings(as.numeric(p[2]))
      if (!is.na(num) && !is.na(den) && den > 0) num / den else NA_real_
    } else {
      NA_real_
    }
  }, numeric(1))
}

#' Discover contrasts from de_sig_*.csv files (or use --contrasts).
discover_contrasts <- function(contrast_arg) {
  if (!is.null(contrast_arg)) {
    cts <- trimws(strsplit(contrast_arg, ",", fixed = TRUE)[[1]])
    cts <- cts[nzchar(cts)]
    if (length(cts) == 0L) {
      stop("--contrasts was provided but empty.", call. = FALSE)
    }
    return(cts)
  }
  if (!dir.exists(de_dir)) return(character(0))
  files <- list.files(de_dir, pattern = "^de_sig_.*\\.csv$")
  sub("^de_sig_", "", sub("\\.csv$", "", files))
}

contrasts <- discover_contrasts(contrast_arg)
if (length(contrasts) == 0L) {
  message("No contrasts found (no de_sig_*.csv under ", de_dir,
          "); DE-dependent figures will be skipped.")
} else {
  message("Contrasts: ", paste(contrasts, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Figure 1: QC overview (library size, detected genes, PCA, correlation)
# ---------------------------------------------------------------------------

fig1_qc_overview <- function() {
  metrics_file <- file.path(qc_dir, "qc_metrics.csv")
  counts_file <- file.path(qc_dir, "counts_filtered.csv")
  if (!file.exists(metrics_file) || !file.exists(counts_file)) {
    message("01_qc results not found; skipping fig1_qc_overview.")
    return(invisible(FALSE))
  }
  metrics <- utils::read.csv(metrics_file, stringsAsFactors = FALSE,
                             check.names = FALSE)
  groups <- unique(as.character(metrics$group))
  cols <- stats::setNames(group_palette(length(groups)), groups)

  counts <- read_counts_matrix(counts_file)
  lib_sizes <- colSums(counts)
  log_cpm <- log2(sweep(counts, 2, lib_sizes, "/") * 1e6 + 1)
  # PCA on the top 500 highly variable genes (same rule as 01_qc_filter.py).
  gene_var <- apply(log_cpm, 1, stats::var)
  hv <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(500, nrow(log_cpm)))]
  pc <- stats::prcomp(t(log_cpm[hv, , drop = FALSE]), center = TRUE, scale. = FALSE)
  var_exp <- summary(pc)$importance[2, 1:2] * 100
  scores <- data.frame(sample_id = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                       stringsAsFactors = FALSE)
  # Merge first (character keys), then fix the sample order for plotting.
  metrics <- merge(scores, metrics, by = "sample_id")
  metrics$sample_id <- factor(metrics$sample_id, levels = metrics$sample_id)
  has_batch <- "batch" %in% colnames(metrics)

  p_lib <- ggplot2::ggplot(metrics,
                           ggplot2::aes(x = sample_id, y = library_size, fill = group)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL, y = "Library size (total counts)",
                  title = "Library size") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p_det <- ggplot2::ggplot(metrics,
                           ggplot2::aes(x = sample_id, y = detected_genes, fill = group)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL, y = "Detected genes (count > 0)",
                  title = "Detected genes") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  shape_aes <- if (has_batch) ggplot2::aes(shape = factor(batch)) else ggplot2::aes()
  p_pca <- ggplot2::ggplot(metrics, ggplot2::aes(x = PC1, y = PC2, colour = group)) +
    ggplot2::geom_point(shape_aes, size = 2) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
                  y = sprintf("PC2 (%.1f%% variance)", var_exp[2]),
                  title = "Sample PCA") +
    theme_nature()
  cor_mat <- stats::cor(log_cpm, method = "spearman")
  cor_df <- data.frame(x = rep(rownames(cor_mat), times = ncol(cor_mat)),
                       y = rep(colnames(cor_mat), each = nrow(cor_mat)),
                       r = as.numeric(cor_mat), stringsAsFactors = FALSE)
  p_cor <- ggplot2::ggplot(cor_df, ggplot2::aes(x = x, y = y, fill = r)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = diverging_palette(), limits = c(-1, 1),
                                  name = "Spearman r") +
    ggplot2::labs(x = NULL, y = NULL, title = "Sample correlation") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Scientific conclusion: balanced library sizes and detected-gene counts
  # indicate comparable sequencing depth; group separation in the PCA and the
  # high within-group sample correlations (r > 0.8) confirm that the treatment
  # effect is visible before modelling.
  combine_or_separate(
    list(p_lib, p_det, p_pca, p_cor),
    c("fig1_qc_overview", "fig1a_library_size", "fig1b_detected_genes",
      "fig1c_sample_pca", "fig1d_sample_correlation"),
    width_mm = 180, height_mm = 160, ncol = 2
  )
  message("Saved fig1 QC overview.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Figure 2: DE overview (volcano + MA + DE counts)
# ---------------------------------------------------------------------------

fig2_de_overview <- function() {
  if (length(contrasts) == 0L) return(invisible(FALSE))
  de_list <- list()
  for (ct in contrasts) {
    f <- file.path(de_dir, paste0("de_all_", ct, ".csv"))
    if (!file.exists(f)) {
      message("de_all_", ct, ".csv not found; skipping the volcano/MA panels ",
              "for contrast '", ct, "'.")
      next
    }
    de <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    de$contrast <- ct
    de_list[[length(de_list) + 1L]] <- de
  }
  if (length(de_list) == 0L) {
    message("No de_all_<contrast>.csv found under ", de_dir,
            "; skipping fig2_de_overview.")
    return(invisible(FALSE))
  }
  de_all <- do.call(rbind, de_list)
  de_all$padj_plot <- ifelse(is.na(de_all$padj), 1, de_all$padj)
  de_all$cls <- ifelse(!is.na(de_all$padj) & de_all$padj <= 0.05 &
                         de_all$log2FoldChange >= 1, "Up",
                       ifelse(!is.na(de_all$padj) & de_all$padj <= 0.05 &
                                de_all$log2FoldChange <= -1, "Down", "NS"))
  de_all$neg_log10_padj <- -log10(pmax(de_all$padj_plot, 1e-300))
  de_all$log10_baseMean <- log10(pmax(de_all$baseMean, 0) + 1)

  # Top 5 genes per contrast by |log2FC| among the significant genes.
  top5 <- do.call(rbind, lapply(split(de_all, de_all$contrast), function(d) {
    d <- d[d$cls != "NS", , drop = FALSE]
    d <- d[order(-abs(d$log2FoldChange)), , drop = FALSE]
    head(d, 5)
  }))

  p_volcano <- ggplot2::ggplot(
    de_all, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, colour = cls)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.75) +
    ggplot2::scale_colour_manual(
      values = c(Up = nature_semantic("up"), Down = nature_semantic("down"),
                 NS = nature_semantic("not_significant")), name = NULL) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed",
                        colour = nature_semantic("reference_line")) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                        colour = nature_semantic("reference_line")) +
    ggrepel::geom_text_repel(data = top5, ggplot2::aes(label = gene_id),
                             size = 2, max.overlaps = 20, seed = 42,
                             show.legend = FALSE) +
    ggplot2::facet_wrap(~ contrast) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 (adjusted p-value)",
                  title = "Volcano plot") +
    theme_nature()

  p_ma <- ggplot2::ggplot(
    de_all, ggplot2::aes(x = log10_baseMean, y = log2FoldChange, colour = cls)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.75) +
    ggplot2::scale_colour_manual(
      values = c(Up = nature_semantic("up"), Down = nature_semantic("down"),
                 NS = nature_semantic("not_significant")), name = NULL) +
    ggplot2::facet_wrap(~ contrast) +
    ggplot2::labs(x = "log10 (base mean)", y = "log2 fold change",
                  title = "MA plot") +
    theme_nature()

  # DE counts: from de_summary.csv when available, else from the de_sig tables.
  summary_file <- file.path(de_dir, "de_summary.csv")
  if (file.exists(summary_file)) {
    s <- utils::read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    rows <- lapply(contrasts, function(ct) {
      f <- file.path(de_dir, paste0("de_sig_", ct, ".csv"))
      if (!file.exists(f)) {
        return(data.frame(contrast = ct, up = 0, down = 0, stringsAsFactors = FALSE))
      }
      d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      data.frame(contrast = ct, up = sum(d$significant == "Up"),
                 down = sum(d$significant == "Down"), stringsAsFactors = FALSE)
    })
    s <- do.call(rbind, rows)
  }
  counts_long <- data.frame(
    contrast = rep(as.character(s$contrast), 2L),
    direction = rep(c("Up", "Down"), each = nrow(s)),
    count = c(s$up, s$down), stringsAsFactors = FALSE
  )
  p_counts <- ggplot2::ggplot(counts_long,
                              ggplot2::aes(x = contrast, y = count, fill = direction)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = c(Up = nature_semantic("up"),
                                          Down = nature_semantic("down"))) +
    ggplot2::labs(x = NULL, y = "Number of DE genes",
                  title = "DE gene counts") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Scientific conclusion: genes in the upper-right (red) and upper-left (blue)
  # volcano wings are significantly up- and down-regulated; the MA plot
  # confirms that significance is not driven by expression level, and the
  # count panel quantifies the transcriptional response per contrast.
  n_ct <- length(unique(de_all$contrast))
  combine_or_separate(
    list(p_volcano, p_ma, p_counts),
    c("fig2_de_overview", "fig2a_volcano", "fig2b_ma", "fig2c_de_counts"),
    width_mm = 180, height_mm = max(89, 45 * n_ct + 50), ncol = 3
  )
  message("Saved fig2 DE overview.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Figure 3: enrichment overview (GO BP + KEGG dot plots)
# ---------------------------------------------------------------------------

fig3_enrichment_overview <- function() {
  if (!dir.exists(enr_dir)) {
    message("03_enrichment results not found; skipping fig3_enrichment_overview.")
    return(invisible(FALSE))
  }
  specs <- list(
    c("GO_BP", "UP", "go_BP_UP.csv"),
    c("GO_BP", "DOWN", "go_BP_DOWN.csv"),
    c("KEGG", "UP", "kegg_UP.csv"),
    c("KEGG", "DOWN", "kegg_DOWN.csv")
  )
  enr_list <- list()
  for (sp in specs) {
    f <- file.path(enr_dir, sp[3])
    if (!file.exists(f)) next
    d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(d) == 0L) next
    d$category <- sp[1]
    d$direction <- sp[2]
    d$ratio <- parse_gene_ratio(d$GeneRatio)
    enr_list[[length(enr_list) + 1L]] <- d
  }
  if (length(enr_list) == 0L) {
    message("No non-empty enrichment tables under ", enr_dir,
            "; skipping fig3_enrichment_overview.")
    return(invisible(FALSE))
  }
  enr <- do.call(rbind, enr_list)
  enr <- enr[!is.na(enr$ratio), , drop = FALSE]
  if (nrow(enr) == 0L) {
    message("Enrichment tables contain no usable GeneRatio values; ",
            "skipping fig3_enrichment_overview.")
    return(invisible(FALSE))
  }
  # Order terms by gene ratio within each category/direction facet.
  enr <- enr[order(enr$category, enr$direction, enr$ratio), , drop = FALSE]
  enr$Description <- factor(enr$Description, levels = unique(enr$Description))

  p <- ggplot2::ggplot(enr,
                       ggplot2::aes(x = ratio, y = Description, size = Count,
                                    colour = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(category ~ direction, scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_gradientn(colours = sequential_palette(), name = "padj") +
    ggplot2::labs(x = "Gene ratio", title = "Enrichment overview") +
    theme_nature()

  # Scientific conclusion: terms in the upper part of each panel (large gene
  # ratio, small padj) are the biological processes/pathways most strongly
  # affected by the treatment in each direction.
  save_nature_figure(p, file.path(outdir, "fig3_enrichment_overview"),
                     width_mm = 180, height_mm = max(89, 14 * 10 + 40))
  message("Saved fig3 enrichment overview.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Figure 4: WGCNA overview (module-trait, module sizes, hub genes)
# ---------------------------------------------------------------------------

fig4_wgcna_overview <- function() {
  mt_file <- file.path(wgcna_dir, "module_trait_correlation.csv")
  gm_file <- file.path(wgcna_dir, "gene_modules.csv")
  hub_file <- file.path(wgcna_dir, "hub_genes.csv")
  vst_file <- file.path(de_dir, "normalized_counts_vst.csv")
  if (!file.exists(mt_file) && !file.exists(gm_file)) {
    message("05_wgcna results not found; skipping fig4_wgcna_overview.")
    return(invisible(FALSE))
  }
  panels <- list()
  stems <- c("fig4_wgcna_overview", "fig4a_module_trait", "fig4b_module_size",
             "fig4c_hub_expression")

  if (file.exists(mt_file)) {
    mt <- utils::read.csv(mt_file, stringsAsFactors = FALSE, check.names = FALSE)
    mt_wide <- with(mt, base::tapply(correlation,
                                     list(module = module, trait = trait),
                                     function(x) x))
    mt_df <- data.frame(
      module = rep(rownames(mt_wide), times = ncol(mt_wide)),
      trait = rep(colnames(mt_wide), each = nrow(mt_wide)),
      correlation = as.numeric(mt_wide), stringsAsFactors = FALSE
    )
    panels[[length(panels) + 1L]] <- ggplot2::ggplot(
      mt_df, ggplot2::aes(x = trait, y = module, fill = correlation)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", correlation)),
                         size = 2) +
      ggplot2::scale_fill_gradientn(colours = diverging_palette(),
                                    limits = c(-1, 1), name = "cor") +
      ggplot2::labs(x = NULL, y = NULL, title = "Module-trait relationships") +
      theme_nature()
  } else {
    message("module_trait_correlation.csv not found; the module-trait panel ",
            "is skipped.")
  }

  if (file.exists(gm_file)) {
    gm <- utils::read.csv(gm_file, stringsAsFactors = FALSE, check.names = FALSE)
    sizes <- as.data.frame(table(gm$module), stringsAsFactors = FALSE)
    panels[[length(panels) + 1L]] <- ggplot2::ggplot(
      sizes, ggplot2::aes(x = stats::reorder(Var1, Freq), y = Freq)) +
      ggplot2::geom_col(fill = rice_palette()$series12[1]) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Module", y = "Number of genes",
                    title = "Module sizes") +
      theme_nature()
  }

  if (file.exists(hub_file) && file.exists(vst_file)) {
    hubs <- utils::read.csv(hub_file, stringsAsFactors = FALSE, check.names = FALSE)
    hubs <- hubs[order(-abs(hubs$kME)), , drop = FALSE]
    hubs <- hubs[!duplicated(hubs$gene_id), , drop = FALSE]
    top20 <- head(hubs, 20)
    vs <- utils::read.csv(vst_file, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(vs) <- vs$gene_id
    vs$gene_id <- NULL
    genes <- intersect(top20$gene_id, rownames(vs))
    if (length(genes) > 0L) {
      mat <- as.matrix(vs[genes, , drop = FALSE])
      z <- t(scale(t(mat)))
      z_df <- data.frame(gene = rep(genes, times = ncol(z)),
                         sample = rep(colnames(z), each = nrow(z)),
                         z = as.numeric(z), stringsAsFactors = FALSE)
      z_df$gene <- factor(z_df$gene, levels = genes)
      m <- max(abs(z_df$z), na.rm = TRUE)
      panels[[length(panels) + 1L]] <- ggplot2::ggplot(
        z_df, ggplot2::aes(x = sample, y = gene, fill = z)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(colours = diverging_palette(),
                                      limits = c(-m, m), name = "Z-score") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Hub gene expression (VST, Z-score)") +
        theme_nature() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    }
  } else {
    message("hub_genes.csv or normalized_counts_vst.csv not found; the hub ",
            "gene panel is skipped.")
  }

  if (length(panels) == 0L) {
    message("No WGCNA panel could be drawn; skipping fig4_wgcna_overview.")
    return(invisible(FALSE))
  }

  # Scientific conclusion: modules correlating with a trait (|r| high) are
  # candidate drivers of the phenotype; the hub genes of the largest modules
  # are the most central regulators of each co-expression programme.
  combine_or_separate(panels, stems, width_mm = 180,
                      height_mm = max(89, 10 * length(panels) * 6 + 40), ncol = 2)
  message("Saved fig4 WGCNA overview.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Figure 5: TF overview (top-20 TF families)
# ---------------------------------------------------------------------------

fig5_tf_overview <- function() {
  stats_file <- file.path(tf_dir, "tf_family_stats.csv")
  if (!file.exists(stats_file)) {
    message("06_tf results not found; skipping fig5_tf_overview.")
    return(invisible(FALSE))
  }
  s <- utils::read.csv(stats_file, stringsAsFactors = FALSE, check.names = FALSE)
  s <- s[order(-s$n_de, -s$n_annotated), , drop = FALSE]
  top20 <- head(s, 20)
  if (nrow(top20) == 0L) {
    message("tf_family_stats.csv is empty; skipping fig5_tf_overview.")
    return(invisible(FALSE))
  }
  long <- data.frame(family = rep(top20$family, 2L),
                     direction = rep(c("Up", "Down"), each = nrow(top20)),
                     count = c(top20$n_up, top20$n_down),
                     stringsAsFactors = FALSE)
  long$family <- factor(long$family, levels = rev(unique(top20$family)))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = family, y = count, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(Up = nature_semantic("up"),
                                          Down = nature_semantic("down"))) +
    ggplot2::labs(x = "TF family", y = "Number of DE transcription factors",
                  title = "Top 20 differentially expressed TF families") +
    theme_nature()

  # Scientific conclusion: TF families with many DE members mark the
  # transcriptional programmes activated or repressed by the treatment.
  save_nature_figure(p, file.path(outdir, "fig5_tf_overview"),
                     width_mm = 120, height_mm = max(89, nrow(top20) * 5 + 30))
  message("Saved fig5 TF overview.")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Run all figures
# ---------------------------------------------------------------------------

fig1_qc_overview()
fig2_de_overview()
fig3_enrichment_overview()
fig4_wgcna_overview()
fig5_tf_overview()

write_run_info_r(
  outdir, "07_plot_transcriptome.R",
  params = list(
    results_dir = results_dir,
    contrasts = paste(contrasts, collapse = ","),
    patchwork = has_patchwork
  ),
  packages = c("ggplot2", "ggrepel")
)
message("Figure assembly complete. Results in ", outdir)
