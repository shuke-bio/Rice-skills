###############################################################################
# 03_plot_proteome.R -- publication figures for differential proteomics
#
# Title:    Rice-skills proteomics figures
# Description:
#   Volcano plot and top-protein expression heatmap from the differential
#   proteomics results, in the shared Nature figure specification (PDF/SVG/
#   TIFF, unified palette, Arial/Helvetica).
#
# Usage:
#   Rscript scripts/proteomics/03_plot_proteome.R \
#       --diff results/proteomics/02_diff/protein_diff_all.csv \
#       --matrix results/proteomics/01_preprocess/proteins_filtered.csv \
#       --metadata demo_data/proteomics/metadata.csv \
#       --outdir results/proteomics/03_plots
#
# Inputs:
#   --diff      protein_diff_all.csv (from 02_diff_protein.R)
#   --matrix    proteins_filtered.csv (log2 normalised; for the heatmap)
#   --metadata  sample metadata CSV (for group annotation)
#   --group-col metadata group column (default "group")
#   --top       number of top proteins in the heatmap (default 50)
#   --lf2fc / --fdr  significance cutoffs (must match 02; default 1 / 0.05)
#   --label-n   number of proteins labelled on the volcano (default 10)
#   --outdir    output directory (default results/proteomics/03_plots)
# Outputs:
#   volcano.*        volcano plot (PDF/SVG/TIFF)
#   top_heatmap.*    z-scored heatmap of the top differential proteins
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
diff_path <- opt(opts, "diff")
matrix_path <- opt(opts, "matrix")
metadata_path <- opt(opts, "metadata")
group_col <- opt(opts, "group-col", "group")
top_n <- suppressWarnings(as.integer(opt(opts, "top", "50")))
lf2fc <- suppressWarnings(as.numeric(opt(opts, "lf2fc", "1")))
fdr_cut <- suppressWarnings(as.numeric(opt(opts, "fdr", "0.05")))
label_n <- suppressWarnings(as.integer(opt(opts, "label-n", "10")))
outdir <- opt(opts, "outdir", "results/proteomics/03_plots")

if (is.null(diff_path)) stop("Required: --diff <protein_diff_all.csv>.", call. = FALSE)
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills proteomics figures ==")

pal <- rice_palette()
col_up <- nature_semantic("up")
col_dn <- nature_semantic("down")
col_ns <- nature_semantic("not_significant")

# ---------------------------------------------------------------------------
# 1) Volcano
# ---------------------------------------------------------------------------
d <- read_csv_safe(diff_path, "differential protein table")
check_columns(d, c("protein_id", "log2FC", "p", "fdr"), "differential protein table")
d$sig <- ifelse(!is.na(d$fdr) & d$fdr <= fdr_cut & d$log2FC >= lf2fc, "up",
         ifelse(!is.na(d$fdr) & d$fdr <= fdr_cut & d$log2FC <= -lf2fc, "down", "ns"))
d$neglog10p <- -log10(pmax(d$p, 1e-300))

n_up <- sum(d$sig == "up"); n_dn <- sum(d$sig == "down")
p_volcano <- ggplot2::ggplot(d, ggplot2::aes(x = log2FC, y = neglog10p)) +
  ggplot2::geom_point(ggplot2::aes(colour = sig), size = 0.8, alpha = 0.75) +
  ggplot2::scale_colour_manual(
    values = c(up = col_up, down = col_dn, ns = col_ns),
    breaks = c("up", "down", "ns"),
    labels = c(sprintf("Up (%d)", n_up), sprintf("Down (%d)", n_dn), "n.s.")) +
  ggplot2::geom_vline(xintercept = c(-lf2fc, lf2fc), linetype = "dashed",
                      colour = "grey55", linewidth = 0.3) +
  ggplot2::labs(x = expression(log[2]~"FC ("*alt~"vs"~ref*")"),
                y = expression(-log[10]~(italic(p))),
                colour = NULL, title = "Differential proteins") +
  theme_nature()

if (label_n > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
  lab <- head(d[order(d$p), , drop = FALSE], label_n)
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = protein_id),
                             size = 1.9, min.segment.length = 0, max.overlaps = 20,
                             segment.size = 0.2, colour = "grey25")
}
save_nature_figure(p_volcano, file.path(outdir, "volcano"), width_mm = 89, height_mm = 80)

# ---------------------------------------------------------------------------
# 2) Top-protein heatmap
# ---------------------------------------------------------------------------
if (!is.null(matrix_path)) {
  mat_df <- read_csv_safe(matrix_path, "protein matrix")
  protein_ids <- as.character(mat_df[[1]])
  X <- as.matrix(mat_df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(X) <- "double")
  rownames(X) <- protein_ids

  top_ids <- head(d$protein_id[order(d$p)], top_n)
  top_ids <- top_ids[top_ids %in% rownames(X)]
  if (length(top_ids) >= 2) {
    Z <- X[top_ids, , drop = FALSE]
    # z-score per protein
    Z <- t(scale(t(Z)))
    Z[is.na(Z)] <- 0

    sample_ids <- colnames(Z)
    group_vec <- NULL
    if (!is.null(metadata_path)) {
      meta <- read_pheno(metadata_path)
      if (group_col %in% colnames(meta)) {
        group_vec <- meta[[group_col]][match(sample_ids, meta$sample_id)]
      }
    }
    # order samples by group then cluster within
    if (!is.null(group_vec) && !all(is.na(group_vec))) {
      ord <- order(group_vec, sample_ids)
    } else {
      hc <- stats::hclust(stats::dist(t(Z)))
      ord <- hc$order
    }
    Z <- Z[, ord, drop = FALSE]
    sample_ids <- colnames(Z)

    df_hm <- data.frame(
      protein = rep(rownames(Z), times = ncol(Z)),
      sample = rep(sample_ids, each = nrow(Z)),
      z = as.vector(Z)
    )
    df_hm$protein <- factor(df_hm$protein, levels = rev(rownames(Z)))
    df_hm$sample <- factor(df_hm$sample, levels = sample_ids)
    df_hm$group <- if (!is.null(group_vec)) rep(group_vec[ord], each = nrow(Z)) else NA

    p_hm <- ggplot2::ggplot(df_hm, ggplot2::aes(x = sample, y = protein, fill = z)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = diverging_palette()[1], mid = "white",
                                    high = diverging_palette()[7], midpoint = 0,
                                    name = "z-score", limits = c(-2.5, 2.5),
                                    oob = scales::squish) +
      ggplot2::labs(x = NULL, y = NULL, title = sprintf("Top %d differential proteins", nrow(Z))) +
      theme_nature() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 4),
                     axis.text.y = ggplot2::element_text(size = 3.6),
                     legend.position = "right")
    if (!all(is.na(df_hm$group))) {
      p_hm <- p_hm + ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free_x")
    }
    save_nature_figure(p_hm, file.path(outdir, "top_heatmap"),
                       width_mm = 140, height_mm = 130)
  }
}

write_run_info_r(outdir, "03_plot_proteome.R", params = list(
  diff = diff_path, matrix = matrix_path, top = top_n, lf2fc = lf2fc, fdr = fdr_cut
))
message("Done. Outputs written to: ", normalizePath(outdir))
