###############################################################################
# 02_oplsda.R -- OPLS-DA differential metabolite analysis (Rice-skills step 02)
#
# Title:    OPLS-DA modelling, permutation testing, VIP and univariate screening
# Description:
#   Fits a two-group OPLS-DA model (ropls package) on a preprocessed metabolome
#   matrix, validates it by permutation testing, ranks metabolites by VIP,
#   and combines VIP with univariate Wilcoxon/t tests (BH-adjusted) to call
#   differential metabolites. Produces publication-ready figures (score plot,
#   S-plot, VIP barplot, permutation plot) via the shared plot_theme.R layer.
#
# Usage:
#   Rscript scripts/metabolomics/02_oplsda.R \
#       --data results/metabolome/01_preprocess/metabolome_preprocessed.csv \
#       --metadata examples/metabolome_data/metadata.csv \
#       --outdir results/metabolome/02_oplsda \
#       --vip 1 --pvalue 0.05 --permutations 200
#
# Inputs:
#   --data        Preprocessed abundance matrix CSV (metabolite_id [, annotations],
#                 samples...), e.g. the output of 01_preprocess.py.
#   --metadata    Sample metadata CSV (sample_id, group).
#
# Outputs (written to --outdir):
#   diff_metabolites.csv  significant metabolites (metabolite_id, VIP, p, padj,
#                         log2FC, significant, plus kegg_id/class/... if present
#                         in --data; consumed by 03_pathway.R / 04_...R)
#   oplsda_scores.csv     sample scores on the predictive/orthogonal components
#   pca_scores.csv        unsupervised PCA scores (reference for the OPLS-DA)
#   vip_scores.csv        VIP of every metabolite
#   oplsda_loadings.csv   predictive loading p(1) and correlation loading p(corr)
#   oplsda_score_plot.*   t1 vs to1 score plot with 95% confidence ellipses
#   s_plot.*              p(1) vs p(corr) S-plot, top-VIP metabolites labelled
#   vip_barplot.*         top-30 VIP barplot coloured by the VIP threshold
#   permutation_test.*    permuted vs observed R2Y/Q2 with p-values
#   run_info.txt          parameters, package versions and the seed
#
# Notes:
#   - ropls (Bioconductor) depends on preprocessCore, which must be installed
#     too: BiocManager::install(c("ropls", "preprocessCore")).
#   - ropls object layout (verified against ropls >= 1.24): slots @vipVn,
#     @scoreMN, @orthoScoreMN, @loadingMN, @modelDF, @summaryDF and
#     @suppLs[["permMN"]] (permutation results); accessors getVipVn(),
#     getScoreMN(orthoL=), getLoadingMN(orthoL=), getSummaryDF(). Field names
#     may differ between ropls versions: run str(opls_obj) to inspect, and the
#     defensive extraction helpers below will report missing fields.
#   - OPLS-DA is binary: with more than two groups, pass --contrast "A,B"
#     (comma-separated) or the first two groups are used with a message.
#   - y is passed as a factor whose levels are (reference, comparison); the
#     sign of the scores follows this level order and log2FC is computed as
#     mean(comparison) - mean(reference), which equals the log2 fold-change
#     when the input matrix was log-transformed (as by 01_preprocess.py).
#   - scaleC="none" is the default because 01_preprocess.py already mean-centres
#     and pareto-scales the data; use --scale center/pareto/standard for raw
#     input (ropls always mean-centres internally).
#   - A metabolite is called differential when VIP >= --vip AND raw p < --pvalue
#     (spec convention); padj is reported as well and is the stricter choice.
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

# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------

need_pkgs("ropls")   # ropls also requires preprocessCore (Bioconductor)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

#' Resolve the two groups to compare.
#'
#' @param groups_vc Character vector of group names (order of appearance).
#' @param contrast Optional "A,B" specification (comma or "_vs_" separated).
#' @return Character vector of length 2 (reference, comparison).
resolve_contrast <- function(groups_vc, contrast = NULL) {
  if (!is.null(contrast) && !isTRUE(contrast)) {
    parts <- trimws(strsplit(contrast, ",|_vs_| vs ", perl = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    if (length(parts) != 2L) {
      stop("--contrast must name exactly two groups, e.g. '--contrast control,salt'. Got: ",
           contrast, call. = FALSE)
    }
    missing <- setdiff(parts, groups_vc)
    if (length(missing)) {
      stop("--contrast group(s) not found in metadata: ", paste(missing, collapse = ", "),
           "\nAvailable groups: ", paste(groups_vc, collapse = ", "), call. = FALSE)
    }
    return(parts)
  }
  if (length(groups_vc) > 2L) {
    message("More than two groups in metadata (", paste(groups_vc, collapse = ", "),
            "); using the first two: ", paste(groups_vc[1:2], collapse = " vs "),
            ". Use --contrast to select another pair.")
  }
  if (length(groups_vc) < 2L) {
    stop("OPLS-DA needs at least two groups; found: ", paste(groups_vc, collapse = ", "),
         call. = FALSE)
  }
  groups_vc[1:2]
}

#' Defensive slot/accessor extraction for ropls objects.
ropls_get <- function(opls_obj, what = c("vip", "score", "ortho_score", "loading",
                                         "ortho_loading", "summary")) {
  what <- match.arg(what)
  tryCatch(
    switch(what,
      vip          = ropls::getVipVn(opls_obj),
      score        = ropls::getScoreMN(opls_obj, orthoL = FALSE),
      ortho_score  = ropls::getScoreMN(opls_obj, orthoL = TRUE),
      loading      = ropls::getLoadingMN(opls_obj, orthoL = FALSE),
      ortho_loading = ropls::getLoadingMN(opls_obj, orthoL = TRUE),
      summary      = ropls::getSummaryDF(opls_obj)),
    error = function(e) {
      # Fall back to direct slot access (field names differ across versions).
      slot_name <- switch(what,
        vip = "vipVn", score = "scoreMN", ortho_score = "orthoScoreMN",
        loading = "loadingMN", ortho_loading = "orthoLoadingMN", summary = "summaryDF")
      if (slot_name %in% slotNames(opls_obj)) slot(opls_obj, slot_name) else {
        stop("Could not extract '", what, "' from the ropls object (slot '", slot_name,
             "' missing). Run str(opls_obj) to inspect this ropls version.",
             call. = FALSE)
      }
    })
}

#' Extract permutation-test results from a ropls object.
#'
#' @return List with permuted R2Y/Q2 vectors, observed values and p-values.
permutation_results <- function(opls_obj, n_perm) {
  perm_mn <- NULL
  if ("suppLs" %in% slotNames(opls_obj)) {
    supp <- slot(opls_obj, "suppLs")
    if (!is.null(supp[["permMN"]])) perm_mn <- supp[["permMN"]]
  }
  sum_df <- ropls_get(opls_obj, "summary")
  stat <- function(col) if (col %in% colnames(sum_df)) as.numeric(sum_df[1, col]) else NA_real_
  out <- list(
    perm_r2y = NULL, perm_q2 = NULL,
    obs_r2y = stat("R2Y(cum)"), obs_q2 = stat("Q2(cum)"),
    p_r2y = stat("pR2Y"), p_q2 = stat("pQ2"),
    r2x = stat("R2X(cum)"), n_pre = stat("pre"), n_ort = stat("ort")
  )
  if (!is.null(perm_mn)) {
    out$perm_r2y <- as.numeric(perm_mn[-1, "R2Y(cum)"])
    out$perm_q2 <- as.numeric(perm_mn[-1, "Q2(cum)"])
    # Manual p-value fallback (same estimator as ropls: (1 + #{perm >= obs}) / (n + 1)).
    if (is.na(out$p_r2y) && length(out$perm_r2y)) {
      out$p_r2y <- (1 + sum(out$perm_r2y >= out$obs_r2y)) / (length(out$perm_r2y) + 1)
    }
    if (is.na(out$p_q2) && length(out$perm_q2)) {
      out$p_q2 <- (1 + sum(out$perm_q2 >= out$obs_q2)) / (length(out$perm_q2) + 1)
    }
  }
  out
}

#' Univariate two-group tests with BH adjustment.
#'
#' @param X Samples x metabolites matrix.
#' @param y Factor with two levels (reference first).
#' @param test "auto" (Wilcoxon when the smallest group has >= 4 samples, else
#'        Welch t-test), "wilcoxon" or "t".
#' @return Data frame with metabolite_id, test, p, padj, log2FC (mean difference).
univariate_tests <- function(X, y, test = "auto") {
  test <- match.arg(test, c("auto", "wilcoxon", "t"))
  lv <- levels(y)
  in_ref <- y == lv[1]
  in_cmp <- y == lv[2]
  n_ref <- sum(in_ref); n_cmp <- sum(in_cmp)
  if (test == "auto") {
    test <- if (min(n_ref, n_cmp) >= 4L) "wilcoxon" else "t"
    message("Univariate test: ", test,
            " (smallest group has ", min(n_ref, n_cmp), " samples; Wilcoxon cannot reach ",
            "p < 0.05 with fewer than 4 replicates per group).")
  }
  p_vn <- vapply(seq_len(ncol(X)), function(j) {
    x_ref <- X[in_ref, j]
    x_cmp <- X[in_cmp, j]
    tryCatch(
      if (test == "wilcoxon") stats::wilcox.test(x_ref, x_cmp)$p.value
      else stats::t.test(x_ref, x_cmp)$p.value,
      error = function(e) NA_real_)
  }, numeric(1))
  names(p_vn) <- colnames(X)
  # Mean difference; equals log2 fold-change for log-transformed input.
  log2fc <- colMeans(X[in_cmp, , drop = FALSE]) - colMeans(X[in_ref, , drop = FALSE])
  data.frame(
    metabolite_id = colnames(X),
    test = test,
    p = as.numeric(p_vn),
    padj = stats::p.adjust(as.numeric(p_vn), method = "BH"),
    log2FC = as.numeric(log2FC),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Figures (ggplot2 + shared Nature theme)
# ---------------------------------------------------------------------------

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

plot_pca_reference <- function(pca_df, evr_vn, outdir) {
  df <- pca_df
  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, colour = group)) +
    ggplot2::geom_point(size = 2)
  p <- maybe_add_ellipses(p)
  p + scale_colour_groups() +
    ggplot2::labs(
      title = "PCA (unsupervised reference)",
      subtitle = paste0("PC1 ", round(evr_vn[1] * 100, 1), "% / PC2 ",
                        round(evr_vn[2] * 100, 1), "% variance; no supervision applied"),
      x = paste0("PC1 (", round(evr_vn[1] * 100, 1), "% variance)"),
      y = paste0("PC2 (", round(evr_vn[2] * 100, 1), "% variance)")) +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "pca_scores"), width_mm = 89, height_mm = 75)
}

plot_oplsda_scores <- function(scores_df, stats_ls, outdir) {
  # Scientific conclusion: separation along the predictive component t1 with
  # orthogonal variation (to1) capturing within-group spread.
  sub_txt <- paste0(
    "R2Y = ", round(stats_ls$r2y, 3), ", Q2 = ", round(stats_ls$q2, 3),
    ", R2X = ", round(stats_ls$r2x, 3), "; groups separate along t1 (predictive), ",
    "to1 = orthogonal (within-group) variation")
  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = t1, y = to1, colour = group)) +
    ggplot2::geom_point(size = 2)
  p <- maybe_add_ellipses(p)
  p + ggplot2::geom_vline(xintercept = 0, colour = nature_semantic("reference_line"),
                          linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, colour = nature_semantic("reference_line"),
                        linewidth = 0.4) +
    scale_colour_groups() +
    ggplot2::labs(title = "OPLS-DA score plot", subtitle = sub_txt,
                  x = "t1 (predictive score)", y = "to1 (orthogonal score)") +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "oplsda_score_plot"), width_mm = 89, height_mm = 80)
}

plot_s_plot <- function(load_df, vip_cut, outdir) {
  # Scientific conclusion: metabolites in the outer arms of the S (high |p(1)|
  # and high |p(corr)|) combine strong covariance with reliable correlation and
  # are the biomarker candidates.
  df <- load_df
  df$cls <- ifelse(!df$significant, "ns",
                   ifelse(df$log2FC >= 0, "up", "down"))
  df$cls <- factor(df$cls, levels = c("up", "down", "ns"))
  col_map <- c(up = nature_semantic("up"), down = nature_semantic("down"),
               ns = nature_semantic("not_significant"))
  top_idx <- head(order(-df$VIP), 8)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = p1, y = p_corr, colour = cls)) +
    ggplot2::geom_point(size = 1.6, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, colour = nature_semantic("reference_line"),
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, colour = nature_semantic("reference_line"),
                        linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = col_map, aesthetics = "colour",
                                 name = "Class",
                                 labels = c(up = paste0("up (log2FC > 0, VIP >= ", vip_cut, ")"),
                                            down = paste0("down (log2FC < 0, VIP >= ", vip_cut, ")"),
                                            ns = paste0("ns (VIP < ", vip_cut, " or p >= cut)"))) +
    ggrepel::geom_label_repel(
      data = df[top_idx, ],
      ggplot2::aes(label = metabolite_id),
      size = 2.2, colour = "grey20", max.overlaps = 20, min.segment.length = 0) +
    ggplot2::labs(
      title = "OPLS-DA S-plot",
      subtitle = "p(1) = predictive loading (covariance with t1); p(corr) = correlation with t1; top VIP labelled",
      x = "p(1) (loading)", y = "p(corr) (correlation loading)") +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "s_plot"), width_mm = 89, height_mm = 80)
}

plot_vip_barplot <- function(vip_df, vip_cut, outdir) {
  # Scientific conclusion: VIP >= 1 metabolites contribute above-average
  # importance to the group separation (red); grey bars are below the cut-off.
  df <- utils::head(vip_df[order(-vip_df$VIP), ], 30)
  df$metabolite_id <- factor(df$metabolite_id, levels = rev(df$metabolite_id))
  df$above <- df$VIP >= vip_cut
  p <- ggplot2::ggplot(df, ggplot2::aes(x = metabolite_id, y = VIP, fill = above)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = vip_cut, colour = nature_semantic("reference_line"),
                        linewidth = 0.5, linetype = "dashed") +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = nature_semantic("up"), `FALSE` = nature_semantic("not_significant")),
      name = NULL, labels = c(`TRUE` = paste0("VIP >= ", vip_cut),
                              `FALSE` = paste0("VIP < ", vip_cut))) +
    ggplot2::labs(title = "Top 30 VIP metabolites",
                  subtitle = paste0("Variable Importance in Projection; dashed line = VIP threshold (",
                                    vip_cut, ")"),
                  x = NULL, y = "VIP") +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "vip_barplot"), width_mm = 89, height_mm = 110)
}

plot_permutation <- function(perm_ls, n_perm, outdir) {
  # Scientific conclusion: permuted R2Y/Q2 far below the observed values (small
  # p) indicate a non-random, well-predictive model.
  if (is.null(perm_ls$perm_r2y) || !length(perm_ls$perm_r2y)) {
    message("No permutation results available (--permutations 0); skipping permutation_test figure.")
    return(invisible(NULL))
  }
  df <- data.frame(
    permutation = rep(seq_along(perm_ls$perm_r2y), 2),
    metric = rep(c("R2Y", "Q2"), each = length(perm_ls$perm_r2y)),
    value = c(perm_ls$perm_r2y, perm_ls$perm_q2)
  )
  fmt_p <- function(x) if (is.na(x)) "NA" else format(x, digits = 2)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = permutation, y = value, colour = metric)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = perm_ls$obs_r2y, colour = "#5B84AE",
                        linewidth = 0.6, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = perm_ls$obs_q2, colour = "#CE6A6B",
                        linewidth = 0.6, linetype = "dashed") +
    scale_colour_nature() +
    ggplot2::labs(
      title = "OPLS-DA permutation test",
      subtitle = paste0(n_perm, " permutations; dashed lines = observed R2Y/Q2; ",
                        "p(R2Y) = ", fmt_p(perm_ls$p_r2y), ", p(Q2) = ", fmt_p(perm_ls$p_q2)),
      x = "Permutation", y = "R2Y / Q2 (cum)") +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "permutation_test"), width_mm = 89, height_mm = 70)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  data_path <- opt(opts, "data")
  if (is.null(data_path)) stop("--data is required (preprocessed abundance matrix).", call. = FALSE)
  meta_path <- opt(opts, "metadata")
  if (is.null(meta_path)) stop("--metadata is required (sample metadata CSV).", call. = FALSE)
  outdir <- opt(opts, "outdir", "results/metabolome/02_oplsda")
  group_col <- opt(opts, "group-col", "group")
  vip_cut <- as.numeric(opt(opts, "vip", 1))
  pvalue <- as.numeric(opt(opts, "pvalue", 0.05))
  n_perm <- as.integer(opt(opts, "permutations", 200))
  contrast <- opt(opts, "contrast", NULL)
  test <- opt(opts, "test", "auto")
  scale <- opt(opts, "scale", "none")

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", normalizePath(outdir))

  # ---- Read and align ------------------------------------------------------
  mat_df <- read_csv_safe(data_path, "preprocessed metabolome matrix")
  meta_df <- read_csv_safe(meta_path, "sample metadata")
  check_columns(meta_df, c("sample_id", group_col), "metadata")

  id_col <- colnames(mat_df)[1]
  shared <- intersect(meta_df$sample_id, colnames(mat_df))
  if (!length(shared)) {
    stop("No overlapping samples between --data and --metadata.\n",
         "matrix samples: ", paste(head(colnames(mat_df)[-1], 5), collapse = ", "),
         "\nmetadata samples: ", paste(head(meta_df$sample_id, 5), collapse = ", "),
         "\nHint: sample IDs must match exactly.", call. = FALSE)
  }
  if (length(shared) < length(unique(meta_df$sample_id)) || length(shared) < ncol(mat_df) - 1) {
    message("Keeping ", length(shared), " shared samples (dropped non-overlapping columns/rows).")
  }
  non_numeric <- shared[!vapply(shared, function(s) is.numeric(mat_df[[s]]), logical(1))]
  if (length(non_numeric)) {
    stop("Non-numeric sample column(s) in --data: ", paste(head(non_numeric, 5), collapse = ", "),
         "\nHint: the matrix should contain only the ID/annotation columns plus numeric abundances.",
         call. = FALSE)
  }
  mat_df <- mat_df[, c(id_col, shared), drop = FALSE]
  meta_df <- meta_df[match(shared, meta_df$sample_id), , drop = FALSE]

  met_ids <- as.character(mat_df[[id_col]])
  X <- t(as.matrix(mat_df[, shared, drop = FALSE]))
  colnames(X) <- met_ids
  rownames(X) <- shared
  storage.mode(X) <- "double"
  message("Model matrix: ", nrow(X), " samples x ", ncol(X), " metabolites.")

  # ---- Groups and design ---------------------------------------------------
  groups <- resolve_contrast(unique(as.character(meta_df[[group_col]])), contrast)
  keep <- as.character(meta_df[[group_col]]) %in% groups
  X <- X[keep, , drop = FALSE]
  y <- factor(as.character(meta_df[[group_col]])[keep], levels = groups)
  message("OPLS-DA contrast: ", groups[1], " (reference) vs ", groups[2], " (comparison); n = ",
          sum(y == groups[1]), " vs ", sum(y == groups[2]), ".")

  # ---- Unsupervised PCA reference ------------------------------------------
  pca <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  evr <- pca$sdev^2 / sum(pca$sdev^2)
  pca_df <- data.frame(sample_id = rownames(X), group = y,
                       PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
  utils::write.csv(pca_df, file.path(outdir, "pca_scores.csv"), row.names = FALSE, quote = FALSE)
  plot_pca_reference(pca_df, evr, outdir)

  # ---- OPLS-DA ---------------------------------------------------------------
  # crossvalI must be smaller than the number of samples (ropls requirement).
  crossval <- max(2L, min(7L, nrow(X) - 1L))
  opls_obj <- ropls::opls(as.matrix(X), y,
                          predI = 1L,      # one predictive component
                          orthoI = 1L,     # one orthogonal component (classic t1/to1 plot)
                          permI = n_perm,  # permutation test (0 disables)
                          scaleC = scale,  # data are already scaled by 01_preprocess.py
                          crossvalI = crossval,
                          fig.pdfC = "none", info.txtC = "none")
  message("Fitted model type: ", ropls::getScoreMN(opls_obj) |> dim() |> paste(collapse = "x"),
          " score matrix; inspect with str(opls_obj) if slot names look unfamiliar.")

  score_mn <- ropls_get(opls_obj, "score")
  ortho_mn <- ropls_get(opls_obj, "ortho_score")
  t1 <- as.numeric(score_mn[, 1])
  to1 <- if (ncol(ortho_mn) >= 1L) as.numeric(ortho_mn[, 1]) else rep(NA_real_, nrow(X))

  # ---- Model statistics and permutation --------------------------------------
  perm_ls <- permutation_results(opls_obj, n_perm)
  stats_ls <- list(r2y = perm_ls$obs_r2y, q2 = perm_ls$obs_q2, r2x = perm_ls$r2x)
  message("OPLS-DA: R2Y = ", round(perm_ls$obs_r2y, 3), ", Q2 = ", round(perm_ls$obs_q2, 3),
          " (permutation p(R2Y) = ", format(perm_ls$p_r2y, digits = 3), ", p(Q2) = ",
          format(perm_ls$p_q2, digits = 3), ").")
  if (!is.na(perm_ls$p_q2) && perm_ls$p_q2 >= 0.05) {
    message("WARNING: permutation Q2 p-value >= 0.05; the model is not significantly ",
            "predictive. Interpret VIP/score plots with caution.")
  }

  # ---- VIP, loadings and S-plot coordinates ----------------------------------
  vip_vn <- ropls_get(opls_obj, "vip")
  loading_mn <- ropls_get(opls_obj, "loading")
  vip_df <- data.frame(metabolite_id = names(vip_vn), VIP = as.numeric(vip_vn),
                       stringsAsFactors = FALSE)
  zero_var <- slot(opls_obj, "xZeroVarVi")
  if (length(zero_var)) {
    message(length(zero_var), " zero-variance metabolite(s) were excluded by ropls: ",
            paste(names(zero_var), collapse = ", "))
  }
  # p(1): predictive loading; p(corr): Pearson correlation of each metabolite
  # with the predictive score t1 (Wiklund et al. 2008 S-plot definition;
  # correlation is invariant to the per-metabolite scaling ropls applied).
  p1_vn <- as.numeric(loading_mn[, 1])
  names(p1_vn) <- rownames(loading_mn)
  pcorr_vn <- apply(X, 2, function(col) stats::cor(col, t1, method = "pearson"))
  pcorr_vn[is.na(pcorr_vn)] <- 0

  # ---- Univariate tests -------------------------------------------------------
  uni <- univariate_tests(X, y, test)

  # ---- Assemble result tables --------------------------------------------------
  # Annotation columns (metabolite_name / class / kegg_id ...) are carried over
  # from --data so that downstream steps (03_pathway.R needs kegg_id,
  # 04_metabolite_classification.R needs class) can consume diff_metabolites.csv
  # directly.
  annot_cols <- setdiff(colnames(mat_df), c(id_col, shared))
  res <- merge(vip_df, uni, by = "metabolite_id", all = TRUE)
  load_df <- data.frame(metabolite_id = names(p1_vn), p1 = as.numeric(p1_vn),
                        p_corr = as.numeric(pcorr_vn[names(p1_vn)]), stringsAsFactors = FALSE)
  res <- merge(res, load_df, by = "metabolite_id", all = TRUE)
  if (length(annot_cols)) {
    annot <- mat_df[, c(id_col, annot_cols), drop = FALSE]
    colnames(annot)[1] <- "metabolite_id"
    annot[[1]] <- as.character(annot[[1]])
    res <- merge(res, annot, by = "metabolite_id", all.x = TRUE)
  }
  res$significant <- !is.na(res$VIP) & res$VIP >= vip_cut & !is.na(res$p) & res$p < pvalue
  res <- res[order(-res$VIP), ]
  core_cols <- c("metabolite_id", "VIP", "p", "padj", "log2FC", "significant")
  res <- res[, c(core_cols, setdiff(colnames(res), core_cols))]

  # All metabolites (full ranking) -> vip_scores.csv
  utils::write.csv(res[, core_cols], file.path(outdir, "vip_scores.csv"),
                   row.names = FALSE, quote = FALSE)
  # Significant metabolites (+ annotations) -> diff_metabolites.csv
  diff_df <- res[which(res$significant), , drop = FALSE]
  utils::write.csv(diff_df, file.path(outdir, "diff_metabolites.csv"),
                   row.names = FALSE, quote = FALSE)
  # S-plot coordinates -> oplsda_loadings.csv
  utils::write.csv(load_df, file.path(outdir, "oplsda_loadings.csv"),
                   row.names = FALSE, quote = FALSE)

  scores_df <- data.frame(sample_id = rownames(X), group = y, t1 = t1, to1 = to1,
                          stringsAsFactors = FALSE)
  utils::write.csv(scores_df, file.path(outdir, "oplsda_scores.csv"),
                   row.names = FALSE, quote = FALSE)

  n_sig <- nrow(diff_df)
  message("Differential metabolites (VIP >= ", vip_cut, " and p < ", pvalue, "): ", n_sig,
          " of ", nrow(res), " (", test, " test, BH-adjusted padj also reported).")
  if (!n_sig) {
    message("NOTE: no differential metabolites passed the thresholds; diff_metabolites.csv ",
            "only contains the header. Loosen --vip/--pvalue or check the data.")
  }

  # ---- Figures -----------------------------------------------------------------
  plot_oplsda_scores(scores_df, stats_ls, outdir)
  plot_s_plot(res, vip_cut, outdir)
  plot_vip_barplot(res, vip_cut, outdir)
  plot_permutation(perm_ls, n_perm, outdir)
  message("Saved figures (PDF/SVG/TIFF): pca_scores, oplsda_score_plot, s_plot, vip_barplot",
          if (!is.null(perm_ls$perm_r2y)) ", permutation_test" else "")

  # ---- Provenance ---------------------------------------------------------------
  write_run_info_r(
    outdir,
    script_name = "02_oplsda.R",
    params = list(
      data = data_path, metadata = meta_path, group_col = group_col,
      contrast = if (is.null(contrast)) paste(groups, collapse = ",") else contrast,
      vip_threshold = vip_cut, pvalue = pvalue, permutations = n_perm,
      univariate_test = test, scaleC = scale,
      n_samples = nrow(X), n_metabolites = ncol(X),
      R2Y = round(perm_ls$obs_r2y, 4), Q2 = round(perm_ls$obs_q2, 4),
      pR2Y = format(perm_ls$p_r2y, digits = 4), pQ2 = format(perm_ls$p_q2, digits = 4),
      n_differential = n_sig
    ),
    packages = c("ropls", "preprocessCore"),
    extra = c(
      "Differential criterion: VIP >= threshold AND raw p < pvalue (padj reported too).",
      "log2FC = mean(comparison) - mean(reference); equals log2 fold-change for log-transformed input.",
      "ropols slot names verified for ropls >= 1.24; run str(opls_obj) if extraction fails on another version."
    )
  )
  message("Done. Results in ", normalizePath(outdir))
  invisible(NULL)
}

main()
