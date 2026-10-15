###############################################################################
# 05_biomarker.R -- biomarker discovery: ROC/AUC + random forest (step 05)
#
# Title:    Single-metabolite biomarker screening (pROC) and multivariate ranking
#           (randomForest) for the two-group metabolome comparison
# Description:
#   For each candidate metabolite a ROC curve is fitted against the group label
#   (pROC), the AUC is reported with a DeLong 95% confidence interval and a
#   significance test against AUC = 0.5. A random-forest classifier (all
#   metabolites as predictors) provides multivariate MeanDecreaseGini /
#   MeanDecreaseAccuracy importance. The final biomarker panel is the
#   intersection of high-AUC metabolites (AUC > 0.9) and the top-20 random
#   forest metabolites.
#
# Usage:
#   Rscript scripts/metabolomics/05_biomarker.R \
#       --data results/metabolome/01_preprocess/metabolome_preprocessed.csv \
#       --metadata examples/metabolome_data/metadata.csv \
#       --outdir results/metabolome/05_biomarker \
#       --n-top 10 --ntree 500
#
# Inputs:
#   --data        Preprocessed abundance matrix (metabolite_id [, annotations],
#                 samples...), e.g. the output of 01_preprocess.py.
#   --metadata    Sample metadata CSV (sample_id, group).
#   --candidates  Optional: a file with one metabolite_id per line (a header
#                 column named metabolite_id is also accepted), or the literal
#                 "top10" (the 10 metabolites most correlated with the group).
#
# Outputs (written to --outdir):
#   biomarker_roc.csv     metabolite_id, AUC, ci_low, ci_high, p, higher_in
#   roc_coordinates.csv   metabolite_id, specificity, sensitivity (for re-plotting)
#   roc_curves.*          top candidate ROC curves (overlaid panel)
#   rf_importance.csv     metabolite_id, MeanDecreaseAccuracy, MeanDecreaseGini
#   rf_importance_plot.*  top-20 importance (horizontal bars)
#   biomarker_panel.csv   intersection of AUC > 0.9 and RF top-20
#   run_info.txt          parameters, package versions and the seed
#
# Notes:
#   - Requires pROC and randomForest: install.packages(c("pROC", "randomForest")).
#   - pROC orients ROC curves with direction="auto"; 'higher_in' records which
#     group has the higher mean abundance so the AUC direction is unambiguous.
#   - The AUC significance test (H0: AUC = 0.5) is attempted with
#     roc.test(roc, 0.5) and, if that fails on the installed pROC version, with
#     an auc-classed null value; a Wald approximation from the DeLong CI is the
#     last resort.
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

need_pkgs(c("pROC", "randomForest"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Read the abundance matrix + metadata and align them (samples x metabolites).
read_inputs <- function(data_path, meta_path, group_col) {
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
  non_numeric <- shared[!vapply(shared, function(s) is.numeric(mat_df[[s]]), logical(1))]
  if (length(non_numeric)) {
    stop("Non-numeric sample column(s) in --data: ", paste(head(non_numeric, 5), collapse = ", "),
         call. = FALSE)
  }
  mat_df <- mat_df[, c(id_col, shared), drop = FALSE]
  meta_df <- meta_df[match(shared, meta_df$sample_id), , drop = FALSE]
  X <- t(as.matrix(mat_df[, shared, drop = FALSE]))
  colnames(X) <- as.character(mat_df[[id_col]])
  rownames(X) <- shared
  storage.mode(X) <- "double"
  y <- factor(as.character(meta_df[[group_col]]))
  list(X = X, y = y)
}

#' Select candidate metabolites.
#'
#' Default: the n_top metabolites with the strongest absolute correlation to
#' the group label (a variance-based alternative is noted in the docs).
select_candidates <- function(X, y, candidates_arg, n_top) {
  all_ids <- colnames(X)
  if (!is.null(candidates_arg) && !isTRUE(candidates_arg)) {
    if (identical(candidates_arg, "top10")) {
      n_top <- 10L
    } else if (file.exists(candidates_arg)) {
      lines <- readLines(candidates_arg, warn = FALSE)
      lines <- trimws(lines)
      lines <- lines[nzchar(lines) & !grepl("^#", lines)]
      if (length(lines) && grepl("metabolite_id", lines[1], fixed = TRUE)) {
        df <- utils::read.csv(candidates_arg, stringsAsFactors = FALSE)
        col <- if ("metabolite_id" %in% colnames(df)) "metabolite_id" else colnames(df)[1]
        ids <- trimws(as.character(df[[col]]))
      } else {
        ids <- lines
      }
      ids <- unique(ids[nzchar(ids)])
      missing <- setdiff(ids, all_ids)
      if (length(missing)) {
        message(length(missing), " candidate ID(s) from ", candidates_arg,
                " are not in the data matrix and are ignored: ",
                paste(head(missing, 5), collapse = ", "))
      }
      ids <- intersect(ids, all_ids)
      if (!length(ids)) {
        stop("None of the candidate IDs in ", candidates_arg, " are present in --data.",
             call. = FALSE)
      }
      message("Using ", length(ids), " candidate metabolites from ", candidates_arg, ".")
      return(ids)
    } else {
      stop("--candidates must be a file path or the literal 'top10'; got: ", candidates_arg,
           call. = FALSE)
    }
  }
  # Default ranking: |correlation with the group| (group-correlated metabolites
  # are the most promising single biomarkers).
  y_num <- as.numeric(y)
  cor_vn <- vapply(all_ids, function(mid) abs(stats::cor(X[, mid], y_num)), numeric(1))
  head(names(sort(cor_vn, decreasing = TRUE)), n_top)
}

#' p-value of the AUC against 0.5 (three fallback strategies).
auc_pvalue <- function(roc_obj) {
  p <- tryCatch(as.numeric(pROC::roc.test(roc_obj, 0.5)$p.value), error = function(e) NA_real_)
  if (is.na(p)) {
    p <- tryCatch({
      auc_half <- structure(0.5, class = "auc")
      attr(auc_half, "roc") <- roc_obj
      as.numeric(pROC::roc.test(roc_obj, auc_half)$p.value)
    }, error = function(e) NA_real_)
  }
  if (is.na(p)) {
    # Last resort: Wald approximation from the DeLong CI (se from CI width).
    ci_vals <- suppressWarnings(as.numeric(unlist(pROC::ci.auc(roc_obj))))
    ci_vals <- sort(stats::na.omit(ci_vals))
    if (length(ci_vals) >= 2L) {
      se <- (ci_vals[length(ci_vals)] - ci_vals[1]) / (2 * stats::qnorm(0.975))
      if (is.finite(se) && se > 0) {
        z <- (as.numeric(roc_obj$auc) - 0.5) / se
        p <- 2 * stats::pnorm(-abs(z))
      }
    }
  }
  p
}

#' DeLong CI bounds from a pROC ci.auc result (robust to length 2 or 3).
ci_bounds <- function(roc_obj) {
  ci_obj <- tryCatch(pROC::ci.auc(roc_obj, conf.level = 0.95), error = function(e) NULL)
  if (is.null(ci_obj)) return(c(NA_real_, NA_real_))
  vals <- sort(stats::na.omit(as.numeric(unlist(ci_obj))))
  if (length(vals) >= 3L) return(c(vals[1], vals[3]))
  if (length(vals) == 2L) return(vals)
  c(NA_real_, NA_real_)
}

# ---------------------------------------------------------------------------
# ROC analysis
# ---------------------------------------------------------------------------

#' Fit ROC curves for all candidate metabolites.
#'
#' @return List with the summary data frame and the coordinate data frame.
run_roc_analysis <- function(X, y, candidates) {
  levels_vc <- levels(y)
  rows <- lapply(candidates, function(mid) {
    values <- X[, mid]
    roc_obj <- pROC::roc(response = y, predictor = values, levels = levels_vc,
                         direction = "auto", quiet = TRUE)
    auc_val <- suppressWarnings(as.numeric(pROC::auc(roc_obj)))
    ci_b <- ci_bounds(roc_obj)
    in_ref <- y == levels_vc[1]
    in_cmp <- y == levels_vc[2]
    higher_in <- if (mean(values[in_cmp]) >= mean(values[in_ref])) levels_vc[2] else levels_vc[1]
    data.frame(
      metabolite_id = mid,
      AUC = auc_val,
      ci_low = ci_b[1],
      ci_high = ci_b[2],
      p = auc_pvalue(roc_obj),
      higher_in = higher_in,
      stringsAsFactors = FALSE
    )
  })
  roc_df <- do.call(rbind, rows)
  rownames(roc_df) <- NULL
  # Coordinates for re-plotting (specificity/sensitivity per threshold).
  coord <- lapply(candidates, function(mid) {
    roc_obj <- pROC::roc(response = y, predictor = X[, mid], levels = levels_vc,
                         direction = "auto", quiet = TRUE)
    data.frame(metabolite_id = mid,
               specificity = as.numeric(roc_obj$specificities),
               sensitivity = as.numeric(roc_obj$sensitivities),
               stringsAsFactors = FALSE)
  })
  coord_df <- do.call(rbind, coord)
  rownames(coord_df) <- NULL
  list(roc_df = roc_df, coord_df = coord_df)
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

plot_roc_curves <- function(roc_df, coord_df, outdir, n_show = 6) {
  # Scientific conclusion: metabolites with AUC close to 1 discriminate the two
  # groups almost perfectly and are the primary biomarker candidates.
  top <- utils::head(roc_df[order(-roc_df$AUC), ], min(n_show, nrow(roc_df)))
  coords <- coord_df[coord_df$metabolite_id %in% top$metabolite_id, , drop = FALSE]
  auc_map <- setNames(round(top$AUC, 3), top$metabolite_id)
  coords$label <- ifelse(coords$metabolite_id %in% names(auc_map),
                         paste0(coords$metabolite_id, " (AUC ", auc_map[coords$metabolite_id], ")"),
                         coords$metabolite_id)
  p <- ggplot2::ggplot(coords, ggplot2::aes(x = 1 - specificity, y = sensitivity,
                                            colour = label)) +
    ggplot2::geom_path(linewidth = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = nature_semantic("reference_line"),
                         linewidth = 0.4, linetype = "dashed") +
    scale_colour_nature(name = "Metabolite") +
    ggplot2::labs(title = "ROC curves of biomarker candidates",
                  subtitle = paste0("Top ", length(unique(coords$label)),
                                    " metabolites by AUC; dashed diagonal = random classifier"),
                  x = "1 - specificity (false positive rate)",
                  y = "Sensitivity (true positive rate)") +
    ggplot2::coord_equal() +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "roc_curves"), width_mm = 89, height_mm = 80)
}

plot_rf_importance <- function(imp_df, outdir, top_n = 20) {
  # Scientific conclusion: MeanDecreaseGini ranks metabolites by their
  # multivariate contribution to group classification; the top metabolites are
  # robust candidate biomarkers.
  df <- utils::head(imp_df, top_n)
  df$metabolite_id <- factor(df$metabolite_id, levels = rev(df$metabolite_id))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = metabolite_id, y = MeanDecreaseGini)) +
    ggplot2::geom_col(fill = rice_palette()$series12[1], width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Random forest importance (top 20)",
                  subtitle = paste0("MeanDecreaseGini from randomForest (all metabolites as ",
                                    "predictors); higher = more important for classification"),
                  x = NULL, y = "MeanDecreaseGini") +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "rf_importance_plot"), width_mm = 89, height_mm = 100)
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
  outdir <- opt(opts, "outdir", "results/metabolome/05_biomarker")
  group_col <- opt(opts, "group-col", "group")
  candidates_arg <- opt(opts, "candidates", NULL)
  n_top <- as.integer(opt(opts, "n-top", 10))
  ntree <- as.integer(opt(opts, "ntree", 500))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", normalizePath(outdir))

  inputs <- read_inputs(data_path, meta_path, group_col)
  X <- inputs$X
  y <- inputs$y
  message("Model matrix: ", nrow(X), " samples x ", ncol(X), " metabolites; groups: ",
          paste(levels(y), collapse = " vs "), " (n = ", paste(table(y), collapse = " / "), ").")
  if (any(table(y) < 2)) {
    stop("ROC/AUC and randomForest need at least 2 samples per group.", call. = FALSE)
  }

  candidates <- select_candidates(X, y, candidates_arg, n_top)
  message("Candidate metabolites: ", length(candidates), " (e.g. ",
          paste(head(candidates, 5), collapse = ", "), ").")

  # ---- ROC / AUC ---------------------------------------------------------------
  roc_res <- run_roc_analysis(X, y, candidates)
  roc_df <- roc_res$roc_df[order(-roc_df$AUC), , drop = FALSE]
  rownames(roc_df) <- NULL
  utils::write.csv(roc_df, file.path(outdir, "biomarker_roc.csv"),
                   row.names = FALSE, quote = FALSE)
  utils::write.csv(roc_res$coord_df, file.path(outdir, "roc_coordinates.csv"),
                   row.names = FALSE, quote = FALSE)
  plot_roc_curves(roc_df, roc_res$coord_df, outdir)
  n_auc09 <- sum(roc_df$AUC > 0.9, na.rm = TRUE)
  message("ROC: ", n_auc09, " of ", nrow(roc_df), " candidates have AUC > 0.9 (best: ",
          roc_df$metabolite_id[1], ", AUC ", round(roc_df$AUC[1], 3), ").")

  # ---- Random forest -------------------------------------------------------------
  message("Fitting randomForest (ntree = ", ntree, ", importance = TRUE) on all ",
          ncol(X), " metabolites...")
  rf <- randomForest::randomForest(x = X, y = y, ntree = ntree, importance = TRUE)
  imp <- rf$importance
  imp_df <- data.frame(
    metabolite_id = rownames(imp),
    MeanDecreaseAccuracy = as.numeric(imp[, "MeanDecreaseAccuracy"]),
    MeanDecreaseGini = as.numeric(imp[, "MeanDecreaseGini"]),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(-imp_df$MeanDecreaseGini), , drop = FALSE]
  rownames(imp_df) <- NULL
  utils::write.csv(imp_df, file.path(outdir, "rf_importance.csv"),
                   row.names = FALSE, quote = FALSE)
  plot_rf_importance(imp_df, outdir, top_n = 20)

  # ---- Biomarker panel -------------------------------------------------------------
  rf_top20 <- utils::head(imp_df$metabolite_id, 20)
  high_auc <- roc_df$metabolite_id[which(roc_df$AUC > 0.9)]
  panel_ids <- intersect(high_auc, rf_top20)
  if (length(panel_ids)) {
    panel <- merge(roc_df, imp_df, by = "metabolite_id")
    panel <- panel[panel$metabolite_id %in% panel_ids, , drop = FALSE]
    panel <- panel[order(-panel$AUC), , drop = FALSE]
    panel <- panel[, c("metabolite_id", "AUC", "ci_low", "ci_high", "p", "higher_in",
                       "MeanDecreaseAccuracy", "MeanDecreaseGini")]
    rownames(panel) <- NULL
    message("Biomarker panel (AUC > 0.9 and RF top-20): ", nrow(panel), " metabolites: ",
            paste(panel$metabolite_id, collapse = ", "))
  } else {
    panel <- data.frame(metabolite_id = character(0), AUC = numeric(0), ci_low = numeric(0),
                        ci_high = numeric(0), p = numeric(0), higher_in = character(0),
                        MeanDecreaseAccuracy = numeric(0), MeanDecreaseGini = numeric(0))
    message("No metabolite satisfies both AUC > 0.9 and RF top-20; biomarker_panel.csv ",
            "contains only the header.")
  }
  utils::write.csv(panel, file.path(outdir, "biomarker_panel.csv"),
                   row.names = FALSE, quote = FALSE)

  write_run_info_r(
    outdir,
    script_name = "05_biomarker.R",
    params = list(
      data = data_path, metadata = meta_path, group_col = group_col,
      candidates = if (is.null(candidates_arg)) paste0("top", n_top, " by |cor| with group")
                   else candidates_arg,
      n_candidates = length(candidates), ntree = ntree,
      n_auc_above_0.9 = n_auc09,
      n_panel = nrow(panel),
      rf_oob_error = round(tail(rf$err.rate[, "OOB"], 1), 4)
    ),
    packages = c("pROC", "randomForest"),
    extra = c(
      "AUC direction: pROC direction='auto'; 'higher_in' gives the group with higher mean.",
      "AUC p-value: H0 = 0.5 (roc.test with fallbacks; see script notes).",
      "Biomarker panel = intersection of AUC > 0.9 and random-forest top-20 (MeanDecreaseGini)."
    )
  )
  message("Done. Results in ", normalizePath(outdir))
  invisible(NULL)
}

main()
