###############################################################################
# 04_metabolite_classification.R -- chemical classification of differential metabolites
#
# Title:    Chemical class statistics for differential metabolites (step 04)
# Description:
#   Counts differential metabolites per chemical class, compares the class
#   distribution with the full-metabolome background (Fisher's exact test per
#   class) and produces a publication-ready barplot. If UpSetR is installed,
#   an UpSet plot of class membership is drawn as well.
#
# Usage:
#   Rscript scripts/metabolomics/04_metabolite_classification.R \
#       --diff-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
#       --all-metabolites results/metabolome/01_preprocess/metabolome_preprocessed.csv \
#       --outdir results/metabolome/04_classification --top-n 15
#
# Inputs:
#   --diff-metabolites  CSV with a 'class' column (diff_metabolites.csv from
#                       02_oplsda.R). If a 'significant' column is present,
#                       only rows flagged TRUE are used.
#   --all-metabolites   Optional background matrix (any file with a 'class'
#                       column, e.g. the preprocessed matrix from
#                       01_preprocess.py). Without it, Fisher p-values are NA.
#
# Outputs (written to --outdir):
#   class_counts.csv   class, n_diff, n_all, ratio_diff, ratio_all, p, padj
#   class_barplot.*    top-N classes (up/down stacked when log2FC is available)
#   class_upset.*      UpSet of class membership (only if UpSetR is installed)
#   run_info.txt       parameters, package versions and the seed
#
# Notes:
#   - Fisher's exact test per class: 2x2 table of (in class / not in class) x
#     (differential / not differential), tested with base R fisher.test.
#   - When log2FC is available the barplot is stacked by regulation direction;
#     scale_fill_nature() is used with levels ordered c("down", "up") so that
#     the first two palette colours coincide with the semantic down (#5B84AE)
#     and up (#CE6A6B) colours.
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
# Input handling
# ---------------------------------------------------------------------------

#' Read the differential-metabolite table and its class annotations.
read_diff_with_class <- function(path) {
  df <- read_csv_safe(path, "differential metabolites table")
  check_columns(df, "class", "differential metabolites table")
  if ("significant" %in% colnames(df)) {
    keep <- as.logical(df$significant)
    keep[is.na(keep)] <- FALSE
    message("Filtering on the 'significant' column: ", sum(keep), " of ", nrow(df),
            " rows are flagged significant.")
    df <- df[keep, , drop = FALSE]
  }
  df$class <- trimws(as.character(df$class))
  df$class[is.na(df$class) | !nzchar(df$class)] <- "Unclassified"
  list(df = df, has_log2fc = "log2FC" %in% colnames(df))
}

#' Read the background metabolite set (any CSV with a 'class' column).
read_background_classes <- function(path) {
  df <- read_csv_safe(path, "background metabolite table")
  check_columns(df, "class", "background metabolite table")
  cls <- trimws(as.character(df$class))
  cls[is.na(cls) | !nzchar(cls)] <- "Unclassified"
  cls
}

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

#' Class counts with Fisher exact tests against the background.
#'
#' @param diff_df Differential metabolite table (must contain 'class').
#' @param bg_classes Optional character vector of background classes.
#' @return Data frame with class, n_diff, n_all, ratio_diff, ratio_all, p, padj.
class_statistics <- function(diff_df, bg_classes = NULL) {
  diff_classes <- diff_df$class
  diff_tab <- table(diff_classes)
  n_diff_total <- length(diff_classes)
  if (is.null(bg_classes)) {
    message("No --all-metabolites background provided: Fisher p-values are reported as NA.")
    classes <- names(diff_tab)
    n_diff <- as.integer(diff_tab)
    n_all <- rep(NA_integer_, length(classes))
    n_all_total <- NA_integer_
  } else {
    bg_tab <- table(bg_classes)
    # Union of classes so that classes present in only one set are kept.
    classes <- sort(union(names(diff_tab), names(bg_tab)))
    n_diff <- as.integer(vapply(classes, function(cl) sum(diff_classes == cl), integer(1)))
    n_all <- as.integer(vapply(classes, function(cl) sum(bg_classes == cl), integer(1)))
    n_all_total <- length(bg_classes)
  }
  names(n_diff) <- classes
  names(n_all) <- classes
  p_vn <- vapply(classes, function(cl) {
    if (is.null(bg_classes)) return(NA_real_)
    a <- n_diff[[cl]]                      # differential in class
    b <- n_diff_total - a                  # differential not in class
    c <- n_all[[cl]] - a                   # background in class
    d <- (n_all_total - n_all[[cl]]) - b   # background not in class
    stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value
  }, numeric(1))
  data.frame(
    class = classes,
    n_diff = n_diff,
    n_all = n_all,
    ratio_diff = n_diff / n_diff_total,
    ratio_all = if (is.null(bg_classes)) NA_real_ else n_all / n_all_total,
    p = p_vn,
    padj = stats::p.adjust(p_vn, method = "BH"),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

plot_class_barplot <- function(counts, diff_df, has_log2fc, top_n, outdir) {
  # Scientific conclusion: the differential pool is dominated by a few chemical
  # classes; over-represented classes (low padj) point at the biochemical
  # processes responding to the treatment.
  df <- utils::head(counts[order(-counts$n_diff), ], top_n)
  df <- df[order(df$n_diff), ]  # ascending for horizontal bars (largest on top)
  df$class <- factor(df$class, levels = df$class)
  if (has_log2fc && !is.null(diff_df$log2FC)) {
    dir_df <- diff_df[diff_df$class %in% df$class, c("class", "log2FC"), drop = FALSE]
    dir_df$direction <- ifelse(dir_df$log2FC >= 0, "up", "down")
    # levels ordered c("down","up"): scale_fill_nature() then maps
    # down -> series12[1] (#5B84AE) and up -> series12[2] (#CE6A6B),
    # which coincide with the semantic down/up colours.
    dir_df$direction <- factor(dir_df$direction, levels = c("down", "up"))
    stack <- as.data.frame(table(class = dir_df$class, direction = dir_df$direction))
    names(stack) <- c("class", "direction", "n")
    stack$class <- factor(stack$class, levels = levels(df$class))
    p <- ggplot2::ggplot(stack, ggplot2::aes(x = class, y = n, fill = direction)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::coord_flip() +
      scale_fill_nature(name = "Direction",
                        labels = c(down = "down (log2FC < 0)", up = "up (log2FC > 0)")) +
      ggplot2::labs(title = "Chemical classes of differential metabolites",
                    subtitle = paste0("Top ", nrow(df), " classes, stacked by regulation ",
                                      "direction; colour ramp = Nature series (down = blue, up = red)"),
                    x = NULL, y = "Number of differential metabolites") +
      theme_nature()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = class, y = n_diff, fill = class)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::coord_flip() +
      scale_fill_nature(name = "Chemical class") +
      ggplot2::labs(title = "Chemical classes of differential metabolites",
                    subtitle = paste0("Top ", nrow(df), " classes by differential count ",
                                      "(no log2FC column: direction stacking unavailable)"),
                    x = NULL, y = "Number of differential metabolites") +
      theme_nature()
  }
  save_nature_figure(p, file.path(outdir, "class_barplot"), width_mm = 89, height_mm = 100)
}

plot_class_upset <- function(diff_df, counts, top_n, outdir) {
  # UpSet of class membership among differential metabolites (top-N classes).
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    message("UpSetR is not installed; skipping the UpSet plot. ",
            "Install with install.packages('UpSetR') to enable it.")
    return(invisible(NULL))
  }
  top_classes <- utils::head(counts[order(-counts$n_diff), "class"], top_n)
  elements <- diff_df$class[diff_df$class %in% top_classes]
  if (length(unique(elements)) < 2L) {
    message("Fewer than two classes among the differential metabolites; skipping the UpSet plot.")
    return(invisible(NULL))
  }
  input_list <- split(seq_along(elements), elements)
  # UpSetR works best with a named list of sets; metabolite IDs label the rows.
  met_ids <- if ("metabolite_id" %in% colnames(diff_df)) diff_df$metabolite_id else
    paste0("met_", seq_along(elements))
  names(input_list) <- top_classes[top_classes %in% names(input_list)]
  input_list <- lapply(input_list, function(idx) met_ids[idx])
  stem <- file.path(outdir, "class_upset")
  for (fmt in c("pdf", "tiff", "png")) {
    out <- paste0(stem, ".", fmt)
    if (fmt == "tiff") {
      grDevices::tiff(out, width = mm_to_inch(89), height = mm_to_inch(100),
                      units = "in", res = 300, compression = "lzw")
    } else if (fmt == "pdf") {
      grDevices::pdf(out, width = mm_to_inch(89), height = mm_to_inch(100))
    } else {
      grDevices::png(out, width = mm_to_inch(89), height = mm_to_inch(100), units = "in", res = 300)
    }
    print(UpSetR::upset(UpSetR::fromList(input_list), order.by = "freq",
                        main.bar.color = rice_palette()$series12[1],
                        sets.bar.color = rice_palette()$series12[1]))
    grDevices::dev.off()
  }
  message("Saved figures: class_upset (PDF/TIFF/PNG)")
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  diff_path <- opt(opts, "diff-metabolites")
  if (is.null(diff_path)) stop("--diff-metabolites is required (CSV with a 'class' column).", call. = FALSE)
  outdir <- opt(opts, "outdir", "results/metabolome/04_classification")
  all_path <- opt(opts, "all-metabolites", NULL)
  top_n <- as.integer(opt(opts, "top-n", 15))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", normalizePath(outdir))

  input <- read_diff_with_class(diff_path)
  diff_df <- input$df
  if (!nrow(diff_df)) {
    message("The differential-metabolite table is empty; writing header-only outputs.")
    utils::write.csv(data.frame(class = character(0), n_diff = integer(0), n_all = integer(0),
                                ratio_diff = numeric(0), ratio_all = numeric(0),
                                p = numeric(0), padj = numeric(0)),
                     file.path(outdir, "class_counts.csv"), row.names = FALSE, quote = FALSE)
    write_run_info_r(outdir, script_name = "04_metabolite_classification.R",
                     params = list(diff_metabolites = diff_path, n_differential = 0),
                     extra = "Empty input: no differential metabolites to classify.")
    return(invisible(NULL))
  }
  message("Differential metabolites: ", nrow(diff_df), " rows, ", length(unique(diff_df$class)),
          " chemical classes.")

  bg_classes <- if (!is.null(all_path)) read_background_classes(all_path) else NULL
  if (!is.null(bg_classes)) {
    message("Background: ", length(bg_classes), " metabolites (", length(unique(bg_classes)),
            " classes) from ", all_path, ".")
  }

  counts <- class_statistics(diff_df, bg_classes)
  counts <- counts[order(-counts$n_diff), , drop = FALSE]
  rownames(counts) <- NULL
  utils::write.csv(counts, file.path(outdir, "class_counts.csv"),
                   row.names = FALSE, quote = FALSE)

  if (!is.null(bg_classes)) {
    sig <- counts[which(counts$padj <= 0.05), , drop = FALSE]
    if (nrow(sig)) {
      message("Over-represented classes (Fisher, padj <= 0.05):")
      for (i in seq_len(nrow(sig))) {
        message("  ", sig$class[i], ": ", sig$n_diff[i], "/", sig$n_all[i],
                " (diff ratio ", round(sig$ratio_diff[i], 3), " vs background ",
                round(sig$ratio_all[i], 3), ", padj ", format(sig$padj[i], digits = 3), ")")
      }
    } else {
      message("No chemical class is over-represented at padj <= 0.05.")
    }
  }

  plot_class_barplot(counts, diff_df, input$has_log2fc, top_n, outdir)
  plot_class_upset(diff_df, counts, top_n, outdir)

  write_run_info_r(
    outdir,
    script_name = "04_metabolite_classification.R",
    params = list(
      diff_metabolites = diff_path,
      all_metabolites = if (is.null(all_path)) "none" else all_path,
      top_n = top_n,
      n_differential = nrow(diff_df),
      n_classes = length(unique(diff_df$class)),
      n_background = if (is.null(bg_classes)) 0 else length(bg_classes),
      fisher_background = if (is.null(bg_classes)) "none (p = NA)" else "provided"
    ),
    extra = c(
      "Fisher exact test per class: 2x2 (class x differential) against the background.",
      "UpSet plot requires UpSetR (not installed -> skipped with a message)."
    )
  )
  message("Done. Results in ", normalizePath(outdir))
  invisible(NULL)
}

main()
