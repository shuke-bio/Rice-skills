###############################################################################
# 05_plot_admixture.R -- plot ADMIXTURE ancestry matrices (.Q) and CV errors
#
# Title:    Rice-skills admixture-plot step
# Description:
#   Draws the classic stacked-bar ancestry plot for one or more K values and,
#   when a CV-error file is supplied, the cross-validation curve used to pick
#   K. Samples are ordered by dominant ancestry component (standard display).
#
# Usage:
#   Rscript scripts/popgen/05_plot_admixture.R \
#       --q-dir results/popgen/04_admixture --cv results/popgen/04_admixture/cv_errors.txt \
#       --outdir results/popgen/05_admixture_figures
#   # or explicitly:
#   Rscript scripts/popgen/05_plot_admixture.R \
#       --q-files "run.2.Q,run.3.Q" --outdir ...
#
# Outputs (per K, PDF + SVG + TIFF):
#   admixture_k<K>.{pdf,svg,tiff}
#   cv_curve.{pdf,svg,tiff}     (only with --cv)
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
q_files_arg <- opt(opts, "q-files")
q_dir <- opt(opts, "q-dir")
cv_path <- opt(opts, "cv")
outdir <- opt(opts, "outdir", "results/popgen/05_admixture_figures")

q_files <- character(0)
if (!is.null(q_files_arg)) {
  q_files <- trimws(strsplit(q_files_arg, ",")[[1]])
} else if (!is.null(q_dir)) {
  q_files <- list.files(q_dir, pattern = "\\.Q$", full.names = TRUE)
}
if (length(q_files) == 0) {
  stop("No .Q files found. Provide --q-files or --q-dir.\n",
       "Hint: .Q matrices come from scripts/popgen/04_admixture.sh (or any ",
       "ADMIXTURE run / web-server export).", call. = FALSE)
}
missing_q <- q_files[!file.exists(q_files)]
if (length(missing_q) > 0) {
  stop("Q file(s) not found: ", paste(missing_q, collapse = ", "), call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pal <- rice_palette()
n_done <- 0L
for (qf in q_files) {
  q <- as.matrix(utils::read.table(qf, header = FALSE))
  k <- ncol(q)
  if (k < 2) {
    message("Skipping ", basename(qf), " (K = ", k, "); nothing to stack.")
    next
  }
  # order samples by dominant component, then by that component's proportion
  ord <- order(max.col(q), -apply(q, 1, max))
  q <- q[ord, , drop = FALSE]
  df <- data.frame(
    sample_idx = rep(seq_len(nrow(q)), times = k),
    cluster = factor(rep(paste0("C", seq_len(k)), each = nrow(q)),
                     levels = paste0("C", seq_len(k))),
    prop = as.vector(q)
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample_idx, y = prop, fill = cluster)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = pal$series12[seq_len(k)]) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = "Individuals (ordered by ancestry)", y = "Ancestry proportion",
                  fill = "Cluster", title = sprintf("ADMIXTURE ancestry (K = %d)", k)) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "right")
  save_nature_figure(p, file.path(outdir, paste0("admixture_k", k)),
                     width_mm = 150, height_mm = 60)
  message(sprintf("K = %d: %d individuals plotted.", k, nrow(q)))
  n_done <- n_done + 1L
}

# --- CV error curve ------------------------------------------------------------
if (!is.null(cv_path)) {
  if (!file.exists(cv_path)) {
    message("CV file not found: ", cv_path, " - skipping the CV curve.")
  } else {
    lines <- readLines(cv_path, warn = FALSE)
    kk <- as.numeric(sub(".*K *= *([0-9]+).*", "\\1", lines[grepl("CV error", lines)]))
    ee <- as.numeric(sub(".*: *([0-9.eE+-]+).*", "\\1", lines[grepl("CV error", lines)]))
    ok <- is.finite(kk) & is.finite(ee)
    if (sum(ok) >= 2) {
      cv_df <- data.frame(K = kk[ok], cv_error = ee[ok])
      best <- cv_df$K[which.min(cv_df$cv_error)]
      message(sprintf("Lowest CV error at K = %d (%.4f).", best, min(cv_df$cv_error)))
      p_cv <- ggplot2::ggplot(cv_df, ggplot2::aes(x = K, y = cv_error)) +
        ggplot2::geom_line(colour = pal$series12[1], linewidth = 0.7) +
        ggplot2::geom_point(colour = pal$series12[1], size = 1.6) +
        ggplot2::geom_point(data = cv_df[cv_df$K == best, , drop = FALSE],
                            colour = nature_semantic("highlight"), size = 2.4) +
        ggplot2::scale_x_continuous(breaks = cv_df$K) +
        ggplot2::labs(x = "K", y = "Cross-validation error",
                      title = "ADMIXTURE CV error",
                      subtitle = sprintf("Best K = %d", best)) +
        theme_nature()
      save_nature_figure(p_cv, file.path(outdir, "cv_curve"), width_mm = 89, height_mm = 70)
    } else {
      message("CV file had fewer than 2 parseable 'CV error (K=..)' lines - skipped.")
    }
  }
}

write_run_info_r(
  outdir, "05_plot_admixture.R",
  params = list(q_files = paste(basename(q_files), collapse = ","),
                cv = if (is.null(cv_path)) NA else cv_path,
                n_q_plotted = n_done, seed = RICE_SEED),
  packages = c("ggplot2")
)
message("Admixture plotting complete: ", n_done, " K value(s). Results in ", outdir)
