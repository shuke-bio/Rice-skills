###############################################################################
# 02_deseq2.R -- differential expression analysis with DESeq2
#
# Title:    Differential expression analysis (DESeq2) for the rice transcriptome
# Description:
#   Reads the (optionally pre-filtered) count matrix and sample metadata,
#   aligns the samples, builds a DESeqDataSet with a configurable design
#   formula (default "~ group"; multi-factor designs such as "~ batch + group"
#   are supported when the metadata contains the referenced columns), runs
#   DESeq() and extracts results for one or more contrasts (format
#   <numerator>_vs_<denominator>). Writes full and significant DE tables per
#   contrast, a VST-normalised count matrix, a dispersion-estimate diagnostic
#   figure and a summary of up/down counts.
#
# Usage:
#   Rscript scripts/transcriptomics/02_deseq2.R \
#       --counts results/transcriptome/01_qc/counts_filtered.csv \
#       --metadata metadata.csv \
#       --contrast salt_vs_control \
#       --outdir results/transcriptome/02_de
#   # multiple contrasts (comma-separated):
#   Rscript scripts/transcriptomics/02_deseq2.R --counts ... \
#       --contrast salt_vs_control,drought_vs_control
#   # multi-factor design (every formula variable must be a metadata column):
#   Rscript scripts/transcriptomics/02_deseq2.R --counts ... --design "~ batch + group"
#   # LFC shrinkage (recommended for reporting and visualisation):
#   Rscript scripts/transcriptomics/02_deseq2.R --counts ... --lfc-shrink
#
# Inputs:
#   --counts    count matrix CSV (first column gene_id, columns = samples)
#   --metadata  sample metadata CSV (first column sample_id, column <group-col>)
#
# Outputs (written to --outdir):
#   de_all_<contrast>.csv                   full DE results per contrast
#   de_sig_<contrast>.csv                   significant DE genes (Up/Down)
#   normalized_counts_vst.csv               VST-normalised counts (blind=FALSE)
#   dispersion_estimates.{pdf,svg,tiff}     DESeq2 dispersion diagnostic
#   de_summary.csv                          up/down/ns counts per contrast
#   run_info.txt                            provenance (parameters, versions, seed)
#
# Notes:
#   - Design: multi-factor or paired designs are configured with --design
#     (e.g. "~ batch + group" or "~ pair + group"); every variable used in
#     the formula must exist as a metadata column, otherwise the script stops
#     with a friendly message.
#   - Time-series experiments are approximated by setting the group column to
#     the sampling time point and using the pairwise contrasts (a dedicated
#     likelihood-ratio workflow is not exposed here).
#   - The first group level appearing in the metadata acts as the reference
#     level of the design; contrasts are evaluated explicitly with
#     results(contrast = ...), so the reference level does not affect the
#     reported log2 fold changes.
#   - DESeq2 needs biological replicates: groups with a single sample produce
#     unreliable dispersion estimates (a warning is printed).
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

counts_path <- opt(opts, "counts")
metadata_path <- opt(opts, "metadata")
outdir <- opt(opts, "outdir", "results/transcriptome/02_de")
group_col <- opt(opts, "group-col", "group")
design <- opt(opts, "design", "~ group")
contrast_arg <- opt(opts, "contrast")
padj_cut <- suppressWarnings(as.numeric(opt(opts, "padj", "0.05")))
lfc_cut <- suppressWarnings(as.numeric(opt(opts, "log2fc", "1")))
# Flags: query with a NULL default so that a bare --flag evaluates to TRUE.
lfc_shrink <- isTRUE(opt(opts, "lfc-shrink", NULL))

if (is.null(counts_path)) {
  stop("--counts is required (count matrix CSV, e.g. results/transcriptome/01_qc/counts_filtered.csv).",
       call. = FALSE)
}
if (is.null(metadata_path)) {
  stop("--metadata is required (sample metadata CSV with a 'group' column).", call. = FALSE)
}
if (is.na(padj_cut) || padj_cut <= 0 || padj_cut >= 1) {
  stop("--padj must be a number in (0, 1); got '", opt(opts, "padj", "0.05"), "'.", call. = FALSE)
}
if (is.na(lfc_cut) || lfc_cut < 0) {
  stop("--log2fc must be a non-negative number; got '", opt(opts, "log2fc", "1"), "'.", call. = FALSE)
}

need_pkgs(c("DESeq2", "ggplot2", "jsonlite"))

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers (script-local equivalents of the Python data_io layer)
# ---------------------------------------------------------------------------

#' Read a count matrix CSV (first column = gene_id) into a numeric matrix.
read_counts_matrix <- function(path) {
  df <- read_csv_safe(path, "count matrix")
  gene_ids <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_ids
  if (any(!is.finite(mat))) {
    stop("The count matrix contains missing or non-finite values: ", path,
         "\nHint: counts must be non-negative integers without gaps.", call. = FALSE)
  }
  if (any(mat < 0)) {
    stop("The count matrix contains negative values: ", path,
         "\nHint: DESeq2 expects raw non-negative integer counts.", call. = FALSE)
  }
  mat
}

#' Read metadata and align it with a count matrix (strict sample matching).
read_metadata_aligned <- function(path, mat) {
  meta <- read_csv_safe(path, "sample metadata")
  colnames(meta)[1] <- "sample_id"
  meta$sample_id <- as.character(meta$sample_id)
  shared <- meta$sample_id[meta$sample_id %in% colnames(mat)]
  only_mat <- setdiff(colnames(mat), meta$sample_id)
  only_meta <- setdiff(meta$sample_id, colnames(mat))
  if (length(shared) == 0L) {
    stop("No overlapping samples between the count matrix and the metadata.\n",
         "matrix samples:   ", paste(head(colnames(mat), 5), collapse = ", "), " ...\n",
         "metadata samples: ", paste(head(meta$sample_id, 5), collapse = ", "), " ...\n",
         "Hint: sample IDs must match exactly (case-sensitive, no extra spaces).",
         call. = FALSE)
  }
  if (length(only_mat) > 0L || length(only_meta) > 0L) {
    stop("Sample mismatch between the count matrix and the metadata.\n",
         "only in matrix:   ", paste(only_mat, collapse = ", "), "\n",
         "only in metadata: ", paste(only_meta, collapse = ", "), "\n",
         "Hint: keep the two files consistent before running the analysis.",
         call. = FALSE)
  }
  meta <- meta[match(shared, meta$sample_id), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  list(matrix = mat[, shared, drop = FALSE], meta = meta)
}

#' Save base-graphics output in the three Nature formats (PDF, SVG, TIFF 300 dpi).
save_base_figure <- function(draw, file_stem, width_mm = 89, height_mm = 89) {
  for (f in c("pdf", "svg", "tiff")) {
    out <- paste0(file_stem, ".", f)
    if (f == "pdf") {
      grDevices::pdf(out, width = mm_to_inch(width_mm), height = mm_to_inch(height_mm))
    } else if (f == "svg") {
      grDevices::svg(out, width = mm_to_inch(width_mm), height = mm_to_inch(height_mm))
    } else {
      grDevices::tiff(out, width = mm_to_inch(width_mm), height = mm_to_inch(height_mm),
                      units = "in", res = 300, compression = "lzw")
    }
    draw()
    grDevices::dev.off()
  }
  invisible(NULL)
}

#' Parse a contrast specification "<numerator>_vs_<denominator>".
parse_contrast <- function(spec, group_col, levels) {
  parts <- trimws(strsplit(spec, "_vs_", fixed = TRUE)[[1]])
  if (length(parts) != 2L || !nzchar(parts[1]) || !nzchar(parts[2])) {
    stop("Invalid contrast '", spec, "'. Expected format: <numerator>_vs_<denominator> ",
         "(e.g. salt_vs_control).", call. = FALSE)
  }
  unknown <- setdiff(parts, levels)
  if (length(unknown) > 0L) {
    stop("Contrast '", spec, "' references unknown group(s): ",
         paste(unknown, collapse = ", "),
         "\nAvailable groups: ", paste(levels, collapse = ", "), call. = FALSE)
  }
  list(name = paste0(parts[1], "_vs_", parts[2]),
       numerator = parts[1], denominator = parts[2])
}

#' File-safe tag for a contrast name.
contrast_tag <- function(name) gsub("[^A-Za-z0-9._-]", "_", name)

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

counts <- read_counts_matrix(counts_path)
aligned <- read_metadata_aligned(metadata_path, counts)
mat <- aligned$matrix
meta <- aligned$meta
message(sprintf("Loaded %d genes x %d samples.", nrow(mat), ncol(mat)))

if (!(group_col %in% colnames(meta))) {
  stop("The group column '", group_col, "' is missing from the metadata.\n",
       "Available columns: ", paste(colnames(meta), collapse = ", "), "\n",
       "Hint: pass the correct name with --group-col.", call. = FALSE)
}
groups <- unique(meta[[group_col]])
if (length(groups) < 2L) {
  stop("At least two groups are required in the '", group_col, "' column; found: ",
       paste(groups, collapse = ", "), call. = FALSE)
}
group_sizes <- table(meta[[group_col]])
if (any(group_sizes < 2L)) {
  message("WARNING: group(s) with fewer than 2 replicates: ",
          paste(names(group_sizes)[group_sizes < 2L], collapse = ", "),
          ". DESeq2 dispersion estimates are unreliable without replicates.")
}

# Design formula: replace the literal token "group" with the actual group column.
design_formula <- stats::as.formula(
  gsub("\\bgroup\\b", group_col, design, perl = TRUE)
)
design_vars <- all.vars(design_formula)
missing_vars <- setdiff(design_vars, colnames(meta))
if (length(missing_vars) > 0L) {
  stop("The design formula '", design, "' references column(s) missing from the ",
       "metadata: ", paste(missing_vars, collapse = ", "), "\n",
       "Available columns: ", paste(colnames(meta), collapse = ", "), "\n",
       "Hint: e.g. '~ batch + group' requires a 'batch' column in the metadata.",
       call. = FALSE)
}

# Contrasts: explicit --contrast (comma-separated) or all ordered group pairs.
if (!is.null(contrast_arg)) {
  contrast_specs <- trimws(strsplit(contrast_arg, ",", fixed = TRUE)[[1]])
  contrast_specs <- contrast_specs[nzchar(contrast_specs)]
  if (length(contrast_specs) == 0L) {
    stop("--contrast was provided but empty. Use e.g. --contrast salt_vs_control.",
         call. = FALSE)
  }
  contrasts <- lapply(contrast_specs, parse_contrast, group_col = group_col, levels = groups)
} else {
  contrasts <- list()
  for (i in seq_along(groups)) {
    for (j in seq_along(groups)) {
      if (i == j) next
      contrasts[[length(contrasts) + 1L]] <- list(
        name = paste0(groups[i], "_vs_", groups[j]),
        numerator = groups[i], denominator = groups[j]
      )
    }
  }
  message("No --contrast given; using all ordered group pairs: ",
          paste(vapply(contrasts, function(x) x$name, character(1)), collapse = ", "))
}

coldata <- meta
coldata[[group_col]] <- factor(coldata[[group_col]], levels = groups)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = coldata,
                                      design = design_formula)
dds <- DESeq2::DESeq(dds, parallel = FALSE)
message("DESeq() finished.")

# VST normalisation (blind = FALSE: uses the experimental design).
vsd <- tryCatch(
  DESeq2::vst(dds, blind = FALSE),
  error = function(e) {
    message("vst() failed (", conditionMessage(e), "); falling back to ",
            "varianceStabilizingTransformation().")
    DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  }
)
norm_mat <- SummarizedExperiment::assay(vsd)
norm_df <- data.frame(gene_id = rownames(norm_mat), as.data.frame(norm_mat),
                      check.names = FALSE, stringsAsFactors = FALSE)
utils::write.csv(norm_df, file.path(outdir, "normalized_counts_vst.csv"),
                 row.names = FALSE, quote = FALSE)
message("Wrote normalized_counts_vst.csv (", nrow(norm_df), " genes x ",
        ncol(norm_df) - 1L, " samples).")

# Diagnostic figure: dispersion estimates (final vs base dispersion).
# Scientific conclusion: gene-wise dispersion estimates should follow the
# fitted trend (red curve); outliers far above the trend indicate noisy genes
# that DESeq2 shrinks towards the curve before testing.
save_base_figure(
  function() DESeq2::plotDispEsts(dds),
  file.path(outdir, "dispersion_estimates"), width_mm = 89, height_mm = 70
)

summary_rows <- list()
for (ct in contrasts) {
  tag <- contrast_tag(ct$name)
  if (lfc_shrink) {
    res <- DESeq2::lfcShrink(dds, contrast = c(group_col, ct$numerator, ct$denominator),
                             type = "ashr")
  } else {
    res <- DESeq2::results(dds, contrast = c(group_col, ct$numerator, ct$denominator),
                           alpha = padj_cut)
  }
  res_df <- data.frame(gene_id = rownames(res), as.data.frame(res),
                       check.names = FALSE, stringsAsFactors = FALSE)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
  rownames(res_df) <- NULL
  utils::write.csv(res_df, file.path(outdir, paste0("de_all_", tag, ".csv")),
                   row.names = FALSE, quote = FALSE)

  is_sig <- !is.na(res_df$padj) & res_df$padj <= padj_cut &
    !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= lfc_cut
  sig_df <- res_df[is_sig, , drop = FALSE]
  sig_df$significant <- ifelse(sig_df$log2FoldChange > 0, "Up", "Down")
  utils::write.csv(sig_df, file.path(outdir, paste0("de_sig_", tag, ".csv")),
                   row.names = FALSE, quote = FALSE)

  n_up <- sum(sig_df$significant == "Up")
  n_down <- sum(sig_df$significant == "Down")
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    contrast = ct$name, up = n_up, down = n_down,
    ns = nrow(res_df) - n_up - n_down, stringsAsFactors = FALSE
  )
  message(sprintf("Contrast %s: %d up, %d down (padj <= %g, |log2FC| >= %g).",
                  ct$name, n_up, n_down, padj_cut, lfc_cut))
}

de_summary <- do.call(rbind, summary_rows)
utils::write.csv(de_summary, file.path(outdir, "de_summary.csv"),
                 row.names = FALSE, quote = FALSE)

write_run_info_r(
  outdir, "02_deseq2.R",
  params = list(
    counts = counts_path,
    metadata = metadata_path,
    group_col = group_col,
    design = design,
    contrasts = paste(vapply(contrasts, function(x) x$name, character(1)),
                      collapse = ","),
    padj = padj_cut,
    log2fc = lfc_cut,
    lfc_shrink = lfc_shrink
  ),
  packages = c("DESeq2")
)
message("Differential expression complete. Results in ", outdir)
