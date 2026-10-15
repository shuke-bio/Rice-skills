###############################################################################
# 01_mofa.R -- Multi-Omics Factor Analysis (MOFA) for rice transcriptome x
#              metabolome integration (Rice-skills skill, module 3, step 01)
#
# Title:    MOFA2-based unsupervised integration of transcriptome and metabolome
# Description:
#   Fits a MOFA (Multi-Omics Factor Analysis) model on two views measured on
#   the same samples: an RNA-seq matrix (VST or counts) and a preprocessed
#   metabolome abundance matrix. MOFA decomposes the aligned multi-omics data
#   into a small number of latent factors; each factor captures one source of
#   variation that may be shared across the two omics layers or specific to one
#   of them, together with the feature weights (loadings) that define it.
#
# Usage:
#   Rscript scripts/integration/01_mofa.R \
#       --rna results/transcriptome/02_de/vst_matrix.csv \
#       --metabolites results/metabolome/01_preprocess/metabolites_preprocessed.csv \
#       --metadata metadata.csv \
#       --outdir results/integration/01_mofa \
#       --n-factors 10 --group-col group --seed 42
#
# Inputs:
#   --rna          gene x sample matrix (VST or normalised counts; header row,
#                  first column = gene_id, remaining columns = sample IDs)
#   --metabolites  metabolite x sample matrix (annotation columns such as
#                  metabolite_name/class/kegg_id are allowed and ignored for
#                  the model; first column = metabolite_id)
#   --metadata     sample metadata (first column = sample_id, must contain the
#                  --group-col column; group labels are only used to annotate
#                  figures, MOFA itself is unsupervised)
#
# Outputs (written to --outdir):
#   mofa_model.hdf5                trained MOFA2 model (HDF5 format)
#   variance_explained.csv         variance explained per factor per view
#   variance_explained.*           plot_variance_explained figure (PDF/SVG/TIFF)
#   factor_scores.csv              sample x factor scores with group annotation
#   factor_scores_heatmap.*        plot_factors figure (factors 1-2, group-coloured)
#   factor_loadings_rna.csv        top 25 weights per factor (RNA view)
#   factor_loadings_metabolites.csv top 25 weights per factor (metabolome view)
#   top_features_per_factor.csv    combined top features across both views
#   run_info.txt                   parameters, package versions, seed
#
# How to choose the number of factors:
#   MOFA is an unsupervised factor model: adding factors always explains more
#   variance, so the variance-explained curve alone cannot pick the number.
#   Recommended workflow:
#     (1) run with --n-factors 10-15 as an upper bound,
#     (2) inspect variance_explained.* and factor_scores_heatmap.*,
#     (3) keep the factors that explain a non-trivial fraction of variance and
#         are interpretable (they separate the groups / correlate with traits),
#     (4) re-run with the selected --n-factors for the final model.
#   Factors that capture no signal are pruned automatically by the ARD prior
#   (their weights shrink to zero), so over-estimating --n-factors is safe.
#
# Data standardisation note:
#   MOFA2 models each view with a Gaussian likelihood and centres each view
#   internally. RNA-seq VST values and metabolite abundances live on very
#   different scales, so this script sets scale_views = TRUE (each view is
#   scaled to unit variance before training). Raw integer counts are NOT a
#   valid input -- pass VST/normalised values (or use a Poisson likelihood,
#   not implemented here). Feature selection (highly variable features) is
#   expected to happen upstream; MOFA2 itself handles centering and scaling
#   once scale_views is set.
#
# Notes:
#   - The MOFA2 R package requires the Python backend mofapy2 (reticulate).
#     Install with pip ("pip install mofapy2") or let MOFA2 create a conda
#     environment via basilisk (run_mofa(..., use_basilisk = TRUE)).
#   - MOFA is exploratory and unsupervised: it does not use the group labels;
#     groups only colour the figures.
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

usage_01 <- paste(
  "Usage: Rscript 01_mofa.R --rna <matrix.csv> --metabolites <matrix.csv>",
  "         --metadata <metadata.csv> [--outdir DIR] [--n-factors N]",
  "         [--group-col COL] [--seed N]",
  sep = "\n"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE))
rna_path <- opt(opts, "rna", NULL)
met_path <- opt(opts, "metabolites", NULL)
metadata_path <- opt(opts, "metadata", NULL)
if (is.null(rna_path) || is.null(met_path) || is.null(metadata_path)) {
  stop("Missing required argument(s): --rna, --metabolites and --metadata are required.\n",
       usage_01, call. = FALSE)
}
outdir <- opt(opts, "outdir", "results/integration/01_mofa")
n_factors <- as.integer(opt(opts, "n-factors", "10"))
group_col <- opt(opts, "group-col", "group")
seed <- as.integer(opt(opts, "seed", "42"))
if (is.na(n_factors) || n_factors < 1L) {
  stop("--n-factors must be a positive integer (recommended: 10-15 for tuning).", call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rice_set_seed(seed)

need_pkgs(c("MOFA2", "HDF5Array"))

# ---------------------------------------------------------------------------
# Shared local helpers (this script is standalone by design)
# ---------------------------------------------------------------------------

#' Read a feature x sample matrix (first column = feature ID, rest = samples).
read_feature_matrix <- function(path, what, id_col) {
  df <- read_csv_safe(path, what)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) {
    stop(what, " contains duplicate IDs in the first column (e.g. ",
         paste(head(unique(ids[duplicated(ids)]), 3), collapse = ", "), ").",
         call. = FALSE)
  }
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  non_numeric <- colnames(df)[!vapply(df, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(what, " has non-numeric sample column(s): ",
         paste(head(non_numeric, 5), collapse = ", "),
         "\nHint: only the first column may hold IDs; all sample columns must be numeric.",
         call. = FALSE)
  }
  as.matrix(df)
}

#' Read a metabolite x sample matrix, splitting annotation columns off.
read_metabolome_matrix <- function(path, what) {
  df <- read_csv_safe(path, what)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) {
    stop(what, " contains duplicate metabolite_id values.", call. = FALSE)
  }
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  is_num <- vapply(df, is.numeric, logical(1))
  annotation <- df[, !is_num, drop = FALSE]   # metabolite_name / class / kegg_id ...
  mat <- df[, is_num, drop = FALSE]
  if (ncol(mat) < 1) {
    stop(what, " has no numeric sample columns.", call. = FALSE)
  }
  list(matrix = as.matrix(mat), annotation = annotation)
}

#' Read sample metadata (first column = sample_id).
read_sample_metadata <- function(path, what) {
  df <- read_csv_safe(path, what)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) stop(what, " contains duplicate sample_id values.", call. = FALSE)
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  df
}

#' Align matrix sample columns with metadata rows (friendly errors).
align_common_samples <- function(mat, meta, mat_name) {
  mat_samples <- colnames(mat)
  meta_samples <- rownames(meta)
  shared <- intersect(meta_samples, mat_samples)
  if (length(shared) == 0) {
    stop("No overlapping samples between ", mat_name, " and metadata.\n",
         "matrix samples: ", paste(head(mat_samples, 5), collapse = ", "), " ...\n",
         "metadata samples: ", paste(head(meta_samples, 5), collapse = ", "), " ...\n",
         "Hint: sample IDs must match exactly (case-sensitive, no extra spaces).",
         call. = FALSE)
  }
  dropped <- setdiff(mat_samples, shared)
  if (length(dropped) > 0) {
    message("Note: ", length(dropped), " sample(s) in ", mat_name,
            " are not present in metadata and were dropped: ",
            paste(head(dropped, 5), collapse = ", "), if (length(dropped) > 5) " ..." else "")
  }
  mat[, shared, drop = FALSE]
}

#' Save a figure that may be a ggplot2 object or a base-graphics draw function.
save_plot_output <- function(x, file_stem, width_mm = 89, height_mm = 89) {
  if (inherits(x, "ggplot")) {
    return(save_nature_figure(x, file_stem, width_mm = width_mm, height_mm = height_mm))
  }
  if (!is.function(x)) {
    stop("save_plot_output() expects a ggplot object or a function that draws a plot.",
         call. = FALSE)
  }
  w_in <- mm_to_inch(width_mm)
  h_in <- mm_to_inch(height_mm)
  for (f in c("pdf", "svg", "tiff")) {
    out <- paste0(file_stem, ".", if (f == "tiff") "tiff" else f)
    if (f == "pdf") {
      grDevices::pdf(out, width = w_in, height = h_in)
    } else if (f == "svg") {
      grDevices::svg(out, width = w_in, height = h_in)
    } else {
      grDevices::tiff(out, width = w_in, height = h_in, units = "in", res = 300,
                      compression = "lzw")
    }
    res <- x()
    if (inherits(res, "ggplot")) print(res)
    grDevices::dev.off()
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Load and align inputs
# ---------------------------------------------------------------------------

rna <- read_feature_matrix(rna_path, "RNA matrix", "gene_id")
met <- read_metabolome_matrix(met_path, "Metabolome matrix")
meta <- read_sample_metadata(metadata_path, "Sample metadata")

if (!group_col %in% colnames(meta)) {
  stop("Sample metadata is missing the group column '", group_col, "'.\n",
       "Available columns: ", paste(colnames(meta), collapse = ", "), "\n",
       "Hint: pass the correct column with --group-col.", call. = FALSE)
}
groups <- factor(meta[[group_col]])
if (nlevels(groups) < 2) {
  message("Note: only one group level ('", levels(groups)[1],
          "'); MOFA is unsupervised so figures will use a single colour.")
}

common <- Reduce(intersect, list(colnames(rna), colnames(met$matrix), rownames(meta)))
if (length(common) < 3) {
  stop("Only ", length(common), " sample(s) are shared by RNA / metabolome / metadata; ",
       "at least 3 are required (MOFA documentation suggests >= 15 samples).\n",
       "RNA samples: ", paste(head(colnames(rna), 5), collapse = ", "), " ...\n",
       "metabolome samples: ", paste(head(colnames(met$matrix), 5), collapse = ", "), " ...\n",
       "metadata samples: ", paste(head(rownames(meta), 5), collapse = ", "), " ...\n",
       "Hint: sample IDs must match exactly (case-sensitive, no extra spaces).",
       call. = FALSE)
}
rna <- align_common_samples(rna, meta, "RNA matrix")
met_mat <- align_common_samples(met$matrix, meta, "Metabolome matrix")
group_vec <- setNames(as.character(groups[common, drop = TRUE]), common)

message("Aligned ", ncol(rna), " samples across RNA (", nrow(rna), " genes) and ",
        "metabolome (", nrow(met_mat), " metabolites).")

# MOFA2 expects a list of matrices: features in rows, samples in columns.
mofa_data <- list(rna = rna, metabolome = met_mat)
if (any(vapply(mofa_data, function(m) any(!is.finite(m)), logical(1)))) {
  stop("Input matrices contain non-finite values (NA/Inf/NaN). ",
       "Please filter or impute them upstream (e.g. via the metabolome ",
       "preprocessing step) before running MOFA.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Build, prepare and train the MOFA model
# ---------------------------------------------------------------------------

mofa <- MOFA2::create_mofa(mofa_data)

data_opts <- MOFA2::get_default_data_options(mofa)
data_opts$scale_views <- TRUE   # RNA and metabolite scales differ -> scale each view
model_opts <- MOFA2::get_default_model_options(mofa)
model_opts$num_factors <- n_factors
train_opts <- MOFA2::get_default_training_options(mofa)
train_opts$seed <- seed

mofa <- MOFA2::prepare_mofa(mofa,
                            data_options = data_opts,
                            model_options = model_opts,
                            training_options = train_opts)

model_file <- file.path(outdir, "mofa_model.hdf5")
message("Training MOFA model with ", n_factors, " factors (this calls the Python mofapy2 backend)...")
mofa <- MOFA2::run_mofa(mofa, outfile = model_file)

# ---------------------------------------------------------------------------
# Variance explained per view (plot_variance_explained + tidy CSV)
# ---------------------------------------------------------------------------

save_plot_output(function() MOFA2::plot_variance_explained(mofa),
                 file.path(outdir, "variance_explained"),
                 width_mm = 120, height_mm = 90)

ve <- tryCatch(MOFA2::calculate_variance_explained(mofa), error = function(e) NULL)
if (!is.null(ve)) {
  ve_rows <- lapply(names(ve), function(v) {
    m <- as.data.frame(ve[[v]])
    m$view <- v
    m$factor <- rownames(ve[[v]])
    m
  })
  ve_df <- do.call(rbind, ve_rows)
  write.csv(ve_df, file.path(outdir, "variance_explained.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Factor scores (samples x factors) + plot_factors coloured by group
# ---------------------------------------------------------------------------

factor_scores <- MOFA2::get_factors(mofa)
if (!is.matrix(factor_scores)) factor_scores <- as.matrix(factor_scores)
factor_scores <- factor_scores[common, , drop = FALSE]
colnames(factor_scores) <- paste0("Factor", seq_len(ncol(factor_scores)))

scores_df <- data.frame(sample_id = common, group = group_vec[common],
                        factor_scores, check.names = FALSE, stringsAsFactors = FALSE)
write.csv(scores_df, file.path(outdir, "factor_scores.csv"), row.names = FALSE)

# plot_factors: colour by the group vector (option 3 of color_by); if MOFA2's
# plotting fails for any reason, fall back to a ggplot2 tile heatmap of the
# factor scores annotated by group.
fs_plot <- tryCatch(
  MOFA2::plot_factors(mofa, factors = c(1, 2), color_by = group_vec[common],
                      color_name = group_col),
  error = function(e) {
    message("Note: MOFA2::plot_factors failed (", conditionMessage(e),
            "); writing a ggplot2 factor-score heatmap instead.")
    NULL
  }
)
if (!is.null(fs_plot)) {
  save_plot_output(fs_plot, file.path(outdir, "factor_scores_heatmap"),
                   width_mm = 89, height_mm = 89)
} else {
  fs_long <- data.frame(
    sample = rep(common, times = ncol(factor_scores)),
    factor = rep(colnames(factor_scores), each = length(common)),
    value = as.vector(factor_scores),
    stringsAsFactors = FALSE
  )
  fs_long$sample <- factor(fs_long$sample, levels = common[order(group_vec[common])])
  fs_heat <- ggplot2::ggplot(fs_long, ggplot2::aes(x = .data$factor, y = .data$sample,
                                                   fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = diverging_palette()) +
    ggplot2::labs(x = "Latent factor", y = paste0("Sample (ordered by ", group_col, ")"),
                  fill = "Factor score",
                  title = "MOFA factor scores",
                  subtitle = "Rows are samples ordered by group; each column is one latent factor") +
    theme_nature()
  save_nature_figure(fs_heat, file.path(outdir, "factor_scores_heatmap"),
                     width_mm = 180, height_mm = max(60, 6 + 2.4 * length(common)))
}

# ---------------------------------------------------------------------------
# Feature weights (loadings): top 25 per factor for each view
# ---------------------------------------------------------------------------

weights <- MOFA2::get_weights(mofa)
view_files <- list()
top_all <- list()
for (view in names(weights)) {
  w <- as.matrix(weights[[view]])
  if (ncol(w) == 0) next
  colnames(w) <- paste0("Factor", seq_len(ncol(w)))
  keep_top <- min(25L, nrow(w))
  top_list <- lapply(seq_len(ncol(w)), function(k) {
    ord <- order(-abs(w[, k]))[seq_len(keep_top)]
    data.frame(factor = colnames(w)[k], rank = seq_len(keep_top),
               feature = rownames(w)[k][ord], weight = w[ord, k],
               stringsAsFactors = FALSE)
  })
  top_df <- do.call(rbind, top_list)
  fname <- if (view == "rna") "factor_loadings_rna.csv" else
           if (view == "metabolome") "factor_loadings_metabolites.csv" else
           paste0("factor_loadings_", view, ".csv")
  write.csv(top_df, file.path(outdir, fname), row.names = FALSE)
  view_files[[view]] <- fname
  top_df$view <- view
  top_all[[view]] <- top_df
}
if (length(top_all) > 0) {
  top_features <- do.call(rbind, top_all)
  top_features <- top_features[, c("view", "factor", "rank", "feature", "weight")]
  write.csv(top_features, file.path(outdir, "top_features_per_factor.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

write_run_info_r(
  outdir = outdir,
  script_name = "01_mofa.R",
  params = list(
    rna = rna_path, metabolites = met_path, metadata = metadata_path,
    n_factors = n_factors, group_col = group_col, seed = seed,
    n_samples = length(common), n_genes = nrow(rna), n_metabolites = nrow(met_mat),
    scale_views = TRUE
  ),
  packages = c("MOFA2", "HDF5Array", "ggplot2"),
  extra = c(
    "MOFA is unsupervised; group labels only annotate figures.",
    "Recommended workflow: run with 10-15 factors, inspect the variance-explained",
    "curve and factor scores, then re-run with the selected number of factors.",
    "MOFA2 requires the Python backend mofapy2 (pip install mofapy2)."
  )
)

message("Done. Outputs written to: ", normalizePath(outdir))
