###############################################################################
# 02_diablo.R -- supervised multi-omics integration of rice transcriptome and
#                metabolome: DIABLO (block.splsda) and O2PLS-style PLS
#                (Rice-skills skill, module 3, step 02)
#
# Title:    DIABLO / two-block PLS integration of transcriptome and metabolome
# Description:
#   Two mixOmics-based strategies for integrating RNA-seq and metabolome data
#   measured on the same samples:
#
#   (A) DIABLO  (--method diablo, default): block.splsda fits a supervised
#       sparse PLS-DA model with one block per omics layer. The design matrix
#       defines how strongly the two blocks are coupled; correlated features
#       are selected across blocks to discriminate the sample groups. This is
#       the right choice for DISCRIMINATIVE questions ("which genes and
#       metabolites jointly separate the treatment groups?"). Model
#       performance is estimated with repeated M-fold cross-validation (perf).
#
#   (B) PLS    (--method o2pls): an unsupervised two-block PLS regression in
#       canonical mode, i.e. the mixOmics-native method for CORRELATION-based
#       N-integration ("which transcript and metabolite modules co-vary?").
#       DEVIATION NOTE: mixOmics does not export an opls() function, so a true
#       O2PLS model (explicit predictive + orthogonal components, as in
#       ropls::opls) cannot be fitted with the installed dependency set. The
#       canonical-mode PLS below is the closest equivalent available in
#       mixOmics: it extracts paired latent components that maximise the
#       cross-block covariance between the RNA and metabolome views. Output
#       file names keep the o2pls_* prefix requested by the pipeline spec.
#
# Usage:
#   Rscript scripts/integration/02_diablo.R \
#       --rna results/transcriptome/02_de/vst_matrix.csv \
#       --metabolites results/metabolome/01_preprocess/metabolites_preprocessed.csv \
#       --metadata metadata.csv \
#       --method diablo --ncomp 2 --keepx 10,10 \
#       --cv-folds 10 --nrepeat 10 --group-col group \
#       --outdir results/integration/02_diablo
#
# Inputs:
#   --rna          gene x sample matrix (VST/normalised values recommended)
#   --metabolites  metabolite x sample matrix (annotation columns ignored)
#   --metadata     sample metadata (first column = sample_id, group column required)
#
# Outputs (DIABLO path, written to --outdir):
#   diablo_perf.*                   cross-validation error-rate plot (PDF/SVG/TIFF)
#   diablo_perf_error_rates.csv     tidy error rates per block/component
#   diablo_scores.csv               per-block and consensus sample scores
#   diablo_loadings.csv             per-block component loadings
#   sample_projection.*             plotIndiv figure (consensus space, group colours)
#   circos_loadings.*               circosPlot of cross-block correlated loadings
#   diablo_selected_features.csv    variables selected per block and component
#   run_info.txt
#
# Outputs (PLS/O2PLS path):
#   o2pls_scores.csv, o2pls_loadings.csv, o2pls_projection.*, run_info.txt
#
# Notes:
#   - keepX guidance: start from the perf() cross-validation results and tune
#     (tune.block.splsda) the number of selected variables per block; too many
#     variables overfit, too few lose signal. The --keepx value is reused for
#     every component.
#   - The design matrix matrix(c(0,1,1,0)) couples the two omics blocks with
#     full strength; reduce the off-diagonal values (e.g. 0.1-0.5) if the two
#     blocks are only weakly related.
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

usage_02 <- paste(
  "Usage: Rscript 02_diablo.R --rna <matrix.csv> --metabolites <matrix.csv>",
  "         --metadata <metadata.csv> [--method diablo|o2pls] [--ncomp N]",
  "         [--keepx 10,10] [--cv-folds N] [--nrepeat N] [--group-col COL]",
  "         [--outdir DIR] [--seed N]",
  sep = "\n"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE))
rna_path <- opt(opts, "rna", NULL)
met_path <- opt(opts, "metabolites", NULL)
metadata_path <- opt(opts, "metadata", NULL)
if (is.null(rna_path) || is.null(met_path) || is.null(metadata_path)) {
  stop("Missing required argument(s): --rna, --metabolites and --metadata are required.\n",
       usage_02, call. = FALSE)
}
outdir <- opt(opts, "outdir", "results/integration/02_diablo")
method <- tolower(opt(opts, "method", "diablo"))
if (!method %in% c("diablo", "o2pls")) {
  stop("--method must be either 'diablo' or 'o2pls', got '", method, "'.", call. = FALSE)
}
ncomp <- as.integer(opt(opts, "ncomp", "2"))
keepx_raw <- opt(opts, "keepx", "10,10")
cv_folds <- as.integer(opt(opts, "cv-folds", "10"))
nrepeat <- as.integer(opt(opts, "nrepeat", "10"))
group_col <- opt(opts, "group-col", "group")
seed <- as.integer(opt(opts, "seed", "42"))
if (is.na(ncomp) || ncomp < 1L) stop("--ncomp must be a positive integer.", call. = FALSE)
if (is.na(cv_folds) || cv_folds < 2L) stop("--cv-folds must be >= 2.", call. = FALSE)
if (is.na(nrepeat) || nrepeat < 1L) stop("--nrepeat must be >= 1.", call. = FALSE)
keepx_vals <- as.integer(strsplit(keepx_raw, ",", fixed = TRUE)[[1]])
keepx_vals <- trimws(keepx_vals)
if (length(keepx_vals) != 2L || any(is.na(keepx_vals)) || any(keepx_vals < 1L)) {
  stop("--keepx must be two positive integers separated by a comma (one per omics ",
       "block), e.g. '10,10'. Got '", keepx_raw, "'.", call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rice_set_seed(seed)

need_pkgs("mixOmics")

# ---------------------------------------------------------------------------
# Shared local helpers (this script is standalone by design)
# ---------------------------------------------------------------------------

read_feature_matrix <- function(path, what) {
  df <- read_csv_safe(path, what)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) stop(what, " contains duplicate IDs in the first column.",
                               call. = FALSE)
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  non_numeric <- colnames(df)[!vapply(df, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(what, " has non-numeric sample column(s): ",
         paste(head(non_numeric, 5), collapse = ", "),
         "\nHint: only the first column may hold IDs; sample columns must be numeric.",
         call. = FALSE)
  }
  as.matrix(df)
}

read_metabolome_matrix <- function(path, what) {
  df <- read_csv_safe(path, what)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) stop(what, " contains duplicate metabolite_id values.", call. = FALSE)
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  is_num <- vapply(df, is.numeric, logical(1))
  if (sum(is_num) < 1) stop(what, " has no numeric sample columns.", call. = FALSE)
  as.matrix(df[, is_num, drop = FALSE])
}

read_sample_metadata <- function(path, what) {
  df <- read_csv_safe(path, what)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) stop(what, " contains duplicate sample_id values.", call. = FALSE)
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  df
}

align_common_samples <- function(mat, meta, mat_name) {
  shared <- intersect(rownames(meta), colnames(mat))
  if (length(shared) == 0) {
    stop("No overlapping samples between ", mat_name, " and metadata.\n",
         "matrix samples: ", paste(head(colnames(mat), 5), collapse = ", "), " ...\n",
         "metadata samples: ", paste(head(rownames(meta), 5), collapse = ", "), " ...\n",
         "Hint: sample IDs must match exactly (case-sensitive, no extra spaces).",
         call. = FALSE)
  }
  dropped <- setdiff(colnames(mat), shared)
  if (length(dropped) > 0) {
    message("Note: ", length(dropped), " sample(s) in ", mat_name,
            " are not present in metadata and were dropped: ",
            paste(head(dropped, 5), collapse = ", "), if (length(dropped) > 5) " ..." else "")
  }
  mat[, shared, drop = FALSE]
}

#' Save a figure that may be a ggplot2 object or a base-graphics draw function
#' (mixOmics circosPlot and plot.perf draw directly to the current device).
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

#' Extract cross-validation error rates from a mixOmics perf object into a
#' tidy data.frame (works whether error.rate is a list per block or a matrix).
tidy_perf_errors <- function(perf_obj) {
  er <- perf_obj$error.rate
  if (is.null(er)) return(NULL)
  if (is.list(er)) {
    rows <- lapply(names(er), function(nm) {
      m <- as.data.frame(er[[nm]])
      m$block <- nm
      m
    })
    out <- do.call(rbind, rows)
  } else {
    out <- as.data.frame(er)
  }
  out
}

# ---------------------------------------------------------------------------
# Load and align inputs
# ---------------------------------------------------------------------------

rna <- read_feature_matrix(rna_path, "RNA matrix")
met <- read_metabolome_matrix(met_path, "Metabolome matrix")
meta <- read_sample_metadata(metadata_path, "Sample metadata")

if (!group_col %in% colnames(meta)) {
  stop("Sample metadata is missing the group column '", group_col, "'.\n",
       "Available columns: ", paste(colnames(meta), collapse = ", "), "\n",
       "Hint: pass the correct column with --group-col.", call. = FALSE)
}

common <- Reduce(intersect, list(colnames(rna), colnames(met), rownames(meta)))
if (length(common) < 4) {
  stop("Only ", length(common), " sample(s) are shared by RNA / metabolome / metadata; ",
       "at least 4 are required for cross-validated multi-omics integration.\n",
       "RNA samples: ", paste(head(colnames(rna), 5), collapse = ", "), " ...\n",
       "metabolome samples: ", paste(head(colnames(met), 5), collapse = ", "), " ...\n",
       "metadata samples: ", paste(head(rownames(meta), 5), collapse = ", "), " ...\n",
       "Hint: sample IDs must match exactly (case-sensitive, no extra spaces).",
       call. = FALSE)
}
rna <- align_common_samples(rna, meta, "RNA matrix")
met <- align_common_samples(met, meta, "Metabolome matrix")
Y <- factor(meta[[group_col]][common, drop = TRUE])
if (nlevels(Y) < 2) {
  stop("The group column '", group_col, "' has a single level; DIABLO/PLS need at ",
       "least two groups (for the unsupervised PLS path a single level is allowed ",
       "but the group annotation is then uninformative).", call. = FALSE)
}
if (min(table(Y)) < 2) {
  stop("Every group needs at least 2 samples for cross-validation; group sizes: ",
       paste(names(table(Y)), table(Y), sep = "=", collapse = ", "), ".", call. = FALSE)
}

message("Aligned ", length(common), " samples; groups: ",
        paste(levels(Y), table(Y), sep = "=", collapse = ", "), ".")

# mixOmics expects blocks as a named list of samples x variables matrices.
blocks <- list(rna = t(rna), metabolome = t(met))
keepx_list <- list(rna = rep(keepx_vals[1], ncomp),
                   metabolome = rep(keepx_vals[2], ncomp))
for (b in names(blocks)) {
  if (any(keepx_list[[b]] >= ncol(blocks[[b]]))) {
    stop("keepX for block '", b, "' (", paste(keepx_list[[b]], collapse = ","),
         ") must be smaller than the number of features (", ncol(blocks[[b]]), ").",
         call. = FALSE)
  }
  if (any(!is.finite(blocks[[b]]))) {
    stop("Block '", b, "' contains non-finite values. Please impute or filter them ",
         "upstream (metabolome preprocessing step).", call. = FALSE)
  }
}

# Design matrix: full coupling between the two omics blocks (0 on the diagonal,
# 1 off-diagonal). Lower off-diagonal values weaken the cross-omics link.
design <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2,
                 dimnames = list(c("rna", "metabolome"), c("rna", "metabolome")))

# ---------------------------------------------------------------------------
# Path A: DIABLO (supervised, discriminative)
# ---------------------------------------------------------------------------

if (method == "diablo") {
  # --- performance assessment with repeated M-fold cross-validation ---------
  folds <- min(cv_folds, min(table(Y)))
  if (folds < cv_folds) {
    message("Note: --cv-folds reduced from ", cv_folds, " to ", folds,
            " (smallest group size).")
  }
  diablo <- mixOmics::block.splsda(X = blocks, Y = Y, ncomp = ncomp,
                                   keepX = keepx_list, design = design)
  message("Running repeated cross-validation (", folds, " folds x ", nrepeat,
          " repeats)...")
  perf_obj <- mixOmics::perf(diablo, validation = "Mfold", folds = folds,
                             nrepeat = nrepeat, progressBar = FALSE)

  save_plot_output(function() plot(perf_obj),
                   file.path(outdir, "diablo_perf"),
                   width_mm = 120, height_mm = 90)
  err_df <- tidy_perf_errors(perf_obj)
  if (!is.null(err_df)) {
    write.csv(err_df, file.path(outdir, "diablo_perf_error_rates.csv"), row.names = FALSE)
  }

  # --- joint sample scores: per-block variates + consensus -----------------
  # The consensus component of DIABLO is the average of the block-specific
  # latent components (same definition as plotIndiv(blocks = "consensus")).
  variates <- diablo$variates
  block_var <- variates[intersect(names(variates), names(blocks))]
  if (length(block_var) < 2) {
    stop("Unexpected DIABLO variates structure: ", paste(names(variates), collapse = ", "),
         call. = FALSE)
  }
  consensus <- Reduce(`+`, block_var) / length(block_var)

  score_frames <- lapply(c(names(block_var), "consensus"), function(b) {
    m <- if (identical(b, "consensus")) consensus else as.matrix(block_var[[b]])
    colnames(m) <- paste0("comp", seq_len(ncol(m)))
    data.frame(sample_id = rownames(m), block = b, m, check.names = FALSE,
               stringsAsFactors = FALSE)
  })
  scores_long <- do.call(rbind, score_frames)
  scores_long <- merge(scores_long,
                       data.frame(sample_id = common, group = as.character(Y)),
                       by = "sample_id", all.x = TRUE, sort = FALSE)
  comp_cols <- grep("^comp[0-9]+$", colnames(scores_long), value = TRUE)
  scores_long <- scores_long[, c("sample_id", "group", "block", comp_cols)]
  write.csv(scores_long, file.path(outdir, "diablo_scores.csv"), row.names = FALSE)

  # --- component loadings per block ----------------------------------------
  loading_rows <- lapply(names(diablo$loadings), function(b) {
    m <- as.matrix(diablo$loadings[[b]])
    colnames(m) <- paste0("comp", seq_len(ncol(m)))
    data.frame(block = b, feature = rownames(m), m, check.names = FALSE,
               stringsAsFactors = FALSE)
  })
  write.csv(do.call(rbind, loading_rows), file.path(outdir, "diablo_loadings.csv"),
            row.names = FALSE)

  # --- sample projection plot (consensus space, coloured by group) ----------
  save_plot_output(
    function() {
      p <- tryCatch(
        mixOmics::plotIndiv(diablo, blocks = "consensus", comp = c(1, 2), group = Y,
                            legend = TRUE, title = "DIABLO consensus projection"),
        error = function(e) {
          message("Note: consensus plotIndiv failed (", conditionMessage(e),
                  "); falling back to the default per-block projection.")
          mixOmics::plotIndiv(diablo, group = Y, legend = TRUE,
                              title = "DIABLO sample projection")
        }
      )
      p
    },
    file.path(outdir, "sample_projection"), width_mm = 120, height_mm = 100
  )

  # --- circos plot of cross-omics correlated loadings ----------------------
  save_plot_output(
    function() mixOmics::circosPlot(
      diablo, cutoff = 0.7, line = TRUE,
      color.blocks = c(rice_palette()$series12[1], rice_palette()$series12[2]),
      color.cor = c(nature_semantic("down"), nature_semantic("up"))
    ),
    file.path(outdir, "circos_loadings"), width_mm = 120, height_mm = 120
  )

  # --- selected features per block/component --------------------------------
  sel <- tryCatch(mixOmics::selectVar(diablo, comp = 1:ncomp), error = function(e) NULL)
  if (!is.null(sel)) {
    sel_rows <- lapply(names(sel), function(b) {
      df <- as.data.frame(sel[[b]])
      df$block <- b
      df
    })
    write.csv(do.call(rbind, sel_rows), file.path(outdir, "diablo_selected_features.csv"),
              row.names = FALSE)
  }

  write_run_info_r(
    outdir = outdir, script_name = "02_diablo.R",
    params = list(method = "diablo", rna = rna_path, metabolites = met_path,
                  metadata = metadata_path, ncomp = ncomp, keepx = keepx_raw,
                  cv_folds = folds, nrepeat = nrepeat, group_col = group_col,
                  seed = seed, design = "c(0,1,1,0)",
                  n_samples = length(common), n_genes = ncol(blocks$rna),
                  n_metabolites = ncol(blocks$metabolome)),
    packages = c("mixOmics", "ggplot2"),
    extra = c(
      "DIABLO is the supervised (discriminative) integration path.",
      "Consensus scores are the average of the block-specific variates.",
      "Tune keepX with tune.block.splsda based on the perf() results."
    )
  )
  message("Done (DIABLO). Outputs written to: ", normalizePath(outdir))
}

# ---------------------------------------------------------------------------
# Path B: two-block PLS (canonical mode) -- correlation-based integration
# ---------------------------------------------------------------------------

if (method == "o2pls") {
  o2 <- mixOmics::pls(blocks$rna, blocks$metabolome, ncomp = ncomp,
                      mode = "canonical", scale = TRUE)

  # Joint scores: X (transcriptome) and Y (metabolome) variates per component.
  x_scores <- as.matrix(o2$variates$X)
  y_scores <- as.matrix(o2$variates$Y)
  colnames(x_scores) <- paste0("rna_comp", seq_len(ncol(x_scores)))
  colnames(y_scores) <- paste0("metabolome_comp", seq_len(ncol(y_scores)))
  score_rows <- list(
    data.frame(sample_id = rownames(x_scores), block = "rna", x_scores,
               check.names = FALSE, stringsAsFactors = FALSE),
    data.frame(sample_id = rownames(y_scores), block = "metabolome", y_scores,
               check.names = FALSE, stringsAsFactors = FALSE)
  )
  scores_long <- do.call(rbind, score_rows)
  scores_long <- merge(scores_long,
                       data.frame(sample_id = common, group = as.character(Y)),
                       by = "sample_id", all.x = TRUE, sort = FALSE)
  write.csv(scores_long, file.path(outdir, "o2pls_scores.csv"), row.names = FALSE)

  # Loadings for both blocks.
  loading_rows <- lapply(c("X", "Y"), function(b) {
    m <- as.matrix(o2$loadings[[b]])
    colnames(m) <- paste0("comp", seq_len(ncol(m)))
    data.frame(block = if (b == "X") "rna" else "metabolite",
               feature = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  })
  write.csv(do.call(rbind, loading_rows), file.path(outdir, "o2pls_loadings.csv"),
            row.names = FALSE)

  # Sample projection in the averaged X-Y space, coloured by group.
  save_plot_output(
    function() mixOmics::plotIndiv(o2, rep.space = "XY-variate", group = Y,
                                   legend = TRUE,
                                   title = "Two-block PLS (canonical) projection"),
    file.path(outdir, "o2pls_projection"), width_mm = 120, height_mm = 100
  )

  write_run_info_r(
    outdir = outdir, script_name = "02_diablo.R",
    params = list(method = "o2pls", rna = rna_path, metabolites = met_path,
                  metadata = metadata_path, ncomp = ncomp, mode = "canonical",
                  scale = TRUE, group_col = group_col, seed = seed,
                  n_samples = length(common), n_genes = ncol(blocks$rna),
                  n_metabolites = ncol(blocks$metabolome)),
    packages = c("mixOmics", "ggplot2"),
    extra = c(
      "DEVIATION: mixOmics exports no opls() function (true O2PLS lives in ropls).",
      "This path uses mixOmics::pls(mode = 'canonical'), the mixOmics-native",
      "two-block method for correlation-based N-integration.",
      "Unsupervised: group labels only annotate the projection plot."
    )
  )
  message("Done (two-block PLS). Outputs written to: ", normalizePath(outdir))
}
