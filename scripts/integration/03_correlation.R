###############################################################################
# 03_correlation.R -- gene x metabolite correlation analysis between the
#                     differential transcriptome and differential metabolome
#                     (Rice-skills skill, module 3, step 03)
#
# Title:    Differential gene x differential metabolite correlation network
# Description:
#   Computes pairwise correlations between the differentially expressed genes
#   (from the transcriptome DE step) and the differential metabolites (from
#   the metabolome diff step), using the aligned VST expression matrix and the
#   preprocessed metabolite abundance matrix. Correlations are computed in
#   blocks of metabolites to keep memory usage bounded for large matrices
#   (thousands of genes x thousands of metabolites), p-values are obtained
#   with the WGCNA Student approximation and adjusted with Benjamini-Hochberg.
#   Significant pairs (|r| >= --cor-cutoff and padj <= --padj) define a
#   gene-metabolite bipartite network that is consumed by the downstream
#   network (05) and figure (06) steps.
#
#   IMPORTANT (biological caveat): a strong gene-metabolite correlation is an
#   association, NOT evidence of causation. Interpret significant pairs in the
#   context of the joint pathway analysis (04) and, where possible, validate
#   with targeted assays.
#
# Usage:
#   Rscript scripts/integration/03_correlation.R \
#       --de-genes results/transcriptome/02_de/de_sig_control_vs_salt.csv \
#       --de-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
#       --rna results/transcriptome/02_de/vst_matrix.csv \
#       --metabolites results/metabolome/01_preprocess/metabolites_preprocessed.csv \
#       --outdir results/integration/03_correlation \
#       --method spearman --cor-cutoff 0.8 --padj 0.05 --block 500
#
# Inputs:
#   --de-genes       DE gene table; required columns: gene_id, log2FoldChange
#   --de-metabolites differential metabolite table; required column: metabolite_id
#   --rna            gene x sample expression matrix (VST/counts)
#   --metabolites    metabolite x sample abundance matrix (annotation columns
#                    are ignored; missing values allowed)
#
# Outputs (written to --outdir):
#   gene_metabolite_pairs_sig.csv   significant pairs (gene_id, metabolite_id,
#                                   r, p, padj, direction)
#   correlation_heatmap.*           heatmap of significant correlations
#   correlation_network_edges.csv   bipartite network edges (for steps 05/06)
#   correlation_network_nodes.csv   bipartite network nodes with degree
#   run_info.txt
#
# Notes:
#   - p-values for Spearman correlations are approximated with the Student
#     t transformation (WGCNA::corPvalueStudent), the standard practice in
#     co-expression network analysis; for a small number of pairs an exact
#     permutation test would be more faithful.
#   - Missing metabolite values: samples with any missing value inside a
#     metabolite block are excluded for that block (listwise deletion), so
#     the effective sample size n can differ slightly between blocks.
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

usage_03 <- paste(
  "Usage: Rscript 03_correlation.R --de-genes <de_sig.csv> --de-metabolites <diff.csv>",
  "         --rna <vst.csv> --metabolites <matrix.csv> [--outdir DIR]",
  "         [--method spearman|pearson] [--cor-cutoff 0.8] [--padj 0.05] [--block 500]",
  sep = "\n"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE))
de_genes_path <- opt(opts, "de-genes", NULL)
de_mets_path <- opt(opts, "de-metabolites", NULL)
rna_path <- opt(opts, "rna", NULL)
met_path <- opt(opts, "metabolites", NULL)
if (is.null(de_genes_path) || is.null(de_mets_path) || is.null(rna_path) || is.null(met_path)) {
  stop("Missing required argument(s): --de-genes, --de-metabolites, --rna and ",
       "--metabolites are all required.\n", usage_03, call. = FALSE)
}
outdir <- opt(opts, "outdir", "results/integration/03_correlation")
method <- tolower(opt(opts, "method", "spearman"))
cor_cutoff <- as.numeric(opt(opts, "cor-cutoff", "0.8"))
padj_cutoff <- as.numeric(opt(opts, "padj", "0.05"))
block_size <- as.integer(opt(opts, "block", "500"))
if (!method %in% c("spearman", "pearson")) {
  stop("--method must be 'spearman' or 'pearson', got '", method, "'.", call. = FALSE)
}
if (is.na(cor_cutoff) || cor_cutoff <= 0 || cor_cutoff > 1) {
  stop("--cor-cutoff must be in (0, 1], e.g. 0.8.", call. = FALSE)
}
if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
  stop("--padj must be in (0, 1], e.g. 0.05.", call. = FALSE)
}
if (is.na(block_size) || block_size < 1L) stop("--block must be a positive integer.", call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rice_set_seed(42)

need_pkgs("WGCNA")

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
  if (anyDuplicated(ids)) stop(what, " contains duplicate metabolite_id values.",
                               call. = FALSE)
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  is_num <- vapply(df, is.numeric, logical(1))
  if (sum(is_num) < 1) stop(what, " has no numeric sample columns.", call. = FALSE)
  as.matrix(df[, is_num, drop = FALSE])
}

#' Save a pheatmap object (class "pheatmap") in PDF + TIFF + PNG (300 dpi).
save_pheatmap <- function(p, file_stem, width_mm, height_mm) {
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("The grid package is required to re-draw pheatmap objects.", call. = FALSE)
  }
  w_in <- mm_to_inch(width_mm)
  h_in <- mm_to_inch(height_mm)
  for (f in c("pdf", "tiff", "png")) {
    out <- paste0(file_stem, ".", if (f == "tiff") "tiff" else f)
    if (f == "pdf") {
      grDevices::pdf(out, width = w_in, height = h_in)
    } else if (f == "png") {
      grDevices::png(out, width = w_in, height = h_in, units = "in", res = 300)
    } else {
      grDevices::tiff(out, width = w_in, height = h_in, units = "in", res = 300,
                      compression = "lzw")
    }
    grid::grid.draw(p$gtable)
    grDevices::dev.off()
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Load inputs and restrict to the differential features
# ---------------------------------------------------------------------------

de_genes <- read_csv_safe(de_genes_path, "DE genes table")
check_columns(de_genes, c("gene_id", "log2FoldChange"), "DE genes table")
de_genes$gene_id <- as.character(de_genes$gene_id)

de_mets <- read_csv_safe(de_mets_path, "Differential metabolites table")
check_columns(de_mets, "metabolite_id", "Differential metabolites table")
de_mets$metabolite_id <- as.character(de_mets$metabolite_id)

rna <- read_feature_matrix(rna_path, "RNA matrix")
met <- read_metabolome_matrix(met_path, "Metabolome matrix")

common <- intersect(colnames(rna), colnames(met))
if (length(common) < 3) {
  stop("Only ", length(common), " sample(s) are shared by the RNA and metabolome ",
       "matrices; at least 3 are required for correlation analysis.\n",
       "RNA samples: ", paste(head(colnames(rna), 5), collapse = ", "), " ...\n",
       "metabolome samples: ", paste(head(colnames(met), 5), collapse = ", "), " ...\n",
       "Hint: sample IDs must match exactly (case-sensitive, no extra spaces).",
       call. = FALSE)
}
if (length(common) < min(ncol(rna), ncol(met))) {
  message("Note: using the ", length(common), " samples shared by both matrices.")
}

genes <- intersect(unique(de_genes$gene_id), rownames(rna))
mets <- intersect(unique(de_mets$metabolite_id), rownames(met))
if (length(genes) == 0) {
  stop("None of the DE genes are present in the RNA matrix.\n",
       "DE genes: ", paste(head(de_genes$gene_id, 3), collapse = ", "), " ...\n",
       "matrix rows: ", paste(head(rownames(rna), 3), collapse = ", "), " ...\n",
       "Hint: gene IDs must match exactly (case-sensitive).", call. = FALSE)
}
if (length(mets) == 0) {
  stop("None of the differential metabolites are present in the metabolome matrix.\n",
       "Hint: metabolite IDs must match exactly (case-sensitive).", call. = FALSE)
}

rna_sub <- rna[genes, common, drop = FALSE]
met_sub <- met[mets, common, drop = FALSE]

# Drop genes with missing expression values (VST matrices should not have any).
na_genes <- rownames(rna_sub)[rowSums(!is.finite(rna_sub)) > 0]
if (length(na_genes) > 0) {
  message("Note: dropping ", length(na_genes), " gene(s) with missing values from the RNA matrix.")
  rna_sub <- rna_sub[setdiff(rownames(rna_sub), na_genes), , drop = FALSE]
  genes <- rownames(rna_sub)
}
# Drop metabolites that cannot be correlated (too few observations).
n_obs <- rowSums(is.finite(met_sub))
bad_mets <- names(n_obs)[n_obs < 3]
if (length(bad_mets) > 0) {
  message("Note: dropping ", length(bad_mets),
          " metabolite(s) with fewer than 3 observed values.")
  met_sub <- met_sub[setdiff(rownames(met_sub), bad_mets), , drop = FALSE]
  mets <- rownames(met_sub)
}
if (nrow(rna_sub) == 0 || nrow(met_sub) == 0) {
  stop("No DE gene / differential metabolite left after filtering; cannot compute correlations.",
       call. = FALSE)
}

message("Correlating ", nrow(rna_sub), " DE genes x ", nrow(met_sub),
        " differential metabolites across ", length(common), " samples ",
        "(method = ", method, ", block = ", block_size, " metabolites).")

# ---------------------------------------------------------------------------
# Blocked correlation: genes (rows) x metabolites (columns) per block
# ---------------------------------------------------------------------------

pair_list <- list()
n_blocks <- ceiling(nrow(met_sub) / block_size)
for (start in seq(1, nrow(met_sub), by = block_size)) {
  idx <- start:min(start + block_size - 1, nrow(met_sub))
  met_b <- met_sub[idx, , drop = FALSE]

  # Listwise deletion within the block: samples complete for all block
  # metabolites (RNA has no missing values at this point).
  ok <- colSums(!is.finite(met_b)) == 0
  samples_b <- common[ok]
  if (length(samples_b) < 3) {
    message("Note: block ", start, "-", max(idx), " skipped (only ",
            length(samples_b), " complete samples).")
    next
  }
  rna_b <- rna_sub[, samples_b, drop = FALSE]
  met_b <- met_b[, samples_b, drop = FALSE]

  r_mat <- stats::cor(t(rna_b), t(met_b), method = method)
  # Spearman p-values via the Student t approximation (WGCNA convention).
  p_mat <- WGCNA::corPvalueStudent(cor = r_mat, nSamples = length(samples_b))

  pair_list[[length(pair_list) + 1L]] <- data.frame(
    gene_id = rep(rownames(r_mat), times = ncol(r_mat)),
    metabolite_id = rep(colnames(r_mat), each = nrow(r_mat)),
    r = as.vector(r_mat),
    p = as.vector(p_mat),
    stringsAsFactors = FALSE
  )
}

if (length(pair_list) == 0) {
  stop("No correlations could be computed (all metabolite blocks lacked complete samples).",
       call. = FALSE)
}
pairs <- do.call(rbind, pair_list)
pairs <- pairs[is.finite(pairs$r) & is.finite(pairs$p), , drop = FALSE]
pairs$padj <- stats::p.adjust(pairs$p, method = "BH")
pairs$significant <- abs(pairs$r) >= cor_cutoff & pairs$padj <= padj_cutoff
pairs$direction <- ifelse(pairs$r >= 0, "positive", "negative")

sig <- pairs[pairs$significant, , drop = FALSE]
sig <- sig[order(-abs(sig$r)), , drop = FALSE]
message("Found ", nrow(sig), " significant gene-metabolite pairs (|r| >= ", cor_cutoff,
        ", padj <= ", padj_cutoff, ") out of ", nrow(pairs), " tested.")

sig_out <- sig[, c("gene_id", "metabolite_id", "r", "p", "padj", "direction")]
write.csv(sig_out, file.path(outdir, "gene_metabolite_pairs_sig.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Correlation heatmap of the significant pairs
# ---------------------------------------------------------------------------

if (nrow(sig) > 0) {
  hm_genes <- unique(sig$gene_id)
  hm_mets <- unique(sig$metabolite_id)
  hm <- matrix(NA_real_, nrow = length(hm_genes), ncol = length(hm_mets),
               dimnames = list(hm_genes, hm_mets))
  hm[cbind(match(sig$gene_id, hm_genes), match(sig$metabolite_id, hm_mets))] <- sig$r
  # Non-significant pairs are encoded as 0 (the neutral mid-colour of the
  # diverging palette) so that hierarchical clustering stays well-defined
  # (NA values would break hclust / ComplexHeatmap clustering).
  hm[is.na(hm)] <- 0
  lim <- max(abs(hm), na.rm = TRUE)

  width_mm <- min(180, max(89, 25 + 2.5 * ncol(hm)))
  height_mm <- min(250, max(60, 25 + 2.4 * nrow(hm)))

  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    if (requireNamespace("circlize", quietly = TRUE)) {
      hm_col <- circlize::colorRamp2(
        c(-lim, 0, lim),
        c(diverging_palette()[1], diverging_palette()[4], diverging_palette()[7])
      )
    } else {
      hm_col <- diverging_palette()
    }
    ch <- ComplexHeatmap::Heatmap(
      hm, name = "correlation", col = hm_col, na_col = "#F5F2EC",
      cluster_rows = nrow(hm) > 1, cluster_columns = ncol(hm) > 1,
      row_title = "DE genes", column_title = "Differential metabolites",
      heatmap_width = grid::unit(mm_to_inch(width_mm) * 0.8, "in"),
      heatmap_height = grid::unit(mm_to_inch(height_mm) * 0.8, "in")
    )
    save_nature_figure(ch, file.path(outdir, "correlation_heatmap"),
                       width_mm = width_mm, height_mm = height_mm)
  } else {
    # pheatmap fallback (ComplexHeatmap not installed): the object is rendered
    # once and re-drawn into the three output devices via its gtable.
    hm_p <- pheatmap::pheatmap(
      hm, color = diverging_palette(), na_col = "#F5F2EC",
      cluster_rows = nrow(hm) > 1, cluster_cols = ncol(hm) > 1,
      silent = TRUE,
      main = "Gene-metabolite correlations"
    )
    save_pheatmap(hm_p, file.path(outdir, "correlation_heatmap"), width_mm, height_mm)
  }
} else {
  message("Note: no significant pairs; correlation heatmap skipped.")
}

# ---------------------------------------------------------------------------
# Bipartite network tables (consumed by steps 05 and 06)
# ---------------------------------------------------------------------------

if (nrow(sig) > 0) {
  edges <- data.frame(
    from = sig$gene_id, to = sig$metabolite_id,
    from_type = "gene", to_type = "metabolite",
    r = sig$r, p = sig$p, padj = sig$padj,
    stringsAsFactors = FALSE
  )
  write.csv(edges, file.path(outdir, "correlation_network_edges.csv"), row.names = FALSE)

  gene_deg <- table(sig$gene_id)
  met_deg <- table(sig$metabolite_id)
  nodes <- rbind(
    data.frame(node = names(gene_deg), type = "gene", degree = as.integer(gene_deg),
               stringsAsFactors = FALSE),
    data.frame(node = names(met_deg), type = "metabolite",
               degree = as.integer(met_deg), stringsAsFactors = FALSE)
  )
  nodes <- nodes[order(-nodes$degree), , drop = FALSE]
  write.csv(nodes, file.path(outdir, "correlation_network_nodes.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(gene_id = character(0), metabolite_id = character(0),
                       r = numeric(0), p = numeric(0), padj = numeric(0),
                       direction = character(0)),
            file.path(outdir, "gene_metabolite_pairs_sig.csv"), row.names = FALSE)
  message("Note: no significant pairs; network tables written empty.")
}

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

write_run_info_r(
  outdir = outdir, script_name = "03_correlation.R",
  params = list(de_genes = de_genes_path, de_metabolites = de_mets_path,
                rna = rna_path, metabolites = met_path, method = method,
                cor_cutoff = cor_cutoff, padj = padj_cutoff, block = block_size,
                n_de_genes = nrow(rna_sub), n_diff_metabolites = nrow(met_sub),
                n_samples = length(common), n_pairs_tested = nrow(pairs),
                n_pairs_significant = nrow(sig)),
  packages = c("WGCNA", "pheatmap"),
  extra = c(
    "Correlation is not causation: interpret significant pairs jointly with the",
    "pathway-level analysis in 04_joint_pathway.R.",
    "Spearman p-values use the WGCNA Student approximation; missing metabolite",
    "values are handled by listwise deletion within each block."
  )
)

message("Done. Outputs written to: ", normalizePath(outdir))
