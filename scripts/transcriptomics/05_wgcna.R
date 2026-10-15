###############################################################################
# 05_wgcna.R -- weighted gene co-expression network analysis (WGCNA)
#
# Title:    Co-expression network construction and module analysis (WGCNA)
# Description:
#   Reads the count matrix and metadata, computes CPM and log2(CPM+1),
#   filters lowly expressed genes, selects the top-N genes by MAD, screens
#   the soft-thresholding power, identifies modules with blockwiseModules,
#   relates modules to sample traits, extracts hub genes and enriches the
#   largest modules with GO.
#
# Usage:
#   Rscript scripts/transcriptomics/05_wgcna.R \
#       --counts results/transcriptome/01_qc/counts_filtered.csv \
#       --metadata metadata.csv \
#       --traits tiller_number,Na_content \
#       --outdir results/transcriptome/05_wgcna
#
# Inputs:
#   --counts    count matrix CSV (first column gene_id, columns = samples)
#   --metadata  sample metadata CSV (first column sample_id, group column,
#               optional numeric trait columns)
#
# Outputs (written to --outdir):
#   gene_modules.csv                     gene -> module colour assignment
#   soft_threshold_selection.*           scale-free fit and connectivity plot
#   module_trait_correlation.csv         module x trait correlations (long)
#   module_trait_heatmap.*               module-trait relationship heatmap
#   module_dendrogram.*                  gene dendrogram + module colours
#   hub_genes.csv                        top 25 hub genes per module (by kME)
#   module_top_hub.csv                   top hub gene of each module
#   module_enrichment.csv                GO (BP) enrichment of the 5 largest modules
#   run_info.txt                         provenance (parameters, versions, seed)
#
# Notes:
#   - Memory: for more than ~20,000 genes increase --top-n or run with more
#     RAM; blockwiseModules already splits the network into blocks of
#     maxBlockSize = 5000 genes. The gene dendrogram saved here covers the
#     first block when several blocks are used.
#   - The grey module collects genes that could not be assigned to any
#     module; it is excluded from the module-trait interpretation and from
#     the "5 largest modules" enrichment (standard WGCNA practice).
#   - Results depend on the random seed; it is fixed at 42 (see run_info.txt).
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
outdir <- opt(opts, "outdir", "results/transcriptome/05_wgcna")
traits_arg <- opt(opts, "traits")
top_n <- suppressWarnings(as.integer(opt(opts, "top-n", "10000")))
min_counts <- suppressWarnings(as.numeric(opt(opts, "min-counts", "10")))
network_type <- opt(opts, "network-type", "signed")
deep_split <- suppressWarnings(as.integer(opt(opts, "deep-split", "4")))
min_module_size <- suppressWarnings(as.integer(opt(opts, "min-module-size", "30")))
merge_cut <- suppressWarnings(as.numeric(opt(opts, "merge-cut", "0.25")))
threads <- suppressWarnings(as.integer(opt(opts, "threads", "0")))

if (is.null(counts_path)) {
  stop("--counts is required (count matrix CSV).", call. = FALSE)
}
if (is.null(metadata_path)) {
  stop("--metadata is required (sample metadata CSV).", call. = FALSE)
}
if (is.na(top_n) || top_n < 100L) {
  stop("--top-n must be an integer >= 100; got '", opt(opts, "top-n", "10000"), "'.",
       call. = FALSE)
}
if (is.na(min_counts) || min_counts < 0) {
  stop("--min-counts must be a non-negative number; got '",
       opt(opts, "min-counts", "10"), "'.", call. = FALSE)
}
if (!network_type %in% c("signed", "unsigned", "signed hybrid")) {
  stop("Invalid --network-type '", network_type,
       "'. Supported: signed, unsigned, signed hybrid.", call. = FALSE)
}
if (is.na(deep_split) || deep_split < 0L) {
  stop("--deep-split must be a non-negative integer; got '",
       opt(opts, "deep-split", "4"), "'.", call. = FALSE)
}
if (is.na(min_module_size) || min_module_size < 2L) {
  stop("--min-module-size must be an integer >= 2; got '",
       opt(opts, "min-module-size", "30"), "'.", call. = FALSE)
}
if (is.na(merge_cut) || merge_cut <= 0) {
  stop("--merge-cut must be a positive number; got '",
       opt(opts, "merge-cut", "0.25"), "'.", call. = FALSE)
}
if (is.na(threads) || threads < 0L) {
  stop("--threads must be 0 (automatic) or a positive integer; got '",
       opt(opts, "threads", "0"), "'.", call. = FALSE)
}

need_pkgs(c("WGCNA", "ggplot2", "jsonlite"))
# WGCNA must be attached (not only loaded): its internal do.call("cor", ...)
# lookups resolve the package's own weighted correlation functions only when
# the package is on the search path.
suppressPackageStartupMessages(library(WGCNA))
# Module GO enrichment is an OPTIONAL add-on: it needs clusterProfiler plus the
# rice annotation database. When either is missing the co-expression network
# itself is still computed and the enrichment section is skipped with a note.
enrich_possible <- requireNamespace("clusterProfiler", quietly = TRUE) &&
  requireNamespace("org.Osativa.eg.db", quietly = TRUE)
if (enrich_possible) {
  suppressPackageStartupMessages(library(org.Osativa.eg.db, character.only = TRUE))
  org_db <- get("org.Osativa.eg.db")
} else {
  message(
    "Note: clusterProfiler and/or org.Osativa.eg.db not installed - module GO ",
    "enrichment will be skipped (the network analysis itself is unaffected).\n",
    "Hint: Rscript -e 'BiocManager::install(c(\"clusterProfiler\", \"org.Osativa.eg.db\"))'"
  )
  org_db <- NULL
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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
  if (any(!is.finite(mat)) || any(mat < 0)) {
    stop("The count matrix contains missing or negative values: ", path,
         "\nHint: counts must be non-negative integers without gaps.", call. = FALSE)
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

#' Save base-graphics / grid output in the three Nature formats (PDF/SVG/TIFF).
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

#' Map SYMBOL gene IDs to ENTREZID; NULL when nothing maps.
map_to_entrez <- function(genes, org_db) {
  genes <- unique(as.character(genes))
  if (length(genes) == 0L) return(NULL)
  mapped <- tryCatch(
    clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org_db),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) == 0L) return(NULL)
  unique(as.character(mapped$ENTREZID))
}

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

counts <- read_counts_matrix(counts_path)
aligned <- read_metadata_aligned(metadata_path, counts)
mat <- aligned$matrix
meta <- aligned$meta
message(sprintf("Loaded %d genes x %d samples.", nrow(mat), ncol(mat)))

# CPM + log2 computed by hand (edgeR::cpm would add a dependency; the simple
# formula is transparent and fully reproducible).
lib_sizes <- colSums(mat)
if (any(lib_sizes == 0)) {
  stop("Sample(s) with zero library size: ",
       paste(colnames(mat)[lib_sizes == 0], collapse = ", "), call. = FALSE)
}
cpm <- sweep(mat, 2, lib_sizes, "/") * 1e6
log_cpm <- log2(cpm + 1)

# Low-expression filter (same rule as 01_qc_filter.py; the minimum number of
# expressing samples follows the pipeline default of 3).
keep <- rowSums(mat >= min_counts) >= 3L
message(sprintf("Low-expression filter: %d -> %d genes (count >= %g in >= 3 samples).",
                nrow(mat), sum(keep), min_counts))
log_cpm <- log_cpm[keep, , drop = FALSE]
if (nrow(log_cpm) == 0L) {
  stop("No genes passed the low-expression filter; lower --min-counts.",
       call. = FALSE)
}

# MAD-based selection of the top-N genes for network construction.
if (nrow(log_cpm) > top_n) {
  mad <- apply(log_cpm, 1, stats::mad)
  ord <- order(mad, decreasing = TRUE)
  log_cpm <- log_cpm[ord[seq_len(top_n)], , drop = FALSE]
  message(sprintf("Selected the top %d genes by MAD.", top_n))
}

datExpr <- t(log_cpm)  # WGCNA expects samples in rows, genes in columns

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  message("Removed ", sum(!gsg$goodSamples), " sample(s) and ", sum(!gsg$goodGenes),
          " gene(s) flagged by goodSamplesGenes().")
}
message(sprintf("Expression matrix for WGCNA: %d samples x %d genes.",
                nrow(datExpr), ncol(datExpr)))

# Soft-threshold screening (powers 1:20).
powers <- 1:20
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
                                networkType = network_type, verbose = 0)
fit <- sft$fitIndices
r2 <- -sign(fit$slope) * fit$SFT.R.sq
ok <- which(r2 >= 0.8)
if (length(ok) > 0L) {
  power <- fit$Power[min(ok)]
} else {
  power <- fit$Power[which.max(r2)]
  message("WARNING: no soft threshold reached signed R^2 >= 0.8; using power ",
          power, " (best available fit).")
}
message("Selected soft-thresholding power: ", power)

# Diagnostic figure: scale-free topology fit and mean connectivity.
# Scientific conclusion: the chosen power is the smallest value whose
# scale-free topology fit (signed R^2) exceeds 0.8, making the network
# approximately scale-free as required by WGCNA; the mean-connectivity panel
# shows how network density decreases with the power.
save_base_figure(
  function() {
    graphics::par(mfrow = c(1, 2))
    graphics::plot(fit$Power, r2, xlab = "Soft threshold (power)",
                   ylab = "Scale-free topology model fit, signed R^2",
                   main = "Scale independence", type = "n")
    graphics::text(fit$Power, r2, labels = fit$Power, cex = 0.7,
                   col = nature_semantic("highlight"))
    graphics::abline(h = 0.8, col = nature_semantic("reference_line"), lty = 2)
    graphics::plot(fit$Power, fit$mean.k., xlab = "Soft threshold (power)",
                   ylab = "Mean connectivity", main = "Mean connectivity",
                   type = "n")
    graphics::text(fit$Power, fit$mean.k., labels = fit$Power, cex = 0.7,
                   col = nature_semantic("highlight"))
  },
  file.path(outdir, "soft_threshold_selection"), width_mm = 180, height_mm = 80
)

# Module identification (blockwise to bound memory usage).
net <- if (threads > 0L) {
  WGCNA::blockwiseModules(datExpr, maxBlockSize = 5000, power = power,
                          networkType = network_type, deepSplit = deep_split,
                          minModuleSize = min_module_size, mergeCutHeight = merge_cut,
                          numericLabels = FALSE, pamStage = TRUE, saveTOMs = FALSE,
                          nThreads = threads, verbose = 0)
} else {
  WGCNA::blockwiseModules(datExpr, maxBlockSize = 5000, power = power,
                          networkType = network_type, deepSplit = deep_split,
                          minModuleSize = min_module_size, mergeCutHeight = merge_cut,
                          numericLabels = FALSE, pamStage = TRUE, saveTOMs = FALSE,
                          verbose = 0)
}

# blockwiseModules runs its own goodSamplesGenes pass; align everything to the
# samples/genes it kept.
kept_samples <- rownames(datExpr)[net$goodSamples]
kept_genes <- colnames(datExpr)[net$goodGenes]
datExpr_used <- datExpr[net$goodSamples, net$goodGenes, drop = FALSE]
moduleColors <- net$colors
names(moduleColors) <- kept_genes
module_sizes <- sort(table(moduleColors), decreasing = TRUE)
message("Detected ", length(module_sizes), " modules; sizes: ",
        paste(paste(names(module_sizes), module_sizes, sep = "="),
              collapse = ", "))

utils::write.csv(
  data.frame(gene_id = names(moduleColors),
             module = as.character(moduleColors), stringsAsFactors = FALSE),
  file.path(outdir, "gene_modules.csv"), row.names = FALSE, quote = FALSE
)

# Module eigengenes (removed samples carry NA rows in net$MEs).
MEs <- net$MEs[net$goodSamples, , drop = FALSE]
rownames(MEs) <- kept_samples
MEs <- WGCNA::orderMEs(MEs)

# ---- Module-trait relationships -------------------------------------------------
trait_cols <- NULL
if (!is.null(traits_arg)) {
  trait_cols <- trimws(strsplit(traits_arg, ",", fixed = TRUE)[[1]])
  trait_cols <- trait_cols[nzchar(trait_cols)]
  missing <- setdiff(trait_cols, colnames(meta))
  if (length(missing) > 0L) {
    stop("Trait column(s) not found in the metadata: ",
         paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(meta), collapse = ", "),
         call. = FALSE)
  }
  for (t in trait_cols) {
    if (!is.numeric(meta[[t]])) {
      v <- suppressWarnings(as.numeric(meta[[t]]))
      if (any(is.na(v))) {
        stop("Trait '", t, "' is not numeric and cannot be converted ",
             "(non-numeric values present).", call. = FALSE)
      }
      meta[[t]] <- v
    }
  }
}

if (!is.null(trait_cols) && length(trait_cols) > 0L) {
  traitData <- as.matrix(meta[kept_samples, trait_cols, drop = FALSE])
  if (nrow(traitData) != nrow(MEs)) {
    stop("Trait rows (", nrow(traitData), ") do not match module eigengene rows (",
         nrow(MEs), ").", call. = FALSE)
  }
  mtc <- stats::cor(MEs, traitData, use = "p")
  mtp <- WGCNA::corPvalueStudent(mtc, nSamples = nrow(MEs))
  mt_long <- data.frame(
    module = rep(rownames(mtc), times = ncol(mtc)),
    trait = rep(colnames(mtc), each = nrow(mtc)),
    correlation = as.numeric(mtc),
    pvalue = as.numeric(mtp),
    stringsAsFactors = FALSE
  )
  utils::write.csv(mt_long, file.path(outdir, "module_trait_correlation.csv"),
                   row.names = FALSE, quote = FALSE)

  # Heatmap: ComplexHeatmap when available, otherwise pheatmap.
  # Scientific conclusion: modules whose eigengenes correlate strongly with a
  # trait (|r| high, p < 0.05) are candidate drivers of that phenotype.
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    ht <- ComplexHeatmap::Heatmap(
      mtc, col = diverging_palette(), name = "cor",
      cluster_rows = TRUE, cluster_cols = TRUE,
      column_title = "Module-trait relationships",
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%.2f", mtc[i, j]), x, y,
                        gp = grid::gpar(fontsize = 7))
      }
    )
    save_nature_figure(ht, file.path(outdir, "module_trait_heatmap"),
                       width_mm = 120, height_mm = 100)
  } else if (requireNamespace("pheatmap", quietly = TRUE)) {
    save_base_figure(
      function() {
        pheatmap::pheatmap(mtc, color = diverging_palette(),
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE, number_format = "%.2f",
                           fontsize = 7, main = "Module-trait relationships")
      },
      file.path(outdir, "module_trait_heatmap"), width_mm = 120, height_mm = 100
    )
  } else {
    stop("Neither ComplexHeatmap nor pheatmap is installed; cannot draw the ",
         "module-trait heatmap.\n",
         "Hint: install.packages('pheatmap') or BiocManager::install('ComplexHeatmap').",
         call. = FALSE)
  }
} else {
  message("No --traits given; skipping the module-trait analysis.")
}

# ---- Gene dendrogram with module colours ----------------------------------------
# blockwiseModules splits genes into blocks; the saved tree covers block 1
# (with maxBlockSize >= the number of genes this is the full gene tree).
if (!is.null(net$dendrograms) && length(net$dendrograms) >= 1L) {
  if (length(net$dendrograms) > 1L) {
    message("Note: genes were split into ", length(net$dendrograms),
            " blocks; the dendrogram figure shows block 1.")
  }
  dendro <- net$dendrograms[[1]]
  # blockwiseModules strips the leaf labels from the stored dendrogram and
  # plotOrderedColors() indexes the colour rows positionally via dendro$order,
  # so the colours must be given in the block's gene order (blockGenes), which
  # is the order in which the tree was built.
  blk_genes <- kept_genes[net$blockGenes[[1]]]
  save_base_figure(
    function() {
      WGCNA::plotDendroAndColors(
        dendro,
        colors = data.frame(module = moduleColors[blk_genes]),
        groupLabels = "module", dendroLabels = FALSE, addGuide = TRUE,
        guideHang = 0.05,
        main = "Gene dendrogram and module colours (block 1)"
      )
    },
    file.path(outdir, "module_dendrogram"), width_mm = 180, height_mm = 100
  )
} else {
  message("No dendrogram available in the blockwiseModules result.")
}

# ---- Hub genes -------------------------------------------------------------------
kME <- stats::cor(datExpr_used, MEs, use = "p")
colnames(kME) <- sub("^ME", "", colnames(kME))
hub_rows <- list()
for (m in unique(moduleColors)) {
  if (!(m %in% colnames(kME))) next
  idx <- which(moduleColors == m)
  if (length(idx) == 0L) next
  km <- kME[idx, m]
  ord <- order(abs(km), decreasing = TRUE)
  top <- head(ord, min(25L, length(ord)))
  for (r in seq_along(top)) {
    hub_rows[[length(hub_rows) + 1L]] <- data.frame(
      module = m, gene_id = names(km)[top[r]], kME = unname(km[top[r]]),
      rank = r, stringsAsFactors = FALSE
    )
  }
}
hub_df <- if (length(hub_rows) > 0L) {
  do.call(rbind, hub_rows)
} else {
  data.frame(module = character(0), gene_id = character(0), kME = numeric(0),
             rank = integer(0), stringsAsFactors = FALSE)
}
utils::write.csv(hub_df, file.path(outdir, "hub_genes.csv"),
                 row.names = FALSE, quote = FALSE)

top_hubs <- WGCNA::chooseTopHubInEachModule(datExpr_used, moduleColors,
                                            omitColors = "")
utils::write.csv(
  data.frame(module = names(top_hubs), gene_id = as.character(top_hubs),
             stringsAsFactors = FALSE),
  file.path(outdir, "module_top_hub.csv"), row.names = FALSE, quote = FALSE
)

# ---- GO enrichment of the 5 largest modules (optional) ---------------------------
# The grey module (unassigned genes) is not biologically interpretable and is
# excluded from the ranking. Skipped (with an empty output table) when the
# optional annotation packages are not installed.
if (!enrich_possible) {
  utils::write.csv(
    data.frame(module = character(0), ID = character(0),
               Description = character(0), GeneRatio = character(0),
               BgRatio = character(0), pvalue = numeric(0),
               p.adjust = numeric(0), geneID = character(0),
               Count = numeric(0), stringsAsFactors = FALSE),
    file.path(outdir, "module_enrichment.csv"), row.names = FALSE, quote = FALSE
  )
  message("Module enrichment skipped: annotation packages not installed.")
} else {
candidate_sizes <- module_sizes[setdiff(names(module_sizes), "grey")]
top_mods <- names(candidate_sizes)[seq_len(min(5L, length(candidate_sizes)))]
universe_entrez <- map_to_entrez(kept_genes, org_db)
enr_rows <- list()
for (m in top_mods) {
  genes <- names(moduleColors)[moduleColors == m]
  entrez <- map_to_entrez(genes, org_db)
  if (length(entrez) == 0L) {
    message("Module ", m, ": no genes could be mapped to ENTREZID; skipping.")
    next
  }
  res <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez, OrgDb = org_db, keyType = "ENTREZID", ont = "BP",
      universe = if (length(universe_entrez) > 0L) universe_entrez else NULL,
      pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
    ),
    error = function(e) NULL
  )
  if (!is.null(res) && nrow(res@result) > 0L) {
    df <- as.data.frame(res@result)
    df$module <- m
    enr_rows[[length(enr_rows) + 1L]] <- df
  } else {
    message("Module ", m, ": no significant BP enrichment.")
  }
}
if (length(enr_rows) > 0L) {
  module_enr <- do.call(rbind, enr_rows)
  utils::write.csv(module_enr, file.path(outdir, "module_enrichment.csv"),
                   row.names = FALSE, quote = FALSE)
} else {
  utils::write.csv(
    data.frame(module = character(0), ID = character(0),
               Description = character(0), GeneRatio = character(0),
               BgRatio = character(0), pvalue = numeric(0),
               p.adjust = numeric(0), geneID = character(0),
               Count = numeric(0), stringsAsFactors = FALSE),
    file.path(outdir, "module_enrichment.csv"), row.names = FALSE, quote = FALSE
  )
}
}

write_run_info_r(
  outdir, "05_wgcna.R",
  params = list(
    counts = counts_path,
    metadata = metadata_path,
    traits = if (is.null(traits_arg)) "none" else traits_arg,
    top_n = top_n,
    min_counts = min_counts,
    network_type = network_type,
    deep_split = deep_split,
    min_module_size = min_module_size,
    merge_cut = merge_cut,
    threads = threads,
    power = power
  ),
  packages = c("WGCNA", "clusterProfiler", "org.Osativa.eg.db")
)
message("WGCNA complete. Results in ", outdir)
