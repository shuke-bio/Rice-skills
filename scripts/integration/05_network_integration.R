###############################################################################
# 05_network_integration.R -- WGCNA module x metabolome integration and the
#                             three-layer multi-omics network
#                             (Rice-skills skill, module 3, step 05)
#
# Title:    Module eigengene-metabolite correlations and a gene-module-
#           metabolite multi-omics network
# Description:
#   Integrates WGCNA co-expression modules (from the transcriptome WGCNA step)
#   with the metabolome:
#     1. Module eigengenes (MEs): the first principal component of the
#        expression matrix of each module's genes (computed with prcomp on the
#        aligned expression matrix). We deliberately do NOT re-run WGCNA's
#        moduleEigenges here: recomputing the ME from the expression matrix is
#        transparent, reproducible and independent of any stored WGCNA object,
#        and it guarantees that the kME values (step 3) are computed against
#        exactly the same ME definition used for the module-metabolite
#        correlations.
#     2. ME x metabolite correlations (Spearman, WGCNA Student p-values, BH).
#     3. Hub genes: the 10 genes with the highest |kME| per module, where
#        kME = correlation between a gene and its module eigengene.
#     4. Hub gene x metabolite correlations.
#     5. A three-layer network: genes -- modules -- metabolites, where edges
#        are significant correlations / memberships. This network is the
#        systems-level summary of the integration and is rendered with igraph
#        (ggraph is used automatically when installed).
#
# Usage:
#   Rscript scripts/integration/05_network_integration.R \
#       --gene-modules results/transcriptome/05_wgcna/gene_modules.csv \
#       --rna results/transcriptome/02_de/vst_matrix.csv \
#       --metabolites results/metabolome/01_preprocess/metabolites_preprocessed.csv \
#       --metadata metadata.csv \
#       --outdir results/integration/05_network --cor-cutoff 0.8 --group-col group
#
# Inputs:
#   --gene-modules  WGCNA module assignment table; required columns:
#                   gene_id, module (the unassigned "grey" module is skipped)
#   --rna           gene x sample expression matrix (VST)
#   --metabolites   metabolite x sample abundance matrix
#   --metadata      sample metadata (first column = sample_id, group column)
#
# Outputs (written to --outdir):
#   module_eigengenes.csv             ME per module per sample
#   module_metabolite_cor.csv         ME x metabolite correlations
#   module_metabolite_heatmap.*       heatmap of significant correlations
#   hub_table.csv                     hub genes with module, kME and degree
#   multiomics_network_edges.csv      network edges (from, to, types, r, padj)
#   multiomics_network_nodes.csv      network nodes with degree
#   network_plot.*                    three-layer network figure (PDF/TIFF/PNG)
#   run_info.txt
#
# Notes:
#   - Significance threshold for network edges is fixed at padj <= 0.05 (BH);
#     --cor-cutoff controls the |r| threshold (default 0.8).
#   - Missing metabolite values are handled per metabolite (complete-case
#     correlation), so the effective n can differ between metabolites.
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

usage_05 <- paste(
  "Usage: Rscript 05_network_integration.R --gene-modules <gene_modules.csv>",
  "         --rna <vst.csv> --metabolites <matrix.csv> --metadata <metadata.csv>",
  "         [--outdir DIR] [--cor-cutoff 0.8] [--group-col COL] [--seed N]",
  sep = "\n"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE))
gene_modules_path <- opt(opts, "gene-modules", NULL)
rna_path <- opt(opts, "rna", NULL)
met_path <- opt(opts, "metabolites", NULL)
metadata_path <- opt(opts, "metadata", NULL)
if (is.null(gene_modules_path) || is.null(rna_path) || is.null(met_path) ||
    is.null(metadata_path)) {
  stop("Missing required argument(s): --gene-modules, --rna, --metabolites and ",
       "--metadata are required.\n", usage_05, call. = FALSE)
}
outdir <- opt(opts, "outdir", "results/integration/05_network")
cor_cutoff <- as.numeric(opt(opts, "cor-cutoff", "0.8"))
group_col <- opt(opts, "group-col", "group")
seed <- as.integer(opt(opts, "seed", "42"))
if (is.na(cor_cutoff) || cor_cutoff <= 0 || cor_cutoff > 1) {
  stop("--cor-cutoff must be in (0, 1], e.g. 0.8.", call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rice_set_seed(seed)

need_pkgs(c("WGCNA", "igraph"))

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
         paste(head(non_numeric, 5), collapse = ", "), call. = FALSE)
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

read_sample_metadata <- function(path, what) {
  df <- read_csv_safe(path, what)
  ids <- as.character(df[[1]])
  if (anyDuplicated(ids)) stop(what, " contains duplicate sample_id values.", call. = FALSE)
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  df
}

#' Save a pheatmap object (class "pheatmap") in PDF + TIFF + PNG (300 dpi).
save_pheatmap <- function(p, file_stem, width_mm, height_mm) {
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
# Load and align inputs
# ---------------------------------------------------------------------------

gene_modules <- read_csv_safe(gene_modules_path, "Gene module table")
check_columns(gene_modules, c("gene_id", "module"), "Gene module table")
gene_modules$gene_id <- as.character(gene_modules$gene_id)
gene_modules$module <- as.character(trimws(gene_modules$module))

rna <- read_feature_matrix(rna_path, "RNA matrix")
met <- read_metabolome_matrix(met_path, "Metabolome matrix")
meta <- read_sample_metadata(metadata_path, "Sample metadata")
if (!group_col %in% colnames(meta)) {
  stop("Sample metadata is missing the group column '", group_col, "'.\n",
       "Available columns: ", paste(colnames(meta), collapse = ", "), call. = FALSE)
}

common <- Reduce(intersect, list(colnames(rna), colnames(met), rownames(meta)))
if (length(common) < 4) {
  stop("Only ", length(common), " sample(s) are shared by RNA / metabolome / metadata; ",
       "at least 4 are required.\n",
       "RNA samples: ", paste(head(colnames(rna), 5), collapse = ", "), " ...\n",
       "metabolome samples: ", paste(head(colnames(met), 5), collapse = ", "), " ...\n",
       "metadata samples: ", paste(head(rownames(meta), 5), collapse = ", "), " ...\n",
       "Hint: sample IDs must match exactly (case-sensitive, no extra spaces).",
       call. = FALSE)
}
rna <- rna[, common, drop = FALSE]
met <- met[, common, drop = FALSE]
group_vec <- setNames(as.character(meta[[group_col]][common, drop = TRUE]), common)

# Restrict the module table to genes present in the expression matrix.
missing_genes <- setdiff(gene_modules$gene_id, rownames(rna))
if (length(missing_genes) > 0) {
  message("Note: ", length(missing_genes), " module gene(s) are not in the RNA matrix ",
          "and were skipped (e.g. ", paste(head(missing_genes, 3), collapse = ", "), ").")
}
gene_modules <- gene_modules[gene_modules$gene_id %in% rownames(rna), , drop = FALSE]

# The WGCNA "grey" module holds unassigned genes; it is not a real module.
grey <- tolower(gene_modules$module) == "grey"
if (any(grey)) {
  message("Note: skipping ", sum(grey), " gene(s) in the unassigned 'grey' module.")
  gene_modules <- gene_modules[!grey, , drop = FALSE]
}
if (nrow(gene_modules) == 0) stop("No module genes left after filtering.", call. = FALSE)

modules <- sort(unique(gene_modules$module))
module_genes <- split(gene_modules$gene_id, gene_modules$module)
module_genes <- module_genes[modules]
module_genes <- lapply(module_genes, function(g) unique(g))

# Drop modules that are too small for a stable PC1.
small <- vapply(module_genes, length, integer(1)) < 3
if (any(small)) {
  message("Note: skipping ", sum(small), " module(s) with fewer than 3 genes: ",
          paste(names(small)[small], collapse = ", "), ".")
  modules <- modules[!small]
  module_genes <- module_genes[!small]
}
if (length(modules) == 0) stop("No module with >= 3 genes remains.", call. = FALSE)

# Drop genes with missing expression values.
na_genes <- rownames(rna)[rowSums(!is.finite(rna)) > 0]
if (length(na_genes) > 0) {
  message("Note: dropping ", length(na_genes), " gene(s) with missing values from the RNA matrix.")
  rna <- rna[setdiff(rownames(rna), na_genes), , drop = FALSE]
  gene_modules <- gene_modules[gene_modules$gene_id %in% rownames(rna), , drop = FALSE]
  module_genes <- lapply(module_genes, function(g) intersect(g, rownames(rna)))
  modules <- modules[vapply(module_genes, length, integer(1)) >= 3]
  module_genes <- module_genes[modules]
}

message("Modules: ", length(modules), " (",
        paste(vapply(module_genes, length, integer(1)), collapse = ", "), " genes); ",
        nrow(met), " metabolites; ", length(common), " samples.")

# ---------------------------------------------------------------------------
# 1) Module eigengenes: PC1 of each module's expression matrix (prcomp)
# ---------------------------------------------------------------------------

me_list <- lapply(modules, function(m) {
  g <- module_genes[[m]]
  expr_m <- t(rna[g, , drop = FALSE])          # samples x genes
  pc <- stats::prcomp(expr_m, center = TRUE, scale. = TRUE)
  scores <- pc$x[, 1]
  # Sign convention: orient the ME like the module's average expression per
  # sample (rowMeans) so that kME values are comparable across modules.
  if (stats::cor(scores, rowMeans(expr_m)) < 0) scores <- -scores
  scores
})
names(me_list) <- modules
me_mat <- do.call(cbind, me_list)
colnames(me_mat) <- modules
me_df <- data.frame(sample_id = common, group = group_vec[common], me_mat,
                    check.names = FALSE, stringsAsFactors = FALSE)
write.csv(me_df, file.path(outdir, "module_eigengenes.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2) Module eigengene x metabolite correlations (complete-case per metabolite)
# ---------------------------------------------------------------------------

cor_with_vector <- function(mat, y, method = "spearman") {
  # mat: samples x k matrix; y: length-samples vector (already complete).
  vapply(seq_len(ncol(mat)), function(k) stats::cor(mat[, k], y, method = method),
         numeric(1))
}

n_met <- nrow(met)
mod_r <- matrix(NA_real_, nrow = n_met, ncol = length(modules),
                dimnames = list(rownames(met), modules))
mod_p <- mod_r
mod_n <- integer(n_met)
for (i in seq_len(n_met)) {
  y <- met[i, ]
  ok <- is.finite(y)
  if (sum(ok) < 3) next
  mod_n[i] <- sum(ok)
  mod_r[i, ] <- cor_with_vector(me_mat[ok, , drop = FALSE], y[ok])
  mod_p[i, ] <- WGCNA::corPvalueStudent(cor = mod_r[i, ], nSamples = sum(ok))
}

mod_long <- data.frame(
  module = rep(colnames(mod_r), each = n_met),
  metabolite_id = rep(rownames(mod_r), times = ncol(mod_r)),
  r = as.vector(mod_r),
  p = as.vector(mod_p),
  n = rep(mod_n, times = ncol(mod_r)),
  stringsAsFactors = FALSE
)
mod_long <- mod_long[is.finite(mod_long$r) & is.finite(mod_long$p), , drop = FALSE]
mod_long$padj <- stats::p.adjust(mod_long$p, method = "BH")
mod_long$significant <- abs(mod_long$r) >= cor_cutoff & mod_long$padj <= 0.05
mod_long <- mod_long[order(mod_long$module, -abs(mod_long$r)), , drop = FALSE]
write.csv(mod_long, file.path(outdir, "module_metabolite_cor.csv"), row.names = FALSE)

# Heatmap of significant module-metabolite correlations.
mod_sig <- mod_long[mod_long$significant, , drop = FALSE]
if (nrow(mod_sig) > 0) {
  keep_mets <- unique(mod_sig$metabolite_id)
  if (length(keep_mets) > 60) {
    rank <- vapply(keep_mets, function(m) max(abs(mod_sig$r[mod_sig$metabolite_id == m])),
                   numeric(1))
    keep_mets <- names(sort(rank, decreasing = TRUE))[seq_len(60)]
    message("Note: heatmap shows the top 60 of ", length(unique(mod_sig$metabolite_id)),
            " significant metabolites (ranked by max |r|).")
  }
  hm <- matrix(NA_real_, nrow = length(modules), ncol = length(keep_mets),
               dimnames = list(modules, keep_mets))
  for (i in seq_len(nrow(mod_sig))) {
    hm[mod_sig$module[i], mod_sig$metabolite_id[i]] <- mod_sig$r[i]
  }
  hm <- hm[, colnames(hm) %in% keep_mets, drop = FALSE]
  # Non-significant module-metabolite pairs are encoded as 0 (neutral
  # mid-colour of the diverging palette) so that clustering is well-defined.
  hm[is.na(hm)] <- 0
  lim <- max(abs(hm), na.rm = TRUE)
  width_mm <- min(180, max(89, 25 + 2.5 * ncol(hm)))
  height_mm <- min(180, max(60, 25 + 5 * nrow(hm)))

  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    hm_col <- if (requireNamespace("circlize", quietly = TRUE)) {
      circlize::colorRamp2(c(-lim, 0, lim),
                           c(diverging_palette()[1], diverging_palette()[4],
                             diverging_palette()[7]))
    } else {
      diverging_palette()
    }
    ch <- ComplexHeatmap::Heatmap(
      hm, name = "correlation", col = hm_col, na_col = "#F5F2EC",
      cluster_rows = nrow(hm) > 1, cluster_columns = ncol(hm) > 1,
      row_title = "WGCNA modules", column_title = "Metabolites",
      heatmap_width = grid::unit(mm_to_inch(width_mm) * 0.8, "in"),
      heatmap_height = grid::unit(mm_to_inch(height_mm) * 0.8, "in")
    )
    save_nature_figure(ch, file.path(outdir, "module_metabolite_heatmap"),
                       width_mm = width_mm, height_mm = height_mm)
  } else {
    hm_p <- pheatmap::pheatmap(
      hm, color = diverging_palette(), na_col = "#F5F2EC",
      cluster_rows = nrow(hm) > 1, cluster_cols = ncol(hm) > 1, silent = TRUE,
      main = "Module-metabolite correlations"
    )
    save_pheatmap(hm_p, file.path(outdir, "module_metabolite_heatmap"), width_mm, height_mm)
  }
} else {
  message("Note: no significant module-metabolite correlations; heatmap skipped.")
}

# ---------------------------------------------------------------------------
# 3) Hub genes: top 10 |kME| per module (kME = cor(gene, ME))
# ---------------------------------------------------------------------------

kme_list <- lapply(modules, function(m) {
  g <- module_genes[[m]]
  r <- cor_with_vector(t(rna[g, , drop = FALSE]), me_mat[, m])
  names(r) <- g
  r
})
names(kme_list) <- modules

hub_rows <- lapply(modules, function(m) {
  k <- kme_list[[m]]
  ord <- order(-abs(k))[seq_len(min(10, length(k)))]
  data.frame(gene_id = names(k)[ord], module = m, kME = k[ord],
             stringsAsFactors = FALSE)
})
hubs <- do.call(rbind, hub_rows)
hub_genes <- unique(hubs$gene_id)

# ---------------------------------------------------------------------------
# 4) Hub gene x metabolite correlations
# ---------------------------------------------------------------------------

hub_mat <- t(rna[hub_genes, , drop = FALSE])   # samples x hub genes
hub_r <- matrix(NA_real_, nrow = n_met, ncol = length(hub_genes),
                dimnames = list(rownames(met), hub_genes))
hub_p <- hub_r
hub_n <- integer(n_met)
for (i in seq_len(n_met)) {
  y <- met[i, ]
  ok <- is.finite(y)
  if (sum(ok) < 3) next
  hub_n[i] <- sum(ok)
  hub_r[i, ] <- cor_with_vector(hub_mat[ok, , drop = FALSE], y[ok])
  hub_p[i, ] <- WGCNA::corPvalueStudent(cor = hub_r[i, ], nSamples = sum(ok))
}

hub_long <- data.frame(
  gene_id = rep(colnames(hub_r), each = n_met),
  metabolite_id = rep(rownames(hub_r), times = ncol(hub_r)),
  r = as.vector(hub_r), p = as.vector(hub_p), n = rep(hub_n, times = ncol(hub_r)),
  stringsAsFactors = FALSE
)
hub_long <- hub_long[is.finite(hub_long$r) & is.finite(hub_long$p), , drop = FALSE]
hub_long$padj <- stats::p.adjust(hub_long$p, method = "BH")
hub_sig <- hub_long[abs(hub_long$r) >= cor_cutoff & hub_long$padj <= 0.05, , drop = FALSE]
message("Significant correlations: ", nrow(mod_sig), " module-metabolite, ",
        nrow(hub_sig), " hub-gene-metabolite.")

# ---------------------------------------------------------------------------
# 5) Three-layer network: genes -- modules -- metabolites
# ---------------------------------------------------------------------------

# Edges: hub gene -> module (membership, weight = kME),
#        module -> metabolite (significant correlation),
#        hub gene -> metabolite (significant correlation).
# Zero-length-safe edge builder: an empty input must yield an empty frame with
# matching columns, otherwise data.frame() fails on 0-length vs scalar columns.
edge_frame <- function(from, to, from_type, to_type, r, padj) {
  if (length(from) == 0L || length(to) == 0L) {
    return(data.frame(from = character(0), to = character(0),
                      from_type = character(0), to_type = character(0),
                      r = numeric(0), padj = numeric(0), stringsAsFactors = FALSE))
  }
  data.frame(from = from, to = to, from_type = from_type, to_type = to_type,
             r = r, padj = padj, stringsAsFactors = FALSE)
}
edge_gene_module <- edge_frame(hubs$gene_id, paste0("ME:", hubs$module),
                               "gene", "module", hubs$kME, NA_real_)
edge_module_met <- edge_frame(paste0("ME:", mod_sig$module), mod_sig$metabolite_id,
                              "module", "metabolite", mod_sig$r, mod_sig$padj)
edge_gene_met <- edge_frame(hub_sig$gene_id, hub_sig$metabolite_id,
                            "gene", "metabolite", hub_sig$r, hub_sig$padj)
edges <- do.call(rbind, list(edge_gene_module, edge_module_met, edge_gene_met))
edges <- edges[!duplicated(edges[, c("from", "to")]), , drop = FALSE]

node_gene <- unique(hubs$gene_id)
node_met <- unique(c(mod_sig$metabolite_id, hub_sig$metabolite_id))
# Zero-length-safe metabolite node frame (same rationale as edge_frame above).
node_met_frame <- if (length(node_met) > 0L) {
  data.frame(node = node_met, type = "metabolite", module = NA_character_,
             stringsAsFactors = FALSE)
} else {
  data.frame(node = character(0), type = character(0), module = character(0),
             stringsAsFactors = FALSE)
}
nodes <- rbind(
  data.frame(node = node_gene, type = "gene",
             module = vapply(node_gene, function(g) paste(hubs$module[hubs$gene_id == g],
                                                          collapse = ";"), character(1)),
             stringsAsFactors = FALSE),
  data.frame(node = paste0("ME:", modules), type = "module", module = modules,
             stringsAsFactors = FALSE),
  node_met_frame
)

if (nrow(edges) == 0) {
  message("Note: no significant network edges; the network plot is skipped.")
  write.csv(edges, file.path(outdir, "multiomics_network_edges.csv"), row.names = FALSE)
  write.csv(nodes, file.path(outdir, "multiomics_network_nodes.csv"), row.names = FALSE)
} else {
  g <- igraph::graph_from_data_frame(edges[, c("from", "to")], directed = FALSE,
                                     vertices = unique(c(edges$from, edges$to)))
  deg <- igraph::degree(g)
  nodes$degree <- deg[nodes$node]
  nodes$degree[is.na(nodes$degree)] <- 0L
  edges$degree_from <- deg[edges$from]
  edges$degree_to <- deg[edges$to]
  write.csv(edges, file.path(outdir, "multiomics_network_edges.csv"), row.names = FALSE)
  nodes <- nodes[order(-nodes$degree), , drop = FALSE]
  write.csv(nodes, file.path(outdir, "multiomics_network_nodes.csv"), row.names = FALSE)

  # Hub table: hub genes (+ module nodes) with module, degree and kME.
  hub_table <- data.frame(
    node = hubs$gene_id, type = "gene", module = hubs$module,
    degree = ifelse(is.na(deg[hubs$gene_id]), 0L, deg[hubs$gene_id]),
    kME = hubs$kME, stringsAsFactors = FALSE
  )
  hub_table <- rbind(hub_table,
                     data.frame(node = paste0("ME:", modules), type = "module",
                                module = modules,
                                degree = ifelse(is.na(deg[paste0("ME:", modules)]), 0L,
                                                deg[paste0("ME:", modules)]),
                                kME = NA_real_, stringsAsFactors = FALSE))
  write.csv(hub_table, file.path(outdir, "hub_table.csv"), row.names = FALSE)

  # --- three-layer layout: genes (left) -- modules (middle) -- metabolites (right)
  pal <- rice_palette()
  type_col <- c(gene = pal$series12[1], module = pal$series12[4],
                metabolite = pal$series12[2])
  vtype <- nodes$type[match(igraph::V(g)$name, nodes$node)]
  vcol <- type_col[vtype]
  vsize <- ifelse(vtype == "module", 8, 3 + 1.2 * sqrt(deg[igraph::V(g)$name]))
  etype <- ifelse(edges$from_type == "gene" & edges$to_type == "module", "gene-module",
                  ifelse(edges$from_type == "gene", "gene-metabolite", "module-metabolite"))
  ecol <- ifelse(etype == "gene-module", nature_semantic("not_significant"),
                 ifelse(edges$r >= 0, nature_semantic("up"), nature_semantic("down")))
  ecol <- ifelse(etype == "gene-metabolite",
                 grDevices::adjustcolor(ecol, alpha.f = 0.45), ecol)
  ewidth <- ifelse(etype == "module-metabolite", 0.4 + 1.6 * abs(edges$r), 0.4)

  # --- three-layer layout: genes (left) -- modules (middle) -- metabolites (right)
  # x positions by layer; y positions: modules evenly spaced, each module's hub
  # genes stacked in a narrow band around the module, metabolites evenly spaced
  # (ordered by correlation strength).
  x <- ifelse(vtype == "gene", 0, ifelse(vtype == "module", 1, 2))
  mod_y <- setNames(seq(0, 1, length.out = length(modules) + 2)[seq_along(modules)], modules)

  gene_node <- igraph::V(g)$name[vtype == "gene"]
  y_gene <- setNames(rep(NA_real_, length(gene_node)), gene_node)
  band <- 0.45 / max(1, length(modules) * 4)
  for (m in modules) {
    gh <- intersect(hubs$gene_id[hubs$module == m], gene_node)
    if (length(gh) == 0) next
    offs <- (seq_along(gh) - (length(gh) + 1) / 2) * band
    y_gene[gh] <- mod_y[[m]] + offs
  }
  # Fallback for genes whose modules could not be placed (should not happen).
  y_gene[is.na(y_gene)] <- 0.5

  met_node <- igraph::V(g)$name[vtype == "metabolite"]
  met_strength <- vapply(met_node, function(m) {
    rr <- abs(c(edges$r[edges$to == m & etype == "module-metabolite"],
                edges$r[edges$to == m & etype == "gene-metabolite"]))
    if (length(rr) == 0) 0 else max(rr)
  }, numeric(1))
  met_order <- met_node[order(-met_strength)]
  met_y <- setNames(seq(0, 1, length.out = length(met_order) + 2)[seq_along(met_order)],
                    met_order)

  y <- numeric(length(igraph::V(g)))
  y[vtype == "module"] <- as.numeric(mod_y[sub("^ME:", "", igraph::V(g)$name[vtype == "module"])])
  y[vtype == "gene"] <- as.numeric(y_gene[gene_node])
  y[vtype == "metabolite"] <- as.numeric(met_y[met_node])
  coords <- cbind(x, y)

  # Labels: modules always; top genes/metabolites by degree only.
  vlabels <- rep(NA_character_, length(igraph::V(g)))
  vlabels[vtype == "module"] <- sub("^ME:", "", igraph::V(g)$name[vtype == "module"])
  top_gene <- gene_node[order(-deg[gene_node])][seq_len(min(10, length(gene_node)))]
  top_met <- met_node[order(-deg[met_node])][seq_len(min(12, length(met_node)))]
  vlabels[igraph::V(g)$name %in% top_gene] <-
    igraph::V(g)$name[igraph::V(g)$name %in% top_gene]
  vlabels[igraph::V(g)$name %in% top_met] <-
    igraph::V(g)$name[igraph::V(g)$name %in% top_met]

  draw_network <- function() {
    op <- graphics::par(mar = c(1, 1, 3, 1))
    on.exit(graphics::par(op), add = TRUE)
    igraph::plot.igraph(g, layout = coords, rescale = FALSE,
                        xlim = c(-0.28, 2.55), ylim = c(-0.06, 1.06),
                        vertex.color = vcol, vertex.size = vsize,
                        vertex.frame.color = "grey30", vertex.label = vlabels,
                        vertex.label.cex = 0.55, vertex.label.dist = 0,
                        vertex.label.degree = ifelse(x < 1, -pi / 2, pi / 2),
                        edge.color = ecol, edge.width = ewidth,
                        main = "Multi-omics network: genes - modules - metabolites")
    graphics::legend("bottomleft", legend = c("Gene", "Module", "Metabolite",
                                              "Module-metabolite (r)",
                                              "Gene-metabolite (r)"),
                     pch = c(16, 16, 16, NA, NA), lty = c(NA, NA, NA, 1, 1),
                     col = c(type_col, nature_semantic("up"),
                             grDevices::adjustcolor(nature_semantic("up"), alpha.f = 0.45)),
                     pt.cex = 1.2, cex = 0.65, bty = "n")
  }
  for (f in c("pdf", "tiff", "png")) {
    out <- paste0(file.path(outdir, "network_plot"), ".", if (f == "tiff") "tiff" else f)
    if (f == "pdf") {
      grDevices::pdf(out, width = mm_to_inch(180), height = mm_to_inch(120))
    } else if (f == "png") {
      grDevices::png(out, width = mm_to_inch(180), height = mm_to_inch(120), units = "in", res = 300)
    } else {
      grDevices::tiff(out, width = mm_to_inch(180), height = mm_to_inch(120),
                      units = "in", res = 300, compression = "lzw")
    }
    draw_network()
    grDevices::dev.off()
  }
}

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

write_run_info_r(
  outdir = outdir, script_name = "05_network_integration.R",
  params = list(gene_modules = gene_modules_path, rna = rna_path,
                metabolites = met_path, metadata = metadata_path,
                cor_cutoff = cor_cutoff, group_col = group_col, seed = seed,
                n_modules = length(modules), n_hub_genes = length(hub_genes),
                n_samples = length(common),
                n_module_metabolite_sig = nrow(mod_sig),
                n_hub_metabolite_sig = nrow(hub_sig)),
  packages = c("WGCNA", "igraph"),
  extra = c(
    "Module eigengenes are PC1 of each module's expression matrix (prcomp),",
    "oriented to correlate positively with the module mean expression.",
    "kME = cor(gene, ME); hubs = top 10 |kME| per module.",
    "Network edges: hub-module membership, module-metabolite and",
    "hub gene-metabolite correlations with |r| >= cutoff and padj <= 0.05."
  )
)

message("Done. Outputs written to: ", normalizePath(outdir))
