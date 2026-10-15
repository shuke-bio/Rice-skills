###############################################################################
# 06_metabolite_network.R -- metabolite-metabolite correlation network (step 06)
#
# Title:    Blockwise correlation network of metabolites with hub detection
# Description:
#   Computes pairwise correlations between all metabolites (Spearman by
#   default) in column blocks to bound memory usage, keeps pairs with |r| >=
#   --cor-cutoff, tests them with WGCNA::corPvalueStudent (avoids the psych
#   dependency), applies BH correction across the tested pairs, and builds an
#   igraph network. Node degree and betweenness centrality identify hub
#   metabolites (top 10% by degree).
#
# Usage:
#   Rscript scripts/metabolome/06_metabolite_network.R \
#       --data results/metabolome/01_preprocess/metabolome_preprocessed.csv \
#       --outdir results/metabolome/06_network \
#       --method spearman --cor-cutoff 0.8 --padj 0.05 --block 500
#
# Inputs:
#   --data  Preprocessed abundance matrix (metabolite_id [, annotations],
#           samples...), e.g. the output of 01_preprocess.py.
#
# Outputs (written to --outdir):
#   metabolite_edges.csv  from, to, r, padj
#   metabolite_nodes.csv  metabolite_id, degree, betweenness, is_hub
#   network_plot.*        network (ggraph when available, else igraph base device)
#   run_info.txt          parameters, package versions and the seed
#
# Notes:
#   - Blockwise computation: metabolites are processed in blocks of --block
#     columns; each block is correlated against all metabolites and only the
#     upper-triangle pairs (j > i) are kept, so each pair is evaluated once and
#     memory stays bounded for large matrices.
#   - p-values are computed only for pairs passing the |r| threshold and BH-
#     corrected across that tested family (pairs below the threshold are not
#     tested and therefore not part of the correction family).
#   - Edge colours follow the semantic convention: positive correlations use
#     the 'up' colour (#CE6A6B), negative correlations the 'down' colour
#     (#5B84AE); node size encodes the degree.
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
# Dependency check (WGCNA provides corPvalueStudent; igraph builds the network)
# ---------------------------------------------------------------------------

need_pkgs(c("WGCNA", "igraph"))

# ---------------------------------------------------------------------------
# Blockwise correlation
# ---------------------------------------------------------------------------

ANNOTATION_COLS <- c("metabolite_id", "metabolite_name", "class", "kegg_id")

#' Identify the numeric sample columns of a metabolite matrix.
sample_columns <- function(df) {
  names(df)[vapply(df, is.numeric, logical(1))]
}

#' Blockwise pairwise correlation with |r| thresholding.
#'
#' @param X Samples x metabolites matrix.
#' @param method Correlation method passed to stats::cor.
#' @param cutoff Absolute correlation threshold.
#' @param block Block size (number of metabolite columns per chunk).
#' @return Data frame with from, to, r (upper-triangle pairs with |r| >= cutoff).
blockwise_correlations <- function(X, method, cutoff, block) {
  met_ids <- colnames(X)
  p <- ncol(X)
  n_samples <- nrow(X)
  pairs <- list()
  for (start in seq(1L, p, by = block)) {
    idx <- start:min(start + block - 1L, p)
    # Correlate the block of metabolites against all metabolites.
    r_block <- stats::cor(X[, idx, drop = FALSE], X, method = method,
                          use = "pairwise.complete.obs")
    if (is.null(dim(r_block))) r_block <- matrix(r_block, ncol = p)  # single-column block
    for (a in seq_along(idx)) {
      i <- idx[a]
      r_row <- r_block[a, ]
      j_all <- seq_len(p)
      keep <- j_all > i  # upper triangle only (each pair once)
      if (!any(keep)) next
      r_cand <- r_row[keep]
      hit <- which(abs(r_cand) >= cutoff)
      if (length(hit)) {
        j_idx <- j_all[keep][hit]
        pairs[[length(pairs) + 1L]] <- data.frame(
          from = met_ids[i], to = met_ids[j_idx], r = as.numeric(r_cand[hit]),
          stringsAsFactors = FALSE)
      }
    }
    message("  correlated metabolites ", min(idx), "-", max(idx), " of ", p, " ...")
  }
  if (!length(pairs)) {
    return(data.frame(from = character(0), to = character(0), r = numeric(0)))
  }
  do.call(rbind, pairs)
}

# ---------------------------------------------------------------------------
# Network construction and plotting
# ---------------------------------------------------------------------------

#' Build the igraph object and centrality table.
build_network <- function(edges_df, met_ids) {
  g <- igraph::graph_from_data_frame(edges_df[, c("from", "to"), drop = FALSE],
                                     directed = FALSE,
                                     vertices = data.frame(name = met_ids, stringsAsFactors = FALSE))
  deg <- igraph::degree(g)
  btw <- igraph::betweenness(g)
  # Hub definition: top 10% by degree (at least one hub when the graph is
  # non-empty).
  n_hub <- max(1L, floor(0.1 * length(met_ids)))
  hub_ids <- names(sort(deg, decreasing = TRUE))[seq_len(min(n_hub, length(met_ids)))]
  nodes <- data.frame(
    metabolite_id = names(deg),
    degree = as.integer(deg),
    betweenness = as.numeric(btw),
    is_hub = names(deg) %in% hub_ids & deg > 0,
    stringsAsFactors = FALSE
  )
  list(graph = g, nodes = nodes)
}

#' Save a base-R igraph plot in PDF + SVG + TIFF (300 dpi).
save_base_network <- function(g, lay, nodes, outdir) {
  deg <- nodes$degree
  deg_range <- max(deg) - min(deg)
  v_size <- if (deg_range > 0) 2 + 8 * (deg - min(deg)) / deg_range else rep(4, length(deg))
  hub <- nodes$is_hub
  v_col <- rice_palette()$series12[1]
  for (fmt in c("pdf", "svg", "tiff")) {
    out <- file.path(outdir, paste0("network_plot.", fmt))
    if (fmt == "tiff") {
      grDevices::tiff(out, width = mm_to_inch(89), height = mm_to_inch(89),
                      units = "in", res = 300, compression = "lzw")
    } else if (fmt == "pdf") {
      grDevices::pdf(out, width = mm_to_inch(89), height = mm_to_inch(89))
    } else {
      grDevices::svg(out, width = mm_to_inch(89), height = mm_to_inch(89))
    }
    igraph::plot.igraph(g, layout = lay,
                 vertex.size = v_size, vertex.color = v_col, vertex.frame.color = NA,
                 vertex.label = ifelse(hub, nodes$metabolite_id, NA),
                 vertex.label.cex = 0.55, vertex.label.dist = 0.6,
                 vertex.label.color = "grey20",
                 edge.color = ifelse(igraph::E(g)$r > 0,
                                     nature_semantic("up"), nature_semantic("down")),
                 edge.width = 0.5, main = "Metabolite correlation network")
    grDevices::dev.off()
  }
  message("ggraph is not installed: saved network_plot with the igraph base device ",
          "(PDF/SVG/TIFF). Install ggraph for a publication-grade layout.")
}

#' ggraph-based network plot (used when ggraph is installed).
plot_ggraph_network <- function(g, nodes, edges_df, outdir) {
  lay <- ggraph::create_layout(g, layout = "fr")
  edge_df <- data.frame(lay$edges, stringsAsFactors = FALSE)
  # Edge colours: positive r -> semantic 'up', negative r -> semantic 'down'.
  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_link(ggplot2::aes(colour = r), alpha = 0.55, linewidth = 0.4) +
    ggraph::geom_node_point(ggplot2::aes(size = degree), colour = rice_palette()$series12[1]) +
    ggraph::geom_node_text(ggplot2::aes(label = ifelse(is_hub, name, "")),
                           size = 2.0, colour = "grey20", repel = TRUE) +
    ggplot2::scale_colour_gradient2(low = nature_semantic("down"),
                                    mid = diverging_palette()[4],
                                    high = nature_semantic("up"),
                                    name = "Correlation (r)") +
    ggplot2::scale_size_continuous(range = c(1.5, 6), name = "Degree") +
    ggplot2::labs(title = "Metabolite correlation network",
                  subtitle = paste0("Spearman |r| >= ", attr(edges_df, "cutoff"),
                                    "; hubs (top 10% degree) labelled; ",
                                    nrow(edges_df), " edges, ",
                                    igraph::vcount(g), " metabolites"),
                  x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  save_nature_figure(p, file.path(outdir, "network_plot"), width_mm = 89, height_mm = 89)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  data_path <- opt(opts, "data")
  if (is.null(data_path)) stop("--data is required (preprocessed abundance matrix).", call. = FALSE)
  outdir <- opt(opts, "outdir", "results/metabolome/06_network")
  method <- opt(opts, "method", "spearman")
  cutoff <- as.numeric(opt(opts, "cor-cutoff", 0.8))
  padj_cut <- as.numeric(opt(opts, "padj", 0.05))
  block <- as.integer(opt(opts, "block", 500))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", normalizePath(outdir))

  # ---- Read matrix -------------------------------------------------------------
  mat_df <- read_csv_safe(data_path, "preprocessed metabolome matrix")
  meta_path <- opt(opts, "metadata", NULL)
  # Sample columns = numeric columns (annotation columns such as metabolite_name,
  # class, kegg_id are character and therefore excluded).
  sample_cols <- setdiff(sample_columns(mat_df), ANNOTATION_COLS)
  if (!is.null(meta_path)) {
    meta_df <- read_csv_safe(meta_path, "sample metadata")
    sample_cols <- intersect(meta_df$sample_id, sample_cols)
    if (!length(sample_cols)) stop("No overlapping samples between --data and --metadata.", call. = FALSE)
  }
  if (!length(sample_cols)) {
    stop("No numeric sample columns found in --data.\n",
         "Hint: the matrix should contain metabolite_id, optional annotation columns, ",
         "and numeric sample columns.", call. = FALSE)
  }
  non_numeric <- sample_cols[!vapply(sample_cols, function(s) is.numeric(mat_df[[s]]), logical(1))]
  if (length(non_numeric)) {
    stop("Non-numeric sample column(s) in --data: ", paste(head(non_numeric, 5), collapse = ", "),
         call. = FALSE)
  }
  met_ids <- as.character(mat_df[[1]])
  X <- t(as.matrix(mat_df[, sample_cols, drop = FALSE]))
  colnames(X) <- met_ids
  rownames(X) <- sample_cols
  storage.mode(X) <- "double"
  message("Correlation matrix: ", nrow(X), " samples x ", ncol(X), " metabolites (method = ",
          method, ", |r| >= ", cutoff, ").")

  # ---- Blockwise correlations -----------------------------------------------------
  cand <- blockwise_correlations(X, method, cutoff, block)
  message(nrow(cand), " candidate pairs with |r| >= ", cutoff, " (before multiple-testing correction).")

  if (!nrow(cand)) {
    message("No metabolite pair reaches |r| >= ", cutoff,
            "; the network is empty (nodes/edges tables are still exported).")
    edges <- data.frame(from = character(0), to = character(0), r = numeric(0), padj = numeric(0))
    nodes <- data.frame(metabolite_id = met_ids, degree = 0L, betweenness = 0,
                        is_hub = FALSE, stringsAsFactors = FALSE)
    utils::write.csv(edges, file.path(outdir, "metabolite_edges.csv"),
                     row.names = FALSE, quote = FALSE)
    utils::write.csv(nodes, file.path(outdir, "metabolite_nodes.csv"),
                     row.names = FALSE, quote = FALSE)
    write_run_info_r(outdir, script_name = "06_metabolite_network.R",
                     params = list(data = data_path, method = method, cor_cutoff = cutoff,
                                   padj_cutoff = padj_cut, block = block, n_edges = 0),
                     packages = c("WGCNA", "igraph"),
                     extra = "No pair passed the correlation threshold; network plot skipped.")
    message("Done (empty network). Results in ", normalizePath(outdir))
    return(invisible(NULL))
  }

  # ---- p-values and BH correction --------------------------------------------------
  p_vn <- as.numeric(WGCNA::corPvalueStudent(cand$r, nrow(X)))
  cand$padj <- stats::p.adjust(p_vn, method = "BH")
  edges <- cand[which(cand$padj <= padj_cut), , drop = FALSE]
  edges <- edges[order(-abs(edges$r)), , drop = FALSE]
  rownames(edges) <- NULL
  attr(edges, "cutoff") <- cutoff
  message(nrow(edges), " edges pass |r| >= ", cutoff, " and padj <= ", padj_cut, ".")

  # ---- Network and centrality --------------------------------------------------------
  net <- build_network(edges, met_ids)
  g <- net$graph
  nodes <- net$nodes
  # Carry node/edge attributes on the graph so that both renderers (ggraph
  # layout and igraph::plot) can use them. graph_from_data_frame keeps the edge
  # order of 'edges', so E(g)$r aligns row-by-row.
  v_match <- match(igraph::V(g)$name, nodes$metabolite_id)
  g <- igraph::set_vertex_attr(g, "degree", value = nodes$degree[v_match])
  g <- igraph::set_vertex_attr(g, "betweenness", value = nodes$betweenness[v_match])
  g <- igraph::set_vertex_attr(g, "is_hub", value = nodes$is_hub[v_match])
  igraph::E(g)$r <- edges$r
  utils::write.csv(edges, file.path(outdir, "metabolite_edges.csv"),
                   row.names = FALSE, quote = FALSE)
  utils::write.csv(nodes, file.path(outdir, "metabolite_nodes.csv"),
                   row.names = FALSE, quote = FALSE)
  hubs <- nodes$metabolite_id[which(nodes$is_hub)]
  message("Hub metabolites (top 10% degree): ", length(hubs), " -> ",
          paste(head(hubs, 10), collapse = ", "), if (length(hubs) > 10) " ..." else "")

  # ---- Figure ---------------------------------------------------------------------------
  lay <- igraph::layout_with_fr(g)
  if (requireNamespace("ggraph", quietly = TRUE)) {
    plot_ggraph_network(g, nodes, edges, outdir)
  } else {
    save_base_network(g, lay, nodes, outdir)
  }

  write_run_info_r(
    outdir,
    script_name = "06_metabolite_network.R",
    params = list(
      data = data_path, method = method, cor_cutoff = cutoff, padj_cutoff = padj_cut,
      block = block, n_samples = nrow(X), n_metabolites = ncol(X),
      n_candidate_pairs = nrow(cand), n_edges = nrow(edges),
      n_hubs = sum(nodes$is_hub),
      renderer = if (requireNamespace("ggraph", quietly = TRUE)) "ggraph" else "igraph base device"
    ),
    packages = c("WGCNA", "igraph"),
    extra = c(
      "p-values: WGCNA::corPvalueStudent on |r| >= cutoff pairs; BH across the tested family.",
      "Edge colours: positive r = semantic 'up', negative r = semantic 'down'; node size = degree."
    )
  )
  message("Done. Results in ", normalizePath(outdir))
  invisible(NULL)
}

main()
