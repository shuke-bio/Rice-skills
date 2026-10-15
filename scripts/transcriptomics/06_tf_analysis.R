###############################################################################
# 06_tf_analysis.R -- transcription factor differential expression and PPI
#
# Title:    Transcription factor (TF) family analysis and optional PPI network
# Description:
#   Merges the differential expression results with a transcription factor
#   annotation (gene_id, family), reports differentially expressed TFs
#   (de_TF.csv), summarises each TF family (number of DE/up/down genes plus a
#   Fisher exact test against the annotated background) and draws a top-20
#   family bar chart. With --ppi, a protein-protein interaction network of the
#   DE genes is fetched from STRING (STRINGdb) or read from an offline edge
#   file; the module is optional and is skipped with exit code 0 when neither
#   source is available.
#
# Usage:
#   Rscript scripts/transcriptomics/06_tf_analysis.R \
#       --de-all results/transcriptome/02_de/de_all_salt_vs_control.csv \
#       --tf-annotation data/rice_TF_annotation.csv \
#       --outdir results/transcriptome/06_tf
#   # with an offline STRING edge file (gene_a, gene_b, score):
#   Rscript scripts/transcriptomics/06_tf_analysis.R --de-all ... \
#       --tf-annotation ... --ppi --string-edges data/string_edges.csv
#
# Inputs:
#   --de-all           full DE table (needs gene_id, log2FoldChange, padj)
#   --tf-annotation    optional CSV with columns gene_id, family
#   --string-edges     optional offline STRING edges (gene_a, gene_b, score)
#
# Outputs (written to --outdir):
#   de_TF.csv                    DE genes with their TF family annotation
#   tf_family_stats.csv          per-family DE counts and Fisher tests
#   tf_family_top20.*            top-20 TF family bar chart (PDF/SVG/TIFF)
#   ppi_edges.csv                PPI edges (only with --ppi)
#   ppi_nodes.csv                PPI nodes with degrees (only with --ppi)
#   ppi_network.*                PPI network figure (only with --ppi)
#   run_info.txt                 provenance (parameters, versions, seed)
#
# Notes:
#   - TF annotation sources for rice:
#       * PlantTFDB  http://planttfdb.gao-lab.org/  -> download the Oryza
#         sativa family annotation (two columns: gene ID and family).
#       * iTAK       http://bioinfo.bti.cornell.edu/tool/itak/  -> run the TF
#         identification step; the output contains a gene_id column and a
#         TF_family column. Keep those two columns (rename to gene_id,family)
#         and save as CSV.
#     Without --tf-annotation the script prints this guidance and exits with
#     code 0 (the TF module is optional and must not break a pipeline).
#   - STRING species for Japonica rice is 39947; the online STRINGdb branch
#     requires the STRINGdb package and internet access.
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

de_all_path <- opt(opts, "de-all")
tf_annotation <- opt(opts, "tf-annotation")
outdir <- opt(opts, "outdir", "results/transcriptome/06_tf")
padj_cut <- suppressWarnings(as.numeric(opt(opts, "padj", "0.05")))
lfc_cut <- suppressWarnings(as.numeric(opt(opts, "log2fc", "1")))
ppi_flag <- isTRUE(opt(opts, "ppi", NULL))
string_edges <- opt(opts, "string-edges")

if (is.null(de_all_path)) {
  stop("--de-all is required (a de_all_<contrast>.csv produced by 02_deseq2.R).",
       call. = FALSE)
}
if (is.na(padj_cut) || padj_cut <= 0 || padj_cut >= 1) {
  stop("--padj must be a number in (0, 1); got '", opt(opts, "padj", "0.05"), "'.",
       call. = FALSE)
}
if (is.na(lfc_cut) || lfc_cut < 0) {
  stop("--log2fc must be a non-negative number; got '", opt(opts, "log2fc", "1"), "'.",
       call. = FALSE)
}

need_pkgs(c("ggplot2", "jsonlite"))

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

de <- read_csv_safe(de_all_path, "DE results table")
check_columns(de, c("gene_id", "log2FoldChange", "padj"), "de-all table")
de$gene_id <- as.character(de$gene_id)
de$log2FoldChange <- as.numeric(de$log2FoldChange)
de$padj <- suppressWarnings(as.numeric(de$padj))
de$significant <- !is.na(de$padj) & de$padj <= padj_cut &
  !is.na(de$log2FoldChange) & abs(de$log2FoldChange) >= lfc_cut
de$direction <- ifelse(de$significant,
                       ifelse(de$log2FoldChange > 0, "Up", "Down"), "NS")
message(sprintf("DE table: %d genes, %d significant (padj <= %g, |log2FC| >= %g).",
                nrow(de), sum(de$significant), padj_cut, lfc_cut))

# The TF module is optional: without an annotation we document how to obtain
# one and exit successfully (exit code 0) instead of failing.
if (is.null(tf_annotation)) {
  message("No --tf-annotation provided; skipping the transcription-factor module.")
  message("Obtain a rice TF annotation (CSV with columns gene_id,family) from:")
  message("  - PlantTFDB: http://planttfdb.gao-lab.org/  (Oryza sativa family")
  message("    annotation: gene ID and family, two columns).")
  message("  - iTAK:      http://bioinfo.bti.cornell.edu/tool/itak/  (run the TF")
  message("    identification step; keep the gene_id and TF_family columns).")
  write_run_info_r(
    outdir, "06_tf_analysis.R",
    params = list(de_all = de_all_path, tf_annotation = "none", padj = padj_cut,
                  log2fc = lfc_cut),
    extra = "TF module skipped: no --tf-annotation provided."
  )
  quit(status = 0)
}

tf <- read_csv_safe(tf_annotation, "TF annotation")
check_columns(tf, c("gene_id", "family"), "TF annotation")
tf$gene_id <- as.character(tf$gene_id)
tf$family <- as.character(tf$family)
tf <- tf[!duplicated(tf$gene_id), , drop = FALSE]
message(sprintf("TF annotation: %d genes in %d families.",
                nrow(tf), length(unique(tf$family))))

# DE genes intersected with the TF annotation.
de_tf <- merge(de, tf, by = "gene_id", all = FALSE, sort = FALSE)
utils::write.csv(de_tf, file.path(outdir, "de_TF.csv"),
                 row.names = FALSE, quote = FALSE)

# Per-family statistics with a Fisher exact test against the background of all
# DE-tested genes (the test needs annotated non-DE genes to be informative).
fam_levels <- unique(tf$family)
fam_rows <- list()
for (f in fam_levels) {
  fam_genes <- tf$gene_id[tf$family == f]
  in_fam <- de$gene_id %in% fam_genes
  is_de <- de$significant
  a <- sum(in_fam & is_de)
  b <- sum(in_fam & !is_de)
  c <- sum(!in_fam & is_de)
  d <- sum(!in_fam & !is_de)
  p_val <- if (a > 0L && b > 0L && c > 0L && d > 0L) {
    stats::fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  } else {
    NA_real_
  }
  fam_dir <- de_tf$direction[de_tf$family == f]
  fam_rows[[length(fam_rows) + 1L]] <- data.frame(
    family = f,
    n_annotated = length(fam_genes),
    n_tested = sum(in_fam),
    n_de = sum(fam_dir != "NS"),
    n_up = sum(fam_dir == "Up"),
    n_down = sum(fam_dir == "Down"),
    fisher_p = p_val,
    stringsAsFactors = FALSE
  )
}
fam_stats <- do.call(rbind, fam_rows)
fam_stats$fisher_padj <- stats::p.adjust(fam_stats$fisher_p, method = "BH")
fam_stats <- fam_stats[order(-fam_stats$n_de, -fam_stats$n_annotated), , drop = FALSE]
rownames(fam_stats) <- NULL
utils::write.csv(fam_stats, file.path(outdir, "tf_family_stats.csv"),
                 row.names = FALSE, quote = FALSE)
message(sprintf("TF families with at least one DE gene: %d.",
                sum(fam_stats$n_de > 0)))

# Top-20 TF family bar chart (stacked up/down).
# Scientific conclusion: families with many DE members (e.g. AP2/ERF, NAC,
# bZIP in salt-stress studies) mark the transcriptional programmes that
# respond to the treatment; the Fisher test highlights families enriched
# among DE genes beyond their overall abundance.
top20 <- head(fam_stats, 20)
if (nrow(top20) > 0L) {
  long <- data.frame(
    family = rep(top20$family, 2L),
    direction = rep(c("Up", "Down"), each = nrow(top20)),
    count = c(top20$n_up, top20$n_down),
    stringsAsFactors = FALSE
  )
  long$family <- factor(long$family,
                        levels = rev(unique(top20$family)))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = family, y = count, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(Up = nature_semantic("up"),
                                          Down = nature_semantic("down"))) +
    ggplot2::labs(x = "TF family", y = "Number of DE transcription factors",
                  title = "Top 20 differentially expressed TF families") +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "tf_family_top20"),
                     width_mm = 120, height_mm = max(89, nrow(top20) * 5 + 30))
} else {
  message("No TF family statistics to plot.")
}

# ---------------------------------------------------------------------------
# Optional PPI module (skipped silently-friendly when no source is available)
# ---------------------------------------------------------------------------

if (ppi_flag) {
  edges <- NULL
  if (requireNamespace("STRINGdb", quietly = TRUE)) {
    message("STRINGdb detected: fetching the PPI of the DE genes from STRING ",
            "(rice, combined score > 400).")
    de_genes <- unique(as.character(de$gene_id))
    string_db <- tryCatch(
      STRINGdb$new(version = "12", species = 39947, score_threshold = 400L,
                   input_directory = ""),
      error = function(e) NULL
    )
    if (!is.null(string_db)) {
      hits <- tryCatch({
        mapped <- string_db$map(data.frame(gene = de_genes), "gene",
                                removeUnmappedRows = TRUE)
        string_db$get_interactions(mapped$STRING_id)
      }, error = function(e) NULL)
      if (!is.null(hits) && nrow(hits) > 0L) {
        score_col <- if ("combined_score" %in% colnames(hits)) {
          "combined_score"
        } else if ("score" %in% colnames(hits)) {
          "score"
        } else {
          NULL
        }
        id2gene <- stats::setNames(mapped$gene, mapped$STRING_id)
        edges <- data.frame(
          gene_a = id2gene[as.character(hits$from)],
          gene_b = id2gene[as.character(hits$to)],
          score = if (!is.null(score_col)) as.numeric(hits[[score_col]]) else 400,
          stringsAsFactors = FALSE
        )
        edges <- edges[!is.na(edges$gene_a) & !is.na(edges$gene_b), , drop = FALSE]
      }
    }
    if (is.null(edges)) {
      message("STRINGdb query returned no interactions; the PPI module is skipped.")
    }
  } else if (!is.null(string_edges)) {
    edges <- read_csv_safe(string_edges, "STRING edges file")
    check_columns(edges, c("gene_a", "gene_b", "score"), "STRING edges file")
    edges$gene_a <- as.character(edges$gene_a)
    edges$gene_b <- as.character(edges$gene_b)
    edges$score <- suppressWarnings(as.numeric(edges$score))
    edges <- edges[!is.na(edges$score) & edges$score > 400, , drop = FALSE]
  } else {
    message("--ppi requested but neither the STRINGdb package nor --string-edges ",
            "is available; skipping the PPI module.\n",
            "  Hint: install.packages('STRINGdb') for the online STRING API, or ",
            "export STRING interactions (score > 400) to a CSV with columns ",
            "gene_a, gene_b, score and pass it via --string-edges.")
  }

  if (!is.null(edges) && nrow(edges) > 0L) {
    edges <- edges[edges$gene_a != edges$gene_b, , drop = FALSE]
    utils::write.csv(edges, file.path(outdir, "ppi_edges.csv"),
                     row.names = FALSE, quote = FALSE)
    degrees <- table(c(edges$gene_a, edges$gene_b))
    nodes <- data.frame(gene_id = names(degrees),
                        degree = as.integer(degrees), stringsAsFactors = FALSE)
    nodes <- merge(nodes, de[, c("gene_id", "log2FoldChange", "direction")],
                   by = "gene_id", all.x = TRUE, sort = FALSE)
    nodes <- nodes[order(-nodes$degree), , drop = FALSE]
    rownames(nodes) <- NULL
    utils::write.csv(nodes, file.path(outdir, "ppi_nodes.csv"),
                     row.names = FALSE, quote = FALSE)

    g <- igraph::graph_from_data_frame(edges[, c("gene_a", "gene_b")],
                                       directed = FALSE, vertices = nodes)
    node_cols <- ifelse(igraph::V(g)$direction == "Up", nature_semantic("up"),
                        ifelse(igraph::V(g)$direction == "Down",
                               nature_semantic("down"),
                               nature_semantic("not_significant")))
    if (requireNamespace("ggraph", quietly = TRUE)) {
      p <- ggraph::ggraph(g, layout = "fr") +
        ggraph::geom_edge_link(colour = nature_semantic("not_significant"),
                               linewidth = 0.3) +
        ggraph::geom_node_point(aes(colour = igraph::V(g)$direction), size = 3) +
        ggraph::geom_node_text(ggplot2::aes(label = igraph::V(g)$name),
                               size = 2.2, repel = TRUE) +
        ggplot2::scale_colour_manual(
          values = c(Up = nature_semantic("up"), Down = nature_semantic("down"),
                     NS = nature_semantic("not_significant")),
          na.value = nature_semantic("not_significant"), name = "DE direction") +
        ggplot2::labs(title = "PPI network of differentially expressed genes") +
        theme_nature()
      save_nature_figure(p, file.path(outdir, "ppi_network"),
                         width_mm = 120, height_mm = 120)
    } else {
      save_base_figure(
        function() {
          igraph::plot.igraph(g, layout = igraph::layout_with_fr(g),
                              vertex.color = node_cols, vertex.label = NA,
                              vertex.size = 4,
                              edge.color = nature_semantic("not_significant"),
                              main = "PPI network of differentially expressed genes")
          graphics::legend("topleft", legend = c("Up", "Down", "NS"), pch = 16,
                           col = c(nature_semantic("up"), nature_semantic("down"),
                                   nature_semantic("not_significant")),
                           title = "DE direction", cex = 0.7, bty = "n")
        },
        file.path(outdir, "ppi_network"), width_mm = 120, height_mm = 120
      )
    }
    message(sprintf("PPI network: %d nodes, %d edges (score > 400).",
                    igraph::vcount(g), igraph::ecount(g)))
  } else if (!is.null(string_edges)) {
    message("No edges with score > 400 in ", string_edges, "; PPI module skipped.")
  }
}

write_run_info_r(
  outdir, "06_tf_analysis.R",
  params = list(
    de_all = de_all_path,
    tf_annotation = tf_annotation,
    padj = padj_cut,
    log2fc = lfc_cut,
    ppi = ppi_flag,
    string_edges = if (is.null(string_edges)) "none" else string_edges
  ),
  packages = c("ggplot2", if (ppi_flag && requireNamespace("STRINGdb", quietly = TRUE))
    "STRINGdb" else character(0))
)
message("TF analysis complete. Results in ", outdir)
