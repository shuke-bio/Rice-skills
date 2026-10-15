###############################################################################
# 04_gsea.R -- gene set enrichment analysis (GSEA) with a custom GMT
#
# Title:    GSEA for the rice transcriptome (fgsea or clusterProfiler::GSEA)
# Description:
#   Builds a ranked gene list from a full DE table
#   (stat = sign(log2FC) * -log10(pvalue)), reads a user-provided GMT file
#   and runs gene set enrichment analysis with fgsea when installed, or with
#   clusterProfiler::GSEA otherwise (--method auto). Writes a unified result
#   table, enrichment curves for the top pathways and an NES summary bubble
#   plot.
#
# Usage:
#   Rscript scripts/transcriptomics/04_gsea.R \
#       --de-all results/transcriptome/02_de/de_all_salt_vs_control.csv \
#       --gmt data/rice_kegg_osa.gmt \
#       --outdir results/transcriptome/04_gsea
#   # force the clusterProfiler implementation:
#   Rscript scripts/transcriptomics/04_gsea.R --de-all ... --gmt ... --method clusterprofiler
#
# Inputs:
#   --de-all    full DE table (needs gene_id, log2FoldChange and pvalue columns)
#   --gmt       gene set file in GMT format (term <TAB> description <TAB> genes)
#
# Outputs (written to --outdir):
#   gsea_results.csv                  unified results (ID, set, NES, p, padj, ...)
#   gsea_curve_<pathway>.*            enrichment curves for the top pathways
#   gsea_nes_bubble.*                 NES vs -log10(padj) bubble plot
#   run_info.txt                      provenance (parameters, versions, seed)
#
# Notes:
#   - Rice is NOT covered by MSigDB. Build a custom GMT for rice, e.g. from
#     the KEGG REST API:
#       curl -s https://rest.kegg.jp/list/pathway/osa            # pathway list
#       curl -s https://rest.kegg.jp/link/osa/<pathway_id>       # pathway genes
#     (one GMT line per pathway: "osa00010<TAB>Glycolysis<TAB>Os01g...<TAB>..."),
#     or export gene sets from PlantTFDB (http://planttfdb.gao-lab.org/).
#     Keep the GMT ID system consistent with the DE table (rice KEGG gene IDs
#     drop the "LOC_" prefix, e.g. Os01g01010).
#   - The ranking statistic is sign(log2FC) * -log10(pvalue) with the p-value
#     floored at 1e-300, the standard GSEA input for two-group comparisons.
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
gmt_path <- opt(opts, "gmt")
id_type <- opt(opts, "id-type", "SYMBOL")
outdir <- opt(opts, "outdir", "results/transcriptome/04_gsea")
min_size <- suppressWarnings(as.integer(opt(opts, "min-size", "10")))
max_size <- suppressWarnings(as.integer(opt(opts, "max-size", "500")))
padj_cut <- suppressWarnings(as.numeric(opt(opts, "padj", "0.05")))
n_top <- suppressWarnings(as.integer(opt(opts, "n-top", "20")))
method <- tolower(opt(opts, "method", "auto"))

if (is.null(de_all_path)) {
  stop("--de-all is required (a de_all_<contrast>.csv produced by 02_deseq2.R).",
       call. = FALSE)
}
if (is.null(gmt_path)) {
  stop("--gmt is required (gene sets in GMT format).\n",
       "Hint: rice has no MSigDB; build a GMT from KEGG REST ",
       "(https://rest.kegg.jp/list/pathway/osa + /link/osa/<id>) or export ",
       "gene sets from PlantTFDB.", call. = FALSE)
}
if (!method %in% c("auto", "fgsea", "clusterprofiler")) {
  stop("Invalid --method '", method, "'. Supported: auto, fgsea, clusterprofiler.",
       call. = FALSE)
}
if (is.na(min_size) || min_size < 1L) {
  stop("--min-size must be a positive integer; got '", opt(opts, "min-size", "10"), "'.",
       call. = FALSE)
}
if (is.na(max_size) || max_size < min_size) {
  stop("--max-size must be an integer >= --min-size; got '",
       opt(opts, "max-size", "500"), "'.", call. = FALSE)
}
if (is.na(padj_cut) || padj_cut <= 0 || padj_cut >= 1) {
  stop("--padj must be a number in (0, 1); got '", opt(opts, "padj", "0.05"), "'.", call. = FALSE)
}
if (is.na(n_top) || n_top < 1L) {
  stop("--n-top must be a positive integer; got '", opt(opts, "n-top", "20"), "'.",
       call. = FALSE)
}

use_fgsea <- method == "fgsea" ||
  (method == "auto" && requireNamespace("fgsea", quietly = TRUE))
if (method == "fgsea" && !requireNamespace("fgsea", quietly = TRUE)) {
  stop("--method fgsea was requested but the fgsea package is not installed.",
       call. = FALSE)
}

need_pkgs(c("clusterProfiler", "ggplot2", "jsonlite"))
if (use_fgsea) need_pkgs("fgsea")

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

#' Pick the first available column from a result data frame.
pick_col <- function(df, candidates, na = NA) {
  for (cn in candidates) {
    if (cn %in% colnames(df)) return(df[[cn]])
  }
  rep(na, nrow(df))
}

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

de <- read_csv_safe(de_all_path, "DE results table")
check_columns(de, c("gene_id", "log2FoldChange", "pvalue"), "de-all table")
de <- de[!is.na(de$pvalue) & !is.na(de$log2FoldChange), , drop = FALSE]
if (nrow(de) == 0L) {
  stop("No rows with both a pvalue and a log2FoldChange remain in ", de_all_path,
       call. = FALSE)
}
# Ranking statistic: sign(log2FC) * -log10(p), p floored at 1e-300.
p_floor <- pmax(as.numeric(de$pvalue), 1e-300)
stat <- sign(as.numeric(de$log2FoldChange)) * (-log10(p_floor))
names(stat) <- as.character(de$gene_id)
stat <- stat[!duplicated(names(stat))]
gene_list <- sort(stat, decreasing = TRUE)
message(sprintf("Ranked gene list: %d genes.", length(gene_list)))

# GMT: clusterProfiler::read.gmt returns a term/gene data.frame.
gmt <- clusterProfiler::read.gmt(gmt_path)
if (!is.data.frame(gmt) || !all(c("term", "gene") %in% colnames(gmt))) {
  stop("Could not parse the GMT file: ", gmt_path,
       "\nHint: expected tab-separated lines 'term <TAB> description <TAB> gene...'.",
       call. = FALSE)
}
pathways_list <- split(as.character(gmt$gene), as.character(gmt$term))
pathways_list <- pathways_list[vapply(pathways_list, length, integer(1)) > 0L]
if (length(pathways_list) == 0L) {
  stop("No gene sets found in the GMT file: ", gmt_path, call. = FALSE)
}
gmt_genes <- unique(unlist(pathways_list))
gene_list_used <- gene_list
n_match <- sum(gmt_genes %in% names(gene_list))
if (n_match == 0L && all(grepl("^LOC_", names(gene_list)))) {
  # Rice KEGG gene IDs drop the "LOC_" prefix; retry with stripped IDs.
  stripped <- sub("^LOC_", "", names(gene_list))
  n_match2 <- sum(gmt_genes %in% stripped)
  if (n_match2 > 0L) {
    message("GMT gene IDs use the KEGG-style naming (no 'LOC_' prefix); ",
            "matching after stripping the prefix from the DE gene IDs.")
    names(gene_list_used) <- stripped
    n_match <- n_match2
  }
}
if (n_match == 0L) {
  stop("No GMT gene ID overlaps the DE table gene IDs.\n",
       "Hint: the GMT and --de-all must use the same ID system (--id-type is ",
       id_type, "); rice KEGG gene IDs drop the 'LOC_' prefix (Os01g01010).",
       call. = FALSE)
}
message(sprintf("GMT: %d gene sets, %d/%d GMT genes present in the DE table.",
                length(pathways_list), n_match, length(gmt_genes)))

gsea_obj <- NULL
if (use_fgsea) {
  message("Running GSEA with fgsea.")
  fgsea_res <- fgsea::fgsea(pathways = pathways_list, stats = gene_list_used,
                            minSize = min_size, maxSize = max_size,
                            nperm = 1000, nproc = 1)
  res_df <- data.frame(
    ID = as.character(fgsea_res$pathway),
    set = as.character(fgsea_res$pathway),
    NES = as.numeric(fgsea_res$NES),
    pvalue = as.numeric(fgsea_res$pval),
    padj = as.numeric(fgsea_res$padj),
    size = as.integer(fgsea_res$size),
    core_enrichment = vapply(fgsea_res$leadingEdge,
                             function(x) paste(x, collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
} else {
  message("Running GSEA with clusterProfiler::GSEA (1000 permutations).")
  # pvalueCutoff = 1 on purpose: clusterProfiler DROPS every gene set above the
  # cutoff from the returned object, which would hide the full ranked table and
  # make this step fail outright on datasets without strong enrichment signals.
  # Significance is flagged in the `significant` column of the output instead.
  gsea_obj <- clusterProfiler::GSEA(
    geneList = gene_list_used, TERM2GENE = gmt, exponent = 1, nPerm = 1000,
    minGSSize = min_size, maxGSSize = max_size, pvalueCutoff = 1,
    pAdjustMethod = "BH", verbose = FALSE
  )
  if (is.null(gsea_obj)) {
    stop("GSEA returned no result.\n",
         "Hint: check that the GMT gene IDs match the DE table IDs ",
         "(rice KEGG gene IDs drop the 'LOC_' prefix) and that gene sets are ",
         "larger than --min-size (", min_size, ").", call. = FALSE)
  }
  r <- as.data.frame(gsea_obj@result)
  res_df <- data.frame(
    ID = as.character(pick_col(r, c("ID", "pathway"))),
    set = as.character(pick_col(r, c("Description", "set", "pathway"))),
    NES = as.numeric(pick_col(r, "NES")),
    pvalue = as.numeric(pick_col(r, c("pvalue", "p", "pval"))),
    padj = as.numeric(pick_col(r, c("p.adjust", "padj"))),
    size = as.integer(pick_col(r, c("setSize", "size"))),
    core_enrichment = as.character(pick_col(r, c("core_enrichment", "leading_edge"))),
    stringsAsFactors = FALSE
  )
}
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), , drop = FALSE]
rownames(res_df) <- NULL
res_df$significant <- !is.na(res_df$padj) & res_df$padj <= padj_cut
utils::write.csv(res_df, file.path(outdir, "gsea_results.csv"),
                 row.names = FALSE, quote = FALSE)
message(sprintf("GSEA finished: %d gene sets tested, %d with padj <= %g.",
                nrow(res_df), sum(res_df$significant), padj_cut))
if (!any(res_df$significant)) {
  message("Note: no gene set reached the significance cutoff; the complete ",
          "ranked table in gsea_results.csv is still written for inspection.")
}

# Enrichment curves for the top pathways (by adjusted p-value).
top_ids <- head(res_df$ID[!is.na(res_df$padj)], n_top)
if (length(top_ids) == 0L) top_ids <- head(res_df$ID, n_top)
for (id in top_ids) {
  safe <- gsub("[^A-Za-z0-9._-]", "_", id)
  stem <- file.path(outdir, paste0("gsea_curve_", safe))
  if (use_fgsea) {
    p <- fgsea::plotEnrichment(pathways_list[[id]], gene_list_used) +
      ggplot2::labs(title = id)
    save_nature_figure(p, stem, width_mm = 89, height_mm = 70)
  } else {
    # enrichplot::gseaplot2 returns a gglist (stacked panels); print() draws it.
    nes_val <- res_df$NES[match(id, res_df$ID)]
    curve_col <- if (!is.na(nes_val) && nes_val < 0) {
      nature_semantic("down")
    } else {
      nature_semantic("up")
    }
    p <- enrichplot::gseaplot2(gsea_obj, geneSetID = id, title = id,
                               color = curve_col)
    save_base_figure(function() print(p), stem, width_mm = 89, height_mm = 70)
  }
}

# NES summary bubble plot.
# Scientific conclusion: pathways in the upper-right quadrant (NES > 0,
# padj <= cutoff) are activated by the treatment, lower-left ones repressed;
# bubble size reflects the gene set size, so small significant sets should be
# interpreted with caution.
plot_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$NES), , drop = FALSE]
if (nrow(plot_df) > 60L) {
  n_all <- nrow(plot_df)
  plot_df <- head(plot_df, 60L)
  message("Note: NES bubble plot shows the top 60 of ", n_all,
          " gene sets (ranked by padj).")
}
if (nrow(plot_df) > 0L) {
  plot_df$direction <- ifelse(plot_df$NES > 0, "Activated", "Repressed")
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = NES, y = -log10(padj), size = size, colour = direction)
  ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_colour_manual(
      values = c(Activated = nature_semantic("up"),
                 Repressed = nature_semantic("down")),
      name = "Direction"
    ) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dashed",
                        colour = nature_semantic("reference_line")) +
    ggplot2::labs(x = "Normalised enrichment score (NES)",
                  y = "-log10 (adjusted p-value)",
                  title = "GSEA pathway summary", size = "Gene set size") +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "gsea_nes_bubble"),
                     width_mm = 120, height_mm = 89)
} else {
  message("No pathway with a valid padj; skipping the NES bubble plot.")
}

write_run_info_r(
  outdir, "04_gsea.R",
  params = list(
    de_all = de_all_path,
    gmt = gmt_path,
    id_type = id_type,
    method = if (use_fgsea) "fgsea" else "clusterprofiler",
    min_size = min_size,
    max_size = max_size,
    padj = padj_cut,
    n_top = n_top
  ),
  packages = c("clusterProfiler", if (use_fgsea) "fgsea" else character(0))
)
message("GSEA complete. Results in ", outdir)
