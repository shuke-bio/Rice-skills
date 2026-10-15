###############################################################################
# 03_enrichment.R -- GO and KEGG enrichment analysis with clusterProfiler
#
# Title:    GO / KEGG over-representation analysis for differential genes
# Description:
#   Reads a significant DE table (produced by 02_deseq2.R), splits the genes
#   into UP and DOWN directions by the sign of log2FoldChange, maps the gene
#   IDs to ENTREZID via the rice OrgDb, and runs GO enrichment (BP/CC/MF)
#   plus KEGG enrichment for each direction. Results are written as CSV
#   tables and dot plots, and a combined overview table (top 10 terms per
#   category/direction) is generated.
#
# Usage:
#   Rscript scripts/transcriptomics/03_enrichment.R \
#       --de-sig results/transcriptome/02_de/de_sig_salt_vs_control.csv \
#       --contrast-name salt_vs_control \
#       --de-all results/transcriptome/02_de/de_all_salt_vs_control.csv \
#       --outdir results/transcriptome/03_enrichment
#
# Inputs:
#   --de-sig     significant DE table (needs gene_id and log2FoldChange columns)
#   --de-all     optional full DE table (gene_id column used as the background)
#   --universe   optional one-column file with background gene IDs
#
# Outputs (written to --outdir):
#   go_<ONT>_<UP|DOWN>.csv          GO results per ontology and direction
#   kegg_<UP|DOWN>.csv              KEGG results per direction
#   go_<ONT>_<UP|DOWN>_dotplot.*    GO dot plots (PDF/SVG/TIFF)
#   kegg_<UP|DOWN>_dotplot.*        KEGG dot plots (PDF/SVG/TIFF)
#   enrichment_overview.csv         top 10 terms per category/direction
#   run_info.txt                    provenance (parameters, versions, seed)
#
# Notes:
#   - ID system: rice LOC_Os IDs (e.g. LOC_Os01g01010) are SYMBOL keys in
#     org.Osativa.eg.db, so the default --id-type is SYMBOL. If the IDs are
#     already ENTREZIDs (e.g. Os01g01010), pass --id-type ENTREZID. KEGG gene
#     IDs drop the "LOC_" prefix (Os01g01010); this conversion is applied
#     automatically for the KEGG branch.
#   - Reactome is NOT available for rice: ReactomePA only covers human, mouse
#     and a few model organisms, so this script deliberately offers GO + KEGG
#     only.
#   - KEGG enrichment requires internet access (KEGG REST API).
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

de_sig_path <- opt(opts, "de-sig")
contrast_name <- opt(opts, "contrast-name")
de_all_path <- opt(opts, "de-all")
universe_path <- opt(opts, "universe")
outdir <- opt(opts, "outdir", "results/transcriptome/03_enrichment")
id_type <- opt(opts, "id-type", "SYMBOL")
org_arg <- tolower(opt(opts, "org", "osa"))
ontologies <- toupper(trimws(strsplit(opt(opts, "ontologies", "BP,CC,MF"),
                                      ",", fixed = TRUE)[[1]]))
ontologies <- ontologies[nzchar(ontologies)]
padj_cut <- suppressWarnings(as.numeric(opt(opts, "padj", "0.05")))
simplify_flag <- isTRUE(opt(opts, "simplify", NULL))

if (is.null(de_sig_path)) {
  stop("--de-sig is required (a de_sig_<contrast>.csv produced by 02_deseq2.R).",
       call. = FALSE)
}
if (is.null(contrast_name)) {
  contrast_name <- sub("^de_sig_", "", sub("\\.csv$", "", basename(de_sig_path)))
  if (!nzchar(contrast_name)) contrast_name <- "contrast"
}
if (is.na(padj_cut) || padj_cut <= 0 || padj_cut >= 1) {
  stop("--padj must be a number in (0, 1); got '", opt(opts, "padj", "0.05"), "'.", call. = FALSE)
}
valid_ont <- c("BP", "CC", "MF")
bad_ont <- setdiff(ontologies, valid_ont)
if (length(bad_ont) > 0L) {
  stop("Unknown ontology/ontologies: ", paste(bad_ont, collapse = ", "),
       "\nSupported values: BP, CC, MF (comma-separated).", call. = FALSE)
}

# Species mapping (edit here to migrate to another organism, see SKILL.md).
org_map <- c(osa = "org.Osativa.eg.db",
             osativa = "org.Osativa.eg.db",
             japonica = "org.Osativa.eg.db")
org_pkg <- unname(org_map[org_arg])
if (is.na(org_pkg)) {
  stop("Unsupported --org '", org_arg, "'. This rice-focused build supports ",
       "'osa' (org.Osativa.eg.db).\n",
       "Hint: for another species install the matching OrgDb and extend org_map.",
       call. = FALSE)
}

need_pkgs(c("clusterProfiler", org_pkg, "ggplot2", "jsonlite"))

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(org_pkg, character.only = TRUE))
org_db <- get(org_pkg)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Map gene IDs to ENTREZID with the species OrgDb; NULL when nothing maps.
map_to_entrez <- function(genes, id_type, org_db) {
  genes <- unique(as.character(genes))
  if (length(genes) == 0L) return(NULL)
  mapped <- tryCatch(
    clusterProfiler::bitr(genes, fromType = id_type, toType = "ENTREZID",
                          OrgDb = org_db),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) == 0L) return(NULL)
  unique(as.character(mapped$ENTREZID))
}

#' GO over-representation test (one ontology, one direction).
run_go <- function(entrez, universe_entrez, ont, org_db, padj_cut) {
  if (length(entrez) == 0L) return(NULL)
  tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez, OrgDb = org_db, keyType = "ENTREZID", ont = ont,
      universe = if (length(universe_entrez) > 0L) universe_entrez else NULL,
      pvalueCutoff = padj_cut, qvalueCutoff = padj_cut, readable = TRUE
    ),
    error = function(e) NULL
  )
}

#' Rice KEGG gene IDs drop the "LOC_" prefix (LOC_Os01g01010 -> Os01g01010).
kegg_gene_id <- function(symbols) sub("^LOC_", "", as.character(symbols))

#' KEGG over-representation test (one direction).
run_kegg <- function(symbols, universe_symbols, org, padj_cut) {
  genes <- kegg_gene_id(symbols)
  if (length(genes) == 0L) return(NULL)
  universe <- if (length(universe_symbols) > 0L) kegg_gene_id(universe_symbols) else NULL
  tryCatch(
    clusterProfiler::enrichKEGG(
      gene = genes, organism = org, keyType = "kegg", universe = universe,
      pvalueCutoff = padj_cut, qvalueCutoff = padj_cut
    ),
    error = function(e) NULL
  )
}

#' Extract the result table of an enrichResult object (NULL when empty).
result_table <- function(res) {
  if (is.null(res)) return(NULL)
  df <- tryCatch(as.data.frame(res@result), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df
}

#' Write an enrichment table; writes a header-only CSV when empty so that the
#' downstream plotting script can detect the (empty) analysis.
write_enrich_csv <- function(df, path) {
  cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
            "geneID", "Count")
  if (is.null(df)) {
    empty <- setNames(
      data.frame(matrix(numeric(0), ncol = length(cols)), stringsAsFactors = FALSE),
      cols
    )
    utils::write.csv(empty, path, row.names = FALSE, quote = FALSE)
  } else {
    utils::write.csv(df, path, row.names = FALSE, quote = FALSE)
  }
}

#' Top-10 rows of a result table for the combined overview.
make_overview <- function(df, category, direction) {
  df <- df[order(df$p.adjust, na.last = TRUE), , drop = FALSE]
  df <- head(df, 10)
  data.frame(category = category, direction = direction,
             rank = seq_len(nrow(df)), ID = df$ID, Description = df$Description,
             GeneRatio = df$GeneRatio, BgRatio = df$BgRatio, pvalue = df$pvalue,
             p.adjust = df$p.adjust, Count = df$Count, geneID = df$geneID,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

de <- read_csv_safe(de_sig_path, "DE significant table")
check_columns(de, c("gene_id", "log2FoldChange"), "de-sig table")
de$gene_id <- as.character(de$gene_id)
de$log2FoldChange <- as.numeric(de$log2FoldChange)
up_genes <- de$gene_id[de$log2FoldChange > 0]
down_genes <- de$gene_id[de$log2FoldChange < 0]
message(sprintf("Contrast '%s': %d UP and %d DOWN genes.", contrast_name,
                length(up_genes), length(down_genes)))

# Background gene universe: --universe file > --de-all table > OrgDb default.
universe_symbols <- NULL
if (!is.null(universe_path)) {
  u <- read_csv_safe(universe_path, "universe file")
  universe_symbols <- as.character(u[[1]])
} else if (!is.null(de_all_path)) {
  de_all <- read_csv_safe(de_all_path, "DE all table")
  check_columns(de_all, "gene_id", "de-all table")
  universe_symbols <- as.character(de_all$gene_id)
}
universe_entrez <- if (!is.null(universe_symbols)) {
  map_to_entrez(universe_symbols, id_type, org_db)
} else {
  NULL
}
if (!is.null(universe_symbols) && length(universe_entrez) == 0L) {
  message("WARNING: none of the background genes could be mapped; ",
          "falling back to the OrgDb-wide universe.")
  universe_entrez <- NULL
}

overview_rows <- list()
directions <- list(UP = up_genes, DOWN = down_genes)

for (dir_tag in names(directions)) {
  genes <- directions[[dir_tag]]
  entrez <- map_to_entrez(genes, id_type, org_db)
  if (length(entrez) == 0L) {
    message("Skipping the ", dir_tag, " direction: no ", id_type,
            " IDs could be mapped to ENTREZID.\n",
            "  Hint: LOC_Os IDs are SYMBOL keys in ", org_pkg,
            "; if your IDs are already ENTREZIDs (e.g. Os01g01010), pass ",
            "--id-type ENTREZID, or convert the IDs beforehand.")
    next
  }
  message(sprintf("Direction %s: %d genes -> %d ENTREZIDs.", dir_tag,
                  length(genes), length(entrez)))

  for (ont in ontologies) {
    res <- run_go(entrez, universe_entrez, ont, org_db, padj_cut)
    if (simplify_flag && !is.null(res) && nrow(res@result) > 1L) {
      res <- tryCatch(
        clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust",
                                  select_fun = min),
        error = function(e) {
          message("simplify() failed for GO ", ont, " ", dir_tag, ": ",
                  conditionMessage(e))
          res
        }
      )
    }
    df <- result_table(res)
    write_enrich_csv(df, file.path(outdir, paste0("go_", ont, "_", dir_tag, ".csv")))
    if (!is.null(df)) {
      p <- clusterProfiler::dotplot(res, showCategory = 10)
      save_nature_figure(
        p, file.path(outdir, paste0("go_", ont, "_", dir_tag, "_dotplot")),
        width_mm = 120, height_mm = max(60, min(nrow(df), 10) * 7 + 30)
      )
      overview_rows[[length(overview_rows) + 1L]] <-
        make_overview(df, paste0("GO_", ont), dir_tag)
    } else {
      message("No GO ", ont, " enrichment for the ", dir_tag,
              " direction at padj <= ", padj_cut, ".")
    }
  }

  # KEGG (online; requires internet access).
  res <- run_kegg(genes, universe_symbols, org_arg, padj_cut)
  df <- result_table(res)
  write_enrich_csv(df, file.path(outdir, paste0("kegg_", dir_tag, ".csv")))
  if (!is.null(df)) {
    p <- clusterProfiler::dotplot(res, showCategory = 10)
    save_nature_figure(
      p, file.path(outdir, paste0("kegg_", dir_tag, "_dotplot")),
      width_mm = 120, height_mm = max(60, min(nrow(df), 10) * 7 + 30)
    )
    overview_rows[[length(overview_rows) + 1L]] <-
      make_overview(df, "KEGG", dir_tag)
  } else {
    message("No KEGG enrichment for the ", dir_tag,
            " direction (check internet access and the KEGG ID mapping).")
  }
}

if (length(overview_rows) > 0L) {
  overview <- do.call(rbind, overview_rows)
  utils::write.csv(overview, file.path(outdir, "enrichment_overview.csv"),
                   row.names = FALSE, quote = FALSE)
} else {
  write_enrich_csv(NULL, file.path(outdir, "enrichment_overview.csv"))
  message("No enriched terms found in any category/direction.")
}

write_run_info_r(
  outdir, "03_enrichment.R",
  params = list(
    de_sig = de_sig_path,
    contrast = contrast_name,
    de_all = if (is.null(de_all_path)) "none" else de_all_path,
    universe = if (is.null(universe_path)) "none" else universe_path,
    id_type = id_type,
    org = org_arg,
    ontologies = paste(ontologies, collapse = ","),
    padj = padj_cut,
    simplify = simplify_flag
  ),
  packages = c("clusterProfiler", org_pkg)
)
message("Enrichment complete. Results in ", outdir)
