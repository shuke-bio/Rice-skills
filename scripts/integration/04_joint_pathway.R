###############################################################################
# 04_joint_pathway.R -- joint transcriptome x metabolome KEGG pathway analysis
#                       (Rice-skills skill, module 3, step 04)
#
# Title:    Joint KEGG pathway enrichment of DE genes and differential metabolites
# Description:
#   Performs KEGG pathway enrichment on the two omics layers separately and
#   then joins them at the pathway level:
#     * Transcriptome: DE genes (split into up- and down-regulated) are
#       converted to ENTREZID with clusterProfiler::bitr and tested with
#       clusterProfiler::enrichKEGG (organism = --org, default osa).
#     * Metabolome: differential metabolites carrying a kegg_id annotation are
#       tested for pathway over-representation with a hypergeometric test.
#       Pathway-compound relationships are retrieved from KEGG via KEGGREST;
#       if KEGGREST is unavailable (or offline), pass the output of
#       metabolomics/03_pathway.R directly with --metabolite-pathway-csv.
#   Pathways significant on both sides (category = "both") are the joint
#   response of the system: e.g. genes and metabolites of the same
#   biosynthetic pathway moving in a coordinated direction.
#
# Usage:
#   Rscript scripts/integration/04_joint_pathway.R \
#       --de-genes results/transcriptome/02_de/de_sig_control_vs_salt.csv \
#       --de-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
#       --outdir results/integration/04_joint_pathway \
#       --org osa --padj 0.05 --id-type SYMBOL
#
#   # offline / KEGGREST-free mode (reuse the metabolome pathway step):
#   Rscript scripts/integration/04_joint_pathway.R ... \
#       --metabolite-pathway-csv results/metabolome/03_pathway/pathway_enrichment.csv
#
# Inputs:
#   --de-genes       DE gene table; required columns: gene_id, log2FoldChange
#   --de-metabolites differential metabolite table; required column: metabolite_id
#                    (a kegg_id column is required for the KEGGREST path)
#   --org            KEGG organism code (default "osa" = Oryza sativa japonica)
#   --padj           adjusted p-value threshold (default 0.05)
#   --id-type        ID type of the gene_id column for bitr (default SYMBOL)
#   --metabolite-pathway-csv  optional pre-computed metabolome pathway table
#
# Outputs (written to --outdir):
#   gene_kegg_enrichment.csv        transcriptome enrichment (up/down merged)
#   metabolite_pathway_enrichment.csv  metabolome enrichment results
#   joint_pathway_table.csv         joined table with category labels
#   joint_pathway_dotplot.*         dot plot of joint pathways (PDF/SVG/TIFF)
#   pathway_links.csv               node-pathway links (input for the sankey in 06)
#   run_info.txt
#
# Notes:
#   - clusterProfiler and org.Osativa.eg.db require internet access to the
#     KEGG/GO databases; the KEGGREST path additionally queries KEGG online.
#   - The gene enrichment background is the full KEGG osa gene set (the DE
#     table does not carry the tested universe); interpret p-values with this
#     in mind.
#   - Category definitions: "both" = significant on both sides;
#     "transcriptome_only" / "metabolome_only" = significant on one side.
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

usage_04 <- paste(
  "Usage: Rscript 04_joint_pathway.R --de-genes <de_sig.csv>",
  "         --de-metabolites <diff.csv> [--outdir DIR] [--org osa] [--padj 0.05]",
  "         [--id-type SYMBOL] [--metabolite-pathway-csv <pathway.csv>]",
  sep = "\n"
)

opts <- parse_args(commandArgs(trailingOnly = TRUE))
de_genes_path <- opt(opts, "de-genes", NULL)
de_mets_path <- opt(opts, "de-metabolites", NULL)
if (is.null(de_genes_path) || is.null(de_mets_path)) {
  stop("Missing required argument(s): --de-genes and --de-metabolites are required.\n",
       usage_04, call. = FALSE)
}
outdir <- opt(opts, "outdir", "results/integration/04_joint_pathway")
org <- opt(opts, "org", "osa")
padj_cutoff <- as.numeric(opt(opts, "padj", "0.05"))
id_type <- opt(opts, "id-type", "SYMBOL")
met_pathway_csv <- opt(opts, "metabolite-pathway-csv", NULL)
if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
  stop("--padj must be in (0, 1], e.g. 0.05.", call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rice_set_seed(42)

need_pkgs(c("clusterProfiler", "org.Osativa.eg.db"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Normalise KEGG pathway identifiers to the bare organism code form
#' ("path:osa00010" / "map00010" -> "osa00010").
norm_pathway_id <- function(x) {
  x <- as.character(x)
  x <- sub("^path:", "", x)
  x <- sub("^map:", "", x)
  tolower(trimws(x))
}

#' Pick the first present column from a list of candidate names.
pick_col <- function(df, candidates, what, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    if (required) {
      stop("Could not find a '", what, "' column in the pathway table.\n",
           "Tried: ", paste(candidates, collapse = ", "), "\n",
           "Available columns: ", paste(colnames(df), collapse = ", "), call. = FALSE)
    }
    return(NULL)
  }
  df[[hit[1]]]
}

# ---------------------------------------------------------------------------
# Load DE tables
# ---------------------------------------------------------------------------

de_genes <- read_csv_safe(de_genes_path, "DE genes table")
check_columns(de_genes, c("gene_id", "log2FoldChange"), "DE genes table")
de_genes$gene_id <- as.character(de_genes$gene_id)
de_genes$log2FoldChange <- suppressWarnings(as.numeric(de_genes$log2FoldChange))

de_mets <- read_csv_safe(de_mets_path, "Differential metabolites table")
check_columns(de_mets, "metabolite_id", "Differential metabolites table")
de_mets$metabolite_id <- as.character(de_mets$metabolite_id)

gene_up <- de_genes$gene_id[de_genes$log2FoldChange > 0]
gene_down <- de_genes$gene_id[de_genes$log2FoldChange < 0]
message("DE genes: ", length(gene_up), " up, ", length(gene_down), " down, ",
        sum(is.na(de_genes$log2FoldChange)), " without log2FoldChange.")

# ---------------------------------------------------------------------------
# Transcriptome side: bitr + enrichKEGG (up/down separately)
# ---------------------------------------------------------------------------

OrgDb <- "org.Osativa.eg.db"
convert_genes <- function(gene_ids, direction) {
  gene_ids <- unique(gene_ids[!is.na(gene_ids)])
  if (length(gene_ids) < 5) {
    message("Note: only ", length(gene_ids), " ", direction,
            "-regulated genes; enrichment needs at least 5, skipping this set.")
    return(NULL)
  }
  bitr_res <- tryCatch(
    clusterProfiler::bitr(gene_ids, fromType = id_type, toType = "ENTREZID",
                          OrgDb = OrgDb, drop = TRUE),
    error = function(e) {
      stop("Gene ID conversion (bitr) failed for the ", direction, " set: ",
           conditionMessage(e),
           "\nHint: check --id-type (currently '", id_type,
           "'); valid types are the keytypes of org.Osativa.eg.db ",
           "(e.g. SYMBOL, ENTREZID, ALIAS).", call. = FALSE)
    }
  )
  if (is.null(bitr_res) || nrow(bitr_res) < 5) {
    message("Note: fewer than 5 ", direction, " genes mapped to ENTREZID; skipping.")
    return(NULL)
  }
  kk <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = bitr_res$ENTREZID, organism = org, keyType = "ncbi-geneid",
      pvalueCutoff = padj_cutoff, pAdjustMethod = "BH", qvalueCutoff = padj_cutoff,
      minGSSize = 5, maxGSSize = 500),
    error = function(e) {
      stop("enrichKEGG failed for the ", direction, " set: ", conditionMessage(e),
           "\nHint: KEGG enrichment needs internet access to rest.kegg.jp; ",
           "check the connection or retry later.", call. = FALSE)
    }
  )
  if (is.null(kk) || nrow(kk) == 0) {
    message("Note: no significant KEGG pathways for ", direction, "-regulated genes.")
    return(NULL)
  }
  res <- as.data.frame(kk)
  res$pathway_id <- norm_pathway_id(res$ID)
  res$direction <- direction
  res
}

gene_enr_up <- convert_genes(gene_up, "up")
gene_enr_down <- convert_genes(gene_down, "down")

gene_enr <- do.call(rbind, Filter(Negate(is.null), list(gene_enr_up, gene_enr_down)))
if (is.null(gene_enr)) {
  message("Note: no transcriptome enrichment results; the joint table will only ",
          "contain metabolome-only pathways (if any).")
}

if (!is.null(gene_enr)) {
  gene_enr$gene_ids <- vapply(gene_enr$geneID, function(s) {
    g <- trimws(unlist(strsplit(s, "/", fixed = TRUE)))
    paste(unique(g), collapse = ";")
  }, character(1))
  gene_enr_out <- gene_enr[, c("pathway_id", "ID", "Description", "direction",
                               "GeneRatio", "BgRatio", "pvalue", "p.adjust",
                               "Count", "gene_ids")]
  names(gene_enr_out)[names(gene_enr_out) == "p.adjust"] <- "padj"
  write.csv(gene_enr_out, file.path(outdir, "gene_kegg_enrichment.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Metabolome side: hypergeometric test against KEGG pathway-compound links
# ---------------------------------------------------------------------------

if (!is.null(met_pathway_csv)) {
  # --- pre-computed metabolome pathway table (e.g. metabolomics/03_pathway.R)
  mp <- read_csv_safe(met_pathway_csv, "Metabolome pathway table")
  pid <- pick_col(mp, c("pathway_id", "ID", "pathway", "KEGG_pathway", "kegg_pathway"),
                  "pathway id")
  pname <- pick_col(mp, c("name", "Description", "description", "pathway_name"),
                    "pathway name", required = FALSE)
  padj_col <- pick_col(mp, c("p.adjust", "padj", "p_adjust", "qvalue", "q.value", "fdr",
                             "P.adjust", "pAdjust"), "adjusted p-value")
  pval_col <- pick_col(mp, c("pvalue", "p.value", "PValue", "Pvalue", "p"),
                       "raw p-value", required = FALSE)
  cpd_col <- pick_col(mp, c("compound_ids", "metabolite_ids", "metabolites", "compounds",
                            "metabolite_id"), "metabolite list", required = FALSE)

  met_enr <- data.frame(
    pathway_id = norm_pathway_id(pid),
    name = if (is.null(pname)) norm_pathway_id(pid) else as.character(pname),
    p = if (is.null(pval_col)) NA_real_ else suppressWarnings(as.numeric(pval_col)),
    padj = suppressWarnings(as.numeric(padj_col)),
    compound_ids = if (is.null(cpd_col)) "" else as.character(cpd_col),
    stringsAsFactors = FALSE
  )
  met_enr <- met_enr[!is.na(met_enr$padj), , drop = FALSE]
  met_enr <- met_enr[order(met_enr$padj), , drop = FALSE]
  message("Loaded ", nrow(met_enr), " metabolite pathway rows from ",
          met_pathway_csv, " (", sum(met_enr$padj <= padj_cutoff), " significant).")
} else {
  # --- KEGGREST path: retrieve osa pathway-compound relationships from KEGG
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("The metabolite pathway analysis needs either KEGGREST (not installed) or ",
         "the pre-computed metabolome pathway table.\n",
         "Hint: install KEGGREST with:\n  install.packages(\"KEGGREST\")\n",
         "or run scripts/metabolomics/03_pathway.R and pass its output with ",
         "--metabolite-pathway-csv <file>.", call. = FALSE)
  }
  if (!"kegg_id" %in% colnames(de_mets)) {
    stop("The differential metabolite table has no 'kegg_id' column, so KEGG ",
         "compound-pathway relationships cannot be used.\n",
         "Hint: add a kegg_id column (e.g. C00031) to ", de_mets_path,
         ", or run scripts/metabolomics/03_pathway.R and pass its output with ",
         "--metabolite-pathway-csv <file>.", call. = FALSE)
  }
  diff_kegg <- unique(trimws(as.character(de_mets$kegg_id)))
  diff_kegg <- diff_kegg[nzchar(diff_kegg) & !is.na(diff_kegg)]
  diff_kegg <- sub("^cpd:", "", diff_kegg)
  if (length(diff_kegg) < 3) {
    stop("Only ", length(diff_kegg), " differential metabolite(s) carry a kegg_id; ",
         "at least 3 are needed for a pathway test.\n",
         "Hint: annotate the metabolites with KEGG compound IDs, or use ",
         "--metabolite-pathway-csv with the output of metabolomics/03_pathway.R.",
         call. = FALSE)
  }

  links <- tryCatch(
    KEGGREST::keggLink("pathway", "compound"),
    error = function(e) {
      stop("KEGGREST::keggLink('pathway', 'compound') failed: ", conditionMessage(e),
           "\nHint: this query needs internet access to rest.kegg.jp. ",
           "Offline? run scripts/metabolomics/03_pathway.R first and pass its ",
           "output with --metabolite-pathway-csv <file>.", call. = FALSE)
    }
  )
  ids <- names(links)
  vals <- unname(links)
  # Orient the pairs robustly (KEGGREST returns source\target; be defensive).
  if (any(grepl("^path:", ids)) && !any(grepl("^path:", vals))) {
    path_ids <- ids; cpd_ids <- vals
  } else if (any(grepl("^path:", vals))) {
    path_ids <- vals; cpd_ids <- ids
  } else {
    stop("Unexpected KEGGREST keggLink output (no 'path:' identifiers found).",
         call. = FALSE)
  }
  keep <- grepl(paste0("^path:", org), path_ids)
  if (!any(keep)) {
    stop("No KEGG pathways found for organism '", org,
         "'. Check the KEGG organism code (rice = 'osa').", call. = FALSE)
  }
  cpd2path <- split(gsub("^cpd:", "", cpd_ids[keep]),
                    gsub("^path:", "", path_ids[keep]))

  # Pathway display names (best effort; IDs are used when unavailable).
  path_names <- tryCatch({
    lst <- KEGGREST::keggList("pathway", org)
    setNames(sub(paste0(" - ", org, ".*$"), "", trimws(unname(lst))),
             norm_pathway_id(names(lst)))
  }, error = function(e) NULL)

  all_cpds <- unique(unlist(cpd2path, use.names = FALSE))
  N <- length(all_cpds)   # universe: all KEGG compounds in osa pathways
  n <- length(diff_kegg)
  met_rows <- lapply(names(cpd2path), function(pid) {
    members <- cpd2path[[pid]]
    k <- sum(diff_kegg %in% members)
    if (k == 0) return(NULL)
    p <- stats::phyper(k - 1, m = length(members), n = N - length(members),
                       k = n, lower.tail = FALSE)
    data.frame(pathway_id = norm_pathway_id(pid), n_metabolites = k,
               p = p, stringsAsFactors = FALSE)
  })
  met_enr <- do.call(rbind, met_rows)
  if (is.null(met_enr)) {
    stop("None of the ", n, " differential metabolites with a kegg_id fall into any ",
         "KEGG ", org, " pathway.", call. = FALSE)
  }
  met_enr$padj <- stats::p.adjust(met_enr$p, method = "BH")
  met_enr$name <- ifelse(met_enr$pathway_id %in% names(path_names),
                         path_names[met_enr$pathway_id], met_enr$pathway_id)
  met_enr$compound_ids <- vapply(met_enr$pathway_id, function(pid) {
    paste(intersect(diff_kegg, cpd2path[[pid]]), collapse = ";")
  }, character(1))
  met_enr <- met_enr[order(met_enr$padj), , drop = FALSE]
  message("Tested ", n, " differential metabolites against ", length(cpd2path),
          " KEGG ", org, " pathways; ", sum(met_enr$padj <= padj_cutoff),
          " pathway(s) significant at padj <= ", padj_cutoff, ".")
}

met_cols <- intersect(c("pathway_id", "name", "p", "padj", "n_metabolites",
                        "compound_ids"), colnames(met_enr))
write.csv(met_enr[, met_cols],
          file.path(outdir, "metabolite_pathway_enrichment.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Join at the pathway level and classify
# ---------------------------------------------------------------------------

gene_sig <- if (!is.null(gene_enr))
  gene_enr[gene_enr$p.adjust <= padj_cutoff, , drop = FALSE] else NULL
met_sig <- met_enr[met_enr$padj <= padj_cutoff, , drop = FALSE]

pathway_ids <- sort(unique(c(
  if (!is.null(gene_sig)) gene_sig$pathway_id else character(0),
  if (!is.null(met_sig)) met_sig$pathway_id else character(0)
)))

if (length(pathway_ids) == 0) {
  message("Note: no significant pathways on either side; writing an empty joint table.")
  joint <- data.frame(pathway_id = character(0), name = character(0),
                      gene_padj = numeric(0), metabolite_padj = numeric(0),
                      n_genes = integer(0), n_metabolites = integer(0),
                      category = character(0), compound_ids = character(0),
                      gene_ids = character(0), gene_direction = character(0),
                      stringsAsFactors = FALSE)
} else {
  gene_by_path <- if (!is.null(gene_sig))
    split(gene_sig, gene_sig$pathway_id) else list()
  met_by_path <- if (!is.null(met_sig))
    split(met_sig, met_sig$pathway_id) else list()

  joint_rows <- lapply(pathway_ids, function(pid) {
    g <- gene_by_path[[pid]]
    m <- met_by_path[[pid]]
    g_ids <- if (!is.null(g))
      unique(trimws(unlist(strsplit(g$gene_ids, ";", fixed = TRUE)))) else character(0)
    g_ids <- g_ids[nzchar(g_ids)]
    c_ids <- if (!is.null(m))
      unique(trimws(unlist(strsplit(m$compound_ids, ";", fixed = TRUE)))) else character(0)
    c_ids <- c_ids[nzchar(c_ids)]
    if (!is.null(g) && !is.null(m)) category <- "both"
    else if (!is.null(g)) category <- "transcriptome_only"
    else category <- "metabolome_only"
    data.frame(
      pathway_id = pid,
      name = if (!is.null(g)) g$Description[1] else if (!is.null(m)) m$name[1] else pid,
      gene_padj = if (!is.null(g)) min(g$p.adjust) else NA_real_,
      metabolite_padj = if (!is.null(m)) min(m$padj) else NA_real_,
      n_genes = length(g_ids),
      n_metabolites = length(c_ids),
      category = category,
      compound_ids = paste(c_ids, collapse = ";"),
      gene_ids = paste(g_ids, collapse = ";"),
      gene_direction = if (!is.null(g))
        paste(unique(g$direction), collapse = "+") else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  joint <- do.call(rbind, joint_rows)
  joint$min_padj <- pmin(joint$gene_padj, joint$metabolite_padj, na.rm = TRUE)
  joint$min_padj[!is.finite(joint$min_padj)] <- 1
  joint <- joint[order(joint$category != "both", joint$category, joint$min_padj), ,
                  drop = FALSE]
  joint$min_padj <- NULL
  rownames(joint) <- NULL
}

write.csv(joint, file.path(outdir, "joint_pathway_table.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Joint pathway dot plot: x = gene count, y = -log10(min padj), colour = category
# ---------------------------------------------------------------------------

if (nrow(joint) > 0) {
  plot_df <- joint
  plot_df$min_padj <- pmin(plot_df$gene_padj, plot_df$metabolite_padj, na.rm = TRUE)
  plot_df$min_padj[!is.finite(plot_df$min_padj)] <- 1
  plot_df$category <- factor(plot_df$category,
                             levels = c("both", "transcriptome_only", "metabolome_only"))
  plot_df$label <- ifelse(nchar(plot_df$name) > 42,
                          paste0(substr(plot_df$name, 1, 40), "..."), plot_df$name)

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data$n_genes,
                                    y = -log10(.data$min_padj),
                                    colour = .data$category)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_metabolites), alpha = 0.9) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed",
                        colour = nature_semantic("reference_line"), linewidth = 0.5) +
    scale_colour_nature(
      labels = c(both = "Both omics", transcriptome_only = "Transcriptome only",
                 metabolome_only = "Metabolome only"),
      drop = FALSE) +
    ggplot2::labs(
      x = "Number of DE genes in pathway",
      y = "-log10(min adjusted p-value)",
      colour = "Pathway category", size = "Differential metabolites",
      title = "Joint transcriptome-metabolome KEGG pathways",
      subtitle = paste0("Dashed line: padj = ", padj_cutoff,
                        "; 'both' = pathway significant on both omics layers")
    ) +
    theme_nature()
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(plot_df) <= 40) {
    top_rows <- head(plot_df[order(-plot_df$min_padj), ], min(10, nrow(plot_df)))
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$label),
      data = top_rows, size = 2.4, colour = "grey30", max.overlaps = 20, seed = 42)
  }
  save_nature_figure(p, file.path(outdir, "joint_pathway_dotplot"),
                     width_mm = 120, height_mm = max(80, 20 + 5 * nrow(plot_df)))

  # --- pathway_links.csv (gene -> pathway -> metabolite, for the sankey in 06)
  link_rows <- list()
  for (i in seq_len(nrow(joint))) {
    g_ids <- trimws(unlist(strsplit(joint$gene_ids[i], ";", fixed = TRUE)))
    c_ids <- trimws(unlist(strsplit(joint$compound_ids[i], ";", fixed = TRUE)))
    lfc <- setNames(de_genes$log2FoldChange, de_genes$gene_id)
    if (length(g_ids) > 0 && nzchar(g_ids[1])) {
      link_rows[[length(link_rows) + 1L]] <- data.frame(
        node = g_ids, pathway = joint$pathway_id[i], node_type = "gene",
        log2fc = ifelse(g_ids %in% names(lfc), lfc[g_ids], NA_real_),
        stringsAsFactors = FALSE)
    }
    if (length(c_ids) > 0 && nzchar(c_ids[1])) {
      link_rows[[length(link_rows) + 1L]] <- data.frame(
        node = c_ids, pathway = joint$pathway_id[i], node_type = "metabolite",
        log2fc = NA_real_, stringsAsFactors = FALSE)
    }
  }
  if (length(link_rows) > 0) {
    write.csv(do.call(rbind, link_rows), file.path(outdir, "pathway_links.csv"),
              row.names = FALSE)
  } else {
    write.csv(data.frame(node = character(0), pathway = character(0),
                         node_type = character(0), log2fc = numeric(0)),
              file.path(outdir, "pathway_links.csv"), row.names = FALSE)
  }
} else {
  write.csv(data.frame(node = character(0), pathway = character(0),
                       node_type = character(0), log2fc = numeric(0)),
            file.path(outdir, "pathway_links.csv"), row.names = FALSE)
  message("Note: joint pathway dot plot skipped (no significant pathways).")
}

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

n_both <- if (nrow(joint) > 0) sum(joint$category == "both") else 0
write_run_info_r(
  outdir = outdir, script_name = "04_joint_pathway.R",
  params = list(de_genes = de_genes_path, de_metabolites = de_mets_path,
                org = org, padj = padj_cutoff, id_type = id_type,
                metabolite_pathway_csv = if (is.null(met_pathway_csv)) "KEGGREST" else met_pathway_csv,
                n_de_genes_up = length(gene_up), n_de_genes_down = length(gene_down),
                n_pathways_joint = nrow(joint), n_pathways_both = n_both),
  packages = c("clusterProfiler", "org.Osativa.eg.db", "KEGGREST"),
  extra = c(
    "Gene enrichment: bitr (SYMBOL->ENTREZID) + enrichKEGG (keyType ncbi-geneid).",
    "Metabolite enrichment: hypergeometric test on KEGG pathway-compound links",
    "(KEGGREST) or a pre-computed table via --metabolite-pathway-csv.",
    "Categories: both / transcriptome_only / metabolome_only at padj <= threshold."
  )
)

message("Done. Outputs written to: ", normalizePath(outdir))
