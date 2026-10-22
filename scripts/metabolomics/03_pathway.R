###############################################################################
# 03_pathway.R -- KEGG pathway enrichment for differential metabolites (step 03)
#
# Title:    Metabolic pathway over-representation analysis (KEGG, Oryza sativa)
# Description:
#   Tests whether the differential metabolites from 02_oplsda.R are enriched in
#   KEGG metabolic pathways. Two interchangeable back-ends are implemented:
#     1. MetaboAnalystR (when installed): its documented KEGG over-representation
#        workflow (InitDataObjects -> Setup.MapData -> CrossReferencing ->
#        SetKEGG.PathLib -> CalculateOraScore).
#     2. Built-in fallback (always available): KEGGREST retrieves the pathway
#        compound lists for the organism and a hypergeometric test is performed
#        with base R stats::phyper + BH correction. This is the same
#        over-representation algorithm MetaboAnalystR uses, without the Java /
#        web-server dependency.
#
# Usage:
#   Rscript scripts/metabolomics/03_pathway.R \
#       --diff-metabolites results/metabolome/02_oplsda/diff_metabolites.csv \
#       --outdir results/metabolome/03_pathway --org osa --padj 0.05
#
# Inputs:
#   --diff-metabolites  CSV with a kegg_id column (e.g. diff_metabolites.csv
#                       from 02_oplsda.R). If a 'significant' column is present,
#                       only rows flagged TRUE are used.
#
# Outputs (written to --outdir):
#   pathway_enrichment.csv  pathway_id, name, hits, background (pathway compound
#                           count), fold_enrichment, p, padj, compound_ids
#   pathway_bubble.*        bubble plot: x = fold enrichment, y = -log10(padj),
#                           size = hits, colour = padj (sequential ramp)
#   run_info.txt            parameters, package versions and the seed
#
# Notes:
#   - --org expects a KEGG organism code ("osa" = Oryza sativa japonica).
#   - The hypergeometric universe (background N) is the set of compounds
#     annotated to the organism's KEGG pathways; k = number of input
#     metabolites carrying a KEGG annotation; hits = input metabolites found
#     in the pathway.
#   - MetaboAnalystR requires internet access (its compound/pathway libraries
#     are downloaded from the MetaboAnalyst web server) and supports ~21 model
#     organisms (rice included). Its API has changed between major versions, so
#     this branch is best-effort: any failure falls through to the built-in
#     KEGGREST implementation, which produces equivalent statistics offline of
#     the KEGG REST API itself.
#   - KEGGREST retrieves at most 10 entries per request; pathways are therefore
#     fetched in batches of 10.
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
# Input handling
# ---------------------------------------------------------------------------

#' Read differential metabolites and their KEGG compound IDs.
#'
#' @param path Path to the diff-metabolites CSV.
#' @return List with the data frame, the unique KEGG IDs, and counts of
#'   metabolites with/without KEGG annotation.
read_diff_metabolites <- function(path) {
  df <- read_csv_safe(path, "differential metabolites table")
  check_columns(df, "kegg_id", "differential metabolites table")
  if ("significant" %in% colnames(df)) {
    keep <- as.logical(df$significant)
    keep[is.na(keep)] <- FALSE
    message("Filtering on the 'significant' column: ", sum(keep), " of ", nrow(df),
            " rows are flagged significant.")
    df <- df[keep, , drop = FALSE]
  }
  kegg_raw <- trimws(as.character(df$kegg_id))
  kegg_raw[is.na(kegg_raw)] <- ""
  has_kegg <- nzchar(kegg_raw)
  n_missing <- sum(!has_kegg)
  if (n_missing) {
    message(n_missing, " of ", nrow(df), " differential metabolites have no KEGG ",
            "annotation (empty kegg_id) and are excluded from the enrichment test.")
  }
  kegg_ids <- unique(kegg_raw[has_kegg])
  if (!length(kegg_ids)) {
    stop("No KEGG-annotated metabolites in ", path,
         "\nHint: annotate kegg_id (KEGG compound IDs such as C00031) before running enrichment.",
         call. = FALSE)
  }
  list(df = df, kegg_ids = kegg_ids, n_missing = n_missing)
}

# ---------------------------------------------------------------------------
# Back-end 1: MetaboAnalystR (best effort, documented workflow)
# ---------------------------------------------------------------------------

#' KEGG over-representation via MetaboAnalystR.
#'
#' Workflow (MetaboAnalystR >= 3.x, Pathway Analysis vignette):
#'   mSet <- MetaboAnalystR::InitDataObjects("conc", "pathora", FALSE)
#'   mSet <- MetaboAnalystR::Setup.MapData(mSet, kegg_ids)
#'   mSet <- MetaboAnalystR::CrossReferencing(mSet, "kegg")
#'   mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)
#'   mSet <- MetaboAnalystR::SetKEGG.PathLib(mSet, org, "current")
#'   mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE)
#'   mSet <- MetaboAnalystR::CalculateOraScore(mSet, "rbc", "hyperg")
#'   ora  <- mSet$analSet$ora.mat
#'
#' Limitations: internet required; organism must be one of the ~21 supported
#' (rice = "osa"); API differs between versions. Any failure is caught by the
#' caller, which falls back to the built-in KEGGREST implementation.
metaboanalyst_enrichment <- function(kegg_ids, org, padj_cut) {
  ma <- "MetaboAnalystR"
  if (!requireNamespace(ma, quietly = TRUE)) {
    stop("MetaboAnalystR is not installed.", call. = FALSE)
  }
  mSet <- get("InitDataObjects", envir = asNamespace(ma))("conc", "pathora", FALSE)
  mSet <- get("Setup.MapData", envir = asNamespace(ma))(mSet, kegg_ids)
  mSet <- get("CrossReferencing", envir = asNamespace(ma))(mSet, "kegg")
  mSet <- get("CreateMappingResultTable", envir = asNamespace(ma))(mSet)
  mSet <- get("SetKEGG.PathLib", envir = asNamespace(ma))(mSet, org, "current")
  mSet <- get("SetMetabolomeFilter", envir = asNamespace(ma))(mSet, FALSE)
  mSet <- get("CalculateOraScore", envir = asNamespace(ma))(mSet, "rbc", "hyperg")
  ora <- mSet$analSet$ora.mat
  if (is.null(ora) || !nrow(ora)) stop("MetaboAnalystR returned no pathway results.", call. = FALSE)
  # Defensive column extraction (field names vary between versions).
  pick <- function(patterns, default = NA) {
    hit <- colnames(ora)[Reduce(`|`, lapply(patterns, function(p) grepl(p, colnames(ora), ignore.case = TRUE)))]
    if (length(hit)) ora[[hit[1]]] else rep(default, nrow(ora))
  }
  hits <- suppressWarnings(as.numeric(pick("^hits$")))
  background <- suppressWarnings(as.numeric(pick("^total$|^background$")))
  p_raw <- suppressWarnings(as.numeric(pick("^raw p$|^p value$|^p$")))
  if (all(is.na(p_raw))) stop("Could not locate the p-value column in MetaboAnalystR ora.mat.", call. = FALSE)
  padj <- suppressWarnings(as.numeric(pick("^fdr$|^holm|^adjusted")))
  if (all(is.na(padj))) padj <- stats::p.adjust(p_raw, method = "BH")
  names_chr <- rownames(ora)
  # Pathway IDs: look for an ID-like column, else extract osaXXXXX from the name.
  id_col <- colnames(ora)[grepl("pathway.*id|^id$", colnames(ora), ignore.case = TRUE)]
  pathway_id <- if (length(id_col)) as.character(ora[[id_col[1]]]) else
    sub(".*(osa[0-9]{5}).*", "\\1", names_chr)
  no_id <- !grepl("osa[0-9]{5}", pathway_id)
  if (any(no_id)) pathway_id[no_id] <- paste0("pathway_", seq_len(sum(no_id)))
  # Fold enrichment needs the library universe N. MetaboAnalystR does not expose
  # it directly, but 'Expected' = k * Total / N, so N can be recovered from the
  # expected-hit counts (median across pathways for robustness).
  expected <- suppressWarnings(as.numeric(pick("^expected$")))
  n_input <- length(kegg_ids)
  fold <- rep(NA_real_, nrow(ora))
  ok <- !is.na(expected) & expected > 0 & !is.na(background) & background > 0
  if (any(ok)) {
    n_universe <- stats::median(n_input * background[ok] / expected[ok])
    if (is.finite(n_universe) && n_universe > 0) {
      fold <- (hits / background) / (n_input / n_universe)
    }
  }
  data.frame(
    pathway_id = pathway_id,
    name = names_chr,
    hits = hits,
    background = background,
    fold_enrichment = fold,
    p = p_raw,
    padj = padj,
    compound_ids = vapply(rownames(ora), function(nm) "", character(1)),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Back-end 2: built-in KEGGREST + hypergeometric fallback
# ---------------------------------------------------------------------------

#' Fetch pathway -> compound mappings for an organism via KEGGREST.
#'
#' @return Named list of compound-ID vectors (names = pathway IDs such as
#'   "osa00010"); the "names" attribute carries pathway display names.
fetch_kegg_pathways <- function(org) {
  pw_list <- KEGGREST::keggList("pathway", org)
  if (!length(pw_list)) stop("No KEGG pathways found for organism '", org, "'.", call. = FALSE)
  pw_ids <- names(pw_list)
  compounds <- list()
  display <- character(0)
  batch <- 10L  # KEGGREST retrieves at most 10 entries per request
  for (i in seq(1L, length(pw_ids), by = batch)) {
    ids <- pw_ids[i:min(i + batch - 1L, length(pw_ids))]
    recs <- KEGGREST::keggGet(ids)
    for (j in seq_along(recs)) {
      rec <- recs[[j]]
      cmps <- rec$COMPOUND
      pid <- sub("^path:", "", ids[j])
      if (!is.null(cmps) && length(cmps)) {
        compounds[[pid]] <- names(cmps)
        display[pid] <- if (!is.null(rec$NAME)) paste(rec$NAME, collapse = " ") else unname(pw_list[[ids[j]]])
      }
    }
  }
  if (!length(compounds)) {
    stop("No compound annotations found for organism '", org,
         "' pathways.\nHint: check the KEGG organism code (rice japonica = 'osa').",
         call. = FALSE)
  }
  attr(compounds, "display") <- display
  compounds
}

#' Built-in KEGG over-representation analysis.
#'
#' Hypergeometric test identical to MetaboAnalystR's "hyperg" method:
#'   P(X >= hits), X ~ Hypergeom(N, K = pathway size, n = input size)
#' with BH correction across the tested pathways.
keggrest_enrichment <- function(kegg_ids, org, padj_cut) {
  need_pkgs("KEGGREST")
  message("Fetching KEGG pathway definitions for '", org, "' (KEGGREST, batched by 10)...")
  pw_compounds <- fetch_kegg_pathways(org)
  display <- attr(pw_compounds, "display")
  universe <- unique(unlist(pw_compounds, use.names = FALSE))
  n_universe <- length(universe)
  k_input <- length(kegg_ids)
  outside <- setdiff(kegg_ids, universe)
  if (length(outside)) {
    message(length(outside), " input KEGG ID(s) are not part of any ", org,
            " pathway and cannot contribute to enrichment: ",
            paste(head(outside, 5), collapse = ", "), if (length(outside) > 5) " ..." else "")
  }
  message("Background: ", n_universe, " compounds in ", length(pw_compounds), " ", org,
          " pathways; ", k_input, " annotated input metabolites.")

  rows <- lapply(names(pw_compounds), function(pid) {
    pw_cmp <- pw_compounds[[pid]]
    hit_ids <- intersect(pw_cmp, kegg_ids)
    m <- length(hit_ids)
    n_pw <- length(pw_cmp)
    if (m == 0L) return(NULL)
    p <- stats::phyper(m - 1L, m = n_pw, n = n_universe - n_pw, k = k_input,
                       lower.tail = FALSE)
    fold <- (m / n_pw) / (k_input / n_universe)
    data.frame(pathway_id = pid, name = display[[pid]], hits = m, background = n_pw,
               fold_enrichment = fold, p = p, padj = NA_real_,
               compound_ids = paste(hit_ids, collapse = ";"), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  if (is.null(res)) {
    stop("None of the ", k_input, " input metabolites mapped to any ", org,
         " KEGG pathway.\nHint: verify the kegg_id column uses KEGG compound IDs (C-numbers).",
         call. = FALSE)
  }
  res$padj <- stats::p.adjust(res$p, method = "BH")
  res <- res[order(res$p), ]
  rownames(res) <- NULL
  res
}

# ---------------------------------------------------------------------------
# Bubble plot
# ---------------------------------------------------------------------------

plot_pathway_bubble <- function(res, padj_cut, outdir, org) {
  # Scientific conclusion: pathways in the upper-right region (high fold
  # enrichment, low adjusted p) are the metabolic processes most strongly
  # perturbed between the compared groups.
  sig <- res[which(res$padj <= padj_cut), , drop = FALSE]
  note <- ""
  if (!nrow(sig)) {
    sig <- utils::head(res[order(-res$fold_enrichment), ], 10)
    note <- paste0("No pathway reached padj <= ", padj_cut,
                   "; showing the top 10 by fold enrichment.")
    message(note)
  } else if (nrow(sig) > 40) {
    sig <- utils::head(sig[order(sig$padj), ], 40)
    note <- "Top 40 significant pathways shown."
  }
  sig <- sig[order(-sig$fold_enrichment), , drop = FALSE]
  sig$label <- ifelse(nchar(sig$name) > 42, paste0(substr(sig$name, 1, 40), "..."), sig$name)
  sig$neg_log10_padj <- -log10(pmax(sig$padj, .Machine$double.xmin))
  p <- ggplot2::ggplot(sig, ggplot2::aes(x = fold_enrichment, y = neg_log10_padj)) +
    ggplot2::geom_point(ggplot2::aes(size = hits, colour = padj), alpha = 0.85) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut),
                        colour = nature_semantic("reference_line"),
                        linewidth = 0.5, linetype = "dashed") +
    ggplot2::scale_colour_gradientn(colours = rev(sequential_palette()),
                                    name = "padj") +
    ggplot2::scale_size_continuous(range = c(1.2, 4.5), name = "Hits") +
    ggrepel::geom_label_repel(ggplot2::aes(label = label), size = 2.0, colour = "grey20",
                              max.overlaps = 25, min.segment.length = 0) +
    ggplot2::labs(
      title = "KEGG pathway enrichment",
      subtitle = paste0("Differential metabolites vs '", org, "' KEGG pathways; ",
                        "dashed line = padj ", padj_cut, if (nzchar(note)) paste0("; ", note) else ""),
      x = "Fold enrichment", y = expression(-log[10](padj))) +
    theme_nature()
  save_nature_figure(p, file.path(outdir, "pathway_bubble"), width_mm = 89, height_mm = 90)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  diff_path <- opt(opts, "diff-metabolites")
  if (is.null(diff_path)) stop("--diff-metabolites is required (CSV with a kegg_id column).", call. = FALSE)
  outdir <- opt(opts, "outdir", "results/metabolome/03_pathway")
  org <- opt(opts, "org", "osa")
  padj_cut <- as.numeric(opt(opts, "padj", 0.05))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", normalizePath(outdir))

  input <- read_diff_metabolites(diff_path)
  kegg_ids <- input$kegg_ids
  message("Enrichment input: ", length(kegg_ids), " unique KEGG-annotated metabolites (",
          nrow(input$df), " rows in the table, ", input$n_missing, " without KEGG annotation).")

  # ---- Back-end selection ----------------------------------------------------
  backend <- "builtin"
  if (requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    message("MetaboAnalystR detected: attempting its KEGG over-representation workflow.")
    res <- tryCatch(metaboanalyst_enrichment(kegg_ids, org, padj_cut),
                    error = function(e) {
                      message("MetaboAnalystR workflow failed (", conditionMessage(e),
                              "); falling back to the built-in KEGGREST implementation.")
                      NULL
                    })
    if (!is.null(res)) backend <- "MetaboAnalystR"
  } else {
    message("MetaboAnalystR not installed; using the built-in KEGGREST hypergeometric fallback.")
  }
  if (backend == "builtin") {
    res <- keggrest_enrichment(kegg_ids, org, padj_cut)
  }
  res <- res[order(res$p), , drop = FALSE]
  rownames(res) <- NULL

  # ---- Export ------------------------------------------------------------------
  utils::write.csv(res, file.path(outdir, "pathway_enrichment.csv"),
                   row.names = FALSE, quote = FALSE)
  n_sig <- sum(res$padj <= padj_cut, na.rm = TRUE)
  message(n_sig, " of ", nrow(res), " tested pathways are significant (padj <= ", padj_cut, ").")
  if (n_sig) {
    top <- utils::head(res[which(res$padj <= padj_cut), ], 5)
    message("Top pathways:")
    for (i in seq_len(nrow(top))) {
      message("  ", top$pathway_id[i], " ", top$name[i], " (hits ", top$hits[i], "/",
              top$background[i], ", fold ", round(top$fold_enrichment[i], 2), ", padj ",
              format(top$padj[i], digits = 3), ")")
    }
  }

  # ---- Figure -------------------------------------------------------------------
  plot_pathway_bubble(res, padj_cut, outdir, org)
  message("Saved figures (PDF/TIFF/PNG): pathway_bubble")

  # ---- Provenance -----------------------------------------------------------------
  write_run_info_r(
    outdir,
    script_name = "03_pathway.R",
    params = list(
      diff_metabolites = diff_path, org = org, padj_cutoff = padj_cut,
      backend = backend,
      n_input_rows = nrow(input$df), n_input_with_kegg = length(kegg_ids),
      n_input_without_kegg = input$n_missing,
      n_pathways_tested = nrow(res), n_significant = n_sig
    ),
    packages = if (backend == "MetaboAnalystR") c("MetaboAnalystR") else "KEGGREST",
    extra = c(
      "Hypergeometric test (phyper, lower.tail = FALSE) with BH correction;",
      "background = number of compounds annotated to the pathway;",
      "universe = compounds annotated to the organism's KEGG pathways."
    )
  )
  message("Done. Results in ", normalizePath(outdir))
  invisible(NULL)
}

main()
