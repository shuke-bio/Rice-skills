###############################################################################
# 01_qtl_mapping.R -- QTL mapping: genotype x molecular trait association
#                     (eQTL / mQTL / pQTL)
#
# Title:    Rice-skills QTL mapping step
# Description:
#   Associates genetic variants with molecular traits (gene expression ->
#   eQTL, metabolites -> mQTL, proteins -> pQTL):
#     - cis mode  : each feature tested against SNPs within +/- window of
#                   its own coordinate (fast, the primary analysis);
#     - trans mode: genome-wide scan (chunked), retaining pairs below a
#                   per-pair p-value cutoff (conservative BH across the
#                   full test count);
#     - both      : union of the two.
#   The engine is self-contained (no MatrixEQTL dependency required) and
#   reuses the shared genotype QC utilities of the GWAS module.
#
# Usage:
#   Rscript scripts/qtl/01_qtl_mapping.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --features expression_matrix.csv \
#       --feature-info gene_info.csv \
#       --outdir results/qtl/01_map --mode cis --window 1000000
#
# Inputs:
#   --genotype      SNP x sample dosage matrix CSV (0/1/2; as GWAS module)
#   --snp-info      snp_id, chr, pos
#   --features      feature x sample matrix CSV (normalized expression /
#                   metabolite intensity / protein abundance); first column
#                   = feature id
#   --feature-info  optional: feature_id, chr, pos  (required for cis mode)
#   --mode          cis | trans | both (default cis)
#   --window        cis window in bp (default 1e6)
#   --min-maf       SNP MAF filter (default 0.05)
#   --max-missing   SNP missingness filter (default 0.1)
#   --trans-p       trans per-pair retention cutoff before BH (default 1e-5)
#   --output-p      cis pairs written only if p <= this (default 0.01)
#   --chunk         features per trans block (default 100; lower = less RAM)
#   --outdir        output directory (default results/qtl/01_map)
# Outputs:
#   qtl_results.csv            cis pairs (p <= output-p) + retained trans
#   qtl_significant.csv        FDR < 0.05 subset
#   qtl_top_per_feature.csv    best SNP per feature (cis priority)
#   qtl_summary.txt            test counts, cis/trans fractions
#   manhattan_qtl.*            Manhattan of the strongest signal per SNP
#   cis_distance_plot.*        cis p-value vs. distance from feature
#   run_info.txt
#
# Notes:
#   - Features are assumed already normalized (e.g. VST/TPM for expression,
#     log+Pareto for metabolomics); the script does NOT re-normalize.
#   - Missing feature values: affected features are skipped in cis mode and
#     excluded per-pair in trans mode.
#   - Population-structure covariates are NOT included by default; for
#     eQTL in structured panels add genotype PCs as covariates (see
#     scripts/gwas/03_association.R for the PC workflow) or restrict to a
#     single subpopulation.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
source(file.path(dirname(script_file), "..", "utils", "genomics_utils.R"))
rice_set_seed(42)

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
  if (isTRUE(v)) return(if (is.logical(default) && length(default) == 1L) default else TRUE); v
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
genotype_path <- opt(opts, "genotype")
snp_info_path <- opt(opts, "snp-info")
features_path <- opt(opts, "features")
feature_info_path <- opt(opts, "feature-info")
mode <- tolower(opt(opts, "mode", "cis"))
window <- suppressWarnings(as.numeric(opt(opts, "window", "1e6")))
maf_cut <- suppressWarnings(as.numeric(opt(opts, "min-maf", "0.05")))
missing_cut <- suppressWarnings(as.numeric(opt(opts, "max-missing", "0.1")))
trans_p <- suppressWarnings(as.numeric(opt(opts, "trans-p", "1e-5")))
output_p <- suppressWarnings(as.numeric(opt(opts, "output-p", "0.01")))
chunk <- suppressWarnings(as.integer(opt(opts, "chunk", "100")))
outdir <- opt(opts, "outdir", "results/qtl/01_map")

if (any(vapply(list(genotype_path, snp_info_path, features_path), is.null, logical(1)))) {
  stop("Required: --genotype, --snp-info, --features.\n",
       "Hint: see the Usage block at the top of this script.", call. = FALSE)
}
if (!mode %in% c("cis", "trans", "both")) stop("--mode must be cis | trans | both.", call. = FALSE)
if (mode %in% c("cis", "both") && is.null(feature_info_path)) {
  stop("cis mode requires --feature-info (feature_id, chr, pos).", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills QTL mapping (", mode, " mode) ==")

# ---------------------------------------------------------------------------
# 1) Load genotypes + features
# ---------------------------------------------------------------------------
gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
fl <- filter_snps(gt$geno, snp_info, min_maf = maf_cut, max_missing = missing_cut)
G <- impute_mean(fl$geno)
snp_info <- fl$snp_info
message(sprintf("  SNPs: %d -> %d after QC (MAF >= %g, missing <= %g)",
                ncol(gt$geno), ncol(G), maf_cut, missing_cut))

feat <- read_csv_safe(features_path, "Feature matrix")
feature_ids <- as.character(feat[[1]])
F <- as.matrix(feat[, -1, drop = FALSE])
suppressWarnings(storage.mode(F) <- "double")
missing_feat <- rowMeans(is.na(F))
F <- F[missing_feat <= 0.5, , drop = FALSE]
feature_ids <- feature_ids[missing_feat <= 0.5]
message(sprintf("  features: %d (dropped %d with >50%% missing)",
                nrow(F), sum(missing_feat > 0.5)))

# sample alignment
shared <- intersect(rownames(G), colnames(F))
if (length(shared) < 5) {
  stop("Fewer than 5 shared samples between genotype and feature matrix. ",
       "Hint: sample IDs must match exactly (case-sensitive).", call. = FALSE)
}
G <- G[shared, , drop = FALSE]
F <- F[, shared, drop = FALSE]
n <- length(shared)
message(sprintf("  samples: %d shared", n))

feat_info <- NULL
if (!is.null(feature_info_path)) {
  feat_info <- read_csv_safe(feature_info_path, "Feature info table")
  check_columns(feat_info, c("feature_id", "chr", "pos"), "feature info table")
  feat_info$chr <- as.character(feat_info$chr)
  feat_info$pos <- as.numeric(feat_info$pos)
  feat_info <- feat_info[match(feature_ids, feat_info$feature_id), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 2) Pairwise regression helper (vectorized over SNPs)
# ---------------------------------------------------------------------------
reg_pairs <- function(y, X) {
  # y: numeric n; X: n x k, no NA (genotypes imputed); returns beta/t/p
  xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  sxx <- colSums(xc^2)
  sxy <- colSums(xc * yc)
  beta <- sxy / sxx
  resid_ss <- sum(yc^2) - beta * sxy
  sigma2 <- pmax(resid_ss, 0) / (n - 2)
  se <- sqrt(sigma2 / sxx)
  tstats <- beta / se
  p <- 2 * stats::pt(-abs(tstats), df = n - 2)
  data.frame(beta = beta, t = tstats, p = pmax(p, 1e-300))
}

all_results <- list()
n_tests <- 0L
cis_p_all <- numeric(0)      # p-values of all cis tests (for conservative BH)
trans_tests_total <- 0L      # total number of trans tests performed

# ---------------------------------------------------------------------------
# 3) cis scan
# ---------------------------------------------------------------------------
if (mode %in% c("cis", "both")) {
  message("--> cis scan (window = ", format(window, scientific = TRUE), " bp)")
  cis_list <- vector("list", 0)
  cis_p_list <- vector("list", 0)
  for (i in seq_len(nrow(F))) {
    fch <- feat_info$chr[i]; fpo <- feat_info$pos[i]
    if (is.na(fch) || is.na(fpo)) next
    idx <- which(snp_info$chr == fch & abs(snp_info$pos - fpo) <= window)
    if (!length(idx)) next
    y <- F[i, ]
    if (anyNA(y)) {
      ok <- !is.na(y)
      if (sum(ok) < n - 5) next
      yv <- y[ok]
      Xv <- G[ok, idx, drop = FALSE]
    } else {
      yv <- y
      Xv <- G[, idx, drop = FALSE]
    }
    if (anyNA(yv) || ncol(Xv) == 0) next
    rp <- reg_pairs(yv, Xv)
    n_tests <- n_tests + nrow(rp)
    cis_p_list[[length(cis_p_list) + 1L]] <- rp$p
    keep_i <- which(rp$p <= output_p)
    if (length(keep_i)) {
      cis_list[[length(cis_list) + 1L]] <- data.frame(
        feature_id = feature_ids[i],
        snp_id = colnames(G)[idx][keep_i],
        snp_chr = snp_info$chr[idx][keep_i],
        snp_pos = snp_info$pos[idx][keep_i],
        feature_chr = fch, feature_pos = fpo,
        distance = abs(snp_info$pos[idx][keep_i] - fpo),
        beta = rp$beta[keep_i], t = rp$t[keep_i], p = rp$p[keep_i],
        type = "cis", stringsAsFactors = FALSE
      )
    }
  }
  cis_df <- if (length(cis_list)) do.call(rbind, cis_list) else data.frame()
  message(sprintf("  cis tests: %d | pairs kept (p <= %g): %d",
                  n_tests, output_p, nrow(cis_df)))
  all_results$cis <- cis_df
  cis_p_all <- unlist(cis_p_list)
}

# ---------------------------------------------------------------------------
# 4) trans scan (chunked genome-wide)
# ---------------------------------------------------------------------------
if (mode %in% c("trans", "both")) {
  message("--> trans scan (chunk = ", chunk, " features)")
  Gs <- scale(G, center = TRUE, scale = FALSE)
  sxx_all <- colSums(Gs^2)
  snp_ids <- colnames(G)
  trans_list <- vector("list", 0)
  trans_tests <- 0L
  starts <- seq(1, nrow(F), by = chunk)
  for (s in starts) {
    e <- min(s + chunk - 1L, nrow(F))
    Y <- F[s:e, , drop = FALSE]
    # per-pair NA handling: features with NA are dropped for this block
    okrow <- rowSums(is.na(Y)) == 0
    if (!any(okrow)) next
    Y <- Y[okrow, , drop = FALSE]
    Ys <- scale(t(Y), center = TRUE, scale = FALSE)  # n x k
    syy <- colSums(Ys^2)
    sxy <- crossprod(Gs, Ys)                          # m x k
    trans_tests <- trans_tests + length(sxy)
    r <- sxy / sqrt(outer(sxx_all, syy))
    r <- pmin(pmax(r, -1 + 1e-12), 1 - 1e-12)
    tstats <- r * sqrt((n - 2) / (1 - r^2))
    pmat <- 2 * stats::pt(-abs(tstats), df = n - 2)
    hit <- which(pmat < trans_p, arr.ind = TRUE)
    if (nrow(hit)) {
      trans_list[[length(trans_list) + 1L]] <- data.frame(
        feature_id = feature_ids[s:e][okrow][hit[, 2]],
        snp_id = snp_ids[hit[, 1]],
        snp_chr = snp_info$chr[hit[, 1]],
        snp_pos = snp_info$pos[hit[, 1]],
        feature_chr = if (!is.null(feat_info)) feat_info$chr[s:e][okrow][hit[, 2]] else NA_character_,
        feature_pos = if (!is.null(feat_info)) feat_info$pos[s:e][okrow][hit[, 2]] else NA_real_,
        distance = NA_real_,
        beta = sxy[hit] / sxx_all[hit[, 1]],
        t = tstats[hit], p = pmat[hit],
        type = "trans", stringsAsFactors = FALSE
      )
    }
  }
  trans_df <- if (length(trans_list)) do.call(rbind, trans_list) else data.frame()
  message(sprintf("  trans tests: %d | pairs retained (p < %g): %d",
                  trans_tests, trans_p, nrow(trans_df)))
  trans_tests_total <- trans_tests
  all_results$trans <- trans_df
}

# ---------------------------------------------------------------------------
# 5) Combine, FDR, top hits
# ---------------------------------------------------------------------------
res <- do.call(rbind, Filter(function(d) nrow(d) > 0, all_results))
if (is.null(res) || nrow(res) == 0) {
  message("No pairs passed the retention cutoffs; nothing to write.")
  write_run_info_r(outdir, "01_qtl_mapping.R", params = list(
    genotype = genotype_path, features = features_path, mode = mode))
  quit(save = "no", status = 0)
}

# Conservative BH across the FULL test counts (retention cutoffs keep only
# the smallest p-values, whose ranks are unchanged, so q = p * m_full / rank
# is exact for everything retained below the cutoffs).
fdr <- rep(NA_real_, nrow(res))
is_cis <- res$type == "cis"

.bh_conservative <- function(p_kept, m_full) {
  # p_kept: retained p-values (all <= retention cutoff)
  ord <- order(p_kept)
  pp <- p_kept[ord]
  qv <- pp * m_full / seq_along(pp)
  qv <- rev(cummin(rev(pmin(qv, 1))))
  out <- rep(NA_real_, length(p_kept))
  out[ord] <- qv
  out
}

if (any(is_cis)) {
  fdr[is_cis] <- .bh_conservative(res$p[is_cis], length(cis_p_all))
}
if (any(!is_cis)) {
  if (trans_tests_total > 0) {
    fdr[!is_cis] <- .bh_conservative(res$p[!is_cis], trans_tests_total)
  }
}
res$fdr <- fdr
res <- res[order(res$p), , drop = FALSE]
utils::write.csv(res, file.path(outdir, "qtl_results.csv"), row.names = FALSE)

sig <- res[!is.na(res$fdr) & res$fdr < 0.05, , drop = FALSE]
utils::write.csv(sig, file.path(outdir, "qtl_significant.csv"), row.names = FALSE)

# best hit per feature (cis first)
res$priority <- ifelse(res$type == "cis", 0, 1)
res <- res[order(res$feature_id, res$priority, res$p), , drop = FALSE]
top <- res[!duplicated(res$feature_id), , drop = FALSE]
utils::write.csv(top, file.path(outdir, "qtl_top_per_feature.csv"), row.names = FALSE)

sum_lines <- c(
  "Rice-skills -- QTL mapping summary",
  sprintf("mode               : %s", mode),
  sprintf("SNPs (QC)          : %d", ncol(G)),
  sprintf("features           : %d", nrow(F)),
  sprintf("samples            : %d", n),
  sprintf("pairs written      : %d (cis %d / trans %d)",
          nrow(res), sum(res$type == "cis"), sum(res$type == "trans")),
  sprintf("significant (FDR<0.05): %d", nrow(sig))
)
writeLines(sum_lines, file.path(outdir, "qtl_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

# ---------------------------------------------------------------------------
# 6) Figures
# ---------------------------------------------------------------------------
pal <- rice_palette()
# strongest signal per SNP -> Manhattan
snp_min <- aggregate(list(p = res$p), by = list(snp_id = res$snp_id), FUN = min)
snp_min <- merge(snp_min, snp_info, by = "snp_id")
if (nrow(snp_min) >= 10) {
  p_man <- manhattan_plot(snp_min, snp_info, title = "QTL Manhattan (strongest per SNP)")
  save_nature_figure(p_man, file.path(outdir, "manhattan_qtl"), width_mm = 180, height_mm = 80)
}
# cis p vs distance
cis_rows <- res[res$type == "cis" & !is.na(res$distance), , drop = FALSE]
if (nrow(cis_rows) >= 10) {
  p_dist <- ggplot2::ggplot(cis_rows, ggplot2::aes(x = distance / 1000, y = -log10(p))) +
    ggplot2::geom_point(alpha = 0.4, size = 0.7, colour = pal$groups5[1]) +
    ggplot2::labs(x = "Distance from feature (kb)", y = expression(-log[10](p)),
                  title = "cis-QTL signal vs. distance") +
    theme_nature()
  save_nature_figure(p_dist, file.path(outdir, "cis_distance_plot"),
                     width_mm = 89, height_mm = 75)
}

write_run_info_r(outdir, "01_qtl_mapping.R", params = list(
  genotype = genotype_path, snp_info = snp_info_path, features = features_path,
  feature_info = feature_info_path, mode = mode, window = window,
  n_pairs_written = nrow(res), n_significant = nrow(sig)
))
message("Done. Outputs written to: ", normalizePath(outdir))
