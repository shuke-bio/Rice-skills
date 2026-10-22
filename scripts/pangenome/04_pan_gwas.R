###############################################################################
# 04_pan_gwas.R -- pan-genome GWAS (PAV / SV markers -> trait association)
#
# Title:    Rice-skills pan-GWAS step
# Description:
#   Association analysis with pangenome-derived markers: presence/absence
#   (PAV) or structural-variant (SV) genotypes are treated as bi-allelic
#   markers (0 = absent, 1 = present) and tested against a trait with the
#   same OLS/GLM engine used by the SNP-GWAS module
#   (scripts/gwas/03_association.R), including PCA covariate correction
#   for population structure.
#   Marker matrices can come from:
#     - the minigraph-cactus VCF (presence/absence coded via bcftools),
#     - APAV output (PAV table), or
#     - scripts/pangenome/03_pav_analysis.R (pav_binary.csv).
#
# Usage:
#   Rscript scripts/pangenome/04_pan_gwas.R \
#       --pav demo_data/pangenome/pav_matrix.csv \
#       --snp-info demo_data/pangenome/marker_info.csv \
#       --phenotype demo_data/pangenome/phenotype.csv \
#       --trait grain_yield --outdir results/pangenome/04_pan_gwas
#
# Inputs:
#   --pav         marker x sample matrix CSV (0/1 or coverage 0-1)
#   --snp-info    marker info table: marker_id, chr, pos
#   --phenotype   phenotype table (first column = sample_id, then traits)
#   --trait       trait column name
#   --values      binary | cov | auto (default auto)
#   --cov-threshold  coverage >= this counts as present (default 0.5)
#   --min-freq    minimum presence frequency to keep a marker (default 0.05)
#   --max-missing maximum marker missingness (default 0.1)
#   --npc         number of PCA covariates (default 5; 0 to disable)
#   --p-threshold significance threshold (default Bonferroni 0.05/n)
#   --outdir      output directory (default results/pangenome/04_pan_gwas)
# Outputs:
#   pan_gwas_results.csv      all tested markers (beta, se, t, p, freq, ...)
#   pan_gwas_significant.csv  markers passing the threshold
#   marker_qc.csv             per-marker presence freq / missingness / kept
#   pca_scores.csv            PCA covariates used
#   manhattan_pan.* / qq_pan.*   figures (PDF/TIFF/PNG)
#   run_info.txt
#
# Notes:
#   - The OLS model on 0/1 markers is the standard linear-probability
#     approach used for PAV-GWAS screening; large-effect hits should be
#     re-checked with mixed models (rMVP supports 0/1 covariates via
#     --geno) or logistic models.
#   - For SV genotypes from vg/cactus VCFs with dosages 0/1/2 the same
#     script works (values pass through unchanged, --values binary only
#     trims >1 values; use --values auto to keep dosage coding).
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
pav_path <- opt(opts, "pav")
snp_info_path <- opt(opts, "snp-info")
phenotype_path <- opt(opts, "phenotype")
trait <- opt(opts, "trait")
outdir <- opt(opts, "outdir", "results/pangenome/04_pan_gwas")
values_mode <- tolower(opt(opts, "values", "auto"))
cov_thr <- suppressWarnings(as.numeric(opt(opts, "cov-threshold", "0.5")))
min_freq <- suppressWarnings(as.numeric(opt(opts, "min-freq", "0.05")))
max_missing <- suppressWarnings(as.numeric(opt(opts, "max-missing", "0.1")))
npc <- suppressWarnings(as.integer(opt(opts, "npc", "5")))
p_threshold <- opt(opts, "p-threshold", NULL)

if (any(vapply(list(pav_path, snp_info_path, phenotype_path, trait), is.null, logical(1)))) {
  stop("Required: --pav, --snp-info, --phenotype, --trait.\n",
       "Hint: see the Usage block at the top of this script.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills pan-GWAS (PAV / SV markers) ==")

# ---------------------------------------------------------------------------
# 1) Load marker matrix (marker x sample) -> binarize -> samples x markers
# ---------------------------------------------------------------------------
pav <- read_csv_safe(pav_path, "PAV/SV marker matrix")
marker_ids <- as.character(pav[[1]])
if (anyDuplicated(marker_ids)) marker_ids <- make.unique(marker_ids)
mat <- as.matrix(pav[, -1, drop = FALSE])
suppressWarnings(storage.mode(mat) <- "double")

vals <- mat[!is.na(mat)]
if (!all(vals %in% c(0, 1, 2))) {
  # coverage-style input: binarize
  message(sprintf("  coverage-style values detected -> binarize at %.2f", cov_thr))
  mat[!is.na(mat)] <- as.numeric(mat[!is.na(mat)] >= cov_thr)
} else if (values_mode == "binary") {
  mat[!is.na(mat) & !(mat %in% c(0, 1))] <- NA_real_
}
geno <- t(mat)
rownames(geno) <- colnames(pav)[-1]
colnames(geno) <- marker_ids
message(sprintf("  markers: %d | samples: %d", ncol(geno), nrow(geno)))

snp_info <- read_csv_safe(snp_info_path, "marker info table")
check_columns(snp_info, c("marker_id", "chr", "pos"), "marker info table")
snp_info$chr <- as.character(snp_info$chr)
snp_info$pos <- as.numeric(snp_info$pos)
pheno <- read_pheno(phenotype_path)
if (!trait %in% colnames(pheno)) {
  stop("Trait '", trait, "' not found in the phenotype table. Available: ",
       paste(setdiff(colnames(pheno), "sample_id"), collapse = ", "), call. = FALSE)
}
al <- align_samples_geno(geno, pheno)
message(sprintf("  %d samples shared between markers and phenotype.", nrow(al$geno)))

# ---------------------------------------------------------------------------
# 2) Marker QC (presence frequency + missingness)
# ---------------------------------------------------------------------------
g <- al$geno
presence <- colMeans(g, na.rm = TRUE)
missing <- colMeans(is.na(g))
keep <- presence >= min_freq & presence <= 1 - min_freq & missing <= max_missing
qc <- data.frame(marker_id = colnames(g), presence_freq = round(presence, 4),
                 missing = round(missing, 4), kept = keep, stringsAsFactors = FALSE)
utils::write.csv(qc, file.path(outdir, "marker_qc.csv"), row.names = FALSE)
message(sprintf("  marker filter (freq %g-%g, missing <= %g): %d -> %d kept",
                min_freq, 1 - min_freq, max_missing, ncol(g), sum(keep)))
if (sum(keep) < 10) {
  stop("Fewer than 10 markers remain after QC. Hint: relax --min-freq / ",
       "--max-missing, or check the input matrix orientation ",
       "(rows = markers, columns = samples).", call. = FALSE)
}
g <- g[, keep, drop = FALSE]
info_kept <- snp_info[match(colnames(g), snp_info$marker_id), , drop = FALSE]
if (anyNA(info_kept$chr)) {
  stop("Some kept markers are missing from the marker info table. ",
       "Hint: marker IDs must match exactly between --pav and --snp-info.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# 3) Trait + PCA covariates
# ---------------------------------------------------------------------------
y <- al$pheno[[trait]]
if (all(is.na(y))) stop("Trait '", trait, "' is entirely missing.", call. = FALSE)
message(sprintf("  trait '%s': n = %d non-missing, range [%s]", trait,
                sum(!is.na(y)), paste(round(range(y, na.rm = TRUE), 3), collapse = ", ")))

covariates <- NULL
if (npc > 0) {
  pca_res <- compute_pca(g, npc = npc)
  covariates <- pca_res$scores[, seq_len(min(npc, ncol(pca_res$scores))), drop = FALSE]
  utils::write.csv(data.frame(sample_id = rownames(pca_res$scores), pca_res$scores),
                   file.path(outdir, "pca_scores.csv"), row.names = FALSE)
  message(sprintf("  PCA covariates: %d PCs (PC1 explains %.1f%%)",
                  ncol(covariates), 100 * pca_res$var_explained[1]))
}

# ---------------------------------------------------------------------------
# 4) Association
# ---------------------------------------------------------------------------
res <- gwas_glm(g, y, covariates)
res$marker_id <- res$snp_id
res$engine <- "simple-GLM (pan)"
res <- merge(res, qc[, c("marker_id", "presence_freq")], by = "marker_id", all.x = TRUE)
res <- res[order(res$p), , drop = FALSE]
utils::write.csv(res, file.path(outdir, "pan_gwas_results.csv"), row.names = FALSE)

if (is.null(p_threshold)) {
  thresh <- 0.05 / nrow(res)
} else {
  thresh <- as.numeric(p_threshold)
}
sig <- res[!is.na(res$p) & res$p <= thresh, , drop = FALSE]
utils::write.csv(sig, file.path(outdir, "pan_gwas_significant.csv"), row.names = FALSE)
message(sprintf("  tested markers: %d | p <= %.2e (Bonferroni): %d",
                nrow(res), thresh, nrow(sig)))

# ---------------------------------------------------------------------------
# 5) Figures
# ---------------------------------------------------------------------------
assoc_plot <- data.frame(snp_id = res$marker_id, p = res$p)
info_plot <- info_kept[, c("marker_id", "chr", "pos")]
names(info_plot) <- c("snp_id", "chr", "pos")

p_man <- manhattan_plot(assoc_plot, info_plot, threshold = thresh,
                        title = paste("Pan-GWAS —", trait))
save_nature_figure(p_man, file.path(outdir, "manhattan_pan"), width_mm = 180, height_mm = 80)

p_qq <- qq_plot(assoc_plot, title = paste("Pan-GWAS QQ —", trait))
save_nature_figure(p_qq, file.path(outdir, "qq_pan"), width_mm = 89, height_mm = 89)

# ---------------------------------------------------------------------------
# 6) Run info
# ---------------------------------------------------------------------------
write_run_info_r(outdir, "04_pan_gwas.R", params = list(
  pav = pav_path, snp_info = snp_info_path, phenotype = phenotype_path,
  trait = trait, n_markers_tested = nrow(res), n_significant = nrow(sig),
  threshold = thresh, npc = npc
))
message("Done. Outputs written to: ", normalizePath(outdir))
