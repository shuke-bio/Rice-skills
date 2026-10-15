###############################################################################
# 02_plink_to_matrix.R -- convert a PLINK .raw dosage export to the
#                        Rice-skills genotype matrix format
#
# Title:    Rice-skills PLINK-to-matrix converter
# Description:
#   Reads the .raw file produced by `plink --recode A` (additive dosages) plus
#   the matching .map file and writes the standard Rice-skills genomics
#   input files:
#       genotype.csv   SNP x sample dosage matrix (0/1/2, missing empty)
#       snp_info.csv   snp_id, chr, pos
#   An optional --prune file (.prune.in) restricts the output to LD-pruned SNPs.
#
# Usage:
#   Rscript scripts/gwas/02_plink_to_matrix.R \
#       --raw results/gwas/01_qc/rice.raw \
#       --map results/gwas/01_qc/rice.map \
#       --prune results/gwas/01_qc/rice.prune.in \
#       --outdir results/gwas/01_qc/matrix
#
# Notes:
#   - .raw column names carry the counted-allele suffix (e.g. s123_A); SNP order
#     matches the .map rows, so positions are matched by order (robust to
#     arbitrary suffixes).
#   - Missing genotypes appear as NA in .raw and stay empty in genotype.csv.
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
raw_path <- opt(opts, "raw")
map_path <- opt(opts, "map")
prune_path <- opt(opts, "prune")
outdir <- opt(opts, "outdir", NULL)
if (is.null(outdir)) {
  outdir <- if (is.null(raw_path)) "." else dirname(raw_path)
}

if (is.null(raw_path) || is.null(map_path)) {
  stop("Required: --raw <plink .raw> and --map <plink .map>.\n",
       "Hint: produce both with scripts/gwas/01_plink_qc.sh.", call. = FALSE)
}
if (!file.exists(raw_path)) stop("File not found: ", raw_path, call. = FALSE)
if (!file.exists(map_path)) stop("File not found: ", map_path, call. = FALSE)

raw <- utils::read.table(raw_path, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
meta_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
if (!all(meta_cols[1:2] %in% colnames(raw)) || ncol(raw) <= 6) {
  stop("Unexpected .raw layout (expected FID/IID/... then SNP dosage columns).",
       call. = FALSE)
}
samples <- as.character(raw$IID)
dosage <- as.matrix(raw[, -(1:6), drop = FALSE])
storage.mode(dosage) <- "double"

map <- utils::read.table(map_path, header = FALSE, stringsAsFactors = FALSE)
if (ncol(map) < 4) stop("Unexpected .map layout (expected chr, snp_id, cm, pos).", call. = FALSE)
map <- data.frame(chr = as.character(map[[1]]), snp_id = as.character(map[[2]]),
                  pos = as.numeric(map[[4]]), stringsAsFactors = FALSE)

if (nrow(map) != ncol(dosage)) {
  stop(sprintf("Column/row mismatch: .raw has %d SNP columns but .map has %d rows.\n",
               ncol(dosage), nrow(map)),
       "Hint: regenerate both with the same PLINK command.", call. = FALSE)
}

# strip the counted-allele suffix from .raw column names for reporting, but use
# the .map IDs as the canonical SNP IDs (order-matched).
if (!is.null(prune_path)) {
  if (!file.exists(prune_path)) stop("Prune file not found: ", prune_path, call. = FALSE)
  keep_ids <- readLines(prune_path, warn = FALSE)
  keep <- map$snp_id %in% keep_ids
  message(sprintf("LD-pruned SNP set: keeping %d of %d SNPs.", sum(keep), nrow(map)))
  map <- map[keep, , drop = FALSE]
  dosage <- dosage[, keep, drop = FALSE]
}

# genotype.csv: SNP x sample
mat <- t(dosage)
colnames(mat) <- samples
geno_df <- data.frame(snp_id = map$snp_id, mat, check.names = FALSE)

save_table <- function(df, path, row_names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = row_names, na = "", quote = FALSE)
  invisible(path)
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save_table(geno_df, file.path(outdir, "genotype.csv"))
save_table(data.frame(snp_id = map$snp_id, chr = map$chr, pos = map$pos),
           file.path(outdir, "snp_info.csv"))

miss_rate <- mean(is.na(mat))
message(sprintf("Converted %d SNPs x %d samples (missing %.2f%%).",
                nrow(map), length(samples), 100 * miss_rate))
message("Files written:")
message("  ", file.path(outdir, "genotype.csv"))
message("  ", file.path(outdir, "snp_info.csv"))
message("Next: Rscript scripts/gwas/03_association.R --genotype ",
        file.path(outdir, "genotype.csv"), " --snp-info ",
        file.path(outdir, "snp_info.csv"), " --phenotype <pheno.csv> --trait <trait>")

write_run_info_r(
  outdir, "02_plink_to_matrix.R",
  params = list(raw = raw_path, map = map_path,
                prune = if (is.null(prune_path)) NA else prune_path,
                n_snps = nrow(map), n_samples = length(samples),
                missing_rate = round(miss_rate, 4), seed = RICE_SEED),
  packages = c("ggplot2")
)
