###############################################################################
# 03_association.R -- GWAS association analysis (three engines)
#
# Title:    Rice-skills GWAS association step
# Description:
#   Runs genome-wide association on a trait using one of three engines:
#     simple : built-in per-SNP OLS ("GLM" family) with PCA covariates; runs
#              anywhere with no extra packages (exploratory baseline);
#     rmvp   : rMVP (Memory-efficient, Visualization-enhanced, Parallel-
#              accelerated; Apache-2.0) with GLM / MLM / FarmCPU models;
#     gapit  : GAPIT (GLM / MLM / CMLM / SUPER / MLMM / FarmCPU / Blink;
#              GPL-2.0+, installed at runtime).
#   In auto mode the best installed engine is chosen (rmvp > gapit > simple).
#
# Usage:
#   Rscript scripts/gwas/03_association.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --phenotype demo_data/genomics/phenotype.csv \
#       --trait salt_tolerance_index --outdir results/gwas/03_association
#       # --engine gapit --model FarmCPU   # GAPIT model choice (gapit engine only)
#
# Inputs:
#   --genotype   SNP x sample dosage matrix CSV (0/1/2, missing = empty)
#   --snp-info   snp_id, chr, pos
#   --phenotype  sample_id + trait columns
# Outputs:
#   association_results.csv      all tested SNPs (snp_id, beta, se, p, engine)
#   association_significant.csv  SNPs passing the significance threshold
#   pca_scores.csv               PCs used as covariates
#   qc_summary.txt / run_info.txt
#
# Notes:
#   - The simple engine imputes missing genotypes to the per-SNP mean dosage
#     and is a legitimate "GLM" baseline; for publication-grade mixed models
#     install rMVP (install.packages("rMVP")) or GAPIT.
#   - Significance threshold defaults to Bonferroni 0.05 / n_SNPs.
#   - --model selects the association model: rMVP supports GLM / MLM / FarmCPU
#     (default FarmCPU); GAPIT additionally supports CMLM / SUPER / MLMM / Blink
#     (default MLM). The simple engine is a per-SNP OLS-GLM baseline and ignores
#     --model. Model choice guidance: docs/references/gapit_gwas_formats.md.
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
phenotype_path <- opt(opts, "phenotype")
trait <- opt(opts, "trait")
outdir <- opt(opts, "outdir", "results/gwas/03_association")
engine <- tolower(opt(opts, "engine", "auto"))
npc <- suppressWarnings(as.integer(opt(opts, "npc", "5")))
maf_cut <- suppressWarnings(as.numeric(opt(opts, "maf", "0.05")))
missing_cut <- suppressWarnings(as.numeric(opt(opts, "missing", "0.1")))
p_threshold <- opt(opts, "p-threshold", NULL)
gapit_models <- c("GLM", "MLM", "CMLM", "SUPER", "MLMM", "FarmCPU", "Blink")
rmvp_models <- c("GLM", "MLM", "FarmCPU")
model_raw <- opt(opts, "model", if (identical(engine, "rmvp")) "FarmCPU" else "MLM")
allowed_models <- if (identical(engine, "rmvp")) rmvp_models else gapit_models
m <- which(tolower(allowed_models) == tolower(model_raw))
if (length(m) == 0) {
  stop("--model must be one of ", paste(allowed_models, collapse = " / "),
       " for engine '", engine, "' (case-insensitive); got '", model_raw, "'.",
       if (identical(engine, "rmvp"))
         "\nHint: CMLM / SUPER / MLMM / Blink require --engine gapit." else "",
       call. = FALSE)
}
model <- allowed_models[m]

if (is.null(genotype_path) || is.null(snp_info_path) || is.null(phenotype_path) || is.null(trait)) {
  stop("Required: --genotype, --snp-info, --phenotype, --trait.\n",
       "Hint: see the Usage block at the top of this script.", call. = FALSE)
}
if (!engine %in% c("auto", "simple", "rmvp", "gapit")) {
  stop("--engine must be one of auto / simple / rmvp / gapit; got '", engine, "'.", call. = FALSE)
}
if (is.na(npc) || npc < 0 || npc > 50) stop("--npc must be between 0 and 50.", call. = FALSE)

need_pkgs(c("ggplot2", "jsonlite"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1) Load and QC
# ---------------------------------------------------------------------------
gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
pheno <- read_pheno(phenotype_path)
if (!trait %in% colnames(pheno)) {
  stop("Trait '", trait, "' not found in the phenotype table. Available: ",
       paste(setdiff(colnames(pheno), "sample_id"), collapse = ", "), call. = FALSE)
}
al <- align_samples_geno(gt$geno, pheno)
message(sprintf("Samples: %d shared between genotype and phenotype.", nrow(al$geno)))

fl <- filter_snps(al$geno, snp_info, min_maf = maf_cut, max_missing = missing_cut)
message(sprintf("SNP filter (MAF >= %g, missing <= %g): %d -> %d SNPs.",
                maf_cut, missing_cut, ncol(al$geno), ncol(fl$geno)))
if (ncol(fl$geno) < 10) {
  stop("Fewer than 10 SNPs remain after QC. Hint: relax --maf / --missing, ",
       "or check that the genotype matrix matches the SNP info table.", call. = FALSE)
}

y <- al$pheno[[trait]]
if (all(is.na(y))) stop("Trait '", trait, "' is entirely missing.", call. = FALSE)
message(sprintf("Trait '%s': n = %d non-missing, range [%s].",
                trait, sum(!is.na(y)), paste(round(range(y, na.rm = TRUE), 3), collapse = ", ")))

# ---------------------------------------------------------------------------
# 2) PCA covariates
# ---------------------------------------------------------------------------
pca_res <- NULL
covariates <- NULL
if (npc > 0) {
  pca_res <- compute_pca(fl$geno, npc = npc)
  covariates <- pca_res$scores[, seq_len(min(npc, ncol(pca_res$scores))), drop = FALSE]
  write.csv(data.frame(sample_id = rownames(pca_res$scores), pca_res$scores),
            file.path(outdir, "pca_scores.csv"), row.names = FALSE)
  message(sprintf("PCA covariates: %d PCs (PC1 explains %.1f%%).",
                  ncol(covariates), 100 * pca_res$var_explained[1]))
}

# ---------------------------------------------------------------------------
# 3) Association
# ---------------------------------------------------------------------------
run_simple <- function() {
  message("Engine: simple (per-SNP OLS / GLM family, PCA covariates; ",
          "install rMVP or GAPIT for MLM/FarmCPU).")
  res <- gwas_glm(fl$geno, y, covariates)
  res$engine <- "simple-GLM"
  res
}

run_rmvp <- function() {
  need_pkgs("rMVP")
  message("Engine: rMVP (model: ", model, ").")
  # rMVP expects samples x SNPs for the numeric input; write a temp file.
  tmp_geno <- file.path(outdir, "_rmvp_geno.txt")
  tmp_phe <- file.path(outdir, "_rmvp_phe.txt")
  g <- impute_mean(fl$geno)
  write.table(data.frame(sample_id = rownames(g), g),
              tmp_geno, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(sample_id = al$pheno$sample_id, trait = y),
              tmp_phe, sep = "\t", row.names = FALSE, quote = FALSE)
  rmvp_data <- rMVP::MVP.Data(fileNum = tmp_geno, filePhe = tmp_phe,
                              sep.num = "\t", sep.phe = "\t", fileKin = FALSE,
                              filePC = FALSE, out = file.path(outdir, "_mvp"))
  rmvp_obj <- rMVP::MVP(rmvp_data$pheno, rmvp_data$geno, rmvp_data$map,
                        K = if (is.null(covariates)) NULL else NULL,
                        CV = covariates, nPC = 0L,
                        method = model, file.output = FALSE)
  res <- as.data.frame(rmvp_obj$map)
  names(res)[names(res) == "SNP"] <- "snp_id"
  if (!"p" %in% names(res)) {
    pcol <- grep("^p", names(res), value = TRUE)[1]
    res$p <- res[[pcol]]
  }
  res$beta <- if ("effect" %in% names(res)) res$effect else NA_real_
  res$se <- NA_real_
  res$engine <- paste0("rMVP-", model)
  res[, c("snp_id", "beta", "se", "p", "engine")]
}

run_gapit <- function() {
  need_pkgs("GAPIT")
  message("Engine: GAPIT (model = ", model, "). GAPIT writes files into the ",
          "working directory; temporary switch to outdir.")
  g <- impute_mean(fl$geno)
  GD <- data.frame(taxa = rownames(g), g, check.names = FALSE)
  GM <- data.frame(SNP = fl$snp_info$snp_id, Chromosome = fl$snp_info$chr,
                   Position = fl$snp_info$pos, check.names = FALSE)
  Y <- data.frame(Taxa = al$pheno$sample_id, trait = y, check.names = FALSE)
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(outdir)
  gapit_res <- GAPIT::GAPIT(Y = Y, GD = GD, GM = GM, PCA.total = npc,
                            model = model, file.output = FALSE)
  res <- as.data.frame(gapit_res$GWAS)
  data.frame(snp_id = res$SNP, beta = if ("effect" %in% names(res)) res$effect else NA_real_,
             se = NA_real_, p = res$P.value, engine = paste0("GAPIT-", model),
             stringsAsFactors = FALSE)
}

chosen <- engine
if (engine == "auto") {
  chosen <- if (requireNamespace("rMVP", quietly = TRUE)) "rmvp"
            else if (requireNamespace("GAPIT", quietly = TRUE)) "gapit"
            else "simple"
}
if (chosen != "gapit" && !is.null(opts$model)) {
  message("Note: --model only applies to the GAPIT engine; ignoring it for '",
          chosen, "'.")
}
assoc <- switch(chosen,
                simple = run_simple(),
                rmvp = run_rmvp(),
                gapit = run_gapit())
assoc <- assoc[order(assoc$p), , drop = FALSE]
rownames(assoc) <- NULL

if (is.null(p_threshold)) {
  p_threshold <- 0.05 / nrow(assoc)
  message(sprintf("Significance threshold: Bonferroni 0.05/%d = %.3g.",
                  nrow(assoc), p_threshold))
} else {
  p_threshold <- as.numeric(p_threshold)
  message("Significance threshold (user): ", p_threshold)
}

sig <- assoc[!is.na(assoc$p) & assoc$p <= p_threshold, , drop = FALSE]
sig <- merge(sig, fl$snp_info, by = "snp_id", all.x = TRUE)
sig <- sig[order(sig$p), , drop = FALSE]

utils::write.csv(assoc, file.path(outdir, "association_results.csv"),
                 row.names = FALSE, quote = FALSE)
utils::write.csv(sig, file.path(outdir, "association_significant.csv"),
                 row.names = FALSE, quote = FALSE)
writeLines(c(
  sprintf("trait: %s", trait),
  sprintf("samples: %d", nrow(fl$geno)),
  sprintf("snps tested: %d (after MAF %g / missing %g filter)", nrow(assoc), maf_cut, missing_cut),
  sprintf("engine: %s", chosen),
  sprintf("model: %s", if (chosen == "gapit") model else "(engine default)"),
  sprintf("covariates: %d PCs", if (is.null(covariates)) 0L else ncol(covariates)),
  sprintf("significance threshold: %.3g", p_threshold),
  sprintf("significant SNPs: %d", nrow(sig)),
  if (nrow(sig) > 0) paste0("top hit: ", sig$snp_id[1],
                            " (chr", sig$chr[1], ":", sig$pos[1],
                            ", p = ", format(sig$p[1], digits = 3), ")") else
    "no SNP passed the threshold; consider a less stringent cutoff or more samples"
), file.path(outdir, "qc_summary.txt"))

write_run_info_r(
  outdir, "03_association.R",
  params = list(genotype = genotype_path, snp_info = snp_info_path,
                phenotype = phenotype_path, trait = trait, engine = chosen,
                model = if (chosen == "gapit") model else NA_character_,
                npc = npc, maf = maf_cut, missing = missing_cut,
                p_threshold = p_threshold, seed = RICE_SEED),
  packages = c("ggplot2", "rMVP", "GAPIT")
)
message(sprintf("Association complete (%s): %d SNPs, %d significant. Results in %s",
                chosen, nrow(assoc), nrow(sig), outdir))
