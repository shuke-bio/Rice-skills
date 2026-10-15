###############################################################################
# 03_gs_gblup.R -- genomic selection via GBLUP with cross-validated accuracy
#
# Title:    Rice-skills genomic selection step
# Description:
#   GBLUP (genomic best linear unbiased prediction) implemented with the
#   VanRaden genomic relationship matrix already available in this package:
#     y = mu + Zg + e,  g ~ N(0, G sigma2_g),  e ~ N(0, sigma2_e I)
#   solved in the equivalent ridge-regression form
#     w = (K_train + lambda I)^-1 y_train,  yhat_test = K_test_train w.
#   The shrinkage lambda is chosen automatically from a small grid by
#   k-fold cross-validation (repeated), and the final GEBVs are produced
#   with the best lambda on all data. Reports prediction accuracy
#   (correlation between GEBV and phenotype in held-out folds).
#
# Usage:
#   Rscript scripts/phenomics/03_gs_gblup.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --phenotype demo_data/genomics/phenotype.csv \
#       --trait salt_tolerance_index \
#       --k-folds 5 --repeats 3 \
#       --outdir results/phenomics/03_gs
#
# Inputs:
#   --genotype   SNP x sample dosage matrix CSV (0/1/2; imputed preferred)
#   --phenotype  phenotype table CSV (sample_id + traits); blups.csv from
#                01_phenotype_blup.R is the recommended input
#   --trait      trait column to predict
#   --k-folds    CV folds (default 5)
#   --repeats    CV repetitions (default 3)
#   --lambdas    comma-separated shrinkage grid (default 0.05,0.1,0.3,1,3,10)
#   --n-markers-cap  optional: random-cap markers for speed (default 20000)
#   --outdir     output directory (default results/phenomics/03_gs)
# Outputs:
#   gebv_all.csv            sample_id, phenotype, GEBV (final model)
#   cv_accuracy.csv         per fold x repeat x lambda: cor, rmse
#   accuracy_summary.txt    best lambda + mean accuracy +/- sd
#   gebv_scatter.*          GEBV vs observed phenotype
#   cv_accuracy_box.*       CV accuracy by lambda
#   run_info.txt
#
# Notes:
#   - Pure base-R implementation (no rrBLUP/BGLR needed); results are
#     numerically equivalent to rrBLUP's kin.blup for the same K and lambda.
#   - Missing phenotypes are dropped; genotypes should be imputed (0/1/2).
#   - Accuracy on the demo data is expected around 0.4-0.6 for the
#     salt_tolerance_index trait (oligogenic architecture).
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
phenotype_path <- opt(opts, "phenotype")
trait <- opt(opts, "trait")
k_folds <- suppressWarnings(as.integer(opt(opts, "k-folds", "5")))
repeats <- suppressWarnings(as.integer(opt(opts, "repeats", "3")))
lambdas_arg <- opt(opts, "lambdas")
n_cap <- suppressWarnings(as.integer(opt(opts, "n-markers-cap", "20000")))
outdir <- opt(opts, "outdir", "results/phenomics/03_gs")

if (any(vapply(list(genotype_path, phenotype_path, trait), is.null, logical(1)))) {
  stop("Required: --genotype, --phenotype, --trait. See --help.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills GBLUP genomic selection ==")

lambdas <- if (!is.null(lambdas_arg)) {
  as.numeric(trimws(strsplit(lambdas_arg, ",")[[1]]))
} else c(0.05, 0.1, 0.3, 1, 3, 10)

# ---------------------------------------------------------------------------
# 1) Load + align
# ---------------------------------------------------------------------------
gt <- read_genotype(genotype_path)
G0 <- gt$geno
pheno <- read_pheno(phenotype_path)
if (!trait %in% colnames(pheno)) {
  stop("Trait '", trait, "' not found. Available: ",
       paste(setdiff(colnames(pheno), "sample_id"), collapse = ", "), call. = FALSE)
}
shared <- intersect(rownames(G0), pheno$sample_id)
if (length(shared) < 20) stop("Fewer than 20 samples with phenotype + genotype.", call. = FALSE)
G0 <- G0[shared, , drop = FALSE]
y_all <- suppressWarnings(as.numeric(pheno[[trait]][match(shared, pheno$sample_id)]))
ok <- !is.na(y_all)
G0 <- G0[ok, , drop = FALSE]; y_all <- y_all[ok]; shared <- shared[ok]
n <- length(shared)
if (ncol(G0) > n_cap) {
  keep_m <- sample.int(ncol(G0), n_cap)
  G0 <- G0[, keep_m, drop = FALSE]
  message(sprintf("  markers capped at %d for speed", n_cap))
}
message(sprintf("  samples: %d | markers: %d | trait: %s", n, ncol(G0), trait))

# ---------------------------------------------------------------------------
# 2) Genomic relationship matrix (VanRaden), diagonal-normalised
# ---------------------------------------------------------------------------
K <- vanraden_kinship(G0)
K <- K / mean(diag(K))
y_scaled <- scale(y_all)

# ---------------------------------------------------------------------------
# 3) Cross-validated lambda selection + accuracy
# ---------------------------------------------------------------------------
cv_rows <- list()
for (rep_i in seq_len(repeats)) {
  folds <- sample(rep(seq_len(k_folds), length.out = n))
  for (fold in seq_len(k_folds)) {
    train <- which(folds != fold); test <- which(folds == fold)
    K_tt <- K[train, train, drop = FALSE]
    K_te <- K[test, train, drop = FALSE]
    y_tr <- y_scaled[train]
    mu_tr <- mean(y_tr)
    for (lam in lambdas) {
      w <- tryCatch(solve(K_tt + lam * diag(length(train)), y_tr - mu_tr),
                    error = function(e) NULL)
      if (is.null(w)) next
      pred <- as.numeric(K_te %*% w) + mu_tr
      r <- suppressWarnings(stats::cor(pred, y_scaled[test]))
      rmse <- sqrt(mean((pred - y_scaled[test])^2))
      cv_rows[[length(cv_rows) + 1L]] <- data.frame(
        rep_id = rep_i, fold = fold, lambda = lam,
        accuracy = round(r, 4), rmse = round(rmse, 4))
    }
  }
}
cv <- do.call(rbind, cv_rows)
utils::write.csv(cv, file.path(outdir, "cv_accuracy.csv"), row.names = FALSE)

acc_summary <- aggregate(accuracy ~ lambda, data = cv, FUN = function(v) {
  c(mean_acc = round(mean(v), 4), sd_acc = round(stats::sd(v), 4))
})
acc_summary <- do.call(data.frame, acc_summary)
colnames(acc_summary) <- c("lambda", "mean_accuracy", "sd_accuracy")
best <- acc_summary$lambda[which.max(acc_summary$mean_accuracy)]
best_acc <- acc_summary$mean_accuracy[which.max(acc_summary$mean_accuracy)]
utils::write.csv(acc_summary, file.path(outdir, "cv_lambda_summary.csv"), row.names = FALSE)
message(sprintf("  CV: best lambda = %g (mean accuracy %.3f +/- %.3f)",
                best, best_acc,
                acc_summary$sd_accuracy[which.max(acc_summary$mean_accuracy)]))

# ---------------------------------------------------------------------------
# 4) Final model on all data
# ---------------------------------------------------------------------------
w <- solve(K + best * diag(n), y_scaled)
gebv <- as.numeric(K %*% w)   # fitted = K * w (standardised scale)
out <- data.frame(sample_id = shared, phenotype = y_all,
                  gebv = round(gebv, 4), stringsAsFactors = FALSE)
utils::write.csv(out, file.path(outdir, "gebv_all.csv"), row.names = FALSE)

pal <- rice_palette()
p_sc <- ggplot2::ggplot(out, ggplot2::aes(x = gebv, y = phenotype)) +
  ggplot2::geom_point(size = 1.3, alpha = 0.75, colour = pal$groups5[1]) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, colour = pal$groups5[2],
                       linewidth = 0.4, fill = "grey85") +
  ggplot2::labs(x = "GEBV (standardised)", y = trait,
                title = sprintf("GBLUP fit — CV accuracy %.3f", best_acc)) +
  theme_nature()
save_nature_figure(p_sc, file.path(outdir, "gebv_scatter"), width_mm = 89, height_mm = 75)

p_cv <- ggplot2::ggplot(cv, ggplot2::aes(x = factor(lambda), y = accuracy)) +
  ggplot2::geom_boxplot(fill = pal$groups5[3], alpha = 0.85, linewidth = 0.3,
                        outlier.size = 0.5, colour = "grey25") +
  ggplot2::labs(x = expression(lambda), y = "CV prediction accuracy (r)",
                title = "Accuracy by shrinkage") +
  theme_nature()
save_nature_figure(p_cv, file.path(outdir, "cv_accuracy_box"), width_mm = 89, height_mm = 70)

sum_lines <- c(
  "Rice-skills -- GBLUP summary",
  sprintf("samples / markers : %d / %d", n, ncol(G0)),
  sprintf("trait             : %s", trait),
  sprintf("CV                : %d-fold x %d repeats", k_folds, repeats),
  sprintf("best lambda       : %g", best),
  sprintf("prediction acc.   : %.3f +/- %.3f (r, held-out)", best_acc,
          acc_summary$sd_accuracy[which.max(acc_summary$mean_accuracy)])
)
writeLines(sum_lines, file.path(outdir, "accuracy_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

write_run_info_r(outdir, "03_gs_gblup.R", params = list(
  genotype = genotype_path, phenotype = phenotype_path, trait = trait,
  k_folds = k_folds, repeats = repeats, best_lambda = best,
  accuracy = best_acc, n_samples = n, n_markers = ncol(G0)
))
message("Done. Outputs written to: ", normalizePath(outdir))
