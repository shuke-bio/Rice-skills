#!/usr/bin/env Rscript
###############################################################################
# 11_multitrait_gs.R -- per-trait GBLUP + multi-trait GEBV correlation
#
# 说明：这不是完整 MTM（多性状混合模型），而是“逐性状 GBLUP + GEBV 相关”
# 的轻量汇总；真正的多性状联合模型需 BGLR/MTM。
#
# 用法：
#   Rscript scripts/phenomics/11_multitrait_gs.R --genotype geno.csv \
#       --phenotype phe.csv --traits trait1,trait2
#   Rscript scripts/phenomics/11_multitrait_gs.R --demo --outdir results/phenomics/11_multitrait_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/phenomics/11_multitrait")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function(n = 300, m = 800) {
  set.seed(42)
  G <- matrix(rbinom(n * m, 2, 0.3), n, m)
  rownames(G) <- sprintf("S%03d", 1:n); colnames(G) <- sprintf("snp%04d", 1:m)
  eff1 <- rnorm(20, 0, 0.25); eff2 <- rnorm(20, 0, 0.25)
  g1 <- as.numeric(G[, 1:20] %*% eff1)
  g2 <- 0.6 * g1 + as.numeric(G[, 21:40] %*% eff2)
  y1 <- g1 + rnorm(n, 0, 0.8)
  y2 <- g2 + rnorm(n, 0, 0.8)
  geno <- data.frame(sample_id = rownames(G), G, check.names = FALSE)
  phe <- data.frame(sample_id = rownames(G), trait1 = y1, trait2 = y2)
  list(geno = geno, phe = phe)
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$geno, file.path(outdir, "demo_genotype.csv"), row.names = FALSE)
  write.csv(d$phe, file.path(outdir, "demo_phenotype.csv"), row.names = FALSE)
  geno_path <- file.path(outdir, "demo_genotype.csv"); phe_path <- file.path(outdir, "demo_phenotype.csv")
  traits <- c("trait1", "trait2")
} else {
  geno_path <- opt("genotype"); phe_path <- opt("phenotype")
  traits <- strsplit(opt("traits", ""), ",")[[1]]
}
if (is.null(geno_path) || is.null(phe_path) || !length(traits)) stop("Need --genotype --phenotype --traits, or --demo")

g <- read.csv(geno_path, check.names = FALSE); p <- read.csv(phe_path, check.names = FALSE)
if (!"sample_id" %in% colnames(g) || !"sample_id" %in% colnames(p)) stop("need sample_id")
d <- merge(g, p, by = "sample_id")
snps <- setdiff(colnames(g), "sample_id")
X <- as.matrix(d[, snps]); storage.mode(X) <- "double"
for (j in seq_len(ncol(X))) { v <- X[, j]; v[is.na(v)] <- mean(v, na.rm = TRUE); X[, j] <- v }
X <- scale(X, center = TRUE, scale = FALSE)
K <- tcrossprod(X) / ncol(X); K <- K / mean(diag(K))
n <- nrow(K)
k_folds <- as.integer(opt("k-folds", 5))
repeats <- as.integer(opt("repeats", 1))
lambdas <- c(0.05, 0.1, 0.3, 1, 3)

acc <- list(); gebv <- list()
for (tr in traits) {
  if (!tr %in% colnames(d)) next
  y <- scale(d[[tr]])
  best_r <- -Inf; best_lam <- NA
  for (rep_i in seq_len(repeats)) {
    folds <- sample(rep(seq_len(k_folds), length.out = n))
    for (lam in lambdas) {
      rs <- c()
      for (fold in seq_len(k_folds)) {
        tr_idx <- which(folds != fold); te_idx <- which(folds == fold)
        mu <- mean(y[tr_idx])
        w <- tryCatch(solve(K[tr_idx, tr_idx] + lam * diag(length(tr_idx)), y[tr_idx] - mu), error = function(e) NULL)
        if (is.null(w)) next
        pred <- as.numeric(K[te_idx, tr_idx] %*% w) + mu
        rs <- c(rs, suppressWarnings(cor(pred, y[te_idx])))
      }
      r <- mean(rs, na.rm = TRUE)
      if (is.finite(r) && r > best_r) { best_r <- r; best_lam <- lam }
    }
  }
  w <- solve(K + best_lam * diag(n), y)
  gebv[[tr]] <- as.numeric(K %*% w)
  acc[[tr]] <- data.frame(trait = tr, best_lambda = best_lam, accuracy = round(best_r, 3))
}
acc_df <- do.call(rbind, acc)
gebv_df <- data.frame(sample_id = d$sample_id, as.data.frame(gebv), check.names = FALSE)
write.csv(acc_df, file.path(outdir, "multitrait_accuracy.csv"), row.names = FALSE)
write.csv(gebv_df, file.path(outdir, "gebv_multitrait.csv"), row.names = FALSE)
if (length(gebv) >= 2) {
  cm <- cor(do.call(cbind, gebv), use = "pairwise.complete.obs")
  write.csv(data.frame(trait = rownames(cm), cm, check.names = FALSE),
            file.path(outdir, "gebv_correlation.csv"), row.names = FALSE)
}
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("traits=", paste(traits, collapse = ",")),
  paste0("samples=", n), "model=per-trait GBLUP; GEBV correlation (not joint MTM)"
), file.path(outdir, "run_info.txt"))
cat("[multitrait] traits:", paste(acc_df$trait, acc_df$accuracy, sep = "=", collapse = ", "), "\n")
cat("[multitrait] outputs ->", outdir, "\n")
