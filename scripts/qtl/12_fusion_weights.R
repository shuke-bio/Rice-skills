#!/usr/bin/env Rscript
###############################################################################
# 12_fusion_weights.R -- train FUSION-like cis-eQTL weights (ridge/elastic net)
#
# 输入：
#   --genotype CSV: sample_id, snp1, ...（0/1/2）
#   --expression CSV: sample_id, gene1, ...（log/TPM）
#   --gene-info 可选 CSV: gene_id, chr, start, end（用于 cis 窗口）
#   --genes 可选逗号分隔；默认全部基因
# 方法：
#   每个基因取 cis SNP（窗口内）-> 标准化 -> ridge/elastic net（glmnet 可选）
#   -> 交叉验证选 lambda -> 输出权重与 R²
#
# 用法：
#   Rscript scripts/qtl/12_fusion_weights.R --genotype geno.csv --expression expr.csv \
#       --gene-info genes.csv --genes G1,G2
#   Rscript scripts/qtl/12_fusion_weights.R --demo --outdir results/qtl/12_fusion_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/qtl/12_fusion")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
window <- as.numeric(opt("window", 1e6))
k_folds <- as.integer(opt("k-folds", 5))

make_demo <- function(n = 200, m = 40) {
  set.seed(42)
  G <- matrix(rbinom(n * m, 2, 0.3), n, m)
  colnames(G) <- sprintf("snp%03d", 1:m); rownames(G) <- sprintf("S%03d", 1:n)
  e1 <- 0.8 * G[, 5] + rnorm(n, 0, 0.5)
  e2 <- rnorm(n)
  info <- data.frame(gene_id = c("G1", "G2"), chr = "chr1",
                     start = c(9e6, 50e6), end = c(9.2e6, 50.2e6))
  list(geno = data.frame(sample_id = rownames(G), G, check.names = FALSE),
       expr = data.frame(sample_id = rownames(G), G1 = e1, G2 = e2),
       info = info)
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$geno, file.path(outdir, "demo_genotype.csv"), row.names = FALSE)
  write.csv(d$expr, file.path(outdir, "demo_expression.csv"), row.names = FALSE)
  write.csv(d$info, file.path(outdir, "demo_gene_info.csv"), row.names = FALSE)
  geno_path <- file.path(outdir, "demo_genotype.csv"); expr_path <- file.path(outdir, "demo_expression.csv")
  info_path <- file.path(outdir, "demo_gene_info.csv"); genes <- c("G1", "G2")
} else {
  geno_path <- opt("genotype"); expr_path <- opt("expression")
  info_path <- opt("gene-info"); genes <- strsplit(opt("genes", ""), ",")[[1]]
}
if (is.null(geno_path) || is.null(expr_path)) stop("Need --genotype and --expression, or --demo")
g <- read.csv(geno_path, check.names = FALSE); e <- read.csv(expr_path, check.names = FALSE)
if (!"sample_id" %in% colnames(g) || !"sample_id" %in% colnames(e)) stop("need sample_id")
d <- merge(g, e, by = "sample_id")
snps <- setdiff(colnames(g), "sample_id")
all_genes <- setdiff(colnames(e), "sample_id")
if (!length(genes) || identical(genes, "")) genes <- all_genes
info <- if (!is.null(info_path) && file.exists(info_path)) read.csv(info_path, check.names = FALSE) else NULL

X <- as.matrix(d[, snps]); storage.mode(X) <- "double"
for (j in seq_len(ncol(X))) { v <- X[, j]; v[is.na(v)] <- mean(v, na.rm = TRUE); X[, j] <- v }

ridge_fit <- function(Xtr, ytr, lam) {
  Xc <- scale(Xtr, center = TRUE, scale = FALSE)
  yc <- ytr - mean(ytr)
  w <- solve(crossprod(Xc) + lam * diag(ncol(Xc)), crossprod(Xc, yc))
  list(w = as.numeric(w), mu = mean(ytr), xm = attr(scale(Xtr), "scaled:center"))
}
ridge_pred <- function(fit, Xte) as.numeric(scale(Xte, center = fit$xm, scale = FALSE) %*% fit$w + fit$mu)

weight_rows <- list(); model_rows <- list()
for (gene in genes) {
  if (!gene %in% colnames(d)) next
  y <- as.numeric(d[[gene]])
  keep <- is.finite(y)
  Xg <- X[keep, , drop = FALSE]; y <- y[keep]
  use <- rep(TRUE, ncol(Xg))
  if (!is.null(info)) {
    gi <- info[info$gene_id == gene, , drop = FALSE]
    if (nrow(gi) && "pos" %in% colnames(info)) {
      # info has no pos; use start/end midpoint
      mid <- (gi$start[1] + gi$end[1]) / 2
      use <- abs(info$pos - mid) <= window
    }
  }
  if (sum(use) < 2) use <- rep(TRUE, ncol(Xg))
  Xg <- Xg[, use, drop = FALSE]
  snp_use <- snps[use]
  lambdas <- c(0.01, 0.1, 1, 10, 100)
  folds <- sample(rep(seq_len(k_folds), length.out = length(y)))
  best <- list(r2 = -Inf, lam = NA)
  for (lam in lambdas) {
    pred <- rep(NA_real_, length(y))
    for (f in seq_len(k_folds)) {
      tr <- which(folds != f); te <- which(folds == f)
      fit <- tryCatch(ridge_fit(Xg[tr, , drop = FALSE], y[tr], lam), error = function(e) NULL)
      if (!is.null(fit)) pred[te] <- ridge_pred(fit, Xg[te, , drop = FALSE])
    }
    r2 <- suppressWarnings(cor(pred, y, use = "complete.obs")^2)
    if (is.finite(r2) && r2 > best$r2) best <- list(r2 = r2, lam = lam)
  }
  fit <- ridge_fit(Xg, y, best$lam)
  weight_rows[[gene]] <- data.frame(gene_id = gene, snp_id = snp_use, weight = round(fit$w, 5))
  model_rows[[gene]] <- data.frame(gene_id = gene, n_snp = length(snp_use),
                                    lambda = best$lam, r2 = round(best$r2, 4))
}
w <- do.call(rbind, weight_rows); m <- do.call(rbind, model_rows)
write.csv(w, file.path(outdir, "fusion_weights.csv"), row.names = FALSE)
write.csv(m, file.path(outdir, "fusion_models.csv"), row.names = FALSE)
writeLines(c(paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
             paste0("genes=", nrow(m)), paste0("window=", window),
             "method=ridge (closed form); glmnet not required"), file.path(outdir, "run_info.txt"))
cat(sprintf("[fusion] genes=%d; mean R2=%.3f\n", nrow(m), mean(m$r2)))
print(m)
cat("[fusion] outputs ->", outdir, "\n")
