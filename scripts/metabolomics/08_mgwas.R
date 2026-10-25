#!/usr/bin/env Rscript
###############################################################################
# 08_mgwas.R -- metabolite GWAS / mQTL (metabolite x genotype association)
#
# 模式：
#   gwas : 自然群体，PCA + VanRaden 亲缘矩阵 + 特征分解 MLM（固定方差比）
#   qtl  : RIL/双亲群体，标记回归，LOD = -log10(p)，可选置换阈值
#
# 输入：
#   --metabolites CSV: sample_id, metabolite1, metabolite2, ...
#   --genotype    CSV: sample_id, snp1, snp2, ...（0/1/2 剂量，缺失留空）
#   --snp-info    CSV: snp_id, chr, pos
#   --covariates  可选 CSV: sample_id + 协变量列
#
# 用法：
#   Rscript scripts/metabolomics/08_mgwas.R --metabolites met.csv --genotype geno.csv \
#       --snp-info snp_info.csv --mode gwas --outdir results/metabolomics/08_mgwas
#   Rscript scripts/metabolomics/08_mgwas.R --demo --outdir results/metabolomics/08_mgwas_demo
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/metabolomics/08_mgwas")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
mode <- tolower(opt("mode", "gwas"))
if (!mode %in% c("gwas", "qtl")) stop("--mode must be gwas or qtl")

vanraden <- function(G) {
  G <- scale(G, center = TRUE, scale = FALSE)
  tcrossprod(G) / ncol(G)
}

mlm_scan <- function(y, G, Z, K, delta = 1e-6) {
  eig <- eigen(K, symmetric = TRUE)
  U <- eig$vectors
  lam <- pmax(eig$values, 0)
  scale_vec <- 1 / sqrt(lam + delta)
  yt <- as.numeric(crossprod(U, y) * scale_vec)
  out <- matrix(NA_real_, ncol(G), 4)
  for (j in seq_len(ncol(G))) {
    D <- cbind(Z, G[, j])
    Dt <- (crossprod(U, D)) * scale_vec
    fit <- lm.fit(Dt, yt)
    if (any(is.na(fit$coefficients))) next
    df <- nrow(Dt) - ncol(Dt)
    if (df <= 0) next
    rss <- sum(fit$residuals^2)
    sigma2 <- rss / df
    XtX_inv <- tryCatch(solve(crossprod(Dt)), error = function(e) NULL)
    if (is.null(XtX_inv)) next
    se <- sqrt(diag(XtX_inv) * sigma2)
    b <- fit$coefficients[ncol(D)]
    s <- se[ncol(D)]
    t <- b / s
    out[j, ] <- c(b, s, t, 2 * pt(-abs(t), df))
  }
  out
}

glm_scan <- function(y, G, Z) {
  out <- matrix(NA_real_, ncol(G), 3)
  for (j in seq_len(ncol(G))) {
    D <- cbind(Z, G[, j])
    fit <- lm.fit(D, y)
    df <- nrow(D) - ncol(D)
    rss <- sum(fit$residuals^2)
    XtX_inv <- tryCatch(solve(crossprod(D)), error = function(e) NULL)
    if (is.null(XtX_inv) || df <= 0) next
    se <- sqrt(diag(XtX_inv) * rss / df)
    b <- fit$coefficients[ncol(D)]
    out[j, ] <- c(b, se[ncol(D)], 2 * pt(-abs(b / se[ncol(D)]), df))
  }
  out
}

make_demo <- function() {
  set.seed(42)
  n <- 80; m <- 400; k <- 30
  G <- matrix(rbinom(n * m, 2, 0.3), n, m)
  colnames(G) <- sprintf("snp%03d", 1:m)
  rownames(G) <- sprintf("S%02d", 1:n)
  # MAF filter
  maf <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)
  G <- G[, maf >= 0.05, drop = FALSE]
  met <- matrix(rnorm(n * k, 0, 1), n, k)
  colnames(met) <- sprintf("M%02d", 1:k)
  rownames(met) <- rownames(G)
  pairs <- list(c(1, 10, 0.90), c(2, 100, 0.85), c(3, 250, 0.95))
  for (p in pairs) {
    mi <- as.integer(p[1]); sj <- as.integer(p[2]); eff <- as.numeric(p[3])
    met[, mi] <- met[, mi] + eff * G[, sj]
  }
  metdf <- data.frame(sample_id = rownames(met), met, check.names = FALSE)
  gdf <- data.frame(sample_id = rownames(G), G, check.names = FALSE)
  info <- data.frame(snp_id = colnames(G), chr = "chr1", pos = seq(1000, by = 1000, length.out = ncol(G)))
  list(met = metdf, geno = gdf, info = info, pairs = pairs)
}

if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$met, file.path(outdir, "demo_metabolites.csv"), row.names = FALSE)
  write.csv(d$geno, file.path(outdir, "demo_genotype.csv"), row.names = FALSE)
  write.csv(d$info, file.path(outdir, "demo_snp_info.csv"), row.names = FALSE)
  met_path <- file.path(outdir, "demo_metabolites.csv")
  geno_path <- file.path(outdir, "demo_genotype.csv")
  info_path <- file.path(outdir, "demo_snp_info.csv")
} else {
  met_path <- opt("metabolites"); geno_path <- opt("genotype"); info_path <- opt("snp-info")
}
if (is.null(met_path) || is.null(geno_path) || is.null(info_path)) {
  stop("Need --metabolites, --genotype, --snp-info, or --demo")
}

met <- read.csv(met_path, check.names = FALSE)
gen <- read.csv(geno_path, check.names = FALSE)
info <- read.csv(info_path, check.names = FALSE)
if (!"sample_id" %in% colnames(met) || !"sample_id" %in% colnames(gen)) {
  stop("metabolites and genotype need a sample_id column")
}
need_info <- c("snp_id", "chr", "pos")
if (length(setdiff(need_info, colnames(info)))) stop("snp_info needs snp_id, chr, pos")

shared <- intersect(met$sample_id, gen$sample_id)
met <- met[match(shared, met$sample_id), , drop = FALSE]
gen <- gen[match(shared, gen$sample_id), , drop = FALSE]
snp_ids <- setdiff(colnames(gen), "sample_id")
G <- as.matrix(gen[, snp_ids, drop = FALSE])
storage.mode(G) <- "double"
# mean impute + MAF
for (j in seq_len(ncol(G))) {
  v <- G[, j]; mu <- mean(v, na.rm = TRUE); v[is.na(v)] <- mu; G[, j] <- v
}
maf <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)
keep <- maf >= as.numeric(opt("maf", 0.05))
G <- G[, keep, drop = FALSE]
snp_ids <- snp_ids[keep]
info <- info[match(snp_ids, info$snp_id), , drop = FALSE]
if (ncol(G) < 10) stop("Fewer than 10 SNPs after MAF filter")

Z <- matrix(1, nrow(G), 1)
if (mode == "gwas") {
  npc <- as.integer(opt("npc", 3))
  pc <- prcomp(G, center = TRUE, scale. = FALSE)$x[, seq_len(min(npc, ncol(G))), drop = FALSE]
  Z <- cbind(Z, pc)
}
cov_path <- opt("covariates")
if (!is.null(cov_path) && !isTRUE(cov_path)) {
  cv <- read.csv(cov_path, check.names = FALSE)
  if (!"sample_id" %in% colnames(cv)) stop("covariates need sample_id")
  cv <- cv[match(shared, cv$sample_id), , drop = FALSE]
  Z <- cbind(Z, as.matrix(cv[, setdiff(colnames(cv), "sample_id"), drop = FALSE]))
}

K <- if (mode == "gwas") vanraden(G) else NULL
met_cols <- setdiff(colnames(met), "sample_id")
all_res <- list()
for (mc in met_cols) {
  y <- as.numeric(met[[mc]])
  ok <- is.finite(y)
  if (sum(ok) < 10) next
  if (mode == "gwas") {
    r <- mlm_scan(y[ok], G[ok, , drop = FALSE], Z[ok, , drop = FALSE], K[ok, ok, drop = FALSE])
    df <- data.frame(metabolite = mc, snp_id = snp_ids, beta = r[, 1], se = r[, 2], t = r[, 3], p = r[, 4])
  } else {
    r <- glm_scan(y[ok], G[ok, , drop = FALSE], Z[ok, , drop = FALSE])
    df <- data.frame(metabolite = mc, snp_id = snp_ids, beta = r[, 1], se = r[, 2], p = r[, 3])
  }
  df <- df[is.finite(df$p), , drop = FALSE]
  if (!nrow(df)) next
  df$padj <- p.adjust(df$p, method = "BH")
  df$chr <- info$chr[match(df$snp_id, info$snp_id)]
  df$pos <- info$pos[match(df$snp_id, info$snp_id)]
  all_res[[mc]] <- df
}
if (!length(all_res)) stop("No valid metabolite results")
res <- do.call(rbind, all_res)
res <- res[order(res$p), ]
write.csv(res, file.path(outdir, if (mode == "gwas") "mgwas_results.csv" else "mqtl_results.csv"), row.names = FALSE)
sig <- res[res$padj < 0.05, , drop = FALSE]
write.csv(sig, file.path(outdir, "significant.csv"), row.names = FALSE)
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("mode=", mode),
  paste0("samples=", nrow(G)), paste0("snps=", ncol(G)), paste0("metabolites=", length(all_res)),
  "method=gwas: eigen-transformed MLM (fixed ratio) + PCA; qtl: marker regression, BH"
), file.path(outdir, "run_info.txt"))

cat(sprintf("[mgwas] mode=%s samples=%d snps=%d metabolites=%d; significant(FDR<0.05)=%d\n",
            mode, nrow(G), ncol(G), length(all_res), nrow(sig)))
print(head(res[, c("metabolite", "snp_id", "beta", "p", "padj")], 5))
cat("[mgwas] outputs ->", outdir, "\n")
