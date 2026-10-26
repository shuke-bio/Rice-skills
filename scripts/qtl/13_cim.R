#!/usr/bin/env Rscript
###############################################################################
# 13_cim.R -- composite interval mapping (CIM) for a single chromosome
#
# 输入：
#   --genotype CSV: sample_id, snp1, ...（0/1/2）
#   --phenotype CSV: sample_id, trait
#   --snp-info CSV: snp_id, chr, pos
#   --chrom 染色体；--cofactor-window bp；--max-cofactors N
# 方法：
#   单标记扫描选 cofactor -> 每个位置 y ~ g + cofactors（排除窗口内）-> LOD
#
# 用法：
#   Rscript scripts/qtl/13_cim.R --genotype geno.csv --phenotype phe.csv --snp-info info.csv --chrom chr1
#   Rscript scripts/qtl/13_cim.R --demo --outdir results/qtl/13_cim_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/qtl/13_cim"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
make_demo <- function(n = 200, m = 80) {
  set.seed(42)
  G <- matrix(rbinom(n * m, 2, 0.3), n, m)
  colnames(G) <- sprintf("snp%03d", 1:m); rownames(G) <- sprintf("S%03d", 1:n)
  y <- 0.7 * G[, 40] + rnorm(n, 0, 0.8)
  info <- data.frame(snp_id = colnames(G), chr = "chr1", pos = seq(1e6, by = 1e5, length.out = m))
  list(geno = data.frame(sample_id = rownames(G), G, check.names = FALSE),
       phe = data.frame(sample_id = rownames(G), trait = y), info = info)
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$geno, file.path(outdir, "demo_genotype.csv"), row.names = FALSE)
  write.csv(d$phe, file.path(outdir, "demo_phenotype.csv"), row.names = FALSE)
  write.csv(d$info, file.path(outdir, "demo_snp_info.csv"), row.names = FALSE)
  geno_path <- file.path(outdir, "demo_genotype.csv"); phe_path <- file.path(outdir, "demo_phenotype.csv")
  info_path <- file.path(outdir, "demo_snp_info.csv"); chrom <- "chr1"
} else {
  geno_path <- opt("genotype"); phe_path <- opt("phenotype"); info_path <- opt("snp-info")
  chrom <- opt("chrom", "chr1")
}
if (is.null(geno_path) || is.null(phe_path) || is.null(info_path)) stop("Need --genotype --phenotype --snp-info, or --demo")
g <- read.csv(geno_path, check.names = FALSE); p <- read.csv(phe_path, check.names = FALSE); info <- read.csv(info_path, check.names = FALSE)
d <- merge(g, p, by = "sample_id")
snps <- setdiff(colnames(g), "sample_id")
X <- as.matrix(d[, snps]); storage.mode(X) <- "double"
for (j in seq_len(ncol(X))) { v <- X[, j]; v[is.na(v)] <- mean(v, na.rm = TRUE); X[, j] <- v }
info <- info[match(snps, info$snp_id), , drop = FALSE]
keep <- info$chr == chrom & is.finite(info$pos)
X <- X[, keep, drop = FALSE]; snps <- snps[keep]; pos <- info$pos[keep]
y <- d$trait
pc <- prcomp(X, center = TRUE)$x[, 1:min(3, ncol(X)), drop = FALSE]
# single-marker p
p1 <- sapply(seq_len(ncol(X)), function(j) {
  fit <- lm(y ~ X[, j] + pc); summary(fit)$coefficients[2, 4]
})
cof_n <- as.integer(opt("max-cofactors", 5)); cof_win <- as.numeric(opt("cofactor-window", 5e6))
ord <- order(p1); cof <- integer(0)
for (i in ord) {
  if (length(cof) >= cof_n) break
  if (all(abs(pos[i] - pos[cof]) > cof_win)) cof <- c(cof, i)
}
lod <- rep(NA_real_, ncol(X)); beta <- rep(NA_real_, ncol(X))
for (j in seq_len(ncol(X))) {
  use <- setdiff(cof, which(abs(pos - pos[j]) <= cof_win))
  Z <- cbind(pc, X[, use, drop = FALSE])
  fit <- lm(y ~ X[, j] + Z)
  co <- summary(fit)$coefficients
  beta[j] <- co[2, 1]; lod[j] <- -log10(co[2, 4])
}
res <- data.frame(snp_id = snps, chr = chrom, pos = pos, beta = beta, LOD = lod)
res <- res[order(-res$LOD), ]
write.csv(res, file.path(outdir, "cim_results.csv"), row.names = FALSE)
writeLines(c(paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
             paste0("chrom=", chrom), paste0("cofactors=", paste(snps[cof], collapse = ",")),
             "model=lm(y ~ SNP + cofactors + PCs)"), file.path(outdir, "run_info.txt"))
cat(sprintf("[cim] SNPs=%d cofactors=%d; top LOD=%.2f at %s\n", nrow(res), length(cof), res$LOD[1], res$snp_id[1]))
print(head(res, 5))
cat("[cim] outputs ->", outdir, "\n")
