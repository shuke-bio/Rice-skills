#!/usr/bin/env Rscript
###############################################################################
# 10_qtl_epistasis.R -- pairwise QTL epistasis scan
#
# 输入：
#   --genotype CSV: sample_id, snp1, ...（0/1/2）
#   --phenotype CSV: sample_id, trait
#   --pairs CSV: snp1, snp2（可选；不给则用 --region 内前 N 对）
#   --region chr:start-end（可选）
# 方法：
#   y ~ g1 + g2 + g1:g2 + PCs；交互项 p -> BH
#
# 用法：
#   Rscript scripts/qtl/10_qtl_epistasis.R --genotype geno.csv --phenotype phe.csv --pairs pairs.csv
#   Rscript scripts/qtl/10_qtl_epistasis.R --demo --outdir results/qtl/10_epistasis_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/qtl/10_epistasis")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function(n = 300, m = 40) {
  set.seed(42)
  G <- matrix(rbinom(n * m, 2, 0.3), n, m)
  colnames(G) <- sprintf("snp%02d", 1:m); rownames(G) <- sprintf("S%03d", 1:n)
  y <- rnorm(n)
  # one epistatic pair
  y <- y + 0.6 * G[, 3] * G[, 7]
  geno <- data.frame(sample_id = rownames(G), G, check.names = FALSE)
  phe <- data.frame(sample_id = rownames(G), trait = y)
  pairs <- data.frame(snp1 = c("snp03", "snp01", "snp10"), snp2 = c("snp07", "snp02", "snp12"))
  list(geno = geno, phe = phe, pairs = pairs)
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$geno, file.path(outdir, "demo_genotype.csv"), row.names = FALSE)
  write.csv(d$phe, file.path(outdir, "demo_phenotype.csv"), row.names = FALSE)
  write.csv(d$pairs, file.path(outdir, "demo_pairs.csv"), row.names = FALSE)
  geno_path <- file.path(outdir, "demo_genotype.csv"); phe_path <- file.path(outdir, "demo_phenotype.csv")
  pairs_path <- file.path(outdir, "demo_pairs.csv")
} else {
  geno_path <- opt("genotype"); phe_path <- opt("phenotype"); pairs_path <- opt("pairs")
}
if (is.null(geno_path) || is.null(phe_path) || is.null(pairs_path)) stop("Need --genotype --phenotype --pairs, or --demo")

g <- read.csv(geno_path, check.names = FALSE); p <- read.csv(phe_path, check.names = FALSE)
pr <- read.csv(pairs_path, check.names = FALSE)
if (!all(c("sample_id") %in% colnames(g)) || !all(c("sample_id","trait") %in% colnames(p))) stop("bad input")
if (!all(c("snp1","snp2") %in% colnames(pr))) stop("pairs needs snp1,snp2")
d <- merge(g, p, by = "sample_id")
snps <- setdiff(colnames(g), "sample_id")
G <- as.matrix(d[, snps]); storage.mode(G) <- "double"
for (j in seq_len(ncol(G))) { v <- G[, j]; v[is.na(v)] <- mean(v, na.rm = TRUE); G[, j] <- v }
y <- d$trait
Z <- matrix(1, nrow(G), 1)
npc <- as.integer(opt("npc", 3))
if (npc > 0) Z <- cbind(Z, prcomp(G, center = TRUE)$x[, seq_len(min(npc, ncol(G))), drop = FALSE])

rows <- list()
for (i in seq_len(nrow(pr))) {
  a <- as.character(pr$snp1[i]); b <- as.character(pr$snp2[i])
  if (!a %in% snps || !b %in% snps) next
  g1 <- G[, a]; g2 <- G[, b]
  fit <- lm(y ~ g1 * g2 + Z)
  co <- summary(fit)$coefficients
  term <- grep(":", rownames(co), value = TRUE)[1]
  if (is.na(term)) next
  rows[[length(rows) + 1]] <- data.frame(
    snp1 = a, snp2 = b, beta_interaction = co[term, "Estimate"],
    se = co[term, "Std. Error"], p = co[term, "Pr(>|t|)"])
}
if (!length(rows)) stop("No valid pairs")
res <- do.call(rbind, rows)
res$padj <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
write.csv(res, file.path(outdir, "epistasis_results.csv"), row.names = FALSE)
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("genotype=", normalizePath(geno_path, mustWork = FALSE)),
  paste0("phenotype=", normalizePath(phe_path, mustWork = FALSE)),
  paste0("pairs=", nrow(res)),
  "model=lm(y ~ g1*g2 + PCs); BH"
), file.path(outdir, "run_info.txt"))
cat(sprintf("[epistasis] pairs tested=%d; FDR<0.05=%d\n", nrow(res), sum(res$padj < 0.05)))
print(head(res, 5))
cat("[epistasis] outputs ->", outdir, "\n")
