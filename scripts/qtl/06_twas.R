#!/usr/bin/env Rscript
###############################################################################
# 06_twas.R -- transcriptome-wide association (weighted-sum, FUSION-like)
#
# 输入：
#   --eqtl   CSV: gene_id, snp_id, weight（可选 beta/se）
#   --gwas   CSV: snp_id, beta, se（或 z）
# 方法：
#   基因 z = sum(w_i * z_i) / sqrt(sum(w_i^2))；BH 校正；可选置换。
#   这是 FUSION 的 weighted-sum 近似，不是完整 FUSION（不含 LD 参考/无权重训练）。
#
# 用法：
#   Rscript scripts/qtl/06_twas.R --eqtl eqtl_weights.csv --gwas gwas.csv \
#       --outdir results/qtl/06_twas
#   Rscript scripts/qtl/06_twas.R --demo --outdir results/qtl/06_twas_demo
###############################################################################

suppressWarnings(suppressMessages({
  # base R only
}))

args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}

outdir <- opt("outdir", "results/qtl/06_twas")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function() {
  set.seed(42)
  genes <- sprintf("Gene%02d", 1:40)
  rows <- list(); gwas <- list()
  causal <- "Gene07"
  for (g in genes) {
    n <- 15
    snps <- sprintf("%s_snp%02d", g, 1:n)
    w <- rnorm(n, 0, 0.35)
    # 只有 causal 基因的遗传成分与性状相关
    signal <- if (g == causal) 0.55 else 0
    b <- signal * w + rnorm(n, 0, 0.10)
    se <- rep(0.08, n)
    rows[[g]] <- data.frame(gene_id = g, snp_id = snps, weight = w)
    gwas[[g]] <- data.frame(snp_id = snps, beta = b, se = se)
  }
  eqtl <- do.call(rbind, rows); gw <- do.call(rbind, gwas)
  list(eqtl = eqtl, gwas = gw, causal = causal)
}

ld_path <- opt("ld")
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$eqtl, file.path(outdir, "demo_eqtl_weights.csv"), row.names = FALSE)
  write.csv(d$gwas, file.path(outdir, "demo_gwas.csv"), row.names = FALSE)
  eqtl_path <- file.path(outdir, "demo_eqtl_weights.csv")
  gwas_path <- file.path(outdir, "demo_gwas.csv")
} else {
  eqtl_path <- opt("eqtl"); gwas_path <- opt("gwas")
}

if (is.null(eqtl_path) || is.null(gwas_path)) {
  stop("Need --eqtl and --gwas, or --demo")
}

eq <- read.csv(eqtl_path, check.names = FALSE)
gw <- read.csv(gwas_path, check.names = FALSE)
need <- function(df, cols, name) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss)) stop(sprintf("%s missing columns: %s", name, paste(miss, collapse = ", ")))
}
need(eq, c("gene_id", "snp_id", "weight"), "eqtl")
need(gw, c("snp_id", "beta", "se"), "gwas")

gw$z <- gw$beta / gw$se
m <- merge(eq, gw[, c("snp_id", "z")], by = "snp_id")
LD <- NULL; ld_ids <- NULL
if (!is.null(ld_path) && file.exists(ld_path)) {
  ldf <- read.csv(ld_path, check.names = FALSE)
  if (!"snp_id" %in% colnames(ldf)) stop("LD file needs snp_id column")
  ld_ids <- ldf$snp_id
  LD <- as.matrix(ldf[, setdiff(colnames(ldf), "snp_id"), drop = FALSE])
  storage.mode(LD) <- "double"
  diag(LD) <- 1
}
if (!nrow(m)) stop("No shared SNPs between eqtl and gwas")

res <- do.call(rbind, lapply(split(m, m$gene_id), function(d) {
  d <- d[is.finite(d$weight) & is.finite(d$z), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  if (!is.null(LD)) {
    idx <- match(d$snp_id, ld_ids)
    ok <- !is.na(idx)
    d <- d[ok, , drop = FALSE]; idx <- idx[ok]
    if (!nrow(d)) return(NULL)
    R <- LD[idx, idx, drop = FALSE]
    R[is.na(R)] <- 0; diag(R) <- 1
    denom <- sqrt(as.numeric(t(d$weight) %*% R %*% d$weight))
  } else {
    denom <- sqrt(sum(d$weight^2))
  }
  if (!is.finite(denom) || denom == 0) return(NULL)
  z <- sum(d$weight * d$z) / denom
  p <- 2 * pnorm(-abs(z))
  data.frame(gene_id = d$gene_id[1], n_snp = nrow(d), z = z, p = p,
             top_snp = d$snp_id[which.max(abs(d$weight))],
             top_weight = d$weight[which.max(abs(d$weight))])
}))
if (is.null(res) || !nrow(res)) stop("TWAS produced no results")
res$padj <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]

write.csv(res, file.path(outdir, "twas_results.csv"), row.names = FALSE)
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("eqtl=", normalizePath(eqtl_path, mustWork = FALSE)),
  paste0("gwas=", normalizePath(gwas_path, mustWork = FALSE)),
  paste0("ld=", if (is.null(LD)) "identity" else normalizePath(ld_path, mustWork = FALSE)),
  paste0("genes=", nrow(res)),
  "method=weighted-sum (FUSION-like); p from normal approximation; BH adjustment"
), file.path(outdir, "run_info.txt"))

cat(sprintf("[twas] genes=%d; min p=%.3g; FDR<0.05: %d\n",
            nrow(res), min(res$p), sum(res$padj < 0.05)))
print(head(res[, c("gene_id", "n_snp", "z", "p", "padj")], 5))
cat("[twas] outputs ->", outdir, "\n")
