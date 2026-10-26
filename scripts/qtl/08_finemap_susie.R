#!/usr/bin/env Rscript
###############################################################################
# 08_finemap_susie.R -- fine-mapping with SuSiE (optional) or ABF fallback
#
# 输入：
#   --z     CSV: snp_id, z（或 beta,se；有 n 时可换算 z）
#   --ld    CSV: 第一列 snp_id，其余为 LD 矩阵（R）
#   --n     样本量（SuSiE 需要；ABF 可省略）
#   --method susie | abf（默认 abf；susie 需安装 susieR）
#
# 输出：
#   finemap_results.csv（每个 SNP 的 ABF/后验/是否在 95% CS）
#   credible_sets.csv（credible set 汇总）
#
# 说明：
#   ABF 为近似 Wakefield 因子，默认 prior r=0.01，不考虑 LD；
#   SuSiE 使用 LD 矩阵，结果更可靠，建议优先。
#
# 用法：
#   Rscript scripts/qtl/08_finemap_susie.R --z z.csv --ld ld.csv --n 500 --method susie
#   Rscript scripts/qtl/08_finemap_susie.R --demo --outdir results/qtl/08_finemap_demo
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/qtl/08_finemap")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function(n_snp = 120, causal = 60) {
  set.seed(42)
  snp <- sprintf("snp%03d", 1:n_snp)
  # AR(1) LD
  rho <- 0.85
  R <- outer(1:n_snp, 1:n_snp, function(i, j) rho^abs(i - j))
  z <- rnorm(n_snp, 0, 0.8)
  z[causal] <- 6.5
  # local correlation spill-over
  z[(causal-2):(causal+2)] <- z[(causal-2):(causal+2)] + c(1.2, 2.0, 0, 2.0, 1.2)
  zdf <- data.frame(snp_id = snp, z = z)
  colnames(R) <- snp; rownames(R) <- snp
  ldf <- data.frame(snp_id = snp, R, check.names = FALSE)
  list(z = zdf, ld = ldf, causal = snp[causal])
}

if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$z, file.path(outdir, "demo_z.csv"), row.names = FALSE)
  write.csv(d$ld, file.path(outdir, "demo_ld.csv"), row.names = FALSE)
  z_path <- file.path(outdir, "demo_z.csv"); ld_path <- file.path(outdir, "demo_ld.csv")
  truth <- d$causal
} else {
  z_path <- opt("z"); ld_path <- opt("ld"); truth <- NULL
}
if (is.null(z_path)) stop("Need --z or --demo")
method <- tolower(opt("method", "abf"))

z <- read.csv(z_path, check.names = FALSE)
if (!all(c("snp_id", "z") %in% colnames(z))) {
  if (all(c("snp_id", "beta", "se") %in% colnames(z))) {
    z$z <- z$beta / z$se
  } else stop("z file needs snp_id,z or snp_id,beta,se")
}
z <- z[is.finite(z$z), , drop = FALSE]

R <- NULL
if (!is.null(ld_path) && file.exists(ld_path)) {
  ld <- read.csv(ld_path, check.names = FALSE)
  if (!"snp_id" %in% colnames(ld)) stop("LD file needs snp_id column")
  common <- intersect(z$snp_id, ld$snp_id)
  if (length(common) < 5) stop("Fewer than 5 shared SNPs between z and LD")
  z <- z[match(common, z$snp_id), , drop = FALSE]
  R <- as.matrix(ld[match(common, ld$snp_id), common, drop = FALSE])
  storage.mode(R) <- "double"
  diag(R) <- 1
}

res <- NULL
if (method == "susie") {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("method=susie requires susieR; install.packages('susieR') or use --method abf")
  }
  if (is.null(R)) stop("SuSiE requires --ld")
  n <- as.numeric(opt("n", 500))
  fit <- susieR::susie_rss(z = z$z, R = R, n = n, L = 10)
  cs <- susieR::susie_get_cs(fit, Xcorr = R)
  res <- data.frame(snp_id = z$snp_id, z = z$z)
  res$posterior <- NA_real_; res$in_cs <- FALSE
  for (k in seq_along(cs$cs)) {
    idx <- cs$cs[[k]]
    res$posterior[idx] <- if (!is.null(cs$purity[[k]])) cs$purity[[k]] else NA_real_
    res$in_cs[idx] <- TRUE
  }
  cs_tab <- data.frame(cs = seq_along(cs$cs),
                       n_snp = vapply(cs$cs, length, 0L),
                       snps = vapply(cs$cs, function(i) paste(z$snp_id[i], collapse = ";"), ""))
} else {
  # approximate Wakefield ABF, prior r = 0.1
  r <- 0.01
  logbf <- 0.5 * (log(1 - r) + z$z^2 * r / (1 - r))
  post <- exp(logbf - max(logbf)); post <- post / sum(post)
  ord <- order(post, decreasing = TRUE)
  cum <- cumsum(post[ord])
  in_cs <- rep(FALSE, nrow(z)); in_cs[ord[cum <= 0.95]] <- TRUE
  if (!any(in_cs)) in_cs[ord[1]] <- TRUE
  res <- data.frame(snp_id = z$snp_id, z = z$z, abf = exp(logbf), posterior = post, in_cs = in_cs)
  cs_tab <- data.frame(cs = 1, n_snp = sum(in_cs),
                       snps = paste(z$snp_id[in_cs], collapse = ";"))
}

res <- res[order(-res$posterior), ]
write.csv(res, file.path(outdir, "finemap_results.csv"), row.names = FALSE)
write.csv(cs_tab, file.path(outdir, "credible_sets.csv"), row.names = FALSE)
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("method=", method), paste0("snps=", nrow(res)),
  "abf=approx Wakefield, prior r=0.01, ignores LD; susie uses LD"
), file.path(outdir, "run_info.txt"))

cat(sprintf("[finemap] method=%s snps=%d; credible SNPs=%d\n", method, nrow(res), sum(res$in_cs)))
print(head(res, 5))
if (!is.null(truth)) cat(sprintf("[finemap] demo causal %s in credible set: %s\n", truth, truth %in% res$snp_id[res$in_cs]))
cat("[finemap] outputs ->", outdir, "\n")
