#!/usr/bin/env Rscript
###############################################################################
# 07_mr.R -- two-sample Mendelian randomization (IVW / MR-Egger / weighted median)
#
# 输入：
#   --exposure CSV: snp_id, beta, se, effect_allele, other_allele, eaf(可选)
#   --outcome  CSV: 同上
# 输出：
#   mr_results.csv（各方法 beta/se/p）、mr_snps.csv（协调后的 SNP 表）
#
# 说明：
#   - 等位基因方向自动协调；回文 SNP（A/T、C/G）且 EAF 接近 0.5 时剔除；
#   - 这是基础 MR，不含 Steiger filtering / 多效性稳健方法（如 MR-PRESSO）。
#
# 用法：
#   Rscript scripts/qtl/07_mr.R --exposure exp.csv --outcome out.csv --outdir results/qtl/07_mr
#   Rscript scripts/qtl/07_mr.R --demo --outdir results/qtl/07_mr_demo
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/qtl/07_mr")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function() {
  set.seed(7)
  n <- 22
  snp <- sprintf("rs%04d", 1:n)
  ea <- sample(c("A", "C", "G", "T"), n, TRUE)
  oa <- vapply(ea, function(a) sample(setdiff(c("A", "C", "G", "T"), a), 1), "")
  eaf <- runif(n, 0.15, 0.85)
  bx <- rnorm(n, 0, 0.28); sx <- runif(n, 0.03, 0.06)
  # 真实因果效应 0.35；加入少量多效性噪声
  by <- 0.35 * bx + rnorm(n, 0, 0.015); sy <- runif(n, 0.04, 0.07)
  exp <- data.frame(snp_id = snp, beta = bx, se = sx,
                    effect_allele = ea, other_allele = oa, eaf = eaf)
  out <- data.frame(snp_id = snp, beta = by, se = sy,
                    effect_allele = ea, other_allele = oa, eaf = eaf)
  # 制造 3 个等位基因互换，测试协调
  idx <- 1:3
  out$effect_allele[idx] <- oa[idx]
  out$other_allele[idx] <- ea[idx]
  out$beta[idx] <- -out$beta[idx]
  out$eaf[idx] <- 1 - out$eaf[idx]
  list(exposure = exp, outcome = out, truth = 0.35)
}

if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$exposure, file.path(outdir, "demo_exposure.csv"), row.names = FALSE)
  write.csv(d$outcome, file.path(outdir, "demo_outcome.csv"), row.names = FALSE)
  exp_path <- file.path(outdir, "demo_exposure.csv")
  out_path <- file.path(outdir, "demo_outcome.csv")
} else {
  exp_path <- opt("exposure"); out_path <- opt("outcome")
}
if (is.null(exp_path) || is.null(out_path)) stop("Need --exposure and --outcome, or --demo")

ex <- read.csv(exp_path, check.names = FALSE)
ou <- read.csv(out_path, check.names = FALSE)
need <- function(df, cols, name) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss)) stop(sprintf("%s missing columns: %s", name, paste(miss, collapse = ", ")))
}
need(ex, c("snp_id", "beta", "se", "effect_allele", "other_allele"), "exposure")
need(ou, c("snp_id", "beta", "se", "effect_allele", "other_allele"), "outcome")

m <- merge(ex, ou, by = "snp_id", suffixes = c("_ex", "_out"))
if (!nrow(m)) stop("No shared SNPs")
for (cc in c("effect_allele_ex", "other_allele_ex", "effect_allele_out", "other_allele_out")) {
  m[[cc]] <- toupper(m[[cc]])
}
keep <- logical(nrow(m))
for (i in seq_len(nrow(m))) {
  a1 <- m$effect_allele_ex[i]; a2 <- m$other_allele_ex[i]
  b1 <- m$effect_allele_out[i]; b2 <- m$other_allele_out[i]
  if (a1 == b1 && a2 == b2) {
    keep[i] <- TRUE
  } else if (a1 == b2 && a2 == b1) {
    m$beta_out[i] <- -m$beta_out[i]
    if ("eaf_out" %in% colnames(m)) m$eaf_out[i] <- 1 - m$eaf_out[i]
    keep[i] <- TRUE
  } else keep[i] <- FALSE
  # 回文 SNP + EAF 接近 0.5 时剔除
  pal <- paste0(sort(c(a1, a2)), collapse = "")
  if (keep[i] && pal %in% c("AT", "CG") && "eaf_ex" %in% colnames(m)) {
    if (!is.na(m$eaf_ex[i]) && m$eaf_ex[i] > 0.42 && m$eaf_ex[i] < 0.58) keep[i] <- FALSE
  }
}
m <- m[keep, , drop = FALSE]
if (nrow(m) < 3) stop("Fewer than 3 usable SNPs after harmonization")

m$ratio <- m$beta_out / m$beta_ex
m$se_ratio <- m$se_out / abs(m$beta_ex)
m$F_stat <- (m$beta_ex / m$se_ex)^2
mean_F <- mean(m$F_stat, na.rm = TRUE)
m <- m[is.finite(m$ratio) & is.finite(m$se_ratio) & m$se_ratio > 0, , drop = FALSE]
if (nrow(m) < 3) stop("Fewer than 3 valid SNPs after ratio calculation")

w <- 1 / m$se_ratio^2
b_ivw <- sum(w * m$ratio) / sum(w)
se_ivw <- sqrt(1 / sum(w))
Q <- sum(w * (m$ratio - b_ivw)^2)
df <- nrow(m) - 1
pQ <- if (df > 0) pchisq(Q, df, lower.tail = FALSE) else NA
se_ivw_r <- if (df > 0) max(se_ivw, sqrt(Q / df / sum(w))) else se_ivw

fit_e <- lm(beta_out ~ beta_ex, data = m, weights = 1 / m$se_out^2)
se_slope <- summary(fit_e)$coefficients[2, 2]
p_slope <- summary(fit_e)$coefficients[2, 4]
se_int <- summary(fit_e)$coefficients[1, 2]
p_int <- summary(fit_e)$coefficients[1, 4]

wmed <- function(x, ww) {
  o <- order(x); x <- x[o]; ww <- ww[o]; ww <- ww / sum(ww)
  x[which(cumsum(ww) >= 0.5)[1]]
}
b_wm <- wmed(m$ratio, w)
set.seed(123)
bs <- replicate(1000, {
  idx <- sample(nrow(m), replace = TRUE)
  wmed(m$ratio[idx], w[idx])
})
se_wm <- sd(bs, na.rm = TRUE)

res <- data.frame(
  method = c("IVW_fixed", "IVW_random", "MR_Egger", "Weighted_median"),
  n_snp = nrow(m),
  beta = c(b_ivw, b_ivw, coef(fit_e)[2], b_wm),
  se = c(se_ivw, se_ivw_r, se_slope, se_wm),
  p = c(2 * pnorm(-abs(b_ivw / se_ivw)), 2 * pnorm(-abs(b_ivw / se_ivw_r)),
        p_slope, 2 * pnorm(-abs(b_wm / se_wm))),
  extra = c(sprintf("Q=%.2f; pQ=%.3g", Q, pQ),
            sprintf("Q=%.2f; pQ=%.3g", Q, pQ),
            sprintf("intercept=%.4f; p_int=%.3g", coef(fit_e)[1], p_int),
            "bootstrap se (1000)")
)

write.csv(res, file.path(outdir, "mr_results.csv"), row.names = FALSE)
if (isTRUE(opt("presso", FALSE))) {
  if (requireNamespace("MRPRESSO", quietly = TRUE)) {
    pr <- tryCatch(MRPRESSO::mr_presso(
      BetaOutcome = m$beta_out, BetaExposure = m$beta_ex,
      SdOutcome = m$se_out, SdExposure = m$se_ex,
      NbDistribution = 1000, SignifThreshold = 0.05,
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE), error = function(e) NULL)
    if (!is.null(pr)) {
      write.csv(pr$`Main MR results`, file.path(outdir, "mr_presso.csv"), row.names = FALSE)
      cat("[mr] MR-PRESSO results written\n")
    }
  } else {
    writeLines(c("# MR-PRESSO plan", "install.packages('MRPRESSO')",
                 "# MRPRESSO::mr_presso(BetaOutcome=..., BetaExposure=..., SdOutcome=..., SdExposure=...)"),
               file.path(outdir, "mr_presso_plan.txt"))
    cat("[mr] MRPRESSO not installed; plan written\n")
  }
}
write.csv(m[, c("snp_id", "beta_ex", "se_ex", "beta_out", "se_out", "ratio", "se_ratio")],
          file.path(outdir, "mr_snps.csv"), row.names = FALSE)
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("exposure=", normalizePath(exp_path, mustWork = FALSE)),
  paste0("outcome=", normalizePath(out_path, mustWork = FALSE)),
  paste0("n_snp=", nrow(m)),
  sprintf("mean_F=%.1f", mean_F),
  "methods=IVW_fixed, IVW_random, MR_Egger, Weighted_median"
), file.path(outdir, "run_info.txt"))

cat(sprintf("[mr] n_snp=%d; IVW beta=%.3f (se %.3f, p=%.3g); Q p=%.3g\n",
            nrow(m), b_ivw, se_ivw, 2 * pnorm(-abs(b_ivw / se_ivw)), pQ))
print(res)
cat("[mr] outputs ->", outdir, "\n")
