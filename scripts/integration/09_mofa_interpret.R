#!/usr/bin/env Rscript
###############################################################################
# 09_mofa_interpret.R -- post-hoc interpretation of MOFA factors
#
# 输入（不需要 MOFA2 对象，直接读导出表）：
#   --factors  CSV: sample_id, Factor1, Factor2, ...
#   --loadings CSV: feature_id, Factor1, Factor2, ...
#   --metadata 可选 CSV: sample_id, trait1, trait2, ...
# 输出：
#   factor_variance.csv（载荷平方和占比）
#   factor_top_features.csv（每因子 top 载荷）
#   factor_trait_correlation.csv（因子-性状 Pearson + BH）
#
# 用法：
#   Rscript scripts/integration/09_mofa_interpret.R --factors factors.csv \
#       --loadings loadings.csv --metadata meta.csv
#   Rscript scripts/integration/09_mofa_interpret.R --demo --outdir results/integration/09_mofa_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/integration/09_mofa")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function(n = 100, f = 4, g = 200) {
  set.seed(42)
  samples <- sprintf("S%03d", 1:n)
  F <- matrix(rnorm(n * f), n, f); colnames(F) <- paste0("Factor", 1:f); rownames(F) <- samples
  L <- matrix(rnorm(g * f, 0, 0.5), g, f); colnames(L) <- colnames(F); rownames(L) <- sprintf("G%03d", 1:g)
  meta <- data.frame(sample_id = samples,
                     yield = 0.8 * F[,1] + rnorm(n, 0, 0.5),
                     stress = -0.7 * F[,2] + rnorm(n, 0, 0.5))
  list(F = data.frame(sample_id = samples, F, check.names = FALSE),
       L = data.frame(feature_id = rownames(L), L, check.names = FALSE), meta = meta)
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$F, file.path(outdir, "demo_factors.csv"), row.names = FALSE)
  write.csv(d$L, file.path(outdir, "demo_loadings.csv"), row.names = FALSE)
  write.csv(d$meta, file.path(outdir, "demo_metadata.csv"), row.names = FALSE)
  f_path <- file.path(outdir, "demo_factors.csv"); l_path <- file.path(outdir, "demo_loadings.csv")
  m_path <- file.path(outdir, "demo_metadata.csv")
} else {
  f_path <- opt("factors"); l_path <- opt("loadings"); m_path <- opt("metadata")
}
if (is.null(f_path) || is.null(l_path)) stop("Need --factors and --loadings, or --demo")

F <- read.csv(f_path, check.names = FALSE); L <- read.csv(l_path, check.names = FALSE)
if (!"sample_id" %in% colnames(F) || !"feature_id" %in% colnames(L)) stop("bad columns")
facs <- intersect(setdiff(colnames(F), "sample_id"), setdiff(colnames(L), "feature_id"))
if (!length(facs)) stop("No shared factors")
var <- data.frame(factor = facs,
                  variance = sapply(facs, function(f) sum(L[[f]]^2, na.rm = TRUE)))
var$variance_pct <- round(100 * var$variance / sum(var$variance), 2)
write.csv(var, file.path(outdir, "factor_variance.csv"), row.names = FALSE)

top_rows <- list()
for (f in facs) {
  d <- L[order(-abs(L[[f]])), c("feature_id", f), drop = FALSE]
  names(d)[2] <- "loading"
  d$factor <- f
  top_rows[[f]] <- head(d, 15)
}
write.csv(do.call(rbind, top_rows), file.path(outdir, "factor_top_features.csv"), row.names = FALSE)

if (!is.null(m_path) && file.exists(m_path)) {
  M <- read.csv(m_path, check.names = FALSE)
  d <- merge(F, M, by = "sample_id")
  traits <- setdiff(colnames(M), "sample_id")
  rows <- list()
  for (f in facs) for (tr in traits) {
    r <- cor(d[[f]], d[[tr]], use = "complete.obs")
    n <- sum(is.finite(d[[f]]) & is.finite(d[[tr]]))
    p <- if (n > 2) 2 * pt(-abs(r * sqrt((n - 2) / (1 - r^2))), n - 2) else NA
    rows[[length(rows) + 1]] <- data.frame(factor = f, trait = tr, r = round(r, 3), p = p)
  }
  ct <- do.call(rbind, rows)
  ct$padj <- p.adjust(ct$p, method = "BH")
  write.csv(ct[order(ct$p), ], file.path(outdir, "factor_trait_correlation.csv"), row.names = FALSE)
}
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("factors=", length(facs)),
  "note=variance from loading squared sums; correlations are not causal"
), file.path(outdir, "run_info.txt"))
cat("[mofa] factors:", paste(facs, collapse = ", "), "\n")
cat("[mofa] variance:", paste(var$factor, var$variance_pct, sep = "=", collapse = ", "), "\n")
cat("[mofa] outputs ->", outdir, "\n")
