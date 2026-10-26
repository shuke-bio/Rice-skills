#!/usr/bin/env Rscript
###############################################################################
# 10_ssgsea.R -- ssGSEA / GSVA pathway activity per sample
#
# 输入：
#   --expression CSV: gene_id, sample1, sample2, ...（log/TPM）
#   --gmt GMT 文件: set_name<TAB>gene1<TAB>gene2...
# 方法：
#   有 GSVA 包 -> GSVA::gsva(method="ssgsea")；否则退化为 set mean z-score
#
# 用法：
#   Rscript scripts/transcriptomics/10_ssgsea.R --expression expr.csv --gmt sets.gmt
#   Rscript scripts/transcriptomics/10_ssgsea.R --demo --outdir results/transcriptomics/10_ssgsea_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/transcriptomics/10_ssgsea")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function() {
  set.seed(42)
  genes <- sprintf("G%03d", 1:200)
  samples <- sprintf("S%02d", 1:30)
  X <- matrix(rnorm(length(genes) * length(samples), 6, 1), length(genes), length(samples))
  rownames(X) <- genes; colnames(X) <- samples
  sets <- list(
    SetA = genes[1:20], SetB = genes[21:40], SetC = genes[41:60],
    SetD = genes[61:80], SetE = genes[81:100], SetF = genes[101:120]
  )
  # SetA higher in second half
  X[sets$SetA, 16:30] <- X[sets$SetA, 16:30] + 1.2
  list(X = X, sets = sets)
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  expr <- data.frame(gene_id = rownames(d$X), d$X, check.names = FALSE)
  write.csv(expr, file.path(outdir, "demo_expression.csv"), row.names = FALSE)
  gmt <- file.path(outdir, "demo_sets.gmt")
  con <- file(gmt, "w")
  for (nm in names(d$sets)) writeLines(paste(c(nm, d$sets[[nm]]), collapse = "\t"), con)
  close(con)
  expr_path <- file.path(outdir, "demo_expression.csv"); gmt_path <- gmt
} else {
  expr_path <- opt("expression"); gmt_path <- opt("gmt")
}
if (is.null(expr_path) || is.null(gmt_path)) stop("Need --expression and --gmt, or --demo")

expr <- read.csv(expr_path, check.names = FALSE)
if (!"gene_id" %in% colnames(expr)) stop("expression needs gene_id")
mat <- as.matrix(expr[, setdiff(colnames(expr), "gene_id"), drop = FALSE])
rownames(mat) <- expr$gene_id
gsets <- list()
for (line in readLines(gmt_path)) {
  f <- strsplit(line, "\t")[[1]]
  if (length(f) >= 2) gsets[[f[1]]] <- f[-1]
}
gsets <- lapply(gsets, function(g) intersect(g, rownames(mat)))
gsets <- gsets[lengths(gsets) >= 2]
if (!length(gsets)) stop("No gene set overlaps expression matrix")

if (requireNamespace("GSVA", quietly = TRUE)) {
  scores <- GSVA::gsva(mat, gsets, method = "ssgsea", verbose = FALSE)
  method <- "GSVA ssgsea"
} else {
  z <- t(scale(t(mat)))
  scores <- do.call(rbind, lapply(gsets, function(g) colMeans(z[g, , drop = FALSE], na.rm = TRUE)))
  method <- "mean z-score fallback (GSVA not installed)"
}
out <- data.frame(gene_set = rownames(scores), scores, check.names = FALSE)
write.csv(out, file.path(outdir, "ssgsea_scores.csv"), row.names = FALSE)
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("expression=", normalizePath(expr_path, mustWork = FALSE)),
  paste0("gmt=", normalizePath(gmt_path, mustWork = FALSE)),
  paste0("method=", method), paste0("sets=", length(gsets))
), file.path(outdir, "run_info.txt"))
cat(sprintf("[ssgsea] sets=%d samples=%d method=%s\n", length(gsets), ncol(mat), method))
print(out[1:min(5, nrow(out)), 1:min(6, ncol(out))])
cat("[ssgsea] outputs ->", outdir, "\n")
