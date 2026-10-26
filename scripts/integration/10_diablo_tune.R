#!/usr/bin/env Rscript
###############################################################################
# 10_diablo_tune.R -- DIABLO (mixOmics) tuning / plan
#
# 输入：--x 基因矩阵, --y 代谢矩阵, --groups 分组
# 有 mixOmics 且 --run -> tune.block.splsda；否则写 diablo_plan.R
#
# 用法：
#   Rscript scripts/integration/10_diablo_tune.R --x genes.csv --y mets.csv --groups groups.csv --run
#   Rscript scripts/integration/10_diablo_tune.R --demo --outdir results/integration/10_diablo_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/integration/10_diablo")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
run <- isTRUE(opt("run", FALSE))

make_demo <- function(n = 80, g = 100, m = 60) {
  set.seed(42)
  grp <- rep(c("A", "B"), each = n / 2)
  X <- matrix(rnorm(n * g), n, g); Y <- matrix(rnorm(n * m), n, m)
  X[grp == "B", 1:10] <- X[grp == "B", 1:10] + 1
  Y[grp == "B", 1:8] <- Y[grp == "B", 1:8] + 1
  colnames(X) <- sprintf("G%03d", 1:g); colnames(Y) <- sprintf("M%03d", 1:m)
  list(X = data.frame(sample_id = sprintf("S%03d", 1:n), X, check.names = FALSE),
       Y = data.frame(sample_id = sprintf("S%03d", 1:n), Y, check.names = FALSE),
       groups = data.frame(sample_id = sprintf("S%03d", 1:n), group = grp))
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$X, file.path(outdir, "demo_genes.csv"), row.names = FALSE)
  write.csv(d$Y, file.path(outdir, "demo_metabolites.csv"), row.names = FALSE)
  write.csv(d$groups, file.path(outdir, "demo_groups.csv"), row.names = FALSE)
  x_path <- file.path(outdir, "demo_genes.csv"); y_path <- file.path(outdir, "demo_metabolites.csv")
  g_path <- file.path(outdir, "demo_groups.csv")
} else {
  x_path <- opt("x"); y_path <- opt("y"); g_path <- opt("groups")
}
if (is.null(x_path) || is.null(y_path) || is.null(g_path)) stop("Need --x --y --groups, or --demo")

plan <- c(
  "library(mixOmics)",
  sprintf("X <- read.csv('%s', check.names=FALSE)", x_path),
  sprintf("Y <- read.csv('%s', check.names=FALSE)", y_path),
  sprintf("G <- read.csv('%s', check.names=FALSE)", g_path),
  "rownames(X) <- X$sample_id; X$sample_id <- NULL",
  "rownames(Y) <- Y$sample_id; Y$sample_id <- NULL",
  "G <- G[match(rownames(X), G$sample_id), ]",
  "data <- list(genes=X, mets=Y)",
  "Ynum <- factor(G$group)",
  "design <- matrix(0.1, 2, 2); diag(design) <- 0",
  "tune <- tune.block.splsda(X=data, Y=Ynum, ncomp=3, design=design,",
  "  test.keepX=list(genes=c(10,20,50), mets=c(10,20,30)), validation='Mfold', folds=5, nrepeat=3)",
  "saveRDS(tune, 'diablo_tune.rds')",
  "write.csv(tune$error.rate, 'diablo_error_rate.csv')"
)
writeLines(plan, file.path(outdir, "diablo_plan.R"))
if (run && requireNamespace("mixOmics", quietly = TRUE)) {
  source(file.path(outdir, "diablo_plan.R"), local = TRUE)
  cat("[diablo] tuning finished\n")
} else {
  msg <- if (!run) "plan only (use --run to execute)" else "mixOmics not installed; plan written"
  cat("[diablo]", msg, "->", file.path(outdir, "diablo_plan.R"), "\n")
}
writeLines(c(paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
             paste0("run=", run), paste0("mixOmics=", requireNamespace("mixOmics", quietly = TRUE))),
           file.path(outdir, "run_info.txt"))
