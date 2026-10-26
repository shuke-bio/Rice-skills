#!/usr/bin/env Rscript
###############################################################################
# 11_batch_timecourse.R -- batch correction (ComBat) + time-course differential model
#
# 输入：
#   --counts   CSV: sample_id, gene1, gene2, ...（counts）
#   --metadata CSV: sample_id, group, time, batch（batch 可选）
# 方法：
#   CPM/log2 -> 可选 ComBat（sva）-> 每基因 lm(y ~ time*group) vs lm(y ~ group)
#   -> F 检验 + BH。时间项/交互项显著即 time-course 差异。
#
# 用法：
#   Rscript scripts/transcriptomics/11_batch_timecourse.R --counts counts.csv \
#       --metadata meta.csv --outdir results/transcriptomics/11_timecourse
#   Rscript scripts/transcriptomics/11_batch_timecourse.R --demo --outdir results/transcriptomics/11_timecourse_demo
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  i <- match(paste0("--", key), args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) return(TRUE)
  args[i + 1]
}
outdir <- opt("outdir", "results/transcriptomics/11_timecourse")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_demo <- function(n_genes = 200, times = c(0, 6, 12, 24), reps = 3, groups = c("ctrl", "salt")) {
  set.seed(42)
  samples <- character(); md <- data.frame()
  for (g in groups) for (t in times) for (r in seq_len(reps)) {
    sid <- sprintf("%s_T%s_R%d", g, t, r)
    samples <- c(samples, sid)
    md <- rbind(md, data.frame(sample_id = sid, group = g, time = t, batch = ifelse(r == 1, "B1", "B2")))
  }
  mu <- matrix(rnbinom(n_genes * length(samples), mu = 80, size = 5), n_genes)
  # 20 time-dependent genes in salt
  for (i in 1:20) {
    idx <- md$group == "salt"
    mu[i, idx] <- mu[i, idx] * (1 + 0.20 * md$time[idx])
  }
  # batch effect
  mu[, md$batch == "B2"] <- mu[, md$batch == "B2"] * 1.3
  counts <- matrix(rpois(length(mu), lambda = mu), nrow = n_genes)
  colnames(counts) <- samples; rownames(counts) <- sprintf("G%03d", 1:n_genes)
  list(counts = data.frame(sample_id = samples, t(counts), check.names = FALSE), meta = md)
}
if (isTRUE(opt("demo", FALSE))) {
  d <- make_demo()
  write.csv(d$counts, file.path(outdir, "demo_counts.csv"), row.names = FALSE)
  write.csv(d$meta, file.path(outdir, "demo_metadata.csv"), row.names = FALSE)
  counts_path <- file.path(outdir, "demo_counts.csv"); meta_path <- file.path(outdir, "demo_metadata.csv")
} else {
  counts_path <- opt("counts"); meta_path <- opt("metadata")
}
if (is.null(counts_path) || is.null(meta_path)) stop("Need --counts and --metadata, or --demo")

cnt <- read.csv(counts_path, check.names = FALSE)
md <- read.csv(meta_path, check.names = FALSE)
if (!"sample_id" %in% colnames(cnt) || !"sample_id" %in% colnames(md)) stop("need sample_id")
md <- md[match(cnt$sample_id, md$sample_id), , drop = FALSE]
if (!all(c("group", "time") %in% colnames(md))) stop("metadata needs group and time")
genes <- setdiff(colnames(cnt), "sample_id")
X <- as.matrix(cnt[, genes, drop = FALSE]); storage.mode(X) <- "double"
X[is.na(X)] <- 0
lib <- colSums(X); lib[lib == 0] <- 1
Y <- t(log2(t(t(X) / lib) * 1e6 + 1))

batch_col <- opt("batch", "batch")
if (batch_col %in% colnames(md) && length(unique(md[[batch_col]])) > 1) {
  if (requireNamespace("sva", quietly = TRUE)) {
    mod <- model.matrix(~ group + time, data = md)
    Y <- sva::ComBat(dat = Y, batch = md[[batch_col]], mod = mod, par.prior = TRUE)
    method <- "ComBat"
  } else {
    method <- "no correction (sva not installed)"
  }
} else method <- "no batch column"

rows <- list()
for (g in genes) {
  y <- as.numeric(Y[g, ])
  fit_full <- lm(y ~ group * time, data = md)
  fit_red <- lm(y ~ group, data = md)
  a <- anova(fit_red, fit_full)
  p <- a$`Pr(>F)`[2]
  if (!is.finite(p)) next
  cf <- coef(fit_full)
  bg <- if (any(grepl("^group", names(cf)))) unname(cf[grep("^group", names(cf))[1]]) else NA_real_
  bt <- if ("time" %in% names(cf)) unname(cf[["time"]]) else NA_real_
  bi <- if (any(grepl(":", names(cf)))) unname(cf[grep(":", names(cf))[1]]) else NA_real_
  rows[[length(rows) + 1]] <- data.frame(
    gene_id = g, p_time = p,
    beta_group = bg, beta_time = bt, beta_interaction = bi)
}
res <- do.call(rbind, rows)
res$padj <- p.adjust(res$p_time, method = "BH")
res <- res[order(res$p_time), ]
write.csv(res, file.path(outdir, "timecourse_results.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = cnt$sample_id, t(Y), check.names = FALSE),
          file.path(outdir, "batch_corrected_log2cpm.csv"), row.names = FALSE)
writeLines(c(
  paste0("date=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  paste0("counts=", normalizePath(counts_path, mustWork = FALSE)),
  paste0("metadata=", normalizePath(meta_path, mustWork = FALSE)),
  paste0("batch=", method),
  "model=lm(y ~ group*time) vs lm(y ~ group); F test; BH"
), file.path(outdir, "run_info.txt"))
cat(sprintf("[timecourse] genes=%d batch=%s; FDR<0.05 time/Interaction: %d\n",
            nrow(res), method, sum(res$padj < 0.05)))
print(head(res, 5))
cat("[timecourse] outputs ->", outdir, "\n")
