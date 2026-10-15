###############################################################################
# 08_2d_association.R -- Pairwise trait-conditioned association (2D-GWAS)
#
# Title:    Rice-skills 2D association mapping (v3.2, method distillation)
# Description:
#   Pairwise trait-conditioned GWAS ("2D association", method distilled from
#   the PATOWAS idea -- see docs/references/patowas-2d-association.md; the
#   GPL-3.0 PATOWAS program itself is NOT bundled):
#
#   For every ordered trait pair (A -> B):
#     1. split the population by trait A into extreme groups
#        (top --tail fraction vs bottom --tail, e.g. 10%);
#     2. inside each extreme group, test B ~ SNP for every SNP
#        (fast closed-form least-squares slope t-test, base R only);
#     3. call a SNP "A-conditioned for B" when B associates in ONE extreme
#        group only at --p-cutoff while the other group does not:
#           direction = "+A"  (B associates only in the HIGH-A group)
#                     = "-A"  (only in the LOW-A group)
#        This finds genetic effects on B that are masked/created by the
#        background of A -- the classic 2D pattern that a standard
#        whole-population GWAS on B misses.
#
# Usage:
#   Rscript scripts/gwas/08_2d_association.R \
#       --geno demo_data/genomics/genotype_matrix.csv \
#       --pheno demo_data/genomics/phenotypes.csv \
#       --traits traitA,traitB --tail 0.1 --p-cutoff 0.01 \
#       --outdir results/gwas/08_2d
#   Rscript scripts/gwas/08_2d_association.R --demo
#
# Inputs:
#   --geno   samples x SNPs genotype matrix CSV (0/1/2 dosages; first col =
#            sample_id, header = SNP ids)
#   --pheno  samples x traits CSV (first col = sample_id)
#   --traits comma list; default = all pairs
#
# Outputs:
#   2d_results.csv     trait_A, trait_B, snp, n_high, n_low, p_high, p_low,
#                      beta_high, beta_low, delta_logp, direction, significant
#   2d_summary.txt     per-pair hit counts + demo self-test result
#   2d_scatter.png     delta_logp (high vs low) per SNP for the best pair
#
# Notes:
#   - Within-group t-test per SNP uses the closed-form slope test
#     (equivalent to lm(B ~ dosage)); no R packages required.
#   - This is a screen: conditional hits should be followed up with the
#     standard GWAS pipeline (module 7) + colocalization (02/03) on the
#     conditioning trait.
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1) return(args[i + 1])
  default
}
has_flag <- function(flag) !is.na(match(flag, args))
outdir <- get_arg("--outdir", if (has_flag("--demo")) "results/gwas/08_2d_demo" else "results/gwas/08_2d")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tail_frac <- as.numeric(get_arg("--tail", "0.1"))
p_cut <- as.numeric(get_arg("--p-cutoff", "0.01"))

# ------------------------------------------------------- fast slope p-value
slope_p <- function(x, y) {
  # closed-form OLS slope test for y ~ x; returns c(p, slope)
  n <- length(x)
  if (n < 5) return(c(1, 0))
  sx <- mean(x); my <- mean(y)
  dx <- x - sx; dy <- y - my
  sxx <- sum(dx^2)
  if (sxx == 0) return(c(1, 0))
  slope <- sum(dx * dy) / sxx
  resid_ss <- sum(dy^2) - slope^2 * sxx
  df <- n - 2
  if (df <= 0 || resid_ss <= 0) return(c(1, slope))
  se <- sqrt(resid_ss / df / sxx)
  if (se == 0) return(c(1, slope))
  tval <- slope / se
  p <- 2 * pt(-abs(tval), df = df)
  c(p, slope)
}

# ------------------------------------------------------------------- demo
demo_data <- function() {
  set.seed(7)
  n <- 120; nsnp <- 80
  snps <- matrix(sample(0:2, n * nsnp, TRUE, prob = c(0.5, 0.35, 0.15)),
                 nrow = n)
  colnames(snps) <- sprintf("SNP%03d", 1:nsnp)
  # latent trait A driven by SNP010
  A <- 0.9 * snps[, "SNP010"] + rnorm(n, 0, 1)
  # trait B: SNP030 has an effect ONLY in the high-A background (epistasis)
  hi <- A > quantile(A, 1 - 0.1)
  lo <- A < quantile(A, 0.1)
  B <- 0.6 * snps[, "SNP055"] + rnorm(n, 0, 1) + 1.5 * snps[, "SNP030"] * hi
  # trait C: pure noise -> no 2D hits expected
  C <- rnorm(n, 0, 1)
  ph <- data.frame(sample_id = sprintf("S%03d", 1:n), A = A, B = B, C = C)
  list(geno = data.frame(sample_id = sprintf("S%03d", 1:n), snps, check.names = FALSE),
       pheno = ph, truth = c(cond_snp = "SNP030", base_snp = "SNP055", A_snp = "SNP010"))
}

report <- character(0)
say <- function(...) { line <- paste0(...); report <<- c(report, line); cat(line, "\n", sep = "") }

if (has_flag("--demo")) {
  d <- demo_data()
  geno <- d$geno; pheno <- d$pheno
  demo_truth <- d$truth
} else {
  gpath <- get_arg("--geno"); ppath <- get_arg("--pheno")
  if (is.null(gpath) || is.null(ppath)) stop("need --geno and --pheno (or --demo)")
  geno <- read.csv(gpath, stringsAsFactors = FALSE, check.names = FALSE)
  pheno <- read.csv(ppath, stringsAsFactors = FALSE, check.names = FALSE)
}

# align samples
common <- intersect(geno[[1]], pheno[[1]])
geno <- geno[match(common, geno[[1]]), ]
pheno <- pheno[match(common, pheno[[1]]), ]
snp_mat <- as.matrix(geno[, -1, drop = FALSE])
snp_ids <- colnames(snp_mat)
trait_names <- if (!is.null(get_arg("--traits"))) {
  strsplit(get_arg("--traits"), ",")[[1]]
} else colnames(pheno)[-1]
trait_names <- intersect(trait_names, colnames(pheno)[-1])
stopifnot(length(trait_names) >= 2)
ph <- as.matrix(pheno[, trait_names, drop = FALSE])
rownames(ph) <- common

pairs <- expand.grid(A = trait_names, B = trait_names,
                     stringsAsFactors = FALSE)
pairs <- pairs[pairs$A != pairs$B, ]

res_list <- list()
for (k in seq_len(nrow(pairs))) {
  A <- ph[, pairs$A[k]]; B <- ph[, pairs$B[k]]
  o <- order(A)
  n_tail <- max(5, round(length(A) * tail_frac))
  hi_idx <- o[(length(A) - n_tail + 1):length(A)]
  lo_idx <- o[1:n_tail]
  for (j in seq_along(snp_ids)) {
    x_all <- snp_mat[, j]
    ph_high <- slope_p(x_all[hi_idx], B[hi_idx])
    ph_low  <- slope_p(x_all[lo_idx], B[lo_idx])
    sig_hi <- ph_high[1] < p_cut
    sig_lo <- ph_low[1] < p_cut
    direction <- if (sig_hi && !sig_lo) "+A" else if (sig_lo && !sig_hi) "-A" else NA
    res_list[[length(res_list) + 1]] <- data.frame(
      trait_A = pairs$A[k], trait_B = pairs$B[k], snp = snp_ids[j],
      n_high = length(hi_idx), n_low = length(lo_idx),
      p_high = ph_high[1], beta_high = ph_high[2],
      p_low = ph_low[1], beta_low = ph_low[2],
      delta_logp = -log10(ph_high[1]) - (-log10(ph_low[1])),
      direction = ifelse(is.na(direction), "", direction),
      significant = sig_hi || sig_lo)
  }
}
res <- do.call(rbind, res_list)
res <- res[order(-abs(res$delta_logp)), ]
write.csv(res, file.path(outdir, "2d_results.csv"), row.names = FALSE)

cond_hits <- res[res$direction != "" & res$significant, ]
say("[08] pairs scanned: ", nrow(pairs), " ; SNP x pair tests: ", nrow(res))
say("[08] A-conditioned hits (one extreme group only, p < ", p_cut, "): ", nrow(cond_hits))
if (nrow(cond_hits) > 0) {
  top <- cond_hits[order(-abs(cond_hits$delta_logp)), ][1:min(10, nrow(cond_hits)), ]
  for (i in seq_len(nrow(top))) {
    say("[08]   ", top$trait_A[i], " -> ", top$trait_B[i], " @ ", top$snp[i],
        " dir=", top$direction[i], " p_high=", sprintf("%.2g", top$p_high[i]),
        " p_low=", sprintf("%.2g", top$p_low[i]))
  }
}

# scatter for the first pair (high vs low -log10p)
if (nrow(res) > 0) {
  first <- res[res$trait_A == pairs$A[1] & res$trait_B == pairs$B[1], ]
  png(file.path(outdir, "2d_scatter.png"), width = 1600, height = 1400, res = 220)
  op <- par(mar = c(4.2, 4.2, 3, 1))
  plot(-log10(first$p_low), -log10(first$p_high),
       pch = 16, cex = 0.8, col = "#4C72B0",
       xlab = paste0("-log10 p in LOW-", pairs$A[1], " group"),
       ylab = paste0("-log10 p in HIGH-", pairs$A[1], " group"),
       main = paste("2D association:", pairs$A[1], "->", pairs$B[1]))
  lim <- max(c(-log10(first$p_low), -log10(first$p_high)))
  abline(0, 1, col = "grey50", lty = 2)
  abline(h = -log10(p_cut), col = "#C44E52", lty = 3)
  abline(v = -log10(p_cut), col = "#C44E52", lty = 3)
  if (nrow(cond_hits) > 0) {
    ch <- cond_hits[cond_hits$trait_A == pairs$A[1] & cond_hits$trait_B == pairs$B[1], ]
    if (nrow(ch) > 0) {
      pts <- first[match(ch$snp, first$snp), ]
      points(-log10(pts$p_low), -log10(pts$p_high), pch = 17, col = "#C44E52", cex = 1.1)
    }
  }
  par(op); dev.off()
  say("[08] wrote ", file.path(outdir, "2d_scatter.png"))
}

if (has_flag("--demo")) {
  t <- demo_truth
  hit <- cond_hits[cond_hits$trait_A == "A" & cond_hits$trait_B == "B" &
                     cond_hits$snp == t[["cond_snp"]], ]
  stopifnot(nrow(hit) == 1, hit$direction == "+A")   # planted epistatic SNP
  stopifnot(!any(cond_hits$snp == t[["base_snp"]] & cond_hits$trait_B == "B" &
                   cond_hits$direction == "-A"))     # base SNP not in low group only
  say("[08] demo self-test PASSED (SNP030 recovered as +A-conditioned for B)")
}

writeLines(report, file.path(outdir, "2d_summary.txt"))
cat("[08] wrote ", file.path(outdir, "2d_summary.txt"), "\n", sep = "")
