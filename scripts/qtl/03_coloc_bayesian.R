###############################################################################
# 03_coloc_bayesian.R -- Bayesian colocalization (coloc.abf self-implemented)
#
# Title:    Rice-skills Bayesian colocalization step (v3.2 upgrade of 02)
# Description:
#   Full Bayesian colocalization (Giambartolomei et al. 2014, PLoS Genet)
#   for one GWAS signal x one molecular QTL, implemented in base R so the
#   pipeline does not depend on the (not-installed) `coloc` package:
#
#     per-SNP approximate Bayes factor (quant trait, sdY=1 parametrisation)
#         ABF_i = sqrt(V_i/(V_i + w^2)) * exp(z_i^2 w^2 / (2 (V_i + w^2)))
#     with V_i = se_i^2, z_i = beta_i/se_i, prior scale w = 0.04 (coloc
#     default when effects are given in standardised units).
#
#     five hypotheses over the region (single causal variant each):
#       H0 no association; H1 GWAS only; H2 QTL only;
#       H3 both, different causal variants; H4 both, SHARED causal variant
#     PP_H4 ∝ p12 * sum_i ABF1_i ABF2_i      (shared)
#     PP_H3 ∝ p1 p2 (sum_i ABF1_i)(sum_j ABF2_j)
#     plus H0/H1/H2 terms; normalised to PP0..PP4 (coloc defaults
#     p1 = p2 = 1e-4, p12 = 1e-5).
#
#   Decision rule (same thresholds as the coloc package guidance):
#       PP4 >= 0.8                      strong colocalization
#       0.5 <= PP4 < 0.8                suggestive
#       PP4 < 0.5 or PP3 > PP4          evidence for distinct signals
#   Also reports the posterior of each SNP being THE shared causal variant
#   (ratio PP4_i / PP4), which the 02_colocalization.R screen lacks.
#
# Usage:
#   Rscript scripts/qtl/03_coloc_bayesian.R \
#       --gwas results/qtl/02_coloc/gwas_region.csv \
#       --qtl  results/qtl/02_coloc/qtl_region.csv \
#       --outdir results/qtl/03_coloc_bayes [--w 0.04] [--p12 1e-5]
#
#   Region CSVs: columns snp_id, beta, se  (one row per SNP in the region;
#   same SNP set and order NOT required -- joined by snp_id).
#   Self-test (no inputs): Rscript scripts/qtl/03_coloc_bayesian.R --demo
#
# Outputs:
#   coloc_pp.csv      PP0..PP4, decision, n_snps, w, priors
#   snp_posteriors.csv per-SNP shared-causal posterior + lead annotation
#   coloc_report.txt  human-readable summary
#
# Notes:
#   - sdY parametrisation: betas/se are assumed to be in units of the trait
#     SD (standardised phenotypes give exactly this). For unstandardised
#     traits scale both traits' se by their sdY, or pass --w scaled by sdY;
#     PP values are robust to moderate misspecification but CHECK with the
#     sensitivity guidance in third_party/kdense-scientific-skills/skills/
#     statistical-analysis/ when PP4 sits in 0.5-0.8.
#   - 02_colocalization.R remains the fast genome-wide screen: feed only the
#     regions it flags (or any locus of interest) to this script.
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1) return(args[i + 1])
  default
}
has_flag <- function(flag) !is.na(match(flag, args))
outdir <- get_arg("--outdir", if (has_flag("--demo")) "results/qtl/03_coloc_bayes_demo" else "results/qtl/03_coloc_bayes")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

w      <- as.numeric(get_arg("--w",   "0.04"))
p1     <- as.numeric(get_arg("--p1",  "1e-4"))
p2     <- as.numeric(get_arg("--p2",  "1e-4"))
p12    <- as.numeric(get_arg("--p12", "1e-5"))

# ------------------------------------------------------------- core maths
abf_vector <- function(beta, se, w) {
  V <- se^2
  tau2 <- w^2
  sqrt(V / (V + tau2)) * exp((beta^2) * tau2 / (2 * V * (V + tau2)))
}

coloc_region <- function(b1, se1, b2, se2, w, p1, p2, p12) {
  abf1 <- abf_vector(b1, se1, w)
  abf2 <- abf_vector(b2, se2, w)
  s1 <- sum(abf1); s2 <- sum(abf2); s12 <- sum(abf1 * abf2)
  h0 <- (1 - p1) * (1 - p2)
  h1 <- p1 * (1 - p2) * s1
  h2 <- p2 * (1 - p1) * s2
  h3 <- p1 * p2 * (s1 * s2 - s12)
  h4 <- p12 * s12
  tot <- h0 + h1 + h2 + h3 + h4
  pp <- c(h0, h1, h2, h3, h4) / tot
  snp_post <- (abf1 * abf2) / s12 * pp[5]   # shared-causal posterior
  list(pp = pp, snp_post = snp_post,
       ncp1 = sum(abf1 > 1), ncp2 = sum(abf2 > 1))
}

decision <- function(pp4, pp3) {
  if (pp4 >= 0.8) "colocalized (PP4 >= 0.8)"
  else if (pp4 >= 0.5 && pp4 > pp3) "suggestive (0.5 <= PP4 < 0.8)"
  else if (pp3 > pp4) "distinct signals likely (PP3 > PP4)"
  else "unclear (PP4 < 0.5)"
}

# ------------------------------------------------------------------- demo
demo_data <- function() {
  set.seed(42)
  n <- 200
  # realistic se scale (N large -> se ~ 0.01); the w=0.04 prior lives in
  # sdY units and needs this scale for the ABF to overcome H0's prior mass
  z1 <- rnorm(n, 0, 1); z2 <- rnorm(n, 0, 1)
  z1[17] <- 5.5; z2[17] <- 5.0                       # same causal variant
  gwas <- data.frame(snp_id = sprintf("SNP%03d", 1:n),
                     beta = z1 * 0.01, se = rep(0.01, n))
  qtl  <- data.frame(snp_id = gwas$snp_id, beta = z2 * 0.03, se = rep(0.03, n))
  # distinct-causal scenario: GWAS peak at 17, QTL peak at 71
  z3 <- rnorm(n, 0, 1); z4 <- rnorm(n, 0, 1)
  z3[17] <- 5.5; z4[71] <- 5.0
  qtl_d <- data.frame(snp_id = gwas$snp_id, beta = z4 * 0.03, se = rep(0.03, n))
  gwas_d <- data.frame(snp_id = gwas$snp_id, beta = z3 * 0.01, se = rep(0.01, n))
  list(gwas = gwas, qtl = qtl, gwas_d = gwas_d, qtl_d = qtl_d)
}

run_and_write <- function(gwas, qtl, label, outdir, tag) {
  m <- merge(gwas, qtl, by = "snp_id", suffixes = c(".g", ".q"))
  stopifnot(nrow(m) > 1)
  res <- coloc_region(m$beta.g, m$se.g, m$beta.q, m$se.q, w, p1, p2, p12)
  pp <- res$pp
  df <- data.frame(hypothesis = c("H0", "H1", "H2", "H3", "H4"),
                   posterior = pp, label = label)
  write.csv(df, file.path(outdir, paste0("coloc_pp_", tag, ".csv")), row.names = FALSE)
  snp_df <- data.frame(snp_id = m$snp_id,
                       pp_shared = res$snp_post)
  snp_df <- snp_df[order(-snp_df$pp_shared), ]
  write.csv(snp_df, file.path(outdir, paste0("snp_posteriors_", tag, ".csv")), row.names = FALSE)
  lead <- snp_df$snp_id[1]
  say("[03b][", label, "] PP0=", sprintf("%.3f", pp[1]),
      " PP1=", sprintf("%.3f", pp[2]), " PP2=", sprintf("%.3f", pp[3]),
      " PP3=", sprintf("%.3f", pp[4]), " PP4=", sprintf("%.3f", pp[5]),
      " -> ", decision(pp[5], pp[4]), " ; top shared-causal SNP: ", lead)
  invisible(list(pp = pp, lead = lead))
}

report <- character(0)
say <- function(...) { line <- paste0(...); report <<- c(report, line); cat(line, "\n", sep = "") }

if (has_flag("--demo")) {
  d <- demo_data()
  r1 <- run_and_write(d$gwas, d$qtl, "demo_shared_causal", outdir, "shared")
  r2 <- run_and_write(d$gwas_d, d$qtl_d, "demo_distinct_causal", outdir, "distinct")
  stopifnot(r1$pp[5] >= 0.8, r1$pp[5] > r1$pp[4])          # shared -> strong PP4
  stopifnot(r1$lead == "SNP017")                            # right SNP
  stopifnot(r2$pp[4] > r2$pp[5])                            # distinct -> PP3 dominates
  say("[03b] demo self-test PASSED (shared: PP4 >= 0.8 lead=SNP017; distinct: PP3 > PP4)")
} else {
  gpath <- get_arg("--gwas"); qpath <- get_arg("--qtl")
  if (is.null(gpath) || is.null(qpath)) stop("need --gwas and --qtl region CSVs (or --demo)")
  gwas <- read.csv(gpath, stringsAsFactors = FALSE)
  qtl  <- read.csv(qpath,  stringsAsFactors = FALSE)
  run_and_write(gwas, qtl, "region", outdir, "region")
}

writeLines(report, file.path(outdir, "coloc_report.txt"))
cat("[03b] wrote ", file.path(outdir, "coloc_report.txt"), "\n", sep = "")
