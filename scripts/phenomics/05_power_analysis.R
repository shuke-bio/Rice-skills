###############################################################################
# 05_power_analysis.R -- Phenotyping / GWAS / GS power & sample-size calc
#
# Title:    Rice-skills experimental power analysis (step 05, phenomics)
# Description:
#   Dependency-free (base stats only) power and sample-size tools for the
#   three questions that come up before a rice phenotyping / GWAS / GS
#   experiment:
#
#   1) Two-group phenotype comparison (e.g. salt stress vs control):
#      power of the two-sample t-test given n per group and Cohen's d,
#      and the n required for a target power (closed form + grid refine).
#
#   2) GWAS single-SNP detection power:
#      given sample size n, MAF, and the SNP's PVE (variance explained),
#      the score/Wald test follows ~ chi2_1 with non-centrality
#          ncp = n * PVE / (1 - PVE)                       (LDAK/GCTA style)
#      power = P(chi2_1'(ncp) > threshold) at a Bonferroni threshold
#      alpha/Nmarkers; also the reverse: PVE needed for 80% power.
#
#   3) Genomic-prediction (GS) training-population size:
#      Daetwyler-type reliability  rho^2 = N h2 / (N h2 + Me)  with effective
#      chromosome segments Me (or r2-based Ne scaling); returns accuracy rho
#      for N, and N needed for a target accuracy.
#
# Usage:
#   Rscript scripts/phenomics/05_power_analysis.R                # demo/self-test
#   Rscript scripts/phenomics/05_power_analysis.R \
#       --ttest --d 1.0 --n-per-group 12 --alpha 0.05 \
#       --gwas --n 200 --maf 0.1 --pve 0.05 --n-markers 40000 \
#       --gs --h2 0.5 --me 5000 --target-acc 0.6 \
#       --outdir results/phenomics/05_power
#
#   (flags combine; each section prints a table and writes CSV)
#
# Outputs (written to --outdir):
#   power_ttest.csv / power_gwas.csv / gs_accuracy_curve.csv
#   power_report.txt    human-readable summary
#
# Notes:
#   - All formulas are the standard closed forms; validate big designs with
#     a simulation when in doubt (see kdense statistical-power skill:
#     third_party/kdense-scientific-skills/skills/statistical-power/).
#   - GWAS power assumes a single-causal-variant score test with population
#     structure controlled; hidden relatedness INFLATES the apparent ncp,
#     so treat results as optimistic when K<0 is neglected.
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1) return(args[i + 1])
  default
}
has_flag <- function(flag) !is.na(match(flag, args))

outdir <- get_arg("--outdir", if (has_flag("--demo")) "results/phenomics/05_power_demo" else "results/phenomics/05_power")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

report <- character(0)
say <- function(...) {
  line <- paste0(...)
  report <<- c(report, line)
  cat(line, "\n", sep = "")
}

# ------------------------------------------------- 1) two-group t-test power
power_ttest <- function(d, n, alpha = 0.05) {
  ncp <- d * sqrt(n / 2)
  tcrit <- qt(1 - alpha / 2, df = 2 * n - 2)
  pt(-tcrit, df = 2 * n - 2, ncp = ncp) + 1 - pt(tcrit, df = 2 * n - 2, ncp = ncp)
}
n_for_power <- function(d, target = 0.8, alpha = 0.05, n_max = 1e4) {
  n <- 2
  while (n < n_max && power_ttest(d, n, alpha) < target) n <- n + 1
  n
}

# --------------------------------------------------- 2) GWAS detection power
gwas_power <- function(n, maf, pve, n_markers = 40000, alpha_genome = 5e-8,
                       threshold = NULL) {
  thr <- if (is.null(threshold)) qchisq(1 - alpha_genome, df = 1) else threshold
  ncp <- n * pve / (1 - pve)
  pchisq(thr, df = 1, ncp = ncp, lower.tail = FALSE)
}
pve_for_power <- function(n, maf, target = 0.8, n_markers = 40000,
                          alpha_genome = 5e-8, tol = 1e-5) {
  thr <- qchisq(1 - alpha_genome, df = 1)
  lo <- 1e-6; hi <- 0.99
  for (i in 1:100) {
    mid <- (lo + hi) / 2
    if (gwas_power(n, maf, mid, n_markers, alpha_genome, thr) < target) lo <- mid else hi <- mid
    if (hi - lo < tol) break
  }
  (lo + hi) / 2
}

# ------------------------------------------------- 3) GS training-set sizing
gs_accuracy <- function(n, h2, me) sqrt(n * h2 / (n * h2 + me))
gs_n_for <- function(target, h2, me) {
  # rho^2 = N h2 / (N h2 + Me)  =>  N = target^2 * Me / (h2 * (1 - target^2))
  target^2 * me / (h2 * (1 - target^2))
}

# --------------------------------------------------------------- assemble
do_ttest <- has_flag("--ttest") || has_flag("--demo")
do_gwas  <- has_flag("--gwas")  || has_flag("--demo")
do_gs    <- has_flag("--gs")    || has_flag("--demo")

if (do_ttest) {
  alpha <- as.numeric(get_arg("--alpha", "0.05"))
  ds <- as.numeric(strsplit(get_arg("--d-list", "0.5,0.8,1.0,1.3,1.5"), ",")[[1]])
  ns <- as.numeric(strsplit(get_arg("--n-list", "6,8,10,12,15,20,30"), ",")[[1]])
  grid <- expand.grid(d = ds, n = ns)
  grid$power <- mapply(power_ttest, grid$d, grid$n, alpha)
  grid$alpha <- alpha
  grid <- grid[order(grid$d, grid$n), ]
  write.csv(grid, file.path(outdir, "power_ttest.csv"), row.names = FALSE)
  if (has_flag("--demo")) {
    # anchors cross-checked against R's built-in power.t.test()
    stopifnot(abs(power_ttest(0.8, 12, 0.05) - 0.4657) < 0.005)  # power.t.test(12,0.8)=0.4657
    stopifnot(abs(n_for_power(0.8, 0.8, 0.05) - 26) <= 1)        # power.t.test n=25.52
    stopifnot(n_for_power(1.2, 0.8, 0.05) <= 12)                 # big d -> small n
    stopifnot(n_for_power(0.5, 0.8, 0.05) >= 60)                 # small d -> big n (63.8)
  }
  say("[06][t-test] power (per group) at alpha = ", alpha,
      " ; n for 80% power: d=0.5 -> ", n_for_power(0.5, 0.8, alpha),
      ", d=0.8 -> ", n_for_power(0.8, 0.8, alpha),
      ", d=1.2 -> ", n_for_power(1.2, 0.8, alpha))
}

if (do_gwas) {
  n <- as.numeric(get_arg("--n", "200"))
  maf <- as.numeric(get_arg("--maf", "0.1"))
  nmk <- as.numeric(get_arg("--n-markers", "40000"))
  ag <- as.numeric(get_arg("--alpha-genome", "5e-8"))
  pves <- as.numeric(strsplit(get_arg("--pve-list", "0.005,0.01,0.02,0.05,0.1,0.2"), ",")[[1]])
  tab <- data.frame(n = n, maf = maf, pve = pves, n_markers = nmk)
  tab$power <- mapply(gwas_power, n, maf, pves, nmk, ag)
  tab$pve_for_80 <- pve_for_power(n, maf, 0.8, nmk, ag)
  write.csv(tab, file.path(outdir, "power_gwas.csv"), row.names = FALSE)
  if (has_flag("--demo")) {
    # anchors computed with the same closed form (cross-checked interactively)
    stopifnot(gwas_power(5000, 0.3, 0.05, 40000) > 0.999)       # big GWAS -> high power
    stopifnot(gwas_power(200, 0.1, 0.05, 40000) < 0.05)         # n=200 sees PVE=5% at 1.4% only
    stopifnot(abs(gwas_power(200, 0.1, 0.20, 40000) - 0.947) < 0.01)
    stopifnot(abs(pve_for_power(1000, 0.25, 0.8, 40000) - 0.0381) < 0.005)
  }
  say("[06][GWAS] n = ", n, ", MAF = ", maf, ", ", nmk, " markers (Bonferroni 5e-8):")
  say("[06][GWAS]   PVE 5% -> power ", sprintf("%.3f", gwas_power(n, maf, 0.05, nmk, ag)),
      " ; PVE 20% -> ", sprintf("%.3f", gwas_power(n, maf, 0.20, nmk, ag)))
  say("[06][GWAS]   PVE needed for 80% power: ", sprintf("%.3f", pve_for_power(n, maf, 0.8, nmk, ag)))
}

if (do_gs) {
  h2 <- as.numeric(get_arg("--h2", "0.5"))
  me <- as.numeric(get_arg("--me", "5000"))
  target <- as.numeric(get_arg("--target-acc", "0.6"))
  ns <- as.numeric(strsplit(get_arg("--n-list-gs", "50,100,200,300,400,500,1000"), ",")[[1]])
  tab <- data.frame(N = ns, h2 = h2, Me = me,
                    accuracy = sapply(ns, gs_accuracy, h2 = h2, me = me))
  write.csv(tab, file.path(outdir, "gs_accuracy_curve.csv"), row.names = FALSE)
  n_need <- ceiling(gs_n_for(target, h2, me))
  if (has_flag("--demo")) {
    # closed-form anchors: acc(500, h2=.5, Me=5000) = sqrt(250/5250) = 0.2182;
    # N for accuracy .6 = .36*5000/(.5*.64) = 5625
    stopifnot(abs(gs_accuracy(500, 0.5, 5000) - 0.2182) < 0.001)
    stopifnot(abs(gs_n_for(0.6, 0.5, 5000) - 5625) < 1)
    stopifnot(gs_accuracy(n_need, h2, me) >= target - 1e-6)
  }
  say("[06][GS] h2 = ", h2, ", Me = ", me, " : accuracy at N=200 -> ",
      sprintf("%.3f", gs_accuracy(200, h2, me)),
      " ; N=500 -> ", sprintf("%.3f", gs_accuracy(500, h2, me)))
  say("[06][GS]   training N needed for accuracy ", target, ": ", n_need)
}

if (length(report) == 0) {
  say("[06] nothing to do: pass --ttest / --gwas / --gs (or --demo)")
}

writeLines(report, file.path(outdir, "power_report.txt"))
cat("[06] wrote ", file.path(outdir, "power_report.txt"), "\n", sep = "")
if (has_flag("--demo")) cat("[06] demo self-test PASSED\n")
