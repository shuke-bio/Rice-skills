###############################################################################
# 07_xpehh.R -- cross-population extended haplotype homozygosity (XP-EHH)
#
# Title:    Rice-skills selection scan (XP-EHH) step
# Description:
#   XP-EHH (Sabeti et al. 2007) detects near-fixation selection by comparing
#   integrated EHH decay between two populations:
#     XP-EHH(site) = ln( iEHH_A / iEHH_B )   over core SNPs,
#   where iEHH integrates the decay of haplotype homozygosity around each
#   core SNP (exact single-haplotype counting on diploid chromosomes with
#   random phase imputation), truncated when EHH < --stop-ehh.
#   Scores are standardised over sliding windows; |z| >= 2 windows are
#   reported as candidates (positive z = selection in population 1).
#   For XP-CLR (composite likelihood) use the XPCLR binary upstream —
#   see docs/genomics_usage.md; this script ships XP-EHH only.
#
# Usage:
#   Rscript scripts/popgen/07_xpehh.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --metadata demo_data/genomics/phenotype.csv \
#       --group-col group --group1 indica --group2 japonica \
#       --chr Chr05 --outdir results/popgen/07_xpehh
#
# Inputs:
#   --genotype   SNP x sample dosage matrix (0/1/2, imputed preferred)
#   --snp-info   snp_id, chr, pos
#   --metadata   sample metadata CSV (sample_id + population column)
#   --group-col  metadata column with the two population labels (default group)
#   --group1/--group2  the two populations (XP-EHH = ln(iEHH_1 / iEHH_2))
#   --chr        restrict to one chromosome (recommended; call per chromosome
#                for genome-wide scans)
#   --max-core   maximum core SNPs scored (default 2000, evenly spaced)
#   --max-dist   EHH decay limit in SNPs each side (default 50)
#   --win-kb     standardisation window (default 100 kb)
#   --maf        core-SNP MAF filter (default 0.05)
#   --stop-ehh   early-stop threshold for the decay integral (default 0.01)
#   --outdir     output directory (default results/popgen/07_xpehh)
# Outputs:
#   xpehh_sites.csv       per core SNP: iEHH_1, iEHH_2, raw log-ratio, window z
#   xpehh_windows.csv     window-level summary (mean/max |z|) + candidate flag
#   candidate_windows.csv candidates sorted by |z| (when any)
#   xpehh_distribution.* / xpehh_manhattan.*    figures
#   xpehh_summary.txt / run_info.txt
#
# Notes:
#   - Diploid dosages are split into two chromosomes with random phase;
#     population-level EHH expectations are insensitive to this (standard
#     practice for genotype-array input).
#   - EHH(d) computed by exact haplotype counting on the cumulative window
#     (core +/- d): EHH = sum_h (n_h / m)^2.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
source(file.path(dirname(script_file), "..", "utils", "genomics_utils.R"))
rice_set_seed(42)

parse_args <- function(args) {
  opts <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) stop("Unexpected argument '", a, "'. Use --key value or --flag.", call. = FALSE)
    key <- sub("^--", "", a)
    if (i == length(args) || grepl("^--", args[i + 1L])) { opts[[key]] <- TRUE; i <- i + 1L }
    else { opts[[key]] <- args[i + 1L]; i <- i + 2L }
  }
  opts
}
opt <- function(opts, key, default = NULL) {
  v <- opts[[key]]; if (is.null(v)) return(default)
  if (isTRUE(v)) return(if (is.logical(default) && length(default) == 1L) default else TRUE); v
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
genotype_path <- opt(opts, "genotype")
snp_info_path <- opt(opts, "snp-info")
metadata_path <- opt(opts, "metadata")
group_col <- opt(opts, "group-col", "group")
g1 <- opt(opts, "group1"); g2 <- opt(opts, "group2")
chr_arg <- opt(opts, "chr")
max_core <- suppressWarnings(as.integer(opt(opts, "max-core", "2000")))
max_dist <- suppressWarnings(as.integer(opt(opts, "max-dist", "50")))
win_kb <- suppressWarnings(as.numeric(opt(opts, "win-kb", "100")))
maf_cut <- suppressWarnings(as.numeric(opt(opts, "maf", "0.05")))
stop_ehh <- suppressWarnings(as.numeric(opt(opts, "stop-ehh", "0.01")))
outdir <- opt(opts, "outdir", "results/popgen/07_xpehh")

if (any(vapply(list(genotype_path, snp_info_path, metadata_path, g1, g2), is.null, logical(1)))) {
  stop("Required: --genotype, --snp-info, --metadata, --group1, --group2. See --help.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills XP-EHH selection scan ==")

# ---------------------------------------------------------------------------
# 1) Load + split by population
# ---------------------------------------------------------------------------
gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
G <- gt$geno
meta <- read_pheno(metadata_path)
if (!group_col %in% colnames(meta)) stop("Metadata has no '", group_col, "' column.", call. = FALSE)
s1 <- intersect(meta$sample_id[meta[[group_col]] == g1], rownames(G))
s2 <- intersect(meta$sample_id[meta[[group_col]] == g2], rownames(G))
if (length(s1) < 10 || length(s2) < 10) {
  stop("Need >= 10 samples per population (got ", length(s1), " / ", length(s2), ").", call. = FALSE)
}
message(sprintf("  populations: %s n=%d | %s n=%d", g1, length(s1), g2, length(s2)))

if (!is.null(chr_arg)) {
  keep_c <- snp_info$chr == chr_arg
  if (!any(keep_c)) stop("No SNPs on chromosome ", chr_arg, ".", call. = FALSE)
  G <- G[, keep_c, drop = FALSE]
  snp_info <- snp_info[keep_c, , drop = FALSE]
  message(sprintf("  chromosome %s: %d SNPs", chr_arg, nrow(snp_info)))
}
if (anyNA(G)) G <- round(impute_mean(G))

# core SNPs: MAF filter + even spacing
info <- snp_info[match(colnames(G), snp_info$snp_id), , drop = FALSE]
maf <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)
core_idx <- which(maf >= maf_cut)
if (length(core_idx) > max_core) {
  core_idx <- core_idx[round(seq(1, length(core_idx), length.out = max_core))]
}
message(sprintf("  core SNPs scored: %d (MAF >= %g)", length(core_idx), maf_cut))

# ---------------------------------------------------------------------------
# 2) Diploid -> two chromosomes (random phase)
# ---------------------------------------------------------------------------
split_diploid <- function(g_sub) {
  n <- nrow(g_sub); m <- ncol(g_sub)
  h1 <- matrix(0L, n, m); h2 <- matrix(0L, n, m)
  ones <- g_sub == 1
  r <- matrix(runif(n * m), n, m)
  h1[ones] <- as.integer(r[ones] < 0.5)
  h2[ones] <- 1L - h1[ones]
  h1[g_sub == 2] <- 1L; h2[g_sub == 2] <- 1L
  rbind(h1, h2)
}
H1 <- split_diploid(G[s1, , drop = FALSE])
H2 <- split_diploid(G[s2, , drop = FALSE])

# ---------------------------------------------------------------------------
# 3) iEHH at a core SNP (exact haplotype counting, cumulative window)
# ---------------------------------------------------------------------------
iehh_core <- function(H, core_j, max_dist, stop_ehh) {
  n_total <- ncol(H)
  carry <- which(H[, core_j] == 1)
  if (length(carry) < 3) return(NA_real_)
  max_l <- min(max_dist, core_j - 1)
  max_r <- min(max_dist, n_total - core_j)
  if (max_l + max_r < 1) return(NA_real_)
  m_carry <- length(carry)
  ehhs <- numeric(max_l + max_r + 1)
  ehhs[1] <- 1
  for (d in seq_len(max_l + max_r)) {
    l_len <- min(d, max_l); r_len <- min(d, max_r)
    cols <- unique(c(core_j - seq_len(l_len), core_j + seq_len(r_len)))
    cols <- cols[cols >= 1 & cols <= n_total]
    code <- as.vector(H[carry, cols, drop = FALSE] %*% 2^(seq_along(cols) - 1))
    tb <- table(code)
    ehhs[d + 1] <- sum(tb^2) / m_carry^2
    if (ehhs[d + 1] <= stop_ehh) {
      ehhs[(d + 2):(max_l + max_r + 1)] <- 0
      break
    }
  }
  # integrate EHH over SNP distance (trapezoid; uniform spacing approximation)
  sum(diff(seq_along(ehhs) - 1) * (head(ehhs, -1) + tail(ehhs, -1)) / 2)
}

# ---------------------------------------------------------------------------
# 4) Score core SNPs
# ---------------------------------------------------------------------------
message("  scoring core SNPs (slow step)...")
res_list <- vector("list", length(core_idx))
pb_every <- max(1, length(core_idx) %/% 10)
k_done <- 0
for (k in seq_along(core_idx)) {
  cj <- core_idx[k]
  if (k %% pb_every == 0) message(sprintf("    %d / %d", k, length(core_idx)))
  i1 <- iehh_core(H1, cj, max_dist, stop_ehh)
  i2 <- iehh_core(H2, cj, max_dist, stop_ehh)
  if (is.na(i1) || is.na(i2) || (i1 + i2) <= 0) next
  k_done <- k_done + 1
  res_list[[k_done]] <- data.frame(
    snp_id = colnames(G)[cj], chr = info$chr[cj], pos = info$pos[cj],
    iehh_1 = i1, iehh_2 = i2,
    xpehh_raw = log((i1 + 1e-9) / (i2 + 1e-9)), stringsAsFactors = FALSE)
}
sites <- do.call(rbind, res_list[seq_len(k_done)])
if (is.null(sites) || nrow(sites) < 50) {
  stop("Too few core SNPs scored; increase --max-dist or check group labels.", call. = FALSE)
}

# window standardisation (MAD-based z). Primary: physical windows of win_kb;
# fallback (sparse data): SNP-count windows of 25 SNPs so that small
# chromosomes/panels still get local standardisation.
win_bp <- win_kb * 1000
sites$window <- floor(sites$pos / win_bp)

standardise_windows <- function(df) {
  out <- lapply(split(df, df$window), function(d) {
    if (nrow(d) < 3) return(NULL)
    med <- stats::median(d$xpehh_raw, na.rm = TRUE)
    madv <- stats::mad(d$xpehh_raw, na.rm = TRUE)
    if (is.na(madv) || madv == 0) return(NULL)
    d$z <- (d$xpehh_raw - med) / (1.4826 * madv)
    d
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

zw_sites <- standardise_windows(sites)
window_scheme <- sprintf("physical %g kb", win_kb)
if (is.null(zw_sites)) {
  sites$window <- ceiling(seq_len(nrow(sites)) / 25)
  zw_sites <- standardise_windows(sites)
  window_scheme <- "SNP-count windows (25 SNPs)"
  message("  (sparse data: falling back to SNP-count windows)")
}
if (!is.null(zw_sites)) sites$z <- zw_sites$z[match(rownames(sites), rownames(zw_sites))]
utils::write.csv(sites, file.path(outdir, "xpehh_sites.csv"), row.names = FALSE)

win_summary <- if (!is.null(zw_sites)) do.call(rbind, lapply(split(zw_sites, zw_sites$window), function(d) {
  data.frame(chr = d$chr[1], start = min(d$pos), end = max(d$pos),
             n_snps = nrow(d), mean_z = round(mean(d$z), 3),
             max_abs_z = round(max(abs(d$z)), 3),
             stringsAsFactors = FALSE)
})) else NULL
if (!is.null(win_summary) && nrow(win_summary)) {
  win_summary$candidate <- win_summary$max_abs_z >= 2
  win_summary <- win_summary[order(win_summary$chr, win_summary$start), , drop = FALSE]
  utils::write.csv(win_summary, file.path(outdir, "xpehh_windows.csv"), row.names = FALSE)
  cand <- win_summary[win_summary$candidate, , drop = FALSE]
  cand <- cand[order(-cand$max_abs_z), , drop = FALSE]
  if (nrow(cand)) utils::write.csv(cand, file.path(outdir, "candidate_windows.csv"), row.names = FALSE)
} else {
  win_summary <- data.frame()
  cand <- data.frame()
}
message(sprintf("  scored: %d core SNPs | windows (%s): %s | candidates (max|z| >= 2): %d",
                nrow(sites), window_scheme,
                if (nrow(win_summary)) nrow(win_summary) else 0, nrow(cand)))

# ---------------------------------------------------------------------------
# 5) Figures
# ---------------------------------------------------------------------------
pal <- rice_palette()
p_dist <- ggplot2::ggplot(sites, ggplot2::aes(x = xpehh_raw)) +
  ggplot2::geom_histogram(bins = 60, fill = pal$groups5[1], colour = "white", linewidth = 0.2) +
  ggplot2::labs(x = sprintf("raw XP-EHH = ln(iEHH_%s / iEHH_%s)", g1, g2),
                y = "Core SNPs", title = "XP-EHH raw score distribution") +
  theme_nature()
save_nature_figure(p_dist, file.path(outdir, "xpehh_distribution"), width_mm = 100, height_mm = 65)

p_man <- ggplot2::ggplot(sites[!is.na(sites$z), ], ggplot2::aes(x = pos / 1e6, y = z)) +  ggplot2::geom_point(size = 0.6, alpha = 0.7, colour = pal$groups5[1]) +
  ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed",
                      colour = "grey50", linewidth = 0.3) +
  ggplot2::labs(x = sprintf("%s position (Mb)", unique(sites$chr)[1]),
                y = sprintf("XP-EHH z (window %g kb; positive = %s)", win_kb, g1),
                title = "XP-EHH selection scan") +
  theme_nature()
save_nature_figure(p_man, file.path(outdir, "xpehh_manhattan"), width_mm = 180, height_mm = 75)

sum_lines <- c(
  "Rice-skills -- XP-EHH summary",
  sprintf("populations        : %s vs %s", g1, g2),
  sprintf("core SNPs scored   : %d", nrow(sites)),
  sprintf("windows            : %d", nrow(win_summary)),
  sprintf("candidates (max|z| >= 2): %d", nrow(cand))
)
if (nrow(cand)) {
  top <- head(sprintf("%s:%d-%d (|z|=%.2f)", cand$chr, cand$start, cand$end, cand$max_abs_z), 5)
  sum_lines <- c(sum_lines, paste0("top candidates     : ", top))
}
writeLines(sum_lines, file.path(outdir, "xpehh_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

write_run_info_r(outdir, "07_xpehh.R", params = list(
  genotype = genotype_path, group1 = g1, group2 = g2,
  chr = if (is.null(chr_arg)) "all" else chr_arg,
  max_core = max_core, max_dist = max_dist, win_kb = win_kb,
  n_core_scored = nrow(sites), n_candidates = nrow(cand)
))
message("Done. Outputs written to: ", normalizePath(outdir))
