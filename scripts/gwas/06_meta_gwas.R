###############################################################################
# 06_meta_gwas.R -- meta-analysis of multiple GWAS studies (inverse variance)
#
# Title:    Rice-skills meta-GWAS step
# Description:
#   Combines association results from independent GWAS runs (cohorts /
#   environments / subpanels) with inverse-variance weighted fixed-effects
#   meta-analysis:
#     beta_meta = sum(w_i * beta_i) / sum(w_i),  w_i = 1 / se_i^2
#     se_meta   = sqrt(1 / sum(w_i))
#     z = beta_meta / se_meta,  p = 2 * pnorm(-|z|)
#     heterogeneity: Cochran's Q, I^2, per-SNP study count
#   Only SNPs present in >= 2 studies are tested. Outputs include
#   forest plots for the top loci.
#   This is the standard route to raise power by pooling samples across
#   studies (e.g. multi-environment / multi-panel rice GWAS).
#
# Usage:
#   Rscript scripts/gwas/06_meta_gwas.R \
#       --studies results/gwas/studyA/association_results.csv,results/gwas/studyB/association_results.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --outdir results/gwas/06_meta
#
# Inputs:
#   --studies    comma-separated list of association result CSVs; each must
#                contain snp_id plus either (beta, se) or a z column.
#                (Scripts in this package emit snp_id, beta, se, t, p.)
#   --study-names optional comma-separated labels (default file basenames)
#   --snp-info   snp_id, chr, pos (for Manhattan/QQ)
#   --min-studies  minimum number of studies per SNP (default 2)
#   --p-threshold  significance threshold (default Bonferroni 0.05 / n)
#   --top-forest   number of top loci shown in the forest plot (default 8)
#   --outdir     output directory (default results/gwas/06_meta)
# Outputs:
#   meta_results.csv       beta, se, z, p, n_studies, Q, Q_p, I2 per SNP
#   meta_significant.csv   SNPs passing the threshold
#   manhattan_meta.* / qq_meta.*   figures (PDF/TIFF/PNG)
#   forest_top.*           forest plot of top loci
#   run_info.txt
#
# Notes:
#   - Fixed-effects assumes the same underlying effect across studies;
#     inspect I^2 (high I^2 -> consider random-effects or stratified views).
#   - For meta across studies of the SAME individuals (overlapping samples)
#     results are anticonservative — keep study panels independent.
#   - Random-effects (DerSimonian-Laird) per SNP is reported when I^2 > 50%:
#     see the `re_beta`/`re_p` columns in meta_results.csv.
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
studies_arg <- opt(opts, "studies")
names_arg <- opt(opts, "study-names")
snp_info_path <- opt(opts, "snp-info")
min_studies <- suppressWarnings(as.integer(opt(opts, "min-studies", "2")))
p_threshold <- opt(opts, "p-threshold", NULL)
top_forest <- suppressWarnings(as.integer(opt(opts, "top-forest", "8")))
outdir <- opt(opts, "outdir", "results/gwas/06_meta")

if (is.null(studies_arg) || is.null(snp_info_path)) {
  stop("Required: --studies a.csv,b.csv and --snp-info <table>.\n",
       "Hint: see the Usage block at the top of this script.", call. = FALSE)
}
study_paths <- trimws(strsplit(studies_arg, ",")[[1]])
study_paths <- study_paths[nzchar(study_paths)]
if (length(study_paths) < 2) stop("Meta-analysis needs at least 2 study files.", call. = FALSE)
for (p in study_paths) {
  if (!file.exists(p)) stop("Study file not found: ", p, call. = FALSE)
}
study_names <- if (!is.null(names_arg)) {
  trimws(strsplit(names_arg, ",")[[1]])
} else {
  # auto names: parent-dir prefix disambiguates same-named files
  # (e.g. studyA/association_results.csv vs studyB/association_results.csv)
  base_names <- tools::file_path_sans_ext(basename(study_paths))
  parent <- basename(dirname(study_paths))
  auto <- ifelse(parent %in% c(".", "", "/"), base_names, paste0(parent, ".", base_names))
  make.unique(auto)
}
if (length(study_names) != length(study_paths)) {
  stop("--study-names count (", length(study_names), ") != --studies count (",
       length(study_paths), ")", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills meta-GWAS (", length(study_paths), " studies) ==")

# ---------------------------------------------------------------------------
# 1) Read studies, harmonize to (snp_id, beta, se)
# ---------------------------------------------------------------------------
studies <- vector("list", length(study_paths))
for (i in seq_along(study_paths)) {
  d <- read_csv_safe(study_paths[i], paste0("study '", study_names[i], "'"))
  check_columns(d, "snp_id", paste0("study '", study_names[i], "'"))
  if (all(c("beta", "se") %in% colnames(d))) {
    d <- d[!is.na(d$beta) & !is.na(d$se) & d$se > 0, c("snp_id", "beta", "se")]
  } else if ("z" %in% colnames(d)) {
    d <- data.frame(snp_id = d$snp_id, beta = d$z, se = 1)
    message("  note: '", study_names[i], "' uses z values (meta on the z scale)")
  } else {
    stop("Study '", study_names[i], "' needs columns (beta, se) or z.", call. = FALSE)
  }
  studies[[i]] <- d
  message(sprintf("  %-18s %8d SNPs", study_names[i], nrow(d)))
}

# ---------------------------------------------------------------------------
# 2) Per-SNP inverse-variance meta
# ---------------------------------------------------------------------------
all_snps <- unique(unlist(lapply(studies, function(d) d$snp_id)))
# SNP count per study
cnt <- table(unlist(lapply(studies, function(d) d$snp_id)))
common <- names(cnt)[cnt >= min_studies]
message(sprintf("  SNPs: %d total | %d in >= %d studies", length(all_snps),
                length(common), min_studies))
if (length(common) < 10) stop("Fewer than 10 SNPs shared across studies.", call. = FALSE)

m <- length(common)
k <- length(studies)
beta_mat <- matrix(NA_real_, k, m, dimnames = list(study_names, common))
se_mat <- beta_mat
for (i in seq_len(k)) {
  idx <- match(common, studies[[i]]$snp_id)
  ok <- !is.na(idx)
  beta_mat[i, ok] <- studies[[i]]$beta[idx[ok]]
  se_mat[i, ok] <- studies[[i]]$se[idx[ok]]
}
w <- 1 / se_mat^2
w[is.na(w)] <- 0
w_sum <- colSums(w)
beta_meta <- colSums(w * beta_mat, na.rm = TRUE) / w_sum
se_meta <- sqrt(1 / w_sum)
z <- beta_meta / se_meta
p <- 2 * stats::pnorm(-abs(z))
n_studies <- colSums(!is.na(beta_mat))

# heterogeneity
Q <- colSums(w * (beta_mat - rep(beta_meta, each = k))^2, na.rm = TRUE)
df_q <- pmax(n_studies - 1, 1)
Q_p <- stats::pchisq(Q, df = df_q, lower.tail = FALSE)
I2 <- pmax(0, (Q - df_q) / pmax(Q, 1e-300)) * 100

# random-effects (DerSimonian-Laird) for high-heterogeneity SNPs
tau2 <- pmax(0, (Q - df_q) / (w_sum - colSums(w^2, na.rm = TRUE) / w_sum))
w_re <- 1 / (se_mat^2 + rep(tau2, each = k))
w_re[is.na(w_re)] <- 0
w_re_sum <- colSums(w_re)
re_beta <- colSums(w_re * beta_mat, na.rm = TRUE) / w_re_sum
re_se <- sqrt(1 / w_re_sum)
re_p <- 2 * stats::pnorm(-abs(re_beta / re_se))

res <- data.frame(
  snp_id = common, beta = beta_meta, se = se_meta, z = z, p = p,
  n_studies = as.integer(n_studies), Q = Q, Q_p = Q_p, I2 = I2,
  re_beta = re_beta, re_p = re_p, stringsAsFactors = FALSE
)
res <- res[order(res$p), , drop = FALSE]

# optional chr/pos join for output table
snp_info <- read_snp_info(snp_info_path)
res_out <- merge(res, snp_info, by = "snp_id", all.x = TRUE)
utils::write.csv(res_out, file.path(outdir, "meta_results.csv"), row.names = FALSE)

thresh <- if (is.null(p_threshold)) 0.05 / nrow(res) else as.numeric(p_threshold)
sig <- res_out[!is.na(res_out$p) & res_out$p <= thresh, , drop = FALSE]
utils::write.csv(sig, file.path(outdir, "meta_significant.csv"), row.names = FALSE)
message(sprintf("  meta-tested SNPs: %d | p <= %.2e: %d | high-I2 (I2>50%%): %d",
                nrow(res), thresh, nrow(sig), sum(res$I2 > 50, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 3) Figures: Manhattan, QQ, forest of top loci
# ---------------------------------------------------------------------------
info_plot <- res_out[!is.na(res_out$chr), c("snp_id", "chr", "pos", "p")]
if (nrow(info_plot) >= 10) {
  p_man <- manhattan_plot(info_plot[, c("snp_id", "p")], info_plot[, c("snp_id", "chr", "pos")],
                          threshold = thresh, title = "Meta-GWAS (inverse variance)")
  save_nature_figure(p_man, file.path(outdir, "manhattan_meta"), width_mm = 180, height_mm = 80)
  p_qq <- qq_plot(info_plot[, c("snp_id", "p")], title = "Meta-GWAS QQ")
  save_nature_figure(p_qq, file.path(outdir, "qq_meta"), width_mm = 89, height_mm = 89)
}

if (nrow(res) >= 1 && top_forest >= 1) {
  top_snps <- head(res$snp_id, top_forest)
  frows <- list()
  for (s in top_snps) {
    for (i in seq_len(k)) {
      idx <- match(s, studies[[i]]$snp_id)
      if (is.na(idx)) next
      frows[[length(frows) + 1L]] <- data.frame(
        snp = s, study = study_names[i], beta = studies[[i]]$beta[idx],
        se = studies[[i]]$se[idx], type = "study", stringsAsFactors = FALSE)
    }
    frows[[length(frows) + 1L]] <- data.frame(
      snp = s, study = "Meta (IVW)", beta = res$beta[match(s, res$snp_id)],
      se = res$se[match(s, res$snp_id)], type = "meta", stringsAsFactors = FALSE)
  }
  fd <- do.call(rbind, frows)
  fd$snp <- factor(fd$snp, levels = top_snps)
  fd$study <- factor(fd$study, levels = rev(c(study_names, "Meta (IVW)")))
  pal <- rice_palette()
  p_forest <- ggplot2::ggplot(fd, ggplot2::aes(x = beta, y = study)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se),
                           orientation = "y", width = 0.18, linewidth = 0.35,
                           colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(fill = type), shape = 21, size = 1.7, colour = "grey20", stroke = 0.2) +
    ggplot2::scale_fill_manual(values = c(study = pal$groups5[1], meta = pal$groups5[2])) +
    ggplot2::facet_wrap(~ snp, scales = "free_x", ncol = 2) +
    ggplot2::labs(x = "Effect size (beta, 95% CI)", y = NULL, title = "Top loci — forest plot") +
    theme_nature() + ggplot2::theme(legend.position = "none")
  save_nature_figure(p_forest, file.path(outdir, "forest_top"), width_mm = 160, height_mm = 150)
}

write_run_info_r(outdir, "06_meta_gwas.R", params = list(
  studies = paste(study_names, collapse = ", "), min_studies = min_studies,
  n_snps_meta = nrow(res), n_significant = nrow(sig), threshold = thresh
))
message("Done. Outputs written to: ", normalizePath(outdir))
