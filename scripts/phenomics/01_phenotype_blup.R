###############################################################################
# 01_phenotype_blup.R -- phenotype QC, BLUP estimation and heritability
#
# Title:    Rice-skills phenotype preprocessing step
# Description:
#   Phenotype pipeline for multi-environment trial data (line x environment):
#     1. descriptive QC and MAD-based outlier flagging (never removes data);
#     2. environment-adjusted, shrinkage-weighted BLUP per line and trait
#        (moment-based variance components; see Notes) — the BLUP table is
#        the recommended input for GWAS (scripts/gwas/) and GS;
#     3. broad-sense heritability h2 per trait;
#     4. figures: trait distributions by environment, trait correlation.
#   If lme4 is installed the script additionally reports REML variance
#   components for cross-checking (the BLUP table itself stays deterministic
#   and dependency-free).
#
# Usage:
#   Rscript scripts/phenomics/01_phenotype_blup.R \
#       --pheno demo_data/phenomics/phenotype_long.csv \
#       --env-col env --traits grain_yield,plant_height \
#       --outdir results/phenomics/01_blup
#
# Inputs:
#   --pheno     phenotype table CSV: first column = line_id, one column of
#               environments (--env-col), then numeric trait columns.
#               Without --env-col a single-environment table still gets QC
#               and BLUP (environment adjustment skipped).
#   --env-col   environment column name (optional)
#   --traits    comma-separated trait names (default: all numeric columns)
#   --mad-k     outlier threshold in MAD units (default 3)
#   --outdir    output directory (default results/phenomics/01_blup)
# Outputs:
#   blups.csv          line_id + per-trait BLUP (shrinkage on the adjusted mean)
#   heritability.csv   per-trait h2, sigma2_g, sigma2_e, n_obs
#   outliers.csv       flagged records (trait, value, z-MAD) — for review only
#   trait_distribution.* / trait_correlation.*   figures
#   qc_summary.txt / run_info.txt
#
# Notes:
#   - Variance components use the moment (ANOVA) estimator on the two-way
#     line + environment model:
#       sigma2_e = MSE,  sigma2_g = max(0, (MS_line - MSE) / r_bar)
#     and BLUP_i = sigma2_g / (sigma2_g + sigma2_e / r_i) * (adjusted_i - mu).
#     This is the classic one-way-shrinkage approximation; for unbalanced
#     designs or pedigrees use lme4/ASReml (cross-check section runs when
#     lme4 is present).
#   - h2 here is broad-sense on a line-mean basis: h2 = sigma2_g / (sigma2_g
#     + sigma2_e / r_bar).
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
pheno_path <- opt(opts, "pheno")
env_col <- opt(opts, "env-col")
traits_arg <- opt(opts, "traits")
mad_k <- suppressWarnings(as.numeric(opt(opts, "mad-k", "3")))
outdir <- opt(opts, "outdir", "results/phenomics/01_blup")

if (is.null(pheno_path)) stop("Required: --pheno <table.csv>. See --help.", call. = FALSE)
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills phenotype BLUP ==")

d <- read_csv_safe(pheno_path, "phenotype table")
line_col <- colnames(d)[1]
d[[line_col]] <- as.character(d[[line_col]])
if (!is.null(env_col)) {
  if (!env_col %in% colnames(d)) stop("--env-col '", env_col, "' not found. Columns: ",
                                      paste(colnames(d), collapse = ", "), call. = FALSE)
  d[[env_col]] <- as.character(d[[env_col]])
}
reserved <- c(line_col, env_col)
cand_traits <- setdiff(colnames(d), reserved)
if (!is.null(traits_arg)) {
  traits <- trimws(strsplit(traits_arg, ",")[[1]])
  miss <- setdiff(traits, colnames(d))
  if (length(miss)) stop("Trait(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
} else {
  traits <- cand_traits[vapply(d[, cand_traits, drop = FALSE], is.numeric, logical(1))]
}
traits <- traits[vapply(d[, traits, drop = FALSE], is.numeric, logical(1))]
if (!length(traits)) stop("No numeric trait columns found.", call. = FALSE)

n_lines <- length(unique(d[[line_col]]))
n_envs <- if (!is.null(env_col)) length(unique(d[[env_col]])) else 1L
message(sprintf("  lines: %d | environments: %d | traits: %s | records: %d",
                n_lines, n_envs, paste(traits, collapse = ", "), nrow(d)))

# ---------------------------------------------------------------------------
# 1) QC + outlier flagging (MAD)
# ---------------------------------------------------------------------------
out_rows <- list()
for (tr in traits) {
  x <- d[[tr]]
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, na.rm = TRUE)
  if (is.na(madv) || madv == 0) next
  z <- abs(x - med) / (1.4826 * madv)
  idx <- which(!is.na(z) & z > mad_k)
  if (length(idx)) {
    out_rows[[length(out_rows) + 1L]] <- data.frame(
      line_id = d[[line_col]][idx],
      env = if (!is.null(env_col)) d[[env_col]][idx] else NA_character_,
      trait = tr, value = x[idx], mad_z = round(z[idx], 2), stringsAsFactors = FALSE)
  }
}
out_df <- if (length(out_rows)) do.call(rbind, out_rows) else
  data.frame(line_id = character(0), env = character(0), trait = character(0),
             value = numeric(0), mad_z = numeric(0))
utils::write.csv(out_df, file.path(outdir, "outliers.csv"), row.names = FALSE)
message(sprintf("  outlier flags: %d record(s) (flagged only, not removed)", nrow(out_df)))

# ---------------------------------------------------------------------------
# 2) BLUP + heritability
# ---------------------------------------------------------------------------
blup_list <- list()
h2_rows <- list()
for (tr in traits) {
  sub <- d[, c(line_col, if (!is.null(env_col)) env_col, tr)]
  colnames(sub) <- c("line", if (!is.null(env_col)) "env", "y")
  sub <- sub[!is.na(sub$y), , drop = FALSE]
  if (!nrow(sub)) next

  mu <- mean(sub$y)
  r_i <- table(sub$line)
  r_bar <- mean(r_i)

  if (!is.null(env_col)) {
    # environment adjustment: subtract env mean, add grand mean
    env_means <- tapply(sub$y, sub$env, mean)
    sub$y_adj <- sub$y - env_means[sub$env] + mu
    ms <- stats::anova(stats::lm(y ~ factor(line) + factor(env), data = sub))
    ms_line <- ms["factor(line)", "Mean Sq"]
    ms_res <- ms["Residuals", "Mean Sq"]
  } else {
    sub$y_adj <- sub$y
    ms <- stats::anova(stats::lm(y ~ factor(line), data = sub))
    ms_line <- ms["factor(line)", "Mean Sq"]
    ms_res <- ms["Residuals", "Mean Sq"]
  }
  sigma2_g <- max(0, (ms_line - ms_res) / r_bar)
  sigma2_e <- ms_res
  h2 <- sigma2_g / (sigma2_g + sigma2_e / r_bar)

  adj_means <- tapply(sub$y_adj, sub$line, mean)
  lambda <- sigma2_g / (sigma2_g + sigma2_e / r_i[names(adj_means)])
  blup <- lambda * (adj_means - mu)

  blup_list[[tr]] <- data.frame(line_id = names(blup), value = as.numeric(blup),
                                stringsAsFactors = FALSE)
  colnames(blup_list[[tr]])[2] <- tr
  h2_rows[[length(h2_rows) + 1L]] <- data.frame(
    trait = tr, h2 = round(h2, 4), sigma2_g = round(sigma2_g, 4),
    sigma2_e = round(sigma2_e, 4), n_obs = nrow(sub), r_bar = round(r_bar, 2),
    stringsAsFactors = FALSE)
  message(sprintf("  [%s] h2(mean-basis) = %.3f (sigma2_g=%.3g, sigma2_e=%.3g)",
                  tr, h2, sigma2_g, sigma2_e))
}

if (length(blup_list)) {
  blups <- Reduce(function(a, b) merge(a, b, by = "line_id", all = TRUE), blup_list)
  utils::write.csv(blups, file.path(outdir, "blups.csv"), row.names = FALSE)
}
h2_df <- do.call(rbind, h2_rows)
utils::write.csv(h2_df, file.path(outdir, "heritability.csv"), row.names = FALSE)

# ---- optional lme4 cross-check (REML variance components) ----
if (requireNamespace("lme4", quietly = TRUE) && !is.null(env_col)) {
  message("  lme4 cross-check (REML):")
  for (tr in traits) {
    sub <- d[, c(line_col, env_col, tr)]
    colnames(sub) <- c("line", "env", "y")
    sub <- sub[!is.na(sub$y), , drop = FALSE]
    fit <- tryCatch(lme4::lmer(y ~ (1 | line) + (1 | env), data = sub), error = function(e) NULL)
    if (!is.null(fit)) {
      vc <- as.data.frame(lme4::VarCorr(fit))
      message(sprintf("    [%s] REML sigma2_line=%.4g sigma2_env=%.4g sigma2_res=%.4g",
                      tr, vc$vcov[vc$grp == "line"], vc$vcov[vc$grp == "env"],
                      vc$vcov[vc$grp == "Residual"]))
    }
  }
}

# ---------------------------------------------------------------------------
# 3) Figures
# ---------------------------------------------------------------------------
long <- do.call(rbind, lapply(traits, function(tr) {
  data.frame(value = d[[tr]], trait = tr,
             env = if (!is.null(env_col)) d[[env_col]] else "single")
}))
long <- long[!is.na(long$value), , drop = FALSE]
pal <- rice_palette()
n_tr <- length(unique(long$trait))
n_env <- length(unique(long$env))
facet_labeller <- if (n_tr > 1 && n_env > 1) {
  ggplot2::facet_grid(trait ~ env, scales = "free_y")
} else if (n_env > 1) {
  ggplot2::facet_wrap(~ env, scales = "free_y")
} else {
  ggplot2::facet_wrap(~ trait, scales = "free_y")
}
p_dist <- ggplot2::ggplot(long, ggplot2::aes(x = env, y = value)) +
  ggplot2::geom_boxplot(outlier.size = 0.4, linewidth = 0.3, fill = pal$groups5[1],
                        colour = "grey25", alpha = 0.85) +
  facet_labeller +
  ggplot2::labs(x = NULL, y = "Phenotype value", title = "Trait distribution by environment") +
  theme_nature()
save_nature_figure(p_dist, file.path(outdir, "trait_distribution"),
                   width_mm = 150, height_mm = 100)

if (length(traits) >= 2) {
  wide <- d[, traits, drop = FALSE]
  cm <- stats::cor(wide, use = "pairwise.complete.obs", method = "pearson")
  df_c <- data.frame(
    t1 = rep(rownames(cm), times = ncol(cm)),
    t2 = rep(colnames(cm), each = nrow(cm)),
    r = as.vector(cm))
  p_cor <- ggplot2::ggplot(df_c, ggplot2::aes(t1, t2, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", r)), size = 2.2) +
    ggplot2::scale_fill_gradient2(low = diverging_palette()[1], mid = "white",
                                  high = diverging_palette()[7], midpoint = 0,
                                  limits = c(-1, 1), name = "r") +
    ggplot2::labs(x = NULL, y = NULL, title = "Trait correlations") +
    theme_nature()
  save_nature_figure(p_cor, file.path(outdir, "trait_correlation"),
                     width_mm = 90, height_mm = 80)
}

# ---------------------------------------------------------------------------
# 4) QC summary
# ---------------------------------------------------------------------------
qc_lines <- c(
  "Rice-skills -- phenotype QC summary",
  sprintf("lines              : %d", n_lines),
  sprintf("environments       : %d", n_envs),
  sprintf("records            : %d", nrow(d)),
  sprintf("outlier flags      : %d (see outliers.csv)", nrow(out_df)),
  sprintf("traits analysed    : %s", paste(traits, collapse = ", "))
)
for (i in seq_len(nrow(h2_df))) {
  qc_lines <- c(qc_lines, sprintf("  %-18s h2 = %.3f", h2_df$trait[i], h2_df$h2[i]))
}
writeLines(qc_lines, file.path(outdir, "qc_summary.txt"))

write_run_info_r(outdir, "01_phenotype_blup.R", params = list(
  pheno = pheno_path,
  env_col = if (is.null(env_col)) NA_character_ else env_col,
  traits = paste(traits, collapse = ", "), mad_k = mad_k,
  n_lines = n_lines, n_envs = n_envs
))
message("Done. BLUP table: ", file.path(outdir, "blups.csv"))
