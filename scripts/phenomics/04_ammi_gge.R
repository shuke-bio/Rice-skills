###############################################################################
# 04_ammi_gge.R -- AMMI / GGE biplot for multi-environment trials
#
# Title:    Rice-skills AMMI-GGE step
# Description:
#   Cultivar x environment interaction analysis on a multi-environment
#   phenotype table (line x env, one trait):
#     - AMMI: additive line + environment effects, SVD of the interaction,
#       AMMI1 biplot (mean vs IPC1) and AMMI2 biplot (IPC1 vs IPC2);
#     - GGE: environment-centred SVD, "which-won-where" polygon view with
#       the average-environment coordination (AEC) circles;
#     - stability table: ASV (AMMI stability value), across-env mean, CV.
#
# Usage:
#   Rscript scripts/phenomics/04_ammi_gge.R \
#       --pheno demo_data/phenomics/phenotype_long.csv --env-col env \
#       --trait grain_yield --outdir results/phenomics/04_ammi
#
# Inputs:
#   --pheno    phenotype table CSV: first column line_id, env column,
#              numeric trait columns
#   --env-col  environment column (default "env")
#   --trait    trait to analyse (default: first numeric column)
#   --n-pc     PCs kept for biplots (default 2)
#   --outdir   output directory (default results/phenomics/04_ammi)
# Outputs:
#   ammi_scores.csv     line/env scores on the retained PCs + means
#   stability.csv       ASV, mean, CV per line
#   ammi1_biplot.*      mean vs IPC1 (PDF/TIFF/PNG)
#   gge_biplot.*        PC1 vs PC2 with which-won-where polygon + AEC
#   ammi_summary.txt    variance explained, top stable/unstable lines
#   run_info.txt
#
# Notes:
#   - Missing line x env combinations are filled with the additive model
#     prediction (row + column effect); the fill count is reported.
#   - AMMI/GGE here are the classic Gauch/Yan constructions on the
#     line-mean matrix; plots follow the shared Nature figure spec.
###############################################################################

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source(file.path(dirname(script_file), "..", "utils", "plot_theme.R"))
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
env_col <- opt(opts, "env-col", "env")
trait <- opt(opts, "trait")
n_pc <- suppressWarnings(as.integer(opt(opts, "n-pc", "2")))
outdir <- opt(opts, "outdir", "results/phenomics/04_ammi")

if (is.null(pheno_path)) stop("Required: --pheno <table.csv>. See --help.", call. = FALSE)
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills AMMI / GGE biplot ==")

d <- read_csv_safe(pheno_path, "phenotype table")
line_col <- colnames(d)[1]
d[[line_col]] <- as.character(d[[line_col]])
if (!env_col %in% colnames(d)) stop("No env column '", env_col, "'.", call. = FALSE)
d[[env_col]] <- as.character(d[[env_col]])
if (is.null(trait)) {
  cand <- setdiff(colnames(d), c(line_col, env_col))
  trait <- cand[vapply(d[, cand, drop = FALSE], is.numeric, logical(1))][1]
}
if (is.na(trait) || !trait %in% colnames(d)) stop("Trait not found; pass --trait.", call. = FALSE)
y_all <- suppressWarnings(as.numeric(d[[trait]]))
d <- d[!is.na(y_all), , drop = FALSE]

# ---------------------------------------------------------------------------
# 1) line x env matrix (mean over reps), additive fill for missing cells
# ---------------------------------------------------------------------------
mat <- tapply(d[[trait]], list(d[[line_col]], d[[env_col]]), mean)
n_missing <- sum(is.na(mat))
if (n_missing) {
  grand0 <- mean(mat, na.rm = TRUE)
  r_eff <- rowMeans(mat, na.rm = TRUE) - grand0
  c_eff <- colMeans(mat, na.rm = TRUE) - grand0
  for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
    if (is.na(mat[i, j])) mat[i, j] <- grand0 + r_eff[i] + c_eff[j]
  }
  message(sprintf("  filled %d missing line x env cells with additive predictions", n_missing))
}
n_lines <- nrow(mat); n_envs <- ncol(mat)
if (n_lines < 5 || n_envs < 2) {
  stop("Need >= 5 lines and >= 2 environments for AMMI/GGE.", call. = FALSE)
}
message(sprintf("  lines: %d | environments: %d | trait: %s", n_lines, n_envs, trait))

# ---------------------------------------------------------------------------
# 2) AMMI: additive effects + SVD of interaction
# ---------------------------------------------------------------------------
grand <- mean(mat)
row_eff <- rowMeans(mat) - grand
col_eff <- colMeans(mat) - grand
resid <- mat - grand - outer(row_eff, rep(1, ncol(mat))) - outer(rep(1, nrow(mat)), col_eff)
sv <- svd(resid)
lam <- sv$d
k <- min(n_pc, length(lam))
ve <- lam^2 / sum(lam^2) * 100
ipc_line <- sv$u[, seq_len(k), drop = FALSE] * rep(lam[seq_len(k)], each = n_lines)
ipc_env <- sv$v[, seq_len(k), drop = FALSE] * rep(lam[seq_len(k)], each = n_envs)

score_lines <- data.frame(
  line_id = rownames(mat), mean = round(rowMeans(mat), 3),
  main_effect = round(row_eff, 3),
  IPC1 = round(ipc_line[, 1], 3),
  IPC2 = if (k >= 2) round(ipc_line[, 2], 3) else NA)
score_env <- data.frame(
  env = colnames(mat), mean = round(colMeans(mat), 3),
  main_effect = round(col_eff, 3),
  IPC1 = round(ipc_env[, 1], 3),
  IPC2 = if (k >= 2) round(ipc_env[, 2], 3) else NA)
utils::write.csv(score_lines, file.path(outdir, "ammi_scores.csv"), row.names = FALSE)

# stability: ASV (AMMI stability value, Purchase 1997 weighting) + CV
pc1_w <- max(ipc_env[, 1]); pc2_w <- if (k >= 2) max(ipc_env[, 2]) else 1
asv <- sqrt((ipc_line[, 1] / pc1_w)^2 + (if (k >= 2) ipc_line[, 2] else 0)^2)
stab <- data.frame(
  line_id = rownames(mat),
  mean = round(rowMeans(mat), 3),
  asv = round(asv, 4),
  cv_env = round(apply(mat, 1, function(v) stats::sd(v) / mean(v) * 100), 2))
stab <- stab[order(stab$asv), , drop = FALSE]
utils::write.csv(stab, file.path(outdir, "stability.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3) Figures
# ---------------------------------------------------------------------------
pal <- rice_palette()
# --- AMMI1: mean vs IPC1 ---
df1 <- rbind(
  data.frame(name = rownames(mat), mean = rowMeans(mat), ipc1 = ipc_line[, 1], type = "line"),
  data.frame(name = colnames(mat), mean = colMeans(mat), ipc1 = ipc_env[, 1], type = "env"))
p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = mean, y = ipc1)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = grand, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  ggplot2::geom_point(ggplot2::aes(shape = type, colour = type), size = 1.9) +
  ggplot2::geom_text(data = df1, ggplot2::aes(label = name), size = 2.1,
                     vjust = -0.7, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = c(line = pal$groups5[1], env = pal$groups5[2])) +
  ggplot2::scale_shape_manual(values = c(line = 16, env = 17)) +
  ggplot2::labs(x = sprintf("Mean %s across environments", trait),
                y = sprintf("IPC1 (%.1f%%)", ve[1]), colour = NULL, shape = NULL,
                title = sprintf("AMMI1 biplot — %s", trait)) +
  theme_nature()
save_nature_figure(p1, file.path(outdir, "ammi1_biplot"), width_mm = 130, height_mm = 100)

# --- GGE: env-centred SVD with which-won-where polygon ---
X <- sweep(mat, 2, colMeans(mat), "-")
sv2 <- svd(X)
ve2 <- sv2$d^2 / sum(sv2$d^2) * 100
g_sc <- sv2$u[, seq_len(min(2, ncol(sv2$u))), drop = FALSE] * rep(sv2$d[seq_len(min(2, ncol(sv2$u)))], each = n_lines)
e_sc <- sv2$v[, seq_len(min(2, ncol(sv2$v))), drop = FALSE] * rep(sv2$d[seq_len(min(2, ncol(sv2$v)))], each = n_envs)
colnames(g_sc) <- c("PC1", "PC2"); colnames(e_sc) <- c("PC1", "PC2")

df_g <- data.frame(name = rownames(mat), PC1 = g_sc[, 1], PC2 = g_sc[, 2], type = "line")
df_e <- data.frame(name = colnames(mat), PC1 = e_sc[, 1], PC2 = e_sc[, 2], type = "env")
hull_idx <- grDevices::chull(df_g$PC1, df_g$PC2)
hull_df <- df_g[c(hull_idx, hull_idx[1]), ]

p2 <- ggplot2::ggplot() +
  ggplot2::geom_polygon(data = hull_df, ggplot2::aes(x = PC1, y = PC2),
                        fill = pal$groups5[1], alpha = 0.08, colour = "grey45",
                        linewidth = 0.3, linetype = "dashed") +
  ggplot2::geom_segment(data = merge(df_g, df_e, by = NULL),
                        ggplot2::aes(x = 0, y = 0, xend = PC1.y, yend = PC2.y),
                        colour = pal$groups5[2], linewidth = 0.25, alpha = 0.6,
                        show.legend = FALSE, inherit.aes = FALSE) +
  ggplot2::geom_point(data = df_e, ggplot2::aes(x = PC1, y = PC2), shape = 17,
                      colour = pal$groups5[2], size = 2.2) +
  ggplot2::geom_point(data = df_g, ggplot2::aes(x = PC1, y = PC2), shape = 16,
                      colour = pal$groups5[1], size = 1.6) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = sprintf("GGE PC1 (%.1f%%)", ve2[1]),
                y = sprintf("GGE PC2 (%.1f%%)", if (length(ve2) >= 2) ve2[2] else 0),
                title = sprintf("GGE biplot — %s (which-won-where)", trait)) +
  theme_nature()
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p2 <- p2 +
    ggrepel::geom_text_repel(data = df_g, ggplot2::aes(x = PC1, y = PC2, label = name),
                             size = 1.9, colour = "grey30", max.overlaps = 20,
                             segment.size = 0.15) +
    ggrepel::geom_text_repel(data = df_e, ggplot2::aes(x = PC1, y = PC2, label = name),
                             size = 2.4, colour = pal$groups5[2],
                             segment.size = 0.15, fontface = "bold")
}
save_nature_figure(p2, file.path(outdir, "gge_biplot"), width_mm = 120, height_mm = 110)

# ---------------------------------------------------------------------------
# 4) Summary
# ---------------------------------------------------------------------------
top_stable <- head(stab$line_id, 3)
top_yield <- head(stab[order(-stab$mean), , drop = FALSE]$line_id, 3)
sum_lines <- c(
  "Rice-skills -- AMMI/GGE summary",
  sprintf("lines x envs     : %d x %d (trait %s)", n_lines, n_envs, trait),
  sprintf("filled cells     : %d", n_missing),
  sprintf("IPC variance     : PC1 %.1f%%, PC2 %.1f%%",
          ve[1], if (length(ve) >= 2) ve[2] else 0),
  sprintf("most stable (ASV): %s", paste(top_stable, collapse = ", ")),
  sprintf("highest mean     : %s", paste(top_yield, collapse = ", "))
)
writeLines(sum_lines, file.path(outdir, "ammi_summary.txt"))
message(paste0("  ", sum_lines[-1], collapse = "\n  "))

write_run_info_r(outdir, "04_ammi_gge.R", params = list(
  pheno = pheno_path, env_col = env_col, trait = trait,
  n_lines = n_lines, n_envs = n_envs, filled = n_missing,
  ve1 = round(ve[1], 1)
))
message("Done. Outputs written to: ", normalizePath(outdir))
