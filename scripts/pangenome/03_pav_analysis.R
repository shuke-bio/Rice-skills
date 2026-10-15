###############################################################################
# 03_pav_analysis.R -- pangenome PAV (presence/absence variation) analysis
#
# Title:    Rice-skills pangenome PAV step
# Description:
#   Presence/absence variation analysis over a pangenome marker matrix:
#     binarize (optional) -> core/variable/private classification -> summary
#     stats -> frequency plot, PAV PCA, marker heatmap, and optional
#     phenotype association (Fisher / Wilcoxon) with a Manhattan plot when
#     marker coordinates are supplied.
#   Input is a marker x sample matrix — build it from APAV output, the
#   minigraph-cactus VCF, or any gene-level presence table
#   (see docs/pangenome_usage.md).
#
# Usage:
#   Rscript scripts/pangenome/03_pav_analysis.R \
#       --pav pav_matrix.csv --outdir results/pangenome/03_pav
#
#   # with phenotype association:
#   Rscript scripts/pangenome/03_pav_analysis.R --pav pav_matrix.csv \
#       --meta marker_info.csv --pheno phenotype.csv --trait grain_yield \
#       --outdir results/pangenome/03_pav
#
# Inputs:
#   --pav        marker x sample matrix CSV; first column = marker id,
#                remaining columns = samples; values 0/1 or coverage 0-1
#   --values     binary | cov | auto (default auto: 0/1 detected automatically)
#   --cov-threshold   coverage >= this counts as present (default 0.5)
#   --meta       optional marker info table: marker_id, chr, pos
#   --pheno      optional phenotype table (sample_id + trait columns)
#   --trait      trait column name (used with --pheno)
#   --n-top      markers drawn in the heatmap (default 120)
#   --outdir     output directory (default results/pangenome/03_pav)
# Outputs:
#   pav_binary.csv        binarized matrix (0/1, NA kept as NA)
#   pav_class.csv         marker_id, class, n_present, n_missing, freq
#   pav_stats.txt         summary counts (core/variable/private/absent)
#   pav_frequency.* pav_class_bar.* pav_pca.* pav_heatmap.*   figures
#   pav_association.csv   (optional) per-marker association + FDR
#   pav_assoc_manhattan.* (optional) Manhattan plot of association
#   run_info.txt
#
# Notes:
#   - NA cells are treated as "not present" in frequency counts and as 0 in
#     the PCA/heatmap; the per-marker n_missing column keeps this auditable.
#   - Classification: core = present in every non-missing sample (marker
#     missingness capped at 20%); private = present in exactly 1 sample;
#     absent = present in none; variable = everything else.
#   - Association: continuous trait -> Wilcoxon rank-sum; categorical trait
#     -> Fisher exact test (RxC). FDR = Benjamini-Hochberg.
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
pav_path <- opt(opts, "pav")
outdir <- opt(opts, "outdir", "results/pangenome/03_pav")
values_mode <- tolower(opt(opts, "values", "auto"))
cov_thr <- suppressWarnings(as.numeric(opt(opts, "cov-threshold", "0.5")))
meta_path <- opt(opts, "meta")
pheno_path <- opt(opts, "pheno")
trait <- opt(opts, "trait")
n_top <- suppressWarnings(as.integer(opt(opts, "n-top", "120")))

if (is.null(pav_path)) {
  stop("Required: --pav <marker x sample matrix.csv>.\n",
       "Hint: see the Usage block at the top of this script.", call. = FALSE)
}
if (is.null(pheno_path) != is.null(trait)) {
  stop("--pheno and --trait must be provided together.", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills PAV analysis ==")

# ---------------------------------------------------------------------------
# 1. Read + binarize
# ---------------------------------------------------------------------------
pav <- read_csv_safe(pav_path, "PAV matrix")
if (ncol(pav) < 3) stop("PAV matrix needs at least 2 samples (columns).", call. = FALSE)
marker_ids <- as.character(pav[[1]])
if (anyDuplicated(marker_ids)) marker_ids <- make.unique(marker_ids)
mat <- as.matrix(pav[, -1, drop = FALSE])
suppressWarnings(storage.mode(mat) <- "double")
rownames(mat) <- marker_ids

n_marker <- nrow(mat); n_sample <- ncol(mat)
n_na_cells <- sum(is.na(mat))
message(sprintf("  markers: %d | samples: %d | NA cells: %d (%.1f%%)",
                n_marker, n_sample, n_na_cells, 100 * n_na_cells / length(mat)))

# auto-detect value scale
if (values_mode == "auto") {
  vals <- mat[!is.na(mat)]
  values_mode <- if (all(vals %in% c(0, 1))) "binary" else "cov"
}
if (values_mode == "cov") {
  message(sprintf("  binarizing coverage values (threshold = %.2f)", cov_thr))
  bin <- matrix(NA_real_, n_marker, n_sample, dimnames = dimnames(mat))
  bin[!is.na(mat)] <- as.numeric(mat[!is.na(mat)] >= cov_thr)
} else {
  bin <- mat
  bin[!is.na(bin) & !(bin %in% c(0, 1))] <- NA_real_
}
message("  value mode: ", values_mode)

# ---------------------------------------------------------------------------
# 2. Classify markers
# ---------------------------------------------------------------------------
n_present <- rowSums(bin == 1, na.rm = TRUE)
n_missing <- rowSums(is.na(bin))
n_obs <- n_sample - n_missing
freq <- n_present / n_sample
# core: present in every NON-missing sample, with limited missingness
cls <- ifelse(n_present >= n_obs & n_missing <= 0.2 * n_sample, "core",
       ifelse(n_present == 1, "private",
       ifelse(n_present == 0, "absent", "variable")))

class_df <- data.frame(
  marker_id = marker_ids, class = cls,
  n_present = n_present, n_missing = n_missing, freq = round(freq, 4),
  stringsAsFactors = FALSE
)
utils::write.csv(class_df, file.path(outdir, "pav_class.csv"), row.names = FALSE)
utils::write.csv(data.frame(marker_id = marker_ids, bin, check.names = FALSE),
                 file.path(outdir, "pav_binary.csv"), row.names = FALSE)

tbl <- table(factor(cls, levels = c("core", "variable", "private", "absent")))
stats_lines <- c(
  "Rice-skills -- PAV summary",
  sprintf("markers            : %d", n_marker),
  sprintf("samples            : %d", n_sample),
  sprintf("core (100%% present) : %d (%.1f%%)", tbl[["core"]], 100 * tbl[["core"]] / n_marker),
  sprintf("variable           : %d (%.1f%%)", tbl[["variable"]], 100 * tbl[["variable"]] / n_marker),
  sprintf("private (1 sample) : %d (%.1f%%)", tbl[["private"]], 100 * tbl[["private"]] / n_marker),
  sprintf("absent             : %d (%.1f%%)", tbl[["absent"]], 100 * tbl[["absent"]] / n_marker)
)
writeLines(stats_lines, file.path(outdir, "pav_stats.txt"))
message(paste0("  ", stats_lines[-1], collapse = "\n  "))

# ---------------------------------------------------------------------------
# 3. Figures
# ---------------------------------------------------------------------------
pal <- rice_palette()
class_cols <- c(core = pal$groups5[1], variable = pal$groups5[2],
                private = pal$groups5[3], absent = "grey70")

# 3a. presence-frequency distribution
df_freq <- data.frame(freq = freq)
p_freq <- ggplot2::ggplot(df_freq, ggplot2::aes(x = freq)) +
  ggplot2::geom_histogram(bins = 50, fill = pal$groups5[1], colour = "white", linewidth = 0.2) +
  ggplot2::labs(x = "Presence frequency across samples", y = "Number of markers",
                title = "PAV frequency distribution") +
  theme_nature()
save_nature_figure(p_freq, file.path(outdir, "pav_frequency"), width_mm = 89, height_mm = 70)

# 3b. class composition bar
df_cls <- as.data.frame(tbl)
colnames(df_cls) <- c("class", "n")
df_cls$pct <- 100 * df_cls$n / sum(df_cls$n)
p_cls <- ggplot2::ggplot(df_cls, ggplot2::aes(x = class, y = n, fill = class)) +
  ggplot2::geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
                     vjust = -0.3, size = 2.2, lineheight = 0.9) +
  ggplot2::scale_fill_manual(values = class_cols) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
  ggplot2::labs(x = NULL, y = "Number of markers", title = "PAV class composition") +
  theme_nature() + ggplot2::theme(legend.position = "none")
save_nature_figure(p_cls, file.path(outdir, "pav_class_bar"), width_mm = 110, height_mm = 75)

# 3c. sample PCA on the PAV matrix
m_fill <- bin; m_fill[is.na(m_fill)] <- 0
if (n_marker >= 3 && n_sample >= 3) {
  pc <- stats::prcomp(t(m_fill), center = TRUE, scale. = FALSE)
  ve <- 100 * (pc$sdev^2 / sum(pc$sdev^2))
  df_pc <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], sample = rownames(pc$x))
  p_pca <- ggplot2::ggplot(df_pc, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 1.6, colour = pal$groups5[1], alpha = 0.85) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", ve[1]),
                  y = sprintf("PC2 (%.1f%%)", ve[2]),
                  title = "PAV-based sample structure") +
    theme_nature()
  save_nature_figure(p_pca, file.path(outdir, "pav_pca"), width_mm = 89, height_mm = 75)
}

# 3d. marker heatmap (core sample + most-variable + private sample)
sel <- integer(0)
core_i <- which(cls == "core"); var_i <- which(cls == "variable"); priv_i <- which(cls == "private")
if (length(core_i)) sel <- c(sel, core_i[seq_len(min(10, length(core_i)))])
if (length(var_i)) {
  vi <- var_i[order(abs(freq[var_i] - 0.5))]
  sel <- c(sel, vi[seq_len(min(n_top, length(vi)))])
}
if (length(priv_i)) sel <- c(sel, priv_i[seq_len(min(10, length(priv_i)))])
if (length(sel) >= 2) {
  m_sub <- bin[sel, , drop = FALSE]
  samp_ord <- if (n_sample >= 3) {
    hc <- stats::hclust(stats::dist(t(m_fill)))
    hc$labels[hc$order]
  } else colnames(m_sub)
  df_hm <- data.frame(
    marker = rep(rownames(m_sub), times = ncol(m_sub)),
    sample = rep(colnames(m_sub), each = nrow(m_sub)),
    present = as.vector(m_sub)
  )
  df_hm$marker <- factor(df_hm$marker, levels = rev(rownames(m_sub)))
  df_hm$sample <- factor(df_hm$sample, levels = samp_ord)
  p_hm <- ggplot2::ggplot(df_hm, ggplot2::aes(x = sample, y = marker, fill = factor(present))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c(`0` = "grey92", `1` = pal$groups5[1]), na.value = "white") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Present",
                  title = sprintf("PAV heatmap (%d markers)", nrow(m_sub))) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 4.5),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "right")
  save_nature_figure(p_hm, file.path(outdir, "pav_heatmap"), width_mm = 140, height_mm = 110)
}

# ---------------------------------------------------------------------------
# 4. Phenotype association (optional)
# ---------------------------------------------------------------------------
if (!is.null(pheno_path)) {
  ph <- read_pheno(pheno_path)
  if (!trait %in% colnames(ph)) {
    stop("Trait '", trait, "' not found in phenotype table. Columns: ",
         paste(colnames(ph), collapse = ", "), call. = FALSE)
  }
  shared <- intersect(colnames(bin), ph$sample_id)
  if (length(shared) < 5) {
    stop("Fewer than 5 shared samples between PAV matrix and phenotype table.",
         call. = FALSE)
  }
  message(sprintf("  association: %d shared samples for trait '%s'", length(shared), trait))
  y <- ph[[trait]][match(shared, ph$sample_id)]
  m_sub <- bin[, shared, drop = FALSE]
  y_num <- suppressWarnings(as.numeric(y))
  is_cont <- is.numeric(y) && length(unique(y[!is.na(y)])) > 5

  pvals <- rep(NA_real_, nrow(m_sub))
  for (i in seq_len(nrow(m_sub))) {
    present <- m_sub[i, ]
    if (all(is.na(present))) next
    ok <- !is.na(present) & !is.na(y); if (sum(ok) < 5) next
    pres <- present[ok] == 1; yy <- y[ok]
    if (!any(pres) || all(pres)) next
    pvals[i] <- tryCatch({
      if (is_cont) {
        stats::wilcox.test(yy[pres], yy[!pres])$p.value
      } else {
        tb <- table(factor(pres, levels = c(TRUE, FALSE)), droplevels(factor(yy)))
        if (any(colSums(tb) == 0)) {
          NA_real_
        } else if (min(dim(tb)) >= 2) {
          suppressWarnings(stats::fisher.test(tb)$p.value)
        } else {
          NA_real_
        }
      }
    }, error = function(e) NA_real_)
  }
  fdr <- stats::p.adjust(pvals, method = "BH")
  assoc <- data.frame(marker_id = marker_ids, class = cls, n_present = n_present,
                      p = pvals, fdr = fdr)
  utils::write.csv(assoc, file.path(outdir, "pav_association.csv"), row.names = FALSE)
  n_sig <- sum(assoc$fdr < 0.05, na.rm = TRUE)
  message(sprintf("  significant markers (FDR < 0.05): %d", n_sig))
  cat(sprintf("association trait: %s | significant (FDR<0.05): %d\n", trait, n_sig),
      file = file.path(outdir, "pav_stats.txt"), append = TRUE)

  # Manhattan if coordinates available
  if (!is.null(meta_path)) {
    info <- read_csv_safe(meta_path, "marker info table")
    check_columns(info, c("marker_id", "chr", "pos"), "marker info table")
    ap <- assoc[!is.na(assoc$p), c("marker_id", "p")]
    names(ap) <- c("snp_id", "p")
    info2 <- info[, c("marker_id", "chr", "pos")]
    names(info2) <- c("snp_id", "chr", "pos")
    keep <- ap$snp_id %in% info2$snp_id
    if (sum(keep) >= 10) {
      p_man <- manhattan_plot(ap[keep, ], info2, title = paste("PAV Manhattan —", trait))
      save_nature_figure(p_man, file.path(outdir, "pav_assoc_manhattan"),
                         width_mm = 180, height_mm = 80)
    } else {
      message("  (skipped Manhattan: fewer than 10 markers with coordinates + p-values)")
    }
  }
}

# ---------------------------------------------------------------------------
# 5. Run info
# ---------------------------------------------------------------------------
write_run_info_r(outdir, "03_pav_analysis.R", params = list(
  pav = pav_path, values_mode = values_mode, cov_threshold = cov_thr,
  n_markers = n_marker, n_samples = n_sample,
  n_core = as.integer(tbl[["core"]]), n_variable = as.integer(tbl[["variable"]]),
  n_private = as.integer(tbl[["private"]]), n_absent = as.integer(tbl[["absent"]]),
  phenotype = if (is.null(pheno_path)) NA_character_ else pheno_path,
  trait = if (is.null(trait)) NA_character_ else trait
))
message("Done. Outputs written to: ", normalizePath(outdir))
