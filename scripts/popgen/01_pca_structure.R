###############################################################################
# 01_pca_structure.R -- population structure: PCA + genomic kinship
#
# Title:    Rice-skills population-structure step
# Description:
#   Computes PCA on the filtered genotype matrix (imputed, scaled), saves the
#   scores for downstream GWAS covariates, computes the VanRaden genomic
#   relationship matrix, and draws the scree + PC1/PC2 scatter (coloured by a
#   group column when a phenotype table is supplied).
#
# Usage:
#   Rscript scripts/popgen/01_pca_structure.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --phenotype demo_data/genomics/phenotype.csv --group-col group \
#       --outdir results/popgen/01_pca
#
# Outputs:
#   pca_scores.csv / pca_variance.csv / kinship_matrix.csv
#   pca_scree.{pdf,svg,tiff} / pca_scatter.{pdf,svg,tiff}
#   fig_structure_overview.{pdf,svg,tiff} (patchwork, when installed)
#   run_info.txt
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
phenotype_path <- opt(opts, "phenotype")
group_col <- opt(opts, "group-col", "group")
outdir <- opt(opts, "outdir", "results/popgen/01_pca")
npc <- suppressWarnings(as.integer(opt(opts, "npc", "10")))
maf_cut <- suppressWarnings(as.numeric(opt(opts, "maf", "0.05")))
missing_cut <- suppressWarnings(as.numeric(opt(opts, "missing", "0.1")))

if (is.null(genotype_path) || is.null(snp_info_path)) {
  stop("Required: --genotype and --snp-info.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
fl <- filter_snps(gt$geno, snp_info, min_maf = maf_cut, max_missing = missing_cut)
message(sprintf("SNP filter (MAF >= %g, missing <= %g): %d -> %d SNPs.",
                maf_cut, missing_cut, ncol(gt$geno), ncol(fl$geno)))

pca_res <- compute_pca(fl$geno, npc = min(npc, nrow(fl$geno) - 1, ncol(fl$geno)))
scores <- data.frame(sample_id = rownames(pca_res$scores), pca_res$scores,
                     check.names = FALSE, row.names = NULL)

group <- NULL
if (!is.null(phenotype_path)) {
  pheno <- read_pheno(phenotype_path)
  if (group_col %in% colnames(pheno)) {
    m <- match(scores$sample_id, pheno$sample_id)
    group <- pheno[[group_col]][m]
    scores[[group_col]] <- group
  } else {
    message("Note: column '", group_col, "' not found in the phenotype table; ",
            "the scatter uses a single colour.")
  }
}
write.csv(scores, file.path(outdir, "pca_scores.csv"), row.names = FALSE)
write.csv(data.frame(PC = seq_along(pca_res$var_explained),
                     var_explained = round(pca_res$var_explained, 6)),
          file.path(outdir, "pca_variance.csv"), row.names = FALSE)
message(sprintf("PCA done: PC1 = %.1f%%, PC2 = %.1f%%.",
                100 * pca_res$var_explained[1], 100 * pca_res$var_explained[2]))

# --- kinship (VanRaden) -------------------------------------------------------
K <- vanraden_kinship(fl$geno)
write.csv(round(K, 5), file.path(outdir, "kinship_matrix.csv"))
message(sprintf("Genomic relationship matrix: %d x %d, mean off-diagonal %.4f.",
                nrow(K), ncol(K), mean(K[upper.tri(K)])))

# --- figures -------------------------------------------------------------------
# Scree: how much structure the leading axes carry (informs how many PCs to
# use as GWAS covariates).
scree_df <- data.frame(PC = factor(seq_along(pca_res$var_explained)),
                       pct = 100 * pca_res$var_explained)
p_scree <- ggplot2::ggplot(scree_df, ggplot2::aes(x = PC, y = pct)) +
  ggplot2::geom_col(fill = rice_palette()$series12[1], width = 0.7) +
  ggplot2::labs(x = "Principal component", y = "Variance explained (%)",
                title = "PCA scree") +
  theme_nature()
save_nature_figure(p_scree, file.path(outdir, "pca_scree"), width_mm = 89, height_mm = 70)

# Scatter: PC1 vs PC2; separated clusters indicate population structure that
# must be controlled in GWAS (the scatter is the visual diagnostic).
scores$group_plot <- if (!is.null(group)) group else "all"
p_scatter <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2)) +
  { if (!is.null(group)) {
      list(
        ggplot2::geom_point(ggplot2::aes(colour = group_plot), size = 1.4, alpha = 0.85),
        scale_colour_groups(name = group_col)
      )
    } else {
      ggplot2::geom_point(colour = rice_palette()$series12[1], size = 1.4, alpha = 0.85)
    } } +
  ggplot2::labs(x = sprintf("PC1 (%.1f%%)", 100 * pca_res$var_explained[1]),
                y = sprintf("PC2 (%.1f%%)", 100 * pca_res$var_explained[2]),
                title = "Population structure (PCA)") +
  theme_nature()
save_nature_figure(p_scatter, file.path(outdir, "pca_scatter"), width_mm = 110, height_mm = 89)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_ov <- patchwork::wrap_plots(p_scatter, p_scree, ncol = 2, widths = c(2, 1)) +
    patchwork::plot_annotation(title = "Population structure overview")
  save_nature_figure(p_ov, file.path(outdir, "fig_structure_overview"),
                     width_mm = 180, height_mm = 90)
} else {
  message("Note: patchwork not installed - overview figure skipped.")
}

write_run_info_r(
  outdir, "01_pca_structure.R",
  params = list(genotype = genotype_path, snp_info = snp_info_path,
                phenotype = if (is.null(phenotype_path)) NA else phenotype_path,
                group_col = group_col, npc = npc, maf = maf_cut,
                missing = missing_cut, seed = RICE_SEED),
  packages = c("ggplot2", "patchwork")
)
message("Population-structure step complete. Results in ", outdir)
