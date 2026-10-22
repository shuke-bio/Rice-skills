###############################################################################
# 06_phylogeny.R -- sample phylogeny from IBS distances (NJ / UPGMA)
#
# Title:    Rice-skills phylogeny step
# Description:
#   Builds a neighbor-joining (or UPGMA) tree from pairwise IBS distances
#   (mean fraction of differing alleles), writes the Newick tree, and draws
#   it with tips coloured by group (when a group column is supplied).
#
# Usage:
#   Rscript scripts/popgen/06_phylogeny.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --phenotype demo_data/genomics/phenotype.csv --group-col group \
#       --outdir results/popgen/06_phylogeny
#
# Outputs:
#   ibs_distance_matrix.csv / tree.nwk
#   phylogeny.{pdf,tiff,png}
#   run_info.txt
#
# Notes:
#   - NJ is unrooted; interpret tip relationships, not a root.
#   - For > ~max-snps SNPs the distance computation subsamples SNPs (seeded);
#     raise --max-snps to use all (slower, O(n^2 * m)).
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
outdir <- opt(opts, "outdir", "results/popgen/06_phylogeny")
method <- tolower(opt(opts, "method", "nj"))
max_snps <- suppressWarnings(as.integer(opt(opts, "max-snps", "20000")))
maf_cut <- suppressWarnings(as.numeric(opt(opts, "maf", "0.05")))

if (is.null(genotype_path) || is.null(snp_info_path)) {
  stop("Required: --genotype and --snp-info.", call. = FALSE)
}
if (!method %in% c("nj", "upgma")) {
  stop("--method must be nj or upgma.", call. = FALSE)
}
need_pkgs(c("ggplot2", "jsonlite", "ape"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
fl <- filter_snps(gt$geno, snp_info, min_maf = maf_cut, max_missing = 0.1)

g <- fl$geno
if (ncol(g) > max_snps) {
  sel <- sort(sample.int(ncol(g), max_snps))
  g <- g[, sel, drop = FALSE]
  message(sprintf("Distance computation uses %d of %d SNPs (--max-snps).", max_snps, ncol(g) + 0L))
}
D <- ibs_distance(g)
write.csv(as.matrix(D), file.path(outdir, "ibs_distance_matrix.csv"))
message(sprintf("IBS distances computed for %d samples (mean %.4f).",
                attr(D, "Size"), mean(as.matrix(D))))

tree <- if (method == "nj") ape::nj(D) else ape::as.phylo(stats::hclust(D, method = "average"))
ape::write.tree(tree, file.path(outdir, "tree.nwk"))

# tip colours by group when available
tip_labels <- tree$tip.label
tip_col <- rep("grey30", length(tip_labels))
if (!is.null(phenotype_path)) {
  pheno <- read_pheno(phenotype_path)
  if (group_col %in% colnames(pheno)) {
    lev <- unique(pheno[[group_col]])
    cols <- rice_palette()$groups5
    cols <- rep(cols, length.out = length(lev))
    lut <- setNames(cols, lev)
    tip_col <- unname(lut[pheno[[group_col]][match(tip_labels, pheno$sample_id)]])
    tip_col[is.na(tip_col)] <- "grey30"
  }
}

save_base_plot_genomics(
  function() {
    ape::plot.phylo(tree, type = "phylogram", cex = 0.35, no.margin = TRUE,
                    tip.color = tip_col, main = paste0("IBS phylogeny (", toupper(method), ")"))
  },
  file.path(outdir, "phylogeny"), width_mm = 140, height_mm = 140
)
message(sprintf("Phylogeny saved (%s, %d tips).", toupper(method), length(tip_labels)))

write_run_info_r(
  outdir, "06_phylogeny.R",
  params = list(genotype = genotype_path, snp_info = snp_info_path,
                phenotype = if (is.null(phenotype_path)) NA else phenotype_path,
                group_col = group_col, method = method, max_snps = max_snps,
                maf = maf_cut, n_tips = length(tip_labels), seed = RICE_SEED),
  packages = c("ape", "ggplot2")
)
message("Phylogeny step complete. Results in ", outdir)
