###############################################################################
# 04_candidate_loci.R -- significant-locus calling and candidate genes
#
# Title:    Rice-skills GWAS candidate-locus step
# Description:
#   Clusters significant SNPs into loci (per chromosome, within --flank bp of
#   each other), reports the lead SNP per locus, and - when a gene annotation
#   file is supplied - overlaps the locus windows with gene models to produce
#   a candidate-gene table.
#
# Usage:
#   Rscript scripts/gwas/04_candidate_loci.R \
#       --assoc results/gwas/03_association/association_results.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --gff demo_data/genomics/annotation.gff3 \
#       --outdir results/gwas/04_candidates
#
# Inputs:
#   --assoc      association_results.csv from 03_association.R
#   --snp-info   snp_id, chr, pos
#   --gff        optional GFF3/GTF gene annotation (MSU / RAP-DB); if omitted,
#                only the loci table is produced (with a download hint)
#   --p-threshold  significance cutoff (default Bonferroni 0.05 / n)
#   --flank      bp added on each side of the significant-SNP span (default 1e5)
# Outputs:
#   gwas_loci.csv            one row per locus (lead SNP, span, min p)
#   significant_snps.csv     significant SNPs annotated with chr/pos
#   candidate_genes.csv      locus x gene overlaps (only with --gff)
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
assoc_path <- opt(opts, "assoc")
snp_info_path <- opt(opts, "snp-info")
gff_path <- opt(opts, "gff")
outdir <- opt(opts, "outdir", "results/gwas/04_candidates")
p_threshold <- opt(opts, "p-threshold", NULL)
flank <- suppressWarnings(as.numeric(opt(opts, "flank", "100000")))

if (is.null(assoc_path) || is.null(snp_info_path)) {
  stop("Required: --assoc and --snp-info.\nHint: run scripts/gwas/03_association.R first.",
       call. = FALSE)
}
if (is.na(flank) || flank < 0) stop("--flank must be a non-negative number of bp.", call. = FALSE)

need_pkgs(c("ggplot2", "jsonlite"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

assoc <- read_csv_safe(assoc_path, "Association results")
check_columns(assoc, c("snp_id", "p"), "Association results")
snp_info <- read_snp_info(snp_info_path)
assoc <- merge(assoc, snp_info, by = "snp_id", all.x = TRUE)

if (is.null(p_threshold)) {
  p_threshold <- 0.05 / nrow(assoc)
} else {
  p_threshold <- as.numeric(p_threshold)
}
sig <- assoc[!is.na(assoc$p) & assoc$p <= p_threshold & !is.na(assoc$chr), , drop = FALSE]
sig <- sig[order(sig$chr, sig$pos), , drop = FALSE]
message(sprintf("Significant SNPs (p <= %.3g): %d.", p_threshold, nrow(sig)))

# --- cluster significant SNPs into loci --------------------------------------
loc_rows <- list()
for (chrom in unique(sig$chr)) {
  d <- sig[sig$chr == chrom, , drop = FALSE]
  if (nrow(d) == 0) next
  cluster <- cumsum(c(TRUE, diff(d$pos) > 2 * flank))
  for (cl in unique(cluster)) {
    dd <- d[cluster == cl, , drop = FALSE]
    lead <- dd[which.min(dd$p), , drop = FALSE]
    loc_rows[[length(loc_rows) + 1L]] <- data.frame(
      locus_id = sprintf("%s_%d-%d", chrom, min(dd$pos), max(dd$pos)),
      chr = chrom,
      start = max(0, min(dd$pos) - flank),
      end = max(dd$pos) + flank,
      n_sig_snps = nrow(dd),
      lead_snp = lead$snp_id,
      lead_pos = lead$pos,
      lead_p = lead$p,
      stringsAsFactors = FALSE
    )
  }
}
loci <- if (length(loc_rows) > 0) do.call(rbind, loc_rows) else
  data.frame(locus_id = character(0), chr = character(0), start = numeric(0),
             end = numeric(0), n_sig_snps = integer(0), lead_snp = character(0),
             lead_pos = numeric(0), lead_p = numeric(0), stringsAsFactors = FALSE)
loci <- loci[order(loci$chr, loci$start), , drop = FALSE]
rownames(loci) <- NULL

utils::write.csv(sig, file.path(outdir, "significant_snps.csv"),
                 row.names = FALSE, quote = FALSE)
utils::write.csv(loci, file.path(outdir, "gwas_loci.csv"),
                 row.names = FALSE, quote = FALSE)
message(sprintf("Loci called: %d (flank %s bp).", nrow(loci), format(flank, scientific = FALSE)))

# --- candidate genes ----------------------------------------------------------
if (!is.null(gff_path)) {
  genes <- read_gff_genes(gff_path)
  message(sprintf("Annotation loaded: %d gene models.", nrow(genes)))
  cand <- if (nrow(loci) > 0) overlap_genes(loci, genes) else
    data.frame(locus_id = character(0), chr = character(0), start = numeric(0),
               end = numeric(0), gene_id = character(0), gene_name = character(0),
               stringsAsFactors = FALSE)
  # annotate each candidate gene with the lead SNP p-value of its locus
  if (nrow(cand) > 0) {
    cand <- merge(cand, loci[, c("locus_id", "lead_snp", "lead_p")], by = "locus_id")
    cand <- cand[order(cand$lead_p), , drop = FALSE]
  }
  utils::write.csv(cand, file.path(outdir, "candidate_genes.csv"),
                   row.names = FALSE, quote = FALSE)
  message(sprintf("Candidate gene overlaps: %d (see candidate_genes.csv).", nrow(cand)))
} else {
  message("Note: no --gff supplied - skipping candidate-gene overlap.\n",
          "Hint: download the MSU/RAP-DB GFF3 (RGAP or RAP-DB download pages) ",
          "and rerun with --gff <annotation.gff3>.")
}

write_run_info_r(
  outdir, "04_candidate_loci.R",
  params = list(assoc = assoc_path, snp_info = snp_info_path, gff = if (is.null(gff_path)) NA else gff_path,
                p_threshold = p_threshold, flank = flank,
                n_sig = nrow(sig), n_loci = nrow(loci), seed = RICE_SEED),
  packages = c("ggplot2")
)
message("Candidate-locus step complete. Results in ", outdir)
