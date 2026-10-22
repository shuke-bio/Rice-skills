###############################################################################
# 07_haplotype.R -- haplotype analysis at a GWAS candidate locus
#
# Title:    Rice-skills haplotype analysis step
# Description:
#   Standard post-GWAS haplotype analysis at a candidate region:
#     mode "haplotype" (default): EM phase inference over the SNPs in the
#       region, samples grouped by their major haplotype, phenotype
#       differences across haplotype groups (Welch/ANOVA + Kruskal),
#       haplotype-phenotype boxplot and frequency barplot;
#     mode "snp": allele-dose grouping at the lead SNP (0/1/2 -> boxes),
#       the quickest allele-effect view.
#
# Usage:
#   Rscript scripts/gwas/07_haplotype.R \
#       --genotype demo_data/genomics/genotype.csv \
#       --snp-info demo_data/genomics/snp_info.csv \
#       --phenotype demo_data/genomics/phenotype.csv --trait salt_tolerance_index \
#       --region Chr01:25000000-26000000 --lead-snp s01704 \
#       --outdir results/gwas/07_haplotype
#
# Inputs:
#   --genotype    SNP x sample dosage matrix CSV (0/1/2)
#   --snp-info    snp_id, chr, pos
#   --phenotype   phenotype table (sample_id + traits)
#   --trait       trait column to test
#   --region      "chr:start-end" (or --chr/--start/--end); SNPs inside are used
#   --lead-snp    lead SNP id (required for mode=snp; recommended for haplotype
#                 ordering: groups are ordered by the lead-SNP minor-allele dose)
#   --mode        haplotype | snp (default haplotype)
#   --max-snp     maximum SNPs used for EM haplotypes (default 5; 2^m states)
#   --min-freq    haplotypes below this frequency are merged into "rare" (default 0.05)
#   --outdir      output directory (default results/gwas/07_haplotype)
# Outputs:
#   haplotypes.csv          haplotype, frequency, n_samples (major-haplotype call)
#   sample_haplotypes.csv   sample -> haplotype pair + group
#   group_stats.csv         per-group n / mean / sd of the trait + test p-values
#   haplotype_phenotype.*   boxplot of the trait by haplotype group (PDF/TIFF/PNG)
#   haplotype_frequency.*   frequency barplot
#   run_info.txt
#
# Notes:
#   - EM phase inference (Excoffier & Slatkin 1995) assumes Hardy-Weinberg
#     equilibrium and no missing genotypes inside the region (imputed input).
#   - Grouping rule: each sample is assigned to the more frequent haplotype
#     it carries (its "major" haplotype); double-heterozygotes therefore may
#     group with either parent haplotype — the standard approximation.
#   - Trait test: ANOVA when group means look normal-ish, Kruskal-Wallis is
#     reported alongside; pairwise Welch tests vs the largest group.
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
trait <- opt(opts, "trait")
region_str <- opt(opts, "region")
chr_arg <- opt(opts, "chr"); start_arg <- opt(opts, "start"); end_arg <- opt(opts, "end")
lead_snp <- opt(opts, "lead-snp")
mode <- tolower(opt(opts, "mode", "haplotype"))
max_snp <- suppressWarnings(as.integer(opt(opts, "max-snp", "5")))
min_freq <- suppressWarnings(as.numeric(opt(opts, "min-freq", "0.05")))
outdir <- opt(opts, "outdir", "results/gwas/07_haplotype")

if (any(vapply(list(genotype_path, snp_info_path, phenotype_path, trait), is.null, logical(1)))) {
  stop("Required: --genotype, --snp-info, --phenotype, --trait. See --help.", call. = FALSE)
}
if (is.null(region_str) && any(vapply(list(chr_arg, start_arg, end_arg), is.null, logical(1)))) {
  stop("Provide --region chr:start-end (or --chr/--start/--end).", call. = FALSE)
}
if (mode == "snp" && is.null(lead_snp)) stop("mode=snp requires --lead-snp.", call. = FALSE)
need_pkgs(c("ggplot2", "jsonlite"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
message("== Rice-skills haplotype analysis (", mode, " mode) ==")

# region parsing
if (!is.null(region_str)) {
  m <- regmatches(region_str, regexec("^(.+):(\\d+)-(\\d+)$", region_str))[[1]]
  if (length(m) != 4) stop("Cannot parse --region '", region_str,
                           "'. Expected chr:start-end.", call. = FALSE)
  chr_r <- m[2]; start_r <- as.numeric(m[3]); end_r <- as.numeric(m[4])
} else {
  chr_r <- chr_arg; start_r <- as.numeric(start_arg); end_r <- as.numeric(end_arg)
}

# ---------------------------------------------------------------------------
# 1) Load + subset region
# ---------------------------------------------------------------------------
gt <- read_genotype(genotype_path)
snp_info <- read_snp_info(snp_info_path)
G <- gt$geno  # samples x SNPs (0/1/2)
info <- snp_info[match(colnames(G), snp_info$snp_id), , drop = FALSE]
in_region <- info$chr == chr_r & info$pos >= start_r & info$pos <= end_r
message(sprintf("  region %s:%d-%d contains %d SNPs (of %d)", chr_r, start_r, end_r,
                sum(in_region), ncol(G)))
if (sum(in_region) < 1) stop("No SNPs inside the region. Check coordinates/chr naming.", call. = FALSE)

G_r <- G[, in_region, drop = FALSE]
info_r <- info[in_region, , drop = FALSE]

pheno <- read_pheno(phenotype_path)
if (!trait %in% colnames(pheno)) {
  stop("Trait '", trait, "' not found. Available: ",
       paste(setdiff(colnames(pheno), "sample_id"), collapse = ", "), call. = FALSE)
}
shared <- intersect(rownames(G_r), pheno$sample_id)
if (length(shared) < 6) stop("Fewer than 6 shared samples.", call. = FALSE)
G_r <- G_r[shared, , drop = FALSE]
y <- pheno[[trait]][match(shared, pheno$sample_id)]
y_num <- suppressWarnings(as.numeric(y))
message(sprintf("  samples: %d | trait '%s' range [%.2f, %.2f]",
                length(shared), trait, min(y_num, na.rm = TRUE), max(y_num, na.rm = TRUE)))
if (anyNA(G_r)) {
  message(sprintf("  imputing %d missing genotype cells to SNP means", sum(is.na(G_r))))
  G_r <- impute_mean(G_r)
  G_r <- round(G_r)
}

# ---------------------------------------------------------------------------
# 2a) SNP mode: lead-SNP allele dose grouping
# ---------------------------------------------------------------------------
if (mode == "snp") {
  if (!lead_snp %in% colnames(G_r)) stop("Lead SNP '", lead_snp, "' not in region/genotype.", call. = FALSE)
  dose <- G_r[, lead_snp]
  group <- factor(c("aa", "Aa", "AA")[dose + 1], levels = c("aa", "Aa", "AA"))
  groups <- data.frame(sample_id = shared, group = as.character(group), stringsAsFactors = FALSE)
  group_label <- "lead SNP genotype"
} else {
  # -----------------------------------------------------------------------
  # 2b) haplotype mode: EM phase inference
  # -----------------------------------------------------------------------
  if (ncol(G_r) > max_snp) {
    # keep the most informative SNPs: highest MAF spread around the lead
    maf <- colMeans(G_r) / 2
    maf <- pmin(maf, 1 - maf)
    if (!is.null(lead_snp) && lead_snp %in% colnames(G_r)) {
      # lead SNP first, then by MAF (most informative)
      ord <- c(match(lead_snp, colnames(G_r)), setdiff(order(-maf), match(lead_snp, colnames(G_r))))
    } else {
      ord <- order(-maf)
    }
    sel <- head(ord, max_snp)
    G_r <- G_r[, sort(sel), drop = FALSE]
    info_r <- info_r[sort(sel), , drop = FALSE]
    message(sprintf("  using %d SNPs for EM (max-snp=%d)", ncol(G_r), max_snp))
  }
  m <- ncol(G_r)
  n <- nrow(G_r)
  # enumerate haplotype states
  states <- as.matrix(expand.grid(rep(list(0:1), m)))
  colnames(states) <- colnames(G_r)
  n_states <- nrow(states)
  # compatible state pairs per sample
  pairs <- vector("list", n)
  for (i in seq_len(n)) {
    g <- G_r[i, ]
    ok <- which(apply(states, 1, function(h1) {
      any(apply(states, 1, function(h2) all(h1 + h2 == g)))
    }))
    comb <- expand.grid(s1 = ok, s2 = ok)
    comb <- comb[comb$s1 <= comb$s2, , drop = FALSE]
    keep <- apply(comb, 1, function(r) all(states[r[1], ] + states[r[2], ] == g))
    pairs[[i]] <- comb[keep, , drop = FALSE]
  }
  # EM
  f <- rep(1 / n_states, n_states)
  for (it in 1:200) {
    f_new <- rep(0, n_states)
    for (i in seq_len(n)) {
      pr <- f[pairs[[i]]$s1] * f[pairs[[i]]$s2]
      if (sum(pr) == 0) next
      pr <- pr / sum(pr)
      f_new[pairs[[i]]$s1] <- f_new[pairs[[i]]$s1] + pr
      f_new[pairs[[i]]$s2] <- f_new[pairs[[i]]$s2] + pr
    }
    f_new <- f_new / (2 * n)
    if (max(abs(f_new - f)) < 1e-8) { f <- f_new; break }
    f <- f_new
  }
  names(f) <- apply(states, 1, paste, collapse = "-")
  # assign each sample its most likely pair -> major haplotype
  group_chr <- character(n)
  hap1 <- character(n); hap2 <- character(n)
  for (i in seq_len(n)) {
    pr <- f[pairs[[i]]$s1] * f[pairs[[i]]$s2]
    k <- which.max(pr)
    h1 <- names(f)[pairs[[i]]$s1[k]]; h2 <- names(f)[pairs[[i]]$s2[k]]
    hap1[i] <- h1; hap2[i] <- h2
    freqs <- c(f[h1], f[h2])
    group_chr[i] <- h1[which.max(freqs)]
  }
  # merge rare haplotypes
  tab <- table(group_chr)
  rare <- names(tab)[tab / n < min_freq]
  group_chr[group_chr %in% rare] <- "rare"
  group <- factor(group_chr)
  groups <- data.frame(sample_id = shared, hap1 = hap1, hap2 = hap2,
                       group = as.character(group), stringsAsFactors = FALSE)
  group_label <- "major haplotype"
  utils::write.csv(data.frame(haplotype = names(f), frequency = round(f, 4),
                              in_rare_pool = names(f) %in% rare),
                   file.path(outdir, "haplotype_states.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 3) Phenotype contrast across groups
# ---------------------------------------------------------------------------
groups$y <- y_num[match(groups$sample_id, shared)]
groups <- groups[!is.na(groups$y), , drop = FALSE]
valid_groups <- names(table(groups$group))[table(groups$group) >= 3]
groups <- groups[groups$group %in% valid_groups, , drop = FALSE]
groups$group <- droplevels(factor(groups$group))
utils::write.csv(groups, file.path(outdir, "sample_haplotypes.csv"), row.names = FALSE)

gstats <- do.call(rbind, lapply(split(groups$y, groups$group), function(v) {
  data.frame(group = names(table(groups$group))[1], n = length(v),
             mean = mean(v), sd = stats::sd(v))
}))
gstats <- do.call(rbind, lapply(split(groups, groups$group), function(d) {
  data.frame(group = d$group[1], n = nrow(d), mean = mean(d$y), sd = stats::sd(d$y))
}))
rownames(gstats) <- NULL

anova_p <- NA; kw_p <- NA
if (nlevels(groups$group) >= 2) {
  anova_p <- tryCatch(summary(stats::aov(y ~ group, data = groups))[[1]][1, "Pr(>F)"],
                      error = function(e) NA_real_)
  kw_p <- tryCatch(stats::kruskal.test(y ~ group, data = groups)$p.value,
                   error = function(e) NA_real_)
}
gstats$anova_p <- anova_p
gstats$kruskal_p <- kw_p
utils::write.csv(gstats, file.path(outdir, "group_stats.csv"), row.names = FALSE)
message(sprintf("  groups: %s | ANOVA p = %s | Kruskal p = %s",
                paste(levels(groups$group), collapse = ", "),
                format(anova_p, digits = 3), format(kw_p, digits = 3)))

# ---------------------------------------------------------------------------
# 4) Figures
# ---------------------------------------------------------------------------
pal <- rice_palette()
cols <- setNames(pal$groups5[seq_len(nlevels(groups$group))], levels(groups$group))
p_box <- ggplot2::ggplot(groups, ggplot2::aes(x = group, y = y, fill = group)) +
  ggplot2::geom_boxplot(outlier.size = 0.6, linewidth = 0.3, colour = "grey25", alpha = 0.9) +
  ggplot2::geom_jitter(width = 0.12, size = 0.7, alpha = 0.5, colour = "grey30") +
  ggplot2::scale_fill_manual(values = cols) +
  ggplot2::labs(x = group_label, y = trait,
                title = sprintf("Phenotype by %s (%s)", group_label,
                                if (!is.na(anova_p)) sprintf("ANOVA p = %.2g", anova_p) else "")) +
  theme_nature() + ggplot2::theme(legend.position = "none")
save_nature_figure(p_box, file.path(outdir, "haplotype_phenotype"), width_mm = 89, height_mm = 75)

if (mode == "haplotype") {
  ft <- sort(table(groups$group), decreasing = TRUE)
  df_f <- data.frame(group = names(ft), n = as.integer(ft))
  p_f <- ggplot2::ggplot(df_f, ggplot2::aes(x = reorder(group, -n), y = n)) +
    ggplot2::geom_col(fill = pal$groups5[2], colour = "grey25", width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4, size = 2.4) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = group_label, y = "Number of samples",
                  title = "Haplotype group sizes") +
    theme_nature()
  save_nature_figure(p_f, file.path(outdir, "haplotype_frequency"), width_mm = 89, height_mm = 65)
}

write_run_info_r(outdir, "07_haplotype.R", params = list(
  genotype = genotype_path, region = region_str %||% sprintf("%s:%d-%d", chr_r, start_r, end_r),
  lead_snp = if (is.null(lead_snp)) NA_character_ else lead_snp,
  mode = mode, trait = trait, n_samples = length(shared),
  anova_p = anova_p, kw_p = kw_p
))
message("Done. Outputs written to: ", normalizePath(outdir))
