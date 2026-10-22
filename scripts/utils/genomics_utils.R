###############################################################################
# genomics_utils.R -- shared helpers for the Rice-skills genomics modules
#                     (scripts/gwas/ and scripts/popgen/)
#
# Title:    Rice genomics shared utilities (genotype I/O, QC, PCA, kinship,
#           association statistics, selection statistics, plotting)
# Description:
#   Self-contained base-R implementations of the core population-genomics and
#   GWAS computations, so the genomics modules run without any external
#   binaries (PLINK / ADMIXTURE / vcftools remain supported as optional
#   enhanced engines via their own script paths).
#
#   Conventions:
#     - Genotype matrix: samples x SNPs, dosages 0/1/2 (alt-allele count),
#       missing values as NA. Missing genotypes are imputed to the per-SNP
#       mean dosage (2*MAF) for PCA / association / distances.
#     - SNP info table: columns snp_id, chr, pos (bp).
#     - Phenotype table: first column sample_id, remaining columns traits.
#
# Usage:
#   source("scripts/utils/plot_theme.R"); source("scripts/utils/genomics_utils.R")
#
# Notes:
#   - Association uses per-SNP OLS (closed-form, fully vectorised); it is an
#     honest baseline ("GLM" model family). For MLM/FarmCPU-grade analyses
#     install rMVP or GAPIT (see scripts/gwas/03_association.R).
#   - Fst uses the Hudson estimator (ratio of averages); pi is computed per
#     site as 2p(1-p)*n/(n-1) and aggregated per window.
###############################################################################

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

#' Read a genotype matrix CSV (first column = snp_id, remaining columns samples).
#'
#' @param path CSV path.
#' @param id_col Name of the first column.
#' @return List with `geno` (samples x SNPs numeric matrix), `samples`, `snps`.
read_genotype <- function(path, id_col = "snp_id") {
  df <- read_csv_safe(path, "Genotype matrix")
  if (colnames(df)[1] != id_col) {
    colnames(df)[1] <- id_col
  }
  snps <- as.character(df[[id_col]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  geno <- t(mat) # samples x SNPs
  rownames(geno) <- colnames(df)[-1]
  colnames(geno) <- snps
  bad <- !is.na(geno) & !(geno %in% c(0, 1, 2))
  if (any(bad)) {
    stop(sum(bad), " genotype value(s) outside {0,1,2} were found. ",
         "Hint: expected alt-allele dosages (0/1/2); check the input matrix.", call. = FALSE)
  }
  list(geno = geno, samples = rownames(geno), snps = snps)
}

#' Read the SNP info table (snp_id, chr, pos).
read_snp_info <- function(path) {
  info <- read_csv_safe(path, "SNP info table")
  check_columns(info, c("snp_id", "chr", "pos"), "SNP info table")
  info$chr <- as.character(info$chr)
  info$pos <- as.numeric(info$pos)
  info
}

#' Read a phenotype table (sample_id + trait columns).
read_pheno <- function(path) {
  ph <- read_csv_safe(path, "Phenotype table")
  colnames(ph)[1] <- "sample_id"
  ph$sample_id <- as.character(ph$sample_id)
  ph
}

#' Align genotype rows (samples) with a phenotype table; returns shared samples.
align_samples_geno <- function(geno, pheno, id_col = "sample_id") {
  shared <- intersect(rownames(geno), pheno[[id_col]])
  if (length(shared) < 3) {
    stop("Fewer than 3 shared samples between genotype and phenotype tables. ",
         "Hint: sample IDs must match exactly (case-sensitive).", call. = FALSE)
  }
  list(geno = geno[shared, , drop = FALSE], pheno = pheno[match(shared, pheno[[id_col]]), , drop = FALSE])
}

# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

#' Filter SNPs by minor allele frequency and missingness.
#'
#' @param geno samples x SNPs dosage matrix.
#' @param snp_info matching SNP info table (same order as columns).
#' @param min_maf Minimum minor allele frequency (default 0.05).
#' @param max_missing Maximum per-SNP missing rate (default 0.1).
#' @return List with filtered geno / snp_info and a statistics data.frame.
filter_snps <- function(geno, snp_info, min_maf = 0.05, max_missing = 0.1) {
  n <- nrow(geno)
  af <- colMeans(geno, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  missing <- colMeans(is.na(geno))
  keep <- maf >= min_maf & missing <= max_missing
  stats_df <- data.frame(
    snp_id = colnames(geno), maf = round(maf, 4), missing = round(missing, 4),
    kept = keep, stringsAsFactors = FALSE
  )
  list(
    geno = geno[, keep, drop = FALSE],
    snp_info = snp_info[match(colnames(geno)[keep], snp_info$snp_id), , drop = FALSE],
    stats = stats_df
  )
}

#' Impute missing genotypes to the per-SNP mean dosage (2*MAF).
impute_mean <- function(geno) {
  for (j in seq_len(ncol(geno))) {
    miss <- is.na(geno[, j])
    if (any(miss)) {
      geno[miss, j] <- mean(geno[, j], na.rm = TRUE)
    }
  }
  geno
}

# ---------------------------------------------------------------------------
# Population structure
# ---------------------------------------------------------------------------

#' PCA on the imputed, scaled genotype matrix.
#'
#' @return List with scores (samples x PCs), variance explained, loadings.
compute_pca <- function(geno, npc = 10) {
  g <- impute_mean(geno)
  g <- scale(g) # centre + unit variance per SNP
  pca <- stats::prcomp(g, center = FALSE, scale. = FALSE, rank. = npc)
  ve <- (pca$sdev^2) / sum(pca$sdev^2)
  list(scores = pca$x, var_explained = ve, rotation = pca$rotation)
}

#' VanRaden genomic relationship matrix (G matrix).
vanraden_kinship <- function(geno) {
  g <- impute_mean(geno)
  p <- colMeans(g) / 2
  Z <- sweep(g, 2, 2 * p, "-")
  denom <- 2 * sum(p * (1 - p))
  K <- tcrossprod(Z) / denom
  rownames(K) <- rownames(g)
  colnames(K) <- rownames(g)
  K
}

# ---------------------------------------------------------------------------
# Association (per-SNP OLS, closed form, vectorised)
# ---------------------------------------------------------------------------

#' Simple GWAS: per-SNP linear regression of a trait on dosage + covariates.
#'
#' Implements the classic "GLM" model family: OLS of y on each SNP's dosage
#' with covariates (e.g. PCs) residualised out of both y and the dosage.
#' Vectorised across SNPs; suitable as an exploratory baseline.
#'
#' @param geno samples x SNPs dosage matrix (imputed internally).
#' @param y numeric trait vector (length = nrow(geno)).
#' @param covariates matrix/data.frame of covariates (e.g. PCs); NULL for none.
#' @return data.frame: snp_id, beta, se, t, p, n.
gwas_glm <- function(geno, y, covariates = NULL) {
  g <- impute_mean(geno)
  n <- nrow(g)
  if (length(y) != n) stop("Trait length does not match genotype row count.", call. = FALSE)
  if (anyNA(y)) {
    keep <- !is.na(y)
    g <- g[keep, , drop = FALSE]
    y <- y[keep]
    if (!is.null(covariates)) covariates <- covariates[keep, , drop = FALSE]
    n <- length(y)
  }
  if (is.null(covariates)) {
    y_res <- y - mean(y)
    g_res <- sweep(g, 2, colMeans(g), "-")
    df_resid <- n - 1L
  } else {
    C <- as.matrix(covariates)
    C <- cbind(1, C) # intercept
    qrC <- qr(C)
    y_res <- qr.resid(qrC, y)
    g_res <- qr.resid(qrC, g)
    df_resid <- n - ncol(C)
  }
  sxx <- colSums(g_res^2)
  sxy <- colSums(g_res * y_res)
  beta <- sxy / sxx
  resid <- y_res - g_res * rep(beta, each = n)
  sigma2 <- colSums(resid^2) / df_resid
  se <- sqrt(sigma2 / sxx)
  tstat <- beta / se
  p <- 2 * stats::pt(-abs(tstat), df = df_resid)
  data.frame(
    snp_id = colnames(g), beta = beta, se = se, t = tstat, p = p, n = n,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Selection statistics
# ---------------------------------------------------------------------------

#' Hudson Fst (ratio of averages) between two sample groups, per SNP.
#'
#' @param geno samples x SNPs dosage matrix.
#' @param group1,group2 character vectors of sample IDs.
#' @return data.frame: snp_id, fst_num, fst_den (per-SNP numerator/denominator).
hudson_fst <- function(geno, group1, group2) {
  g1 <- geno[group1, , drop = FALSE]
  g2 <- geno[group2, , drop = FALSE]
  p1 <- colMeans(g1, na.rm = TRUE) / 2
  p2 <- colMeans(g2, na.rm = TRUE) / 2
  n1 <- 2 * colSums(!is.na(g1))
  n2 <- 2 * colSums(!is.na(g2))
  num <- (p1 - p2)^2 - p1 * (1 - p1) / pmax(n1 - 1, 1) - p2 * (1 - p2) / pmax(n2 - 1, 1)
  den <- p1 * (1 - p2) + p2 * (1 - p1)
  data.frame(snp_id = colnames(geno), fst_num = num, fst_den = den,
             p1 = p1, p2 = p2, stringsAsFactors = FALSE)
}

#' Per-SNP nucleotide diversity pi = 2p(1-p) * n/(n-1); NA when n < 2.
site_pi <- function(geno) {
  n <- colSums(!is.na(geno))
  p <- colMeans(geno, na.rm = TRUE) / 2
  out <- 2 * p * (1 - p) * n / pmax(n - 1, 1)
  out[n < 2] <- NA_real_
  out
}

#' Aggregate per-SNP values into genomic windows (mean per window).
#'
#' @param snp_info SNP info table (aligned with `values`).
#' @param values numeric vector (e.g. per-SNP pi).
#' @param window_bp Window size in bp (default 100 kb).
#' @return data.frame: chr, start, end, mid, n_snp, value (mean).
window_summarise <- function(snp_info, values, window_bp = 1e5) {
  df <- data.frame(chr = snp_info$chr, pos = snp_info$pos, value = values,
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$value), , drop = FALSE]
  df$window <- floor(df$pos / window_bp)
  out <- stats::aggregate(value ~ chr + window, data = df, FUN = mean)
  counts <- stats::aggregate(value ~ chr + window, data = df, FUN = length)
  out$chr <- as.character(out$chr)
  out$n_snp <- counts$value
  out$start <- out$window * window_bp
  out$end <- out$start + window_bp - 1
  out$mid <- (out$start + out$end) / 2
  out <- out[order(out$chr, out$start), , drop = FALSE]
  rownames(out) <- NULL
  out[, c("chr", "start", "end", "mid", "n_snp", "value")]
}

# ---------------------------------------------------------------------------
# Linkage disequilibrium
# ---------------------------------------------------------------------------

#' LD decay curve: r^2 between SNP pairs binned by physical distance.
#'
#' @param geno samples x SNPs dosage matrix.
#' @param snp_info matching SNP info table.
#' @param max_dist_bp Maximum distance (default 500 kb).
#' @param max_pairs Cap on evaluated pairs (random subsample beyond it).
#' @return data.frame with distance bins and mean r^2.
ld_decay <- function(geno, snp_info, max_dist_bp = 5e5, max_pairs = 2e5) {
  g <- impute_mean(geno)
  chr <- snp_info$chr[match(colnames(g), snp_info$snp_id)]
  pos <- snp_info$pos[match(colnames(g), snp_info$snp_id)]
  n <- nrow(g)
  # standardise once; r^2 of a pair = (mean(x*y))^2 with centred/scaled columns
  gs <- scale(g)
  dists <- numeric(0)
  r2s <- numeric(0)
  for (chrom in unique(chr)) {
    idx <- which(chr == chrom)
    if (length(idx) < 2) next
    idx <- idx[order(pos[idx])]
    p <- pos[idx]
    # candidate pairs within max_dist
    for (a in seq_along(idx)) {
      b_max <- which(p - p[a] > 0 & p - p[a] <= max_dist_bp)
      if (length(b_max) == 0) next
      bsel <- idx[b_max]
      aa <- idx[a]
      r <- colMeans(gs[, aa] * gs[, bsel, drop = FALSE]) # cor for scaled columns
      dists <- c(dists, p[b_max] - p[a])
      r2s <- c(r2s, r^2)
    }
  }
  if (length(dists) == 0) {
    return(data.frame(bin = numeric(0), dist = numeric(0), mean_r2 = numeric(0),
                      n_pairs = integer(0)))
  }
  if (length(dists) > max_pairs) {
    sel <- sort(sample.int(length(dists), max_pairs))
    dists <- dists[sel]
    r2s <- r2s[sel]
  }
  bins <- seq(0, max_dist_bp, by = 2.5e4)
  bin_id <- cut(dists, breaks = bins, include.lowest = TRUE, labels = FALSE)
  out <- data.frame(
    bin = seq_along(bins[-1]),
    dist = bins[-1] - 1.25e4,
    mean_r2 = tapply(r2s, bin_id, mean),
    n_pairs = tapply(r2s, bin_id, length)
  )
  out <- out[!is.na(out$mean_r2), , drop = FALSE]
  out$mean_r2 <- round(out$mean_r2, 5)
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Phylogeny
# ---------------------------------------------------------------------------

#' IBS distance matrix between samples (mean fraction of differing alleles).
ibs_distance <- function(geno) {
  g <- impute_mean(geno)
  n <- nrow(g)
  D <- matrix(0, n, n, dimnames = list(rownames(g), rownames(g)))
  for (i in seq_len(n)) {
    diff <- abs(sweep(g, 2, g[i, ], "-"))
    D[i, ] <- rowMeans(diff) / 2
  }
  stats::as.dist(D)
}

# ---------------------------------------------------------------------------
# Annotation (GFF3 gene models)
# ---------------------------------------------------------------------------

#' Read gene records from a GFF3/GTF file (only 'gene' features).
#'
#' @param path GFF3/GTF path.
#' @return data.frame: chr, start, end, gene_id, gene_name.
read_gff_genes <- function(path) {
  if (!file.exists(path)) {
    stop("Annotation file not found: ", path,
         "\nHint: download the MSU/RAP-DB GFF3 (e.g. from the RGAP or RAP-DB ",
         "download pages) or skip candidate annotation.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[vapply(fields, length, integer(1)) >= 9]
  is_gene <- vapply(fields, function(f) f[3] %in% c("gene", "mRNA"), logical(1))
  # Prefer 'gene' features; fall back to 'mRNA' when no gene lines exist.
  want <- if (any(vapply(fields, function(f) f[3] == "gene", logical(1)))) "gene" else "mRNA"
  fields <- fields[vapply(fields, function(f) f[3] == want, logical(1))]
  out <- data.frame(
    chr = vapply(fields, `[`, character(1), 1),
    start = as.numeric(vapply(fields, `[`, character(1), 4)),
    end = as.numeric(vapply(fields, `[`, character(1), 5)),
    stringsAsFactors = FALSE
  )
  attrs <- vapply(fields, `[`, character(1), 9)
  out$gene_id <- sub(".*(?:ID|gene_id)=([^;]+).*", "\\1", attrs)
  out$gene_name <- ifelse(grepl("Name=", attrs),
                          sub(".*Name=([^;]+).*", "\\1", attrs), NA_character_)
  out
}

#' Overlap candidate windows with gene models.
#'
#' @param loci data.frame with chr, start, end (candidate windows).
#' @param genes gene table from read_gff_genes().
#' @return data.frame: one row per window-gene overlap.
overlap_genes <- function(loci, genes) {
  hits <- list()
  for (i in seq_len(nrow(loci))) {
    g <- genes[genes$chr == loci$chr[i] &
                 genes$start <= loci$end[i] & genes$end >= loci$start[i], , drop = FALSE]
    if (nrow(g) > 0) {
      hits[[length(hits) + 1L]] <- data.frame(
        locus_id = loci$locus_id[i], chr = loci$chr[i],
        start = loci$start[i], end = loci$end[i],
        gene_id = g$gene_id, gene_name = g$gene_name,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(hits) == 0) {
    return(data.frame(locus_id = character(0), chr = character(0),
                      start = numeric(0), end = numeric(0),
                      gene_id = character(0), gene_name = character(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, hits)
}

# ---------------------------------------------------------------------------
# Plotting helpers (Nature theme via plot_theme.R)
# ---------------------------------------------------------------------------

#' Manhattan plot (ggplot2, Nature theme).
#'
#' @param assoc result table with snp_id, p (and optional beta).
#' @param snp_info SNP info table.
#' @param threshold Significance threshold in p (default Bonferroni 0.05/m).
#' @param highlight_colour Semantic colour for significant points.
#' @return ggplot object.
manhattan_plot <- function(assoc, snp_info, threshold = NULL, title = "Manhattan plot") {
  df <- merge(assoc[, c("snp_id", "p")], snp_info, by = "snp_id")
  df <- df[order(df$chr, df$pos), , drop = FALSE]
  if (is.null(threshold)) threshold <- 0.05 / nrow(df)
  df$logp <- -log10(pmax(df$p, 1e-300))
  # cumulative position across chromosomes
  chr_levels <- unique(df$chr)
  chr_len <- vapply(chr_levels, function(c) max(df$pos[df$chr == c]), numeric(1))
  offset <- c(0, cumsum(chr_len + max(chr_len) * 0.02))[seq_along(chr_levels)]
  names(offset) <- chr_levels
  df$cum_pos <- df$pos + offset[df$chr]
  df$chr_f <- factor(df$chr, levels = chr_levels)
  df$sig <- df$p <= threshold
  pal <- rice_palette()
  alt_cols <- rep(pal$series12[c(1, 5)], length.out = length(chr_levels))
  names(alt_cols) <- chr_levels
  ticks <- vapply(chr_levels, function(c) (min(df$cum_pos[df$chr == c]) +
                                             max(df$cum_pos[df$chr == c])) / 2, numeric(1))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cum_pos, y = logp)) +
    ggplot2::geom_point(ggplot2::aes(colour = chr_f), size = 0.7, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = alt_cols, guide = "none") +
    ggplot2::geom_point(data = df[df$sig, , drop = FALSE],
                        colour = nature_semantic("up"), size = 1.0) +
    ggplot2::geom_hline(yintercept = -log10(threshold), linetype = "dashed",
                        colour = nature_semantic("reference_line"), linewidth = 0.4) +
    ggplot2::scale_x_continuous(breaks = ticks, labels = chr_levels, expand = c(0.01, 0)) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](italic(p))), title = title) +
    theme_nature()
  p
}

#' QQ plot with 95% confidence band.
qq_plot <- function(assoc, title = "QQ plot") {
  p <- sort(assoc$p)
  n <- length(p)
  if (n == 0) stop("QQ plot: no p-values supplied.", call. = FALSE)
  obs <- -log10(pmax(p, 1e-300))
  exp_mid <- -log10(stats::ppoints(n))
  ci_lo <- -log10(stats::qbeta(0.975, seq_len(n), n - seq_len(n) + 1))
  ci_hi <- -log10(stats::qbeta(0.025, seq_len(n), n - seq_len(n) + 1))
  df <- data.frame(exp = exp_mid, obs = obs, lo = ci_lo, hi = ci_hi)
  lim <- max(c(df$obs, df$exp), na.rm = TRUE) * 1.05
  ggplot2::ggplot(df, ggplot2::aes(x = exp, y = obs)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey85", alpha = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = nature_semantic("reference_line"), linewidth = 0.4) +
    ggplot2::geom_point(colour = rice_palette()$series12[1], size = 0.9) +
    ggplot2::coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
    ggplot2::labs(x = expression(Expected~-log[10](italic(p))),
                  y = expression(Observed~-log[10](italic(p))), title = title) +
    theme_nature()
}

#' Save a base-R plot (e.g. ape tree) in the three Nature formats.
save_base_plot_genomics <- function(draw_fun, file_stem, width_mm = 120, height_mm = 100) {
  pal <- rice_palette()
  for (fmt in pal$figure$formats) {
    fmt <- tolower(fmt)
    out <- paste0(file_stem, ".", if (fmt %in% c("tif", "tiff")) "tiff" else fmt)
    if (fmt %in% c("tif", "tiff")) {
      grDevices::tiff(out, width = mm_to_inch(width_mm), height = mm_to_inch(height_mm),
                      units = "in", res = pal$figure$dpi_raster, compression = "lzw")
    } else if (fmt == "pdf") {
      grDevices::pdf(out, width = mm_to_inch(width_mm), height = mm_to_inch(height_mm))
    } else {
      grDevices::png(out, width = mm_to_inch(width_mm), height = mm_to_inch(height_mm), units = "in", res = 300)
    }
    draw_fun()
    grDevices::dev.off()
  }
  invisible(NULL)
}
