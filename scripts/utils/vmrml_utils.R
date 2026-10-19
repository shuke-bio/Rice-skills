###############################################################################
# vmrml_utils.R -- helpers for the 3VmrMLM-family GWAS wrappers
#
# Shared helpers for scripts/gwas/09_multilocus_qtn.R (3VmrMLM family:
# IIIVmrMLM / Fast3VmrMLM), scripts/gwas/10_case_control_qtn.R
# (FastBiCmrMLM) and scripts/gwas/11_multiparent_qtl.R (mppQTL).
#
# Nothing here executes on source(); functions are called explicitly.
# All upstream engines are installed at RUNTIME (never vendored):
#   - IIIVmrMLM (GPL-3, CRAN)        : CRAN binary, loads on any R 4.x
#   - Fast3VmrMLM (GPL-3 per package DESCRIPTION; repo LICENSE says
#     BSD-3 - treated conservatively as GPL-3) : GitHub zip, R 4.3-built
#     binary; older R versions load it, R >= 4.4 usually fails at DLL load
#   - FastBiCmrMLM (no LICENSE, Linux-only)     : GitHub zip, Linux only
#   - mppQTL (GPL-3)                            : GitHub zip, loads on R 4.x
###############################################################################

# ---------------------------------------------------------------------------
# Engine load check with actionable guidance
# ---------------------------------------------------------------------------
engine_ready <- function(pkg, extra_msg = NULL) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  msg <- NULL
  if (!ok) {
    msg <- paste0("R package '", pkg, "' is not installed or fails to load.")
    if (!is.null(extra_msg)) msg <- paste(msg, extra_msg)
  }
  list(ok = ok, message = msg)
}

print_install_hint <- function(lines) {
  for (l in lines) message("  ", l)
}

# ---------------------------------------------------------------------------
# PLINK binary writer (BED v1f, 0/1/2 dosage)
#
# Writes <prefix>.bed / .bim / .fam so the synthetic demo data (or a dosage
# matrix exported by 02_plink_to_matrix.R) can be fed directly to the
# 3VmrMLM engines, which read PLINK binary sets natively.
# ---------------------------------------------------------------------------
write_plink_set <- function(prefix, geno, snp_info) {
  # geno:      samples (rows) x SNPs (cols), values 0/1/2/NA
  # snp_info:  data.frame with snp_id, chr, pos
  n <- nrow(geno); m <- ncol(geno)
  snp_info <- snp_info[snp_info$snp_id %in% colnames(geno), ]
  snp_info <- snp_info[match(colnames(geno), snp_info$snp_id), ]
  if (any(is.na(snp_info$pos))) stop("snp_info is missing rows for some SNPs.", call. = FALSE)

  fam <- data.frame(
    FID = rownames(geno), IID = rownames(geno),
    PAT = 0, MAT = 0, SEX = 0, PHENO = -9
  )
  utils::write.table(fam, paste0(prefix, ".fam"), sep = " ", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  bim <- data.frame(
    chr = snp_info$chr, snp = snp_info$snp_id, gd = 0,
    pos = as.integer(snp_info$pos), a1 = "A", a2 = "G"
  )
  utils::write.table(bim, paste0(prefix, ".bim"), sep = " ", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)

  # BED: magic 6c 1b 01, then per-SNP 2-bit codes (LSB-first, 4 samples/byte)
  con <- file(paste0(prefix, ".bed"), open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  bytes_per_snp <- ceiling(n / 4)
  # integer code matrix (SNP x sample): 0/1/2 observed, 3 = missing
  g <- matrix(3L, nrow = n, ncol = m)
  valid <- !is.na(geno)
  g[valid] <- as.integer(round(geno[valid]))
  g[g < 0 | g > 2] <- 3L
  gi <- t(g)                                            # SNP x sample
  rawmat <- matrix(0L, nrow = m, ncol = bytes_per_snp)
  for (k in seq_len(bytes_per_snp)) {
    s0 <- ((k - 1L) * 4L + 1L)
    cols <- seq(s0, min(s0 + 3L, n))
    for (j in seq_along(cols)) {
      # NOTE: bitwOr(), not `|` -- `|` is logical OR and returns TRUE/FALSE.
      rawmat[, k] <- bitwOr(rawmat[, k], bitwShiftL(gi[, cols[j]], (j - 1L) * 2L))
    }
  }
  writeBin(as.raw(as.vector(t(rawmat))), con)
  invisible(prefix)
}

# ---------------------------------------------------------------------------
# Synthetic PLINK dataset with planted QTNs (demo data, seed = 42)
# ---------------------------------------------------------------------------
make_demo_plink <- function(n_ind = 220, n_snp = 600, n_chr = 3,
                            n_qtl = 2, h2 = 0.5, seed = 42) {
  set.seed(seed)
  chr <- rep(seq_len(n_chr), length.out = n_snp)
  per_chr <- table(chr)
  pos <- unlist(lapply(names(per_chr), function(cc) {
    k <- per_chr[[cc]]
    sort(sample(1:2e6, k, replace = FALSE))
  }))
  snp_id <- paste0("S", seq_len(n_snp))
  snp_info <- data.frame(snp_id = snp_id, chr = chr, pos = pos,
                         stringsAsFactors = FALSE)

  af <- runif(n_snp, 0.15, 0.5)
  geno <- matrix(rbinom(n_ind * n_snp, 2, rep(af, each = n_ind)),
                 nrow = n_ind, ncol = n_snp,
                 dimnames = list(paste0("Ind", seq_len(n_ind)), snp_id))

  qtl_idx <- sort(sample(n_snp, n_qtl))
  qtl_beta <- runif(n_qtl, 1.2, 2.0) * sample(c(-1, 1), n_qtl, TRUE)
  y <- matrix(0, nrow = n_ind, ncol = 2)
  for (j in seq_len(n_qtl)) {
    y[, 1] <- y[, 1] + geno[, qtl_idx[j]] * qtl_beta[j]
  }
  # polygenic background: 30 tiny effects
  bg_idx <- sample(setdiff(seq_len(n_snp), qtl_idx), 30)
  bg_beta <- rnorm(30, 0, 0.25)
  for (j in seq_along(bg_idx)) y[, 1] <- y[, 1] + geno[, bg_idx[j]] * bg_beta[j]
  y[, 1] <- y[, 1] + rnorm(n_ind, 0, sqrt(1 - h2) / sqrt(h2) * sd(y[, 1]))
  y[, 1] <- scale(y[, 1]) * 10 + 100
  # trait 2: uncorrelated background + mild QTL sharing
  y[, 2] <- scale(rnorm(n_ind, 100, 10))

  pheno <- data.frame(sample_id = rownames(geno), trait1 = y[, 1], trait2 = y[, 2],
                      stringsAsFactors = FALSE)
  list(geno = geno, snp_info = snp_info, pheno = pheno, qtl_idx = qtl_idx,
       qtl_snp = snp_id[qtl_idx])
}

# Phenotype in engine format: ID column named "<Phenotype>"
write_engine_pheno <- function(pheno, path, env_labels = NULL) {
  # pheno: data.frame, first column sample_id, remaining numeric traits.
  # env_labels: NULL = single environment; otherwise e.g. c("E1","E1","E2","E2")
  #             naming each trait column's environment (multi-env wide format).
  df <- pheno
  colnames(df)[1] <- "<Phenotype>"
  if (!is.null(env_labels)) {
    colnames(df)[-1] <- paste(names(pheno)[-1], env_labels, sep = "__")
  }
  utils::write.csv(df, path, row.names = FALSE, quote = FALSE)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Normalise engine result workbooks into one tidy table
# ---------------------------------------------------------------------------
# openxlsx renamed sheetNames -> getSheetNames (exported) in 4.2.5.3+; support both.
xlsx_sheet_names <- function(path) {
  ns <- asNamespace("openxlsx")
  if (exists("getSheetNames", envir = ns, mode = "function")) {
    get("getSheetNames", envir = ns)(path)
  } else {
    get("sheetNames", envir = ns)(path)
  }
}

read_engine_result <- function(xlsx_path, sheet = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("openxlsx is required to read engine result workbooks.", call. = FALSE)
  }
  sheets <- xlsx_sheet_names(xlsx_path)
  if (is.null(sheet)) {
    # first sheet that looks like a result table
    sheet <- sheets[1]
  }
  df <- openxlsx::read.xlsx(xlsx_path, sheet = sheet)
  names(df) <- trimws(names(df))
  df
}

# Canonical QTN summary table: trait, marker, chr, pos, lod, additive, pve, p
# Covers the actual column spellings of the engines:
#   IIIVmrMLM : Marker/Chromosome/Position.(bp)/LOD/add/dom/variance/r2(%)/P-value
#   Fast3VmrMLM: same family with 加性效应/显性效应/P-value style headers in some sheets
normalise_qtns <- function(df, trait_label = NULL) {
  pick <- function(...) {
    cand <- c(...)
    hit <- cand[cand %in% names(df)]
    if (length(hit) == 0) return(rep(NA, nrow(df)))
    df[[hit[1]]]
  }
  out <- data.frame(
    trait = if (is.null(trait_label)) pick("Trait", "trait") else trait_label,
    marker = pick("Marker", "marker", "标记名"),
    chr = pick("Chromosome", "chr", "chromosome", "染色体"),
    pos = pick("Position.(bp)", "Position", "pos", "position", "标记位置"),
    lod = pick("LOD.score", "LOD得分", "LOD", "LOD score"),
    additive = pick("Additive.effect", "add", "加性效应", "Additive"),
    dominance = pick("Dominance.effect", "dom", "显性效应", "Dominance"),
    pve = pick("PVE", "pve", "variance", "r2(%)", "QTN遗传方差贡献率", "贡献率"),
    p = pick("P-value", "P.value", "pvalue", "p值", "P值"),
    stringsAsFactors = FALSE
  )
  out
}

# ---------------------------------------------------------------------------
# run_info.txt writer (project convention)
# ---------------------------------------------------------------------------
write_run_info <- function(path, entries) {
  lines <- c("# Rice-skills run info",
             paste0("# generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
             paste0("script: ", entries$script), "")
  for (nm in setdiff(names(entries), "script")) {
    lines <- c(lines, paste0(nm, ": ", entries[[nm]]))
  }
  writeLines(lines, path)
  invisible(path)
}
