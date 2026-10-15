#!/usr/bin/env bash
#
# 04_admixture.sh -- ADMIXTURE ancestry estimation (template runner)
#
# Description:
#   Runs ADMIXTURE for K = k_min..k_max on a PLINK fileset with cross
#   validation, and collects the CV error curve (the standard way to choose K).
#   If the ADMIXTURE binary is not installed the script explains how to get it
#   and exits without touching any files.
#
#   No-ADMIXTURE alternative: use scripts/popgen/01_pca_structure.R (PCA) and
#   scripts/popgen/05_plot_admixture.R can still plot any .Q matrix you obtain
#   elsewhere (e.g. from a web server run).
#
# Usage:
#   bash scripts/popgen/04_admixture.sh --bed plink/rice.bed \
#        --k-min 2 --k-max 8 --outdir results/popgen/04_admixture
#
# Inputs:
#   --bed   PLINK fileset prefix or .bed path (needs .bed/.bim/.fam together);
#           build it from a VCF with: plink --vcf in.vcf --make-bed --out rice
# Outputs:
#   <outdir>/<prefix>.<K>.Q / .P       ancestry proportions / frequencies
#   <outdir>/cv_errors.txt             "CV error (K=x): y" lines
#   <outdir>/admixture_log.txt         full log
#
# Notes:
#   - Install ADMIXTURE (Linux / WSL / HPC): conda install -c bioconda admixture
#     or download the official binaries: https://dalexander.github.io/admixture/
#   - ADMIXTURE is not available for native Windows; run this step on Linux/WSL.
#   - Then plot with: Rscript scripts/popgen/05_plot_admixture.R --q-dir ...

set -u

BED=""
K_MIN=2
K_MAX=8
THREADS=4
OUTDIR="results/popgen/04_admixture"
ADMIXTURE_BIN="admixture"

while [ $# -gt 0 ]; do
  case "$1" in
    --bed) BED="$2"; shift 2 ;;
    --k-min) K_MIN="$2"; shift 2 ;;
    --k-max) K_MAX="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --admixture) ADMIXTURE_BIN="$2"; shift 2 ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

if [ -z "$BED" ]; then
  echo "ERROR: --bed is required (PLINK fileset prefix or .bed path)."
  echo "Hint: build it from a VCF with:"
  echo "  plink --vcf genotypes.vcf --maf 0.05 --geno 0.1 --make-bed --out rice"
  exit 1
fi
PREFIX="${BED%.bed}"
if [ ! -f "${PREFIX}.bed" ] || [ ! -f "${PREFIX}.bim" ] || [ ! -f "${PREFIX}.fam" ]; then
  echo "ERROR: PLINK fileset not found for prefix '${PREFIX}'"
  echo "Hint: expected ${PREFIX}.bed / .bim / .fam next to each other."
  echo "      Build it from a VCF with: plink --vcf in.vcf --make-bed --out ${PREFIX}"
  exit 1
fi

if ! command -v "$ADMIXTURE_BIN" >/dev/null 2>&1; then
  echo "ADMIXTURE binary not found: '$ADMIXTURE_BIN'"
  echo ""
  echo "ADMIXTURE runs on Linux / WSL / HPC (no native Windows build). Options:"
  echo "  1. conda install -c bioconda admixture"
  echo "  2. Official binaries: https://dalexander.github.io/admixture/"
  echo ""
  echo "No-ADMIXTURE alternatives available in this package:"
  echo "  - Rscript scripts/popgen/01_pca_structure.R   (PCA structure)"
  echo "  - Rscript scripts/popgen/05_plot_admixture.R  (plots any .Q matrix)"
  exit 0
fi

mkdir -p "$OUTDIR"
LOG="$OUTDIR/admixture_log.txt"
CV="$OUTDIR/cv_errors.txt"
: > "$LOG"
: > "$CV"

echo "[$(date '+%F %T')] ADMIXTURE K = ${K_MIN}..${K_MAX} on ${PREFIX}" | tee -a "$LOG"
for K in $(seq "$K_MIN" "$K_MAX"); do
  echo "--> K = $K" | tee -a "$LOG"
  "$ADMIXTURE_BIN" --cv "${PREFIX}.bed" "$K" -j"$THREADS" 2>&1 | tee -a "$LOG" | \
    grep -E "CV error" | tee -a "$CV"
  # ADMIXTURE writes <prefix>.<K>.Q/.P next to the .bed file; move into outdir
  if [ -f "${PREFIX}.${K}.Q" ]; then
    mv "${PREFIX}.${K}.Q" "${PREFIX}.${K}.P" "$OUTDIR"/ 2>/dev/null || true
  fi
done

echo "" | tee -a "$LOG"
echo "Done. Outputs:"
echo "  Q/P matrices : $OUTDIR"
echo "  CV errors    : $CV"
echo ""
echo "Choose the K with the lowest CV error, then plot:"
echo "  Rscript scripts/popgen/05_plot_admixture.R --q-dir $OUTDIR --cv $CV \\"
echo "      --outdir $OUTDIR/figures"
