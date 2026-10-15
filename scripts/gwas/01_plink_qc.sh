#!/usr/bin/env bash
#
# 01_plink_qc.sh -- VCF -> PLINK conversion + genotype QC (template runner)
#
# Description:
#   Standard PLINK front-end for the GWAS module: converts a VCF to a PLINK
#   fileset, applies MAF / missingness / HWE filters, LD-prunes, and writes a
#   .raw dosage export for scripts/gwas/02_plink_to_matrix.R.
#   If PLINK is not installed the script prints install instructions and exits.
#
# Usage:
#   bash scripts/gwas/01_plink_qc.sh --vcf genotypes.vcf.gz \
#        --outdir results/gwas/01_qc --prefix rice
#
# Inputs:
#   --vcf      input VCF (can be .vcf.gz) OR --bed <existing plink prefix>
#   --outdir   output directory (default results/gwas/01_qc)
#   --prefix   output fileset name (default "rice")
#   --maf      minor allele frequency filter (default 0.05)
#   --geno     per-SNP missingness filter (default 0.1)
#   --hwe      HWE filter (default 1e-6; use 0 to skip)
#   --no-prune skip LD pruning (--indep-pairwise 50 5 0.2)
# Outputs:
#   <outdir>/<prefix>.{bed,bim,fam}          filtered PLINK fileset
#   <outdir>/<prefix>.prune.in               LD-pruned SNP list
#   <outdir>/<prefix>.raw + .map             additive dosage export (--recode A)
#   <outdir>/plink_log.txt                   full log
#
# Notes:
#   - PLINK install: conda install -c bioconda plink  (or
#     https://www.cog-genomics.org/plink/ official static binaries).
#   - Then run: Rscript scripts/gwas/02_plink_to_matrix.R \
#         --raw <outdir>/<prefix>.raw --map <outdir>/<prefix>.map \
#         --prune <outdir>/<prefix>.prune.in --outdir <outdir>/matrix

set -u

VCF=""
BED=""
OUTDIR="results/gwas/01_qc"
PREFIX="rice"
MAF=0.05
GENO=0.1
HWE=1e-6
PRUNE=1
PLINK_BIN="${PLINK_BIN:-}"

while [ $# -gt 0 ]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --bed) BED="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --maf) MAF="$2"; shift 2 ;;
    --geno) GENO="$2"; shift 2 ;;
    --hwe) HWE="$2"; shift 2 ;;
    --no-prune) PRUNE=0; shift ;;
    --plink) PLINK_BIN="$2"; shift 2 ;;
    -h|--help) sed -n '2,38p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

if [ -z "$VCF" ] && [ -z "$BED" ]; then
  echo "ERROR: provide --vcf <genotypes.vcf[.gz]> (or --bed <plink prefix>)."
  exit 1
fi

# --- locate PLINK -------------------------------------------------------------
if [ -z "$PLINK_BIN" ]; then
  if command -v plink2 >/dev/null 2>&1; then PLINK_BIN="plink2"
  elif command -v plink >/dev/null 2>&1; then PLINK_BIN="plink"
  fi
fi
if [ -z "$PLINK_BIN" ]; then
  echo "PLINK not found (plink / plink2)."
  echo ""
  echo "Install options:"
  echo "  1. conda install -c bioconda plink        (Linux / WSL / macOS)"
  echo "  2. Official binaries: https://www.cog-genomics.org/plink/"
  echo ""
  echo "Without PLINK you can still run the whole GWAS/popgen pipeline from a"
  echo "dosage matrix CSV: see scripts/utils/make_demo_genomics.py for the"
  echo "expected format (genotype.csv + snp_info.csv), then start at"
  echo "scripts/gwas/03_association.R."
  exit 0
fi
echo "Using PLINK: $PLINK_BIN ($(${PLINK_BIN} --version 2>/dev/null | head -1))"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/plink_log.txt"
: > "$LOG"

run() {
  echo "+ $*" | tee -a "$LOG"
  "$@" >>"$LOG" 2>&1
}

# --- convert + QC -------------------------------------------------------------
if [ -n "$VCF" ]; then
  [ -f "$VCF" ] || { echo "ERROR: VCF not found: $VCF"; exit 1; }
  echo "--> converting VCF to PLINK fileset" | tee -a "$LOG"
  run "$PLINK_BIN" --vcf "$VCF" --allow-extra-chr --make-bed --out "$OUTDIR/$PREFIX"
fi
if [ -n "$BED" ]; then
  BASE="${BED%.bed}"
  echo "--> using existing PLINK fileset: $BASE" | tee -a "$LOG"
  cp "$BASE".bed "$BASE".bim "$BASE".fam "$OUTDIR/" 2>/dev/null || true
fi

echo "--> QC filters: --maf $MAF --geno $GENO --hwe $HWE" | tee -a "$LOG"
HWE_ARG="--hwe $HWE"
[ "$HWE" = "0" ] && HWE_ARG=""
run "$PLINK_BIN" --bfile "$OUTDIR/$PREFIX" --allow-extra-chr \
    --maf "$MAF" --geno "$GENO" $HWE_ARG --make-bed --out "$OUTDIR/${PREFIX}_qc"
mv "$OUTDIR/${PREFIX}_qc".bed "$OUTDIR/$PREFIX".bed 2>/dev/null || true
mv "$OUTDIR/${PREFIX}_qc".bim "$OUTDIR/$PREFIX".bim 2>/dev/null || true
mv "$OUTDIR/${PREFIX}_qc".fam "$OUTDIR/$PREFIX".fam 2>/dev/null || true

# --- LD pruning + dosage export -------------------------------------------------
if [ "$PRUNE" -eq 1 ]; then
  echo "--> LD pruning (indep-pairwise 50 5 0.2)" | tee -a "$LOG"
  run "$PLINK_BIN" --bfile "$OUTDIR/$PREFIX" --allow-extra-chr \
      --indep-pairwise 50 5 0.2 --out "$OUTDIR/$PREFIX"
fi

REMAIN_ARG=""
[ "$PRUNE" -eq 1 ] && REMAIN_ARG="--extract $OUTDIR/$PREFIX.prune.in"
echo "--> exporting additive dosages (.raw) for 02_plink_to_matrix.R" | tee -a "$LOG"
run "$PLINK_BIN" --bfile "$OUTDIR/$PREFIX" --allow-extra-chr $REMAIN_ARG \
    --recode A --out "$OUTDIR/$PREFIX"

echo ""
echo "QC complete. Files in $OUTDIR :"
ls "$OUTDIR" | sed 's/^/  /'
echo ""
echo "Next step:"
echo "  Rscript scripts/gwas/02_plink_to_matrix.R \\"
echo "      --raw $OUTDIR/$PREFIX.raw --map $OUTDIR/$PREFIX.map \\"
echo "      --prune $OUTDIR/$PREFIX.prune.in --outdir $OUTDIR/matrix"
