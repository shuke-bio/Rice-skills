#!/usr/bin/env bash
#
# 03_vcf_filter.sh -- raw VCF -> filtered VCF (GATK hard filters / bcftools)
#
# Description:
#   Standard hard-filtering of raw variant calls before downstream analyses:
#     engine "gatk":     VariantFiltration (soft-flag) + SelectVariants
#                        --exclude-filtered (hard-remove flagged sites)
#     engine "bcftools": bcftools filter with equivalent expression
#   Defaults follow the common GATK best-practice hard-filter cutoffs,
#   adapted for rice WGS/GBS data:
#     QUAL < 30 | QD < 2.0 | FS > 60.0 | MQ < 40.0 | SOR > 3.0
#     MQRankSum < -12.5 | ReadPosRankSum < -8.0
#   If the selected toolchain is not installed the script prints install
#   instructions and exits.
#
# Usage:
#   bash scripts/variants/03_vcf_filter.sh --vcf cohort.vcf.gz \
#        --ref reference.fa --outdir results/variants/03_filter
#
#   # custom thresholds
#   bash scripts/variants/03_vcf_filter.sh --vcf cohort.vcf.gz --ref ref.fa \
#        --min-qual 40 --min-qd 2.5 --max-fs 50 --outdir results/variants/03_filter
#
# Inputs:
#   --vcf       raw VCF (.vcf.gz; indexed or auto-indexed)
#   --ref       reference FASTA (.fai needed for GATK mode)
#   --outdir    output directory (default results/variants/03_filter)
#   --engine    gatk | bcftools (default: gatk if present, else bcftools)
#   --min-qual  QUAL cutoff        (default 30)
#   --min-qd    QD cutoff          (default 2.0)
#   --max-fs    FS cutoff          (default 60.0)
#   --min-mq    MQ cutoff          (default 40.0)
#   --max-sor   SOR cutoff         (default 3.0)
#   --java-mem  Java heap for GATK (default 8g)
# Outputs:
#   <outdir>/filtered.vcf.gz (+ .tbi)      analysis-ready variants
#   <outdir>/flagged.vcf.gz (gatk mode)    soft-flagged calls (audit trail)
#   <outdir>/filter_stats.txt              before/after stats
#   <outdir>/filter_log.txt                full command log
#
# Notes:
#   - Hard filters suit small/moderate cohorts. For >100 deep-WGS samples
#     prefer GATK VQSR (rice is not a supported VQSR training species —
#     hard filtering or machine-learning filters like GATK CNNScoreVariants
#     are the practical options).
#   - Next step: bash scripts/gwas/01_plink_qc.sh --vcf <outdir>/filtered.vcf.gz

set -u

VCF=""
REF=""
OUTDIR="results/variants/03_filter"
ENGINE=""
MIN_QUAL=30
MIN_QD=2.0
MAX_FS=60.0
MIN_MQ=40.0
MAX_SOR=3.0
JAVA_MEM="8g"

while [ $# -gt 0 ]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --engine) ENGINE="$2"; shift 2 ;;
    --min-qual) MIN_QUAL="$2"; shift 2 ;;
    --min-qd) MIN_QD="$2"; shift 2 ;;
    --max-fs) MAX_FS="$2"; shift 2 ;;
    --min-mq) MIN_MQ="$2"; shift 2 ;;
    --max-sor) MAX_SOR="$2"; shift 2 ;;
    --java-mem) JAVA_MEM="$2"; shift 2 ;;
    -h|--help) sed -n '2,47p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

[ -n "$VCF" ] || { echo "ERROR: --vcf is required. See --help."; exit 1; }
[ -f "$VCF" ] || { echo "ERROR: VCF not found: $VCF"; exit 1; }

# engine auto-select
if [ -z "$ENGINE" ]; then
  if command -v gatk >/dev/null 2>&1; then ENGINE="gatk"
  elif command -v bcftools >/dev/null 2>&1; then ENGINE="bcftools"
  fi
fi

if [ "$ENGINE" = "gatk" ]; then
  if ! command -v gatk >/dev/null 2>&1; then
    echo "GATK (gatk4) not found. Install: conda install -c bioconda gatk4"
    echo "Or re-run with --engine bcftools (needs bcftools)."
    exit 0
  fi
  [ -n "$REF" ] || { echo "ERROR: --ref is required for GATK mode."; exit 1; }
  [ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
elif [ "$ENGINE" = "bcftools" ]; then
  command -v bcftools >/dev/null 2>&1 || {
    echo "bcftools not found. Install: conda install -c bioconda bcftools"; exit 0; }
else
  echo "ERROR: --engine must be 'gatk' or 'bcftools'"; exit 1
fi
echo "Engine: $ENGINE"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/filter_log.txt"
: > "$LOG"
run() { echo "+ $*" | tee -a "$LOG"; "$@" >>"$LOG" 2>&1; }

# auto-index input if needed
if command -v bcftools >/dev/null 2>&1 && [ ! -f "${VCF}.tbi" ] && [ ! -f "${VCF}.csi" ]; then
  echo "--> indexing input VCF"
  run bcftools index -t "$VCF"
fi

echo "--> stats before filtering"
if command -v bcftools >/dev/null 2>&1; then
  bcftools stats "$VCF" 2>/dev/null | grep -E "^(SN|number of records)" \
    | sed 's/^/  /' >> "$LOG"
fi

if [ "$ENGINE" = "gatk" ]; then
  FLAGGED="$OUTDIR/flagged.vcf.gz"
  echo "--> VariantFiltration (soft flags)"
  run gatk --java-options "-Xmx$JAVA_MEM" VariantFiltration \
      -R "$REF" -V "$VCF" -O "$FLAGGED" \
      --filter-expression "QD < $MIN_QD"                 --filter-name "LowQD" \
      --filter-expression "FS > $MAX_FS"                 --filter-name "HighFS" \
      --filter-expression "MQ < $MIN_MQ"                 --filter-name "LowMQ" \
      --filter-expression "SOR > $MAX_SOR"               --filter-name "HighSOR" \
      --filter-expression "MQRankSum < -12.5"            --filter-name "LowMQRankSum" \
      --filter-expression "ReadPosRankSum < -8.0"        --filter-name "LowReadPosRankSum" \
      --filter-expression "QUAL < $MIN_QUAL"             --filter-name "LowQual"

  echo "--> removing flagged sites (SelectVariants --exclude-filtered)"
  run gatk --java-options "-Xmx$JAVA_MEM" SelectVariants \
      -R "$REF" -V "$FLAGGED" -O "$OUTDIR/filtered.vcf.gz" \
      --exclude-filtered
else
  # bcftools engine: expression mirrors the GATK hard filters
  EXPR="QUAL<$MIN_QUAL || INFO/QD<$MIN_QD || INFO/FS>$MAX_FS || INFO/MQ<$MIN_MQ || INFO/SOR>$MAX_SOR"
  echo "--> bcftools filter: $EXPR"
  run bcftools filter -e "$EXPR" -s "HardFilter" -Oz -o "$OUTDIR/filtered.vcf.gz" "$VCF"
  run bcftools index -t -f "$OUTDIR/filtered.vcf.gz"
  # keep only PASS sites as a convenience copy
  run bash -c "bcftools view -f PASS -Oz -o '$OUTDIR/filtered.PASS.vcf.gz' '$OUTDIR/filtered.vcf.gz'"
  run bash -c "bcftools index -t -f '$OUTDIR/filtered.PASS.vcf.gz'"
fi

echo "--> stats after filtering"
if command -v bcftools >/dev/null 2>&1; then
  bcftools stats "$OUTDIR/filtered.vcf.gz" 2>/dev/null | grep -E "^SN" \
    | tee -a "$LOG" | tee "$OUTDIR/filter_stats.txt" | sed 's/^/  /'
fi

echo ""
echo "Filtering complete. Analysis-ready VCF: $OUTDIR/filtered.vcf.gz"
echo ""
echo "Next steps:"
echo "  1. (optional) functional annotation:"
echo "     bash scripts/variants/04_snpeff_annotate.sh --vcf $OUTDIR/filtered.vcf.gz \\"
echo "         --outdir results/variants/04_annot"
echo "  2. GWAS / population genetics:"
echo "     bash scripts/gwas/01_plink_qc.sh --vcf $OUTDIR/filtered.vcf.gz \\"
echo "         --outdir results/gwas/01_qc"
