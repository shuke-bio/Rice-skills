#!/usr/bin/env bash
###############################################################################
# 08_phasing.sh -- read-backed / statistical phasing (plan-first)
#
# 用法：
#   bash scripts/variants/08_phasing.sh --vcf sample.vcf.gz --bam sample.bam --sample S01 --tool whatshap
#   bash scripts/variants/08_phasing.sh --vcf cohort.vcf.gz --tool shapeit
#   bash scripts/variants/08_phasing.sh --demo --outdir results/variants/08_phasing_demo
###############################################################################
set -euo pipefail
VCF=""; BAM=""; SAMPLE="sample"; OUT="results/variants/08_phasing"; TOOL="whatshap"; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --tool) TOOL="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,30p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then VCF="$OUT/demo.vcf.gz"; BAM="$OUT/demo.bam"; : > "$VCF"; : > "$BAM"; fi
[[ -n "$VCF" ]] || { echo "Need --vcf, or --demo" >&2; exit 1; }
PLAN="$OUT/tool_plan.txt"
{
  echo "# Phasing plan ($(date -Iseconds)) tool=$TOOL"
  if [[ "$TOOL" == "whatshap" ]]; then
    [[ -n "$BAM" ]] || { echo "# whatshap needs --bam" >&2; }
    echo "whatshap polyphase $VCF $BAM -o $OUT/${SAMPLE}.phased.vcf.gz --reference=ref.fa"
  else
    echo "shapeit5 phase --input $VCF --output $OUT/${SAMPLE}.phased.vcf.gz --thread 8"
  fi
  echo "bcftools index -t $OUT/${SAMPLE}.phased.vcf.gz"
} > "$PLAN"
cat "$PLAN"
if [[ $RUN -eq 1 ]]; then
  if [[ "$TOOL" == "whatshap" ]]; then
    command -v whatshap >/dev/null 2>&1 || { echo "[phasing] whatshap not found" >&2; exit 1; }
    whatshap polyphase "$VCF" "$BAM" -o "$OUT/${SAMPLE}.phased.vcf.gz" > "$OUT/whatshap.log" 2>&1
  else
    command -v shapeit5 >/dev/null 2>&1 || { echo "[phasing] shapeit5 not found" >&2; exit 1; }
    shapeit5 phase --input "$VCF" --output "$OUT/${SAMPLE}.phased.vcf.gz" --thread 8 > "$OUT/shapeit.log" 2>&1
  fi
  echo "[phasing] finished -> $OUT/${SAMPLE}.phased.vcf.gz"
fi
