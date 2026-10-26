#!/usr/bin/env bash
###############################################################################
# 07_deepvariant.sh -- DeepVariant calling (plan-first)
#
# 用法：
#   bash scripts/variants/07_deepvariant.sh --reference ref.fa --bam sample.bam --sample S01
#   bash scripts/variants/07_deepvariant.sh --demo --outdir results/variants/07_deepvariant_demo
###############################################################################
set -euo pipefail
REF=""; BAM=""; SAMPLE="sample"; OUT="results/variants/07_deepvariant"; MODEL="WGS"; SHARDS=8; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference) REF="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --model) MODEL="$2"; shift 2;;
    --shards) SHARDS="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,30p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then REF="$OUT/ref.fa"; BAM="$OUT/demo.bam"; : > "$BAM"; printf '>chr1\nACGT\n' > "$REF"; fi
[[ -n "$REF" && -n "$BAM" ]] || { echo "Need --reference and --bam, or --demo" >&2; exit 1; }
PLAN="$OUT/tool_plan.txt"
{
  echo "# DeepVariant plan ($(date -Iseconds))"
  echo "run_deepvariant --model_type=$MODEL --ref=$REF --reads=$BAM --output_vcf=$OUT/$SAMPLE.vcf.gz --output_gvcf=$OUT/$SAMPLE.g.vcf.gz --sample_name=$SAMPLE --num_shards=$SHARDS"
  echo "bcftools index -t $OUT/$SAMPLE.vcf.gz"
} > "$PLAN"
cat "$PLAN"
if [[ $RUN -eq 1 ]]; then
  command -v run_deepvariant >/dev/null 2>&1 || { echo "[deepvariant] not found" >&2; exit 1; }
  run_deepvariant --model_type="$MODEL" --ref="$REF" --reads="$BAM" \
    --output_vcf="$OUT/$SAMPLE.vcf.gz" --output_gvcf="$OUT/$SAMPLE.g.vcf.gz" \
    --sample_name="$SAMPLE" --num_shards="$SHARDS" > "$OUT/deepvariant.log" 2>&1
  bcftools index -t "$OUT/$SAMPLE.vcf.gz" || true
  echo "[deepvariant] finished -> $OUT/$SAMPLE.vcf.gz"
fi
