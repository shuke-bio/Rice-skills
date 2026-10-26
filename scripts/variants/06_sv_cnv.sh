#!/usr/bin/env bash
###############################################################################
# 06_sv_cnv.sh -- long/short-read SV and CNV calling (plan-first)
#
# 路线：
#   ont   : Sniffles2 + cuteSV
#   short : Delly + Manta（可选）
#   both  : 全部
# 合并：SURVIVOR（可选）或 bcftools concat
#
# 默认只打印计划；--run 才执行。缺工具会提示。
#
# 用法：
#   bash scripts/variants/06_sv_cnv.sh --bam sample.bam --reference ref.fa --mode ont \
#       --outdir results/variants/06_sv --threads 8
#   bash scripts/variants/06_sv_cnv.sh --demo --outdir results/variants/06_sv_demo
###############################################################################
set -euo pipefail

BAM=""; REF=""; MODE="both"; OUT="results/variants/06_sv"; THREADS=8; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2;;
    --reference) REF="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then
  BAM="$OUT/demo.bam"; REF="$OUT/demo_ref.fa"
  : > "$BAM"; printf '>chr1\nACGTACGTACGT\n' > "$REF"
  echo "[sv] demo placeholders -> $OUT"
fi
[[ -n "$BAM" && -n "$REF" ]] || { echo "Need --bam and --reference, or --demo" >&2; exit 1; }

PLAN="$OUT/tool_plan.txt"
{
  echo "# SV/CNV calling plan ($(date -Iseconds))"
  echo "# bam=$BAM ref=$REF mode=$MODE"
  if [[ "$MODE" == "ont" || "$MODE" == "both" ]]; then
    echo "sniffles --input $BAM --vcf $OUT/sniffles.vcf --threads $THREADS --reference $REF"
    echo "cuteSV $BAM $REF $OUT/cutesv.vcf $OUT/cutesv_work --threads $THREADS"
  fi
  if [[ "$MODE" == "short" || "$MODE" == "both" ]]; then
    echo "delly call -g $REF -o $OUT/delly.bcf $BAM"
    echo "bcftools view $OUT/delly.bcf -Oz -o $OUT/delly.vcf.gz"
    echo "# optional: manta configManta.py --bam $BAM --referenceFasta $REF --runDir $OUT/manta"
  fi
  echo "# merge (optional): SURVIVOR merge $OUT/sv_list.txt 1000 1 1 1 0 50 $OUT/merged.vcf"
} > "$PLAN"
cat "$PLAN"

if [[ $RUN -eq 1 ]]; then
  for t in sniffles cuteSV delly bcftools; do
    command -v "$t" >/dev/null 2>&1 || echo "[sv] missing: $t" >&2
  done
  if [[ "$MODE" == "ont" || "$MODE" == "both" ]] && command -v sniffles >/dev/null 2>&1; then
    sniffles --input "$BAM" --vcf "$OUT/sniffles.vcf" --threads "$THREADS" --reference "$REF" \
      > "$OUT/sniffles.log" 2>&1 || echo "[sv] sniffles failed" >&2
  fi
  if [[ "$MODE" == "ont" || "$MODE" == "both" ]] && command -v cuteSV >/dev/null 2>&1; then
    mkdir -p "$OUT/cutesv_work"
    cuteSV "$BAM" "$REF" "$OUT/cutesv.vcf" "$OUT/cutesv_work" --threads "$THREADS" \
      > "$OUT/cutesv.log" 2>&1 || echo "[sv] cuteSV failed" >&2
  fi
  if [[ "$MODE" == "short" || "$MODE" == "both" ]] && command -v delly >/dev/null 2>&1; then
    delly call -g "$REF" -o "$OUT/delly.bcf" "$BAM" > "$OUT/delly.log" 2>&1 || echo "[sv] delly failed" >&2
  fi
  echo "[sv] finished; see $OUT"
fi
