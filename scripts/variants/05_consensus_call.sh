#!/usr/bin/env bash
###############################################################################
# 05_consensus_call.sh -- consensus variant calling (bcftools + FreeBayes + GATK)
#
# 流程：
#   1) 三个呼叫器分别 call（GATK 需 BAM；无 BAM 时可用 --bam 提供，或用 bwa/samtools 先比对）
#   2) bcftools norm 规范化
#   3) bcftools isec -n=2 保留至少 2 个呼叫器支持的位点
#
# 默认只打印计划；--run 才执行。缺工具会明确提示。
#
# 用法：
#   bash scripts/variants/05_consensus_call.sh --reference ref.fa --reads1 R1.fq --reads2 R2.fq \
#       --sample S01 --outdir results/variants/05_consensus
#   bash scripts/variants/05_consensus_call.sh --demo --outdir results/variants/05_consensus_demo
###############################################################################
set -euo pipefail

REF=""; R1=""; R2=""; BAM=""; SAMPLE="sample"; OUT="results/variants/05_consensus"
THREADS=8; RUN=0; DEMO=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference) REF="$2"; shift 2;;
    --reads1) R1="$2"; shift 2;;
    --reads2) R2="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
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
  REF="$OUT/demo_ref.fa"; R1="$OUT/demo_R1.fq"
  printf '>chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' > "$REF"
  printf '@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n' > "$R1"
  echo "[consensus] demo reference/reads -> $OUT"
fi

[[ -n "$REF" && -f "$REF" ]] || { echo "Need --reference or --demo" >&2; exit 1; }

PLAN="$OUT/tool_plan.txt"
{
  echo "# Consensus variant calling plan ($(date -Iseconds))"
  echo "# reference: $REF"
  echo "# sample: $SAMPLE"
  if [[ -z "$BAM" ]]; then
    echo "# No --bam: align first (bwa + samtools)"
    echo "bwa mem -t $THREADS $REF $R1 ${R2:-} | samtools sort -@ $THREADS -o $OUT/$SAMPLE.bam -"
    echo "samtools index $OUT/$SAMPLE.bam"
    BAM="$OUT/$SAMPLE.bam"
  fi
  echo "bcftools mpileup -f $REF -Ou -a AD,DP -q 20 -Q 20 $BAM | bcftools call -mv -Oz -o $OUT/bcftools.vcf.gz"
  echo "freebayes -f $REF -p 1 $BAM > $OUT/freebayes.vcf"
  echo "gatk HaplotypeCaller -R $REF -I $BAM -O $OUT/gatk.vcf"
  echo "bcftools norm -f $REF -m -any $OUT/bcftools.vcf.gz -Oz -o $OUT/bcftools.norm.vcf.gz"
  echo "bcftools norm -f $REF -m -any $OUT/freebayes.vcf -Oz -o $OUT/freebayes.norm.vcf.gz"
  echo "bcftools norm -f $REF -m -any $OUT/gatk.vcf -Oz -o $OUT/gatk.norm.vcf.gz"
  echo "bcftools isec -p $OUT/isec -n=2 -w1 $OUT/bcftools.norm.vcf.gz $OUT/freebayes.norm.vcf.gz $OUT/gatk.norm.vcf.gz"
  echo "# consensus calls: $OUT/isec/0000.vcf (present in >=2 callers)"
} > "$PLAN"
cat "$PLAN"

if [[ $RUN -eq 1 ]]; then
  for t in bcftools freebayes gatk bwa samtools; do
    command -v "$t" >/dev/null 2>&1 || echo "[consensus] missing tool: $t" >&2
  done
  if [[ -z "$BAM" && -n "$R1" ]]; then
    bwa mem -t "$THREADS" "$REF" "$R1" ${R2:+$R2} 2>/dev/null | samtools sort -@ "$THREADS" -o "$OUT/$SAMPLE.bam" -
    samtools index "$OUT/$SAMPLE.bam"
    BAM="$OUT/$SAMPLE.bam"
  fi
  bcftools mpileup -f "$REF" -Ou -a AD,DP -q 20 -Q 20 "$BAM" | bcftools call -mv -Oz -o "$OUT/bcftools.vcf.gz"
  freebayes -f "$REF" -p 1 "$BAM" > "$OUT/freebayes.vcf"
  if command -v gatk >/dev/null 2>&1; then
    gatk HaplotypeCaller -R "$REF" -I "$BAM" -O "$OUT/gatk.vcf" > "$OUT/gatk.log" 2>&1 || echo "[consensus] GATK failed" >&2
  fi
  for v in bcftools freebayes gatk; do
    [[ -f "$OUT/$v.vcf" || -f "$OUT/$v.vcf.gz" ]] || continue
    if [[ "$v" == "bcftools" ]]; then
      bcftools norm -f "$REF" -m -any "$OUT/$v.vcf.gz" -Oz -o "$OUT/$v.norm.vcf.gz"
    else
      bcftools norm -f "$REF" -m -any "$OUT/$v.vcf" -Oz -o "$OUT/$v.norm.vcf.gz"
    fi
  done
  bcftools isec -p "$OUT/isec" -n=2 -w1 "$OUT/bcftools.norm.vcf.gz" "$OUT/freebayes.norm.vcf.gz" "$OUT/gatk.norm.vcf.gz" || true
  echo "[consensus] done -> $OUT/isec/0000.vcf"
fi
