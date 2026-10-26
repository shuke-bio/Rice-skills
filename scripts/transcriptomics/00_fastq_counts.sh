#!/usr/bin/env bash
###############################################################################
# 00_fastq_counts.sh -- FASTQ -> gene counts (fastp + HISAT2 + featureCounts)
#
# 默认只打印计划；--run 才执行。
#
# 用法：
#   bash scripts/transcriptomics/00_fastq_counts.sh --r1 R1.fq.gz --r2 R2.fq.gz \
#       --reference genome.fa --annotation genes.gtf --sample S01 --outdir results/transcriptomics/00_counts
#   bash scripts/transcriptomics/00_fastq_counts.sh --demo --outdir results/transcriptomics/00_counts_demo
###############################################################################
set -euo pipefail

R1=""; R2=""; REF=""; GTF=""; IDX=""; SAMPLE="sample"; OUT="results/transcriptomics/00_counts"
THREADS=8; STRAND=0; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --r1) R1="$2"; shift 2;;
    --r2) R2="$2"; shift 2;;
    --reference) REF="$2"; shift 2;;
    --annotation) GTF="$2"; shift 2;;
    --hisat2-index) IDX="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --strand) STRAND="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then
  R1="$OUT/demo_R1.fq"; R2="$OUT/demo_R2.fq"; REF="$OUT/demo_ref.fa"; GTF="$OUT/demo.gtf"
  printf '@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n' > "$R1"; cp "$R1" "$R2"
  printf '>chr1\nACGTACGTACGTACGTACGTACGT\n' > "$REF"
  printf 'chr1\ttest\tgene\t1\t12\t.\t+\t.\tgene_id "G1";\n' > "$GTF"
  echo "[counts] demo inputs -> $OUT"
fi
[[ -n "$R1" && -n "$REF" && -n "$GTF" ]] || { echo "Need --r1 --reference --annotation, or --demo" >&2; exit 1; }
[[ -n "$IDX" ]] || IDX="$OUT/hisat2_index/genome"

PLAN="$OUT/tool_plan.txt"
{
  echo "# FASTQ -> counts plan ($(date -Iseconds))"
  echo "fastp -i $R1 ${R2:+-I $R2} -o $OUT/${SAMPLE}_R1.clean.fq.gz ${R2:+-O $OUT/${SAMPLE}_R2.clean.fq.gz} --thread $THREADS -j $OUT/fastp.json -h $OUT/fastp.html"
  if [[ ! -f "${IDX}.1.ht2" ]]; then
    echo "hisat2-build -p $THREADS $REF $IDX"
  fi
  echo "hisat2 -p $THREADS -x $IDX -1 $OUT/${SAMPLE}_R1.clean.fq.gz ${R2:+-2 $OUT/${SAMPLE}_R2.clean.fq.gz} | samtools sort -@ $THREADS -o $OUT/${SAMPLE}.bam -"
  echo "samtools index $OUT/${SAMPLE}.bam"
  echo "featureCounts -T $THREADS -a $GTF -o $OUT/${SAMPLE}_counts.txt -s $STRAND $OUT/${SAMPLE}.bam"
} > "$PLAN"
cat "$PLAN"

if [[ $RUN -eq 1 ]]; then
  for t in fastp hisat2 samtools featureCounts; do
    command -v "$t" >/dev/null 2>&1 || { echo "[counts] missing tool: $t" >&2; exit 1; }
  done
  fastp -i "$R1" ${R2:+-I "$R2"} -o "$OUT/${SAMPLE}_R1.clean.fq.gz" ${R2:+-O "$OUT/${SAMPLE}_R2.clean.fq.gz"} \
    --thread "$THREADS" -j "$OUT/fastp.json" -h "$OUT/fastp.html" > "$OUT/fastp.log" 2>&1
  if [[ ! -f "${IDX}.1.ht2" ]]; then hisat2-build -p "$THREADS" "$REF" "$IDX" > "$OUT/hisat2_build.log" 2>&1; fi
  hisat2 -p "$THREADS" -x "$IDX" -1 "$OUT/${SAMPLE}_R1.clean.fq.gz" ${R2:+-2 "$OUT/${SAMPLE}_R2.clean.fq.gz"} \
    | samtools sort -@ "$THREADS" -o "$OUT/${SAMPLE}.bam" -
  samtools index "$OUT/${SAMPLE}.bam"
  featureCounts -T "$THREADS" -a "$GTF" -o "$OUT/${SAMPLE}_counts.txt" -s "$STRAND" "$OUT/${SAMPLE}.bam" \
    > "$OUT/featurecounts.log" 2>&1
  echo "[counts] finished -> $OUT/${SAMPLE}_counts.txt"
fi
