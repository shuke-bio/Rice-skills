#!/usr/bin/env bash
#
# 02_gatk_call.sh -- BAM -> raw VCF (GATK HaplotypeCaller / bcftools)
#
# Description:
#   Variant discovery from analysis-ready BAM(s):
#     mode "gatk" (default, recommended for rice):
#       per-sample HaplotypeCaller -ERC GVCF
#       -> GenomicsDBImport (multi-sample) -> GenotypeGVCFs
#       -> single cohort VCF. Single-sample runs call GVCF then GT directly.
#     mode "bcftools" (lightweight alternative):
#       bcftools mpileup | bcftools call -mv  (fast, lower memory)
#   If the selected toolchain is not installed the script prints install
#   instructions and exits.
#
# Usage:
#   # single sample (GATK)
#   bash scripts/variants/02_gatk_call.sh --bam sample.sorted.bam \
#        --ref reference.fa --sample IR24 --outdir results/variants/02_call
#
#   # multi-sample cohort (repeat --bam, or use --bamdir)
#   bash scripts/variants/02_gatk_call.sh \
#        --bam a.sorted.bam --bam b.sorted.bam --bam c.sorted.bam \
#        --ref reference.fa --outdir results/variants/02_call
#
#   # lightweight engine
#   bash scripts/variants/02_gatk_call.sh --engine bcftools --bam a.sorted.bam \
#        --ref reference.fa --outdir results/variants/02_call
#
# Inputs:
#   --bam      BAM file (repeatable); or --bamdir <dir> (all *.sorted.bam)
#   --ref      reference genome FASTA (+ .fai; created if missing)
#   --sample   sample name (single-sample mode; default from BAM filename)
#   --outdir   output directory (default results/variants/02_call)
#   --engine   gatk | bcftools (default gatk)
#   --threads  threads for GATK/bcftools (default 8)
#   --java-mem Java heap for GATK (default 8g; raise for big cohorts)
#   --intervals <chr[,chr]>  restrict calling to region(s) (optional)
# Outputs:
#   <outdir>/cohort.vcf.gz (+ .tbi)   raw variant calls (all samples)
#   <outdir>/<sample>.g.vcf.gz        per-sample GVCFs (gatk mode)
#   <outdir>/call_log.txt             full command log
#
# Notes:
#   - Tools install: conda install -c bioconda gatk4  (or bcftools).
#     GATK needs Java 17+; the conda package bundles it.
#   - BAMs must carry read groups (@RG) — 01_bwa_align.sh adds them.
#   - Next step: bash scripts/variants/03_vcf_filter.sh

set -u

BAMS=()
BAMDIR=""
REF=""
SAMPLE=""
OUTDIR="results/variants/02_call"
ENGINE="gatk"
THREADS=8
JAVA_MEM="8g"
INTERVALS=""

while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAMS+=("$2"); shift 2 ;;
    --bamdir) BAMDIR="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --engine) ENGINE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --java-mem) JAVA_MEM="$2"; shift 2 ;;
    --intervals) INTERVALS="$2"; shift 2 ;;
    -h|--help) sed -n '2,50p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

if [ -n "$BAMDIR" ]; then
  for b in "$BAMDIR"/*.sorted.bam; do [ -f "$b" ] && BAMS+=("$b"); done
fi
if [ "${#BAMS[@]}" -eq 0 ] || [ -z "$REF" ]; then
  echo "ERROR: provide --bam (repeatable) or --bamdir, plus --ref. See --help."
  exit 1
fi
for b in "${BAMS[@]}"; do [ -f "$b" ] || { echo "ERROR: BAM not found: $b"; exit 1; }; done
[ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }

# single-sample name fallback
if [ -z "$SAMPLE" ] && [ "${#BAMS[@]}" -eq 1 ]; then
  SAMPLE=$(basename "${BAMS[0]}" .sorted.bam)
fi

# --- locate engine ------------------------------------------------------------
if [ "$ENGINE" = "gatk" ]; then
  if ! command -v gatk >/dev/null 2>&1; then
    echo "GATK (gatk4) not found."
    echo ""
    echo "Install options:"
    echo "  1. conda install -c bioconda gatk4        (recommended; bundles Java)"
    echo "  2. https://gatk.broadinstitute.org  (manual, needs Java 17+)"
    echo ""
    echo "Or re-run with the lightweight alternative engine:"
    echo "  bash scripts/variants/02_gatk_call.sh --engine bcftools ..."
    exit 0
  fi
  echo "Using: $(gatk --version 2>&1 | sed -n '2p')"
elif [ "$ENGINE" = "bcftools" ]; then
  if ! command -v bcftools >/dev/null 2>&1; then
    echo "bcftools not found."
    echo ""
    echo "Install: conda install -c bioconda bcftools samtools"
    exit 0
  fi
  echo "Using: $(bcftools --version | head -1)"
else
  echo "ERROR: --engine must be 'gatk' or 'bcftools'"; exit 1
fi

mkdir -p "$OUTDIR"
LOG="$OUTDIR/call_log.txt"
: > "$LOG"
run() { echo "+ $*" | tee -a "$LOG"; "$@" >>"$LOG" 2>&1; }

# --- reference .fai -----------------------------------------------------------
if [ ! -f "${REF}.fai" ]; then
  if command -v samtools >/dev/null 2>&1; then
    echo "--> creating reference index (.fai)"
    run samtools faidx "$REF"
  else
    echo "ERROR: ${REF}.fai missing and samtools unavailable to build it."; exit 1
  fi
fi

# region list for joint genotyping
if [ -n "$INTERVALS" ]; then
  REGION_ARG="$INTERVALS"
else
  REGION_ARG=$(awk '{print $1}' "${REF}.fai" | paste -sd, -)
fi

# --- calling ------------------------------------------------------------------
if [ "$ENGINE" = "gatk" ]; then
  if [ "${#BAMS[@]}" -eq 1 ]; then
    GVCF="$OUTDIR/${SAMPLE}.g.vcf.gz"
    echo "--> single-sample HaplotypeCaller: $SAMPLE"
    if [ -n "$INTERVALS" ]; then
      run gatk --java-options "-Xmx$JAVA_MEM" HaplotypeCaller \
          -R "$REF" -I "${BAMS[0]}" -O "$GVCF" -ERC GVCF -L "$INTERVALS"
    else
      run gatk --java-options "-Xmx$JAVA_MEM" HaplotypeCaller \
          -R "$REF" -I "${BAMS[0]}" -O "$GVCF" -ERC GVCF
    fi
    run gatk --java-options "-Xmx$JAVA_MEM" GenotypeGVCFs \
        -R "$REF" -V "$GVCF" -O "$OUTDIR/cohort.vcf.gz"
  else
    echo "--> multi-sample joint calling (${#BAMS[@]} samples)"
    GVCFS=()
    for b in "${BAMS[@]}"; do
      s=$(basename "$b" .sorted.bam)
      g="$OUTDIR/${s}.g.vcf.gz"
      echo "    HaplotypeCaller: $s"
      if [ -n "$INTERVALS" ]; then
        run gatk --java-options "-Xmx$JAVA_MEM" HaplotypeCaller \
            -R "$REF" -I "$b" -O "$g" -ERC GVCF -L "$INTERVALS"
      else
        run gatk --java-options "-Xmx$JAVA_MEM" HaplotypeCaller \
            -R "$REF" -I "$b" -O "$g" -ERC GVCF
      fi
      GVCFS+=("$g")
    done
    DB="$OUTDIR/genomicsdb"
    MAP_ARGS=()
    for g in "${GVCFS[@]}"; do
      MAP_ARGS+=(--sample-name "$(basename "$g" .g.vcf.gz)" -V "$g")
    done
    run gatk --java-options "-Xmx$JAVA_MEM" GenomicsDBImport \
        -R "$REF" "${MAP_ARGS[@]}" --genomicsdb-workspace-path "$DB" -L "$REGION_ARG"
    run gatk --java-options "-Xmx$JAVA_MEM" GenotypeGVCFs \
        -R "$REF" -V "gendb://$DB" -O "$OUTDIR/cohort.vcf.gz"
  fi
else
  # ---- bcftools engine: per-sample call, then merge ----
  OUT_VCF="$OUTDIR/cohort.vcf.gz"
  RAWS=()
  for b in "${BAMS[@]}"; do
    s=$(basename "$b" .sorted.bam)
    echo "    bcftools mpileup/call: $s"
    run bash -c "bcftools mpileup -Ou -f '$REF' -a AD,DP -q 20 -Q 20 '$b' | bcftools call -mv -Oz -o '$OUTDIR/${s}.raw.vcf.gz'"
    RAWS+=("$OUTDIR/${s}.raw.vcf.gz")
  done
  if [ "${#RAWS[@]}" -eq 1 ]; then
    run cp "${RAWS[0]}" "$OUT_VCF"
    run bcftools index -t -f "$OUT_VCF"
  else
    echo "    merging samples (bcftools merge)"
    for r in "${RAWS[@]}"; do run bcftools index -t "$r"; done
    run bcftools merge -Oz -o "$OUT_VCF" "${RAWS[@]}"
    run bcftools index -t -f "$OUT_VCF"
  fi
fi

echo "--> variant counts"
if command -v bcftools >/dev/null 2>&1; then
  bcftools stats "$OUTDIR/cohort.vcf.gz" 2>/dev/null \
    | grep -E "^SN" | sed 's/^/  /' | tee -a "$LOG"
fi

echo ""
echo "Calling complete. Raw VCF: $OUTDIR/cohort.vcf.gz"
echo ""
echo "Next step:"
echo "  bash scripts/variants/03_vcf_filter.sh --vcf $OUTDIR/cohort.vcf.gz \\"
echo "      --ref $REF --outdir results/variants/03_filter"
