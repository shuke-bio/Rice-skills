#!/usr/bin/env bash
#
# 01_bwa_align.sh -- FASTQ -> analysis-ready BAM (BWA-MEM + samtools)
#
# Description:
#   Read alignment front-end for the variant-calling module:
#     bwa index (if needed) -> bwa mem -> name-sort -> fixmate -> coord-sort
#     -> markdup -> index -> flagstat/idxstats QC report.
#   If BWA / samtools are not installed the script prints install
#   instructions and exits.
#
# Usage:
#   bash scripts/variants/01_bwa_align.sh \
#        --fastq1 R1.fq.gz --fastq2 R2.fq.gz \
#        --ref reference.fa --sample IR24 \
#        --outdir results/variants/01_align --threads 16
#
# Inputs:
#   --fastq1   R1 FASTQ (.fq / .fq.gz); with --fastq2 for paired-end
#   --fastq2   R2 FASTQ (optional; omit for single-end)
#   --ref      reference genome FASTA (index built automatically if missing)
#   --sample   sample name (used for RG, output file names)  [required]
#   --outdir   output directory (default results/variants/01_align)
#   --threads  CPU threads (default 8)
#   --skip-markdup  skip duplicate marking (faster, not recommended for WGS)
# Outputs:
#   <outdir>/<sample>.sorted.bam (+ .bai)   analysis-ready alignments
#   <outdir>/<sample>.flagstat.txt          mapping summary
#   <outdir>/<sample>.idxstats.txt          per-chromosome counts
#   <outdir>/<sample>.stats.txt             samtools stats (optional detail)
#   <outdir>/align_log.txt                  full command log
#
# Notes:
#   - Tools install (one option): conda install -c bioconda bwa samtools
#   - The reference FASTA index (bwa .bwt/.pac/...) is REUSED across samples;
#     indexing runs once and takes minutes for a ~400 Mb rice genome.
#   - Next step: bash scripts/variants/02_gatk_call.sh

set -u

FASTQ1=""
FASTQ2=""
REF=""
SAMPLE=""
OUTDIR="results/variants/01_align"
THREADS=8
MARKDUP=1

while [ $# -gt 0 ]; do
  case "$1" in
    --fastq1) FASTQ1="$2"; shift 2 ;;
    --fastq2) FASTQ2="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --skip-markdup) MARKDUP=0; shift ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

if [ -z "$FASTQ1" ] || [ -z "$REF" ] || [ -z "$SAMPLE" ]; then
  echo "ERROR: --fastq1, --ref and --sample are required. See --help."
  exit 1
fi
for f in "$FASTQ1" "$REF"; do
  [ -f "$f" ] || { echo "ERROR: file not found: $f"; exit 1; }
done
[ -n "$FASTQ2" ] && [ ! -f "$FASTQ2" ] && { echo "ERROR: file not found: $FASTQ2"; exit 1; }

# --- locate tools -------------------------------------------------------------
MISSING=""
command -v bwa >/dev/null 2>&1 || MISSING="$MISSING bwa"
command -v samtools >/dev/null 2>&1 || MISSING="$MISSING samtools"
if [ -n "$MISSING" ]; then
  echo "Required tool(s) not found:$MISSING"
  echo ""
  echo "Install options (Linux / WSL / server is recommended for alignment):"
  echo "  1. conda install -c bioconda bwa samtools"
  echo "  2. Prebuilt binaries: http://bio-bwa.sourceforge.net  /  https://www.htslib.org"
  echo ""
  echo "Note: alignment is compute-heavy (tens of GB RAM-disk); run it on a"
  echo "Linux server or WSL2. On Windows-native environments use conda +"
  echo "bioconda (WSL recommended)."
  exit 0
fi
echo "Using bwa     : $(bwa 2>&1 | awk '/^Version/{print; exit}')"
echo "Using samtools: $(samtools --version | head -1)"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/align_log.txt"
: > "$LOG"
run() { echo "+ $*" | tee -a "$LOG"; "$@" >>"$LOG" 2>&1; }

# --- reference index (reused across samples) ----------------------------------
if [ ! -f "${REF}.bwt" ]; then
  echo "--> building BWA index for $REF (one-off, takes minutes)"
  run bwa index "$REF"
else
  echo "--> reusing existing BWA index for $REF"
fi

# --- alignment (streamed) -----------------------------------------------------
RG=$(printf '@RG\tID:%s\tSM:%s\tPL:ILLUMINA' "$SAMPLE" "$SAMPLE")
NAMEBAM="$OUTDIR/${SAMPLE}.name.bam"
FIXBAM="$OUTDIR/${SAMPLE}.fixmate.bam"
POSBAM="$OUTDIR/${SAMPLE}.pos.bam"
FINALBAM="$OUTDIR/${SAMPLE}.sorted.bam"

echo "--> aligning $SAMPLE with bwa mem ($THREADS threads)"
if [ -n "$FASTQ2" ]; then
  run bash -c "bwa mem -t $THREADS -R '$RG' '$REF' '$FASTQ1' '$FASTQ2' | samtools sort -@ $THREADS -n -o '$NAMEBAM' -"
else
  run bash -c "bwa mem -t $THREADS -R '$RG' '$REF' '$FASTQ1' | samtools sort -@ $THREADS -n -o '$NAMEBAM' -"
fi

if [ "$MARKDUP" -eq 1 ]; then
  echo "--> duplicate marking (fixmate + markdup)"
  run samtools fixmate -@ "$THREADS" -m "$NAMEBAM" "$FIXBAM"
  rm -f "$NAMEBAM"
  run samtools sort -@ "$THREADS" -o "$POSBAM" "$FIXBAM"
  rm -f "$FIXBAM"
  run samtools markdup -@ "$THREADS" "$POSBAM" "$FINALBAM"
  rm -f "$POSBAM"
else
  run samtools sort -@ "$THREADS" -o "$FINALBAM" "$NAMEBAM"
  rm -f "$NAMEBAM"
fi

run samtools index "$FINALBAM"

# --- QC reports ---------------------------------------------------------------
echo "--> QC reports (flagstat / idxstats / stats)"
run bash -c "samtools flagstat '$FINALBAM' > '$OUTDIR/${SAMPLE}.flagstat.txt'"
run bash -c "samtools idxstats '$FINALBAM' > '$OUTDIR/${SAMPLE}.idxstats.txt'"
run bash -c "samtools stats '$FINALBAM' | grep ^SN > '$OUTDIR/${SAMPLE}.stats.txt'"

echo ""
echo "Alignment complete for $SAMPLE"
grep -E "mapped \(|properly paired \(" "$OUTDIR/${SAMPLE}.flagstat.txt" 2>/dev/null | sed 's/^/  /'
echo ""
echo "Next step:"
echo "  bash scripts/variants/02_gatk_call.sh --bam $FINALBAM --ref $REF \\"
echo "      --sample $SAMPLE --outdir results/variants/02_call"
