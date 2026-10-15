#!/usr/bin/env bash
#
# 02_sv_detect.sh -- structural variant (SV) detection from assemblies
#
# Description:
#   Assembly-based SV detection with the minimap2 + SYRI standard route:
#     minimap2 (asm5/asm10) -> samtools sort -> syri -c
#     -> per-sample SV table (INS/DEL/INV/TRA/CPG/CPL/...)
#   Also documents the graph-based alternative (graph genotyping via vg
#   giraffe / minigraph-cactus VCF) for cohorts already carrying a graph.
#   If the required tools are not installed the script prints install
#   instructions and exits.
#
# Usage:
#   # single query assembly vs reference
#   bash scripts/pangenome/02_sv_detect.sh --ref IRGSP-1.0.fa \
#        --query sample_asm.fa --sample sampleA \
#        --outdir results/pangenome/02_sv --threads 16
#
#   # batch: a directory / list of query assemblies
#   bash scripts/pangenome/02_sv_detect.sh --ref IRGSP-1.0.fa \
#        --queries assemblies/ --outdir results/pangenome/02_sv
#
# Inputs:
#   --ref      reference genome FASTA
#   --query    single query assembly FASTA  (or use --queries)
#   --queries  directory of FASTA files or a file list
#   --sample   sample name (single-query mode; default from filename)
#   --outdir   output directory (default results/pangenome/02_sv)
#   --threads  threads for minimap2/samtools (default 16)
#   --preset   minimap2 preset: asm5 | asm10 | asm20 (default asm10)
# Outputs:
#   <outdir>/<sample>.sv.syri.out       SV table (per SV: type, coords, size)
#   <outdir>/<sample>.syri.summary      SV counts by type
#   <outdir>/<sample>.bam (+ .bai)      whole-genome alignment
#   <outdir>/sv_all.tsv                 merged SV table (all samples)
#   <outdir>/sv_log.txt                 full command log
#
# Notes:
#   - Install: conda install -c bioconda minimap2 syri samtools
#   - SYRI requires the alignment BAM sorted and indexed; both done here.
#   - asm5 = <1% divergence, asm10 = ~1-5%, asm20 = >5% (rice intraspecies
#     assemblies: asm10 is a safe default; same-species polished: asm5).
#   - Graph route (cohort with pangenome graph): genotype SVs directly from
#     the cactus VCF (already produced by 01_pangenome_graph.sh
#     --tool minigraph-cactus), then run scripts/pangenome/04_pan_gwas.R.
#   - Next step: bash scripts/pangenome/04_pan_gwas.R (SV/PAV association)

set -u

REF=""
QUERY=""
QUERIES=""
SAMPLE=""
OUTDIR="results/pangenome/02_sv"
THREADS=16
PRESET="asm10"

while [ $# -gt 0 ]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --query) QUERY="$2"; shift 2 ;;
    --queries) QUERIES="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --preset) PRESET="$2"; shift 2 ;;
    -h|--help) sed -n '2,44p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

[ -n "$REF" ] || { echo "ERROR: --ref is required. See --help."; exit 1; }
[ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
if [ -z "$QUERY" ] && [ -z "$QUERIES" ]; then
  echo "ERROR: provide --query <assembly.fa> or --queries <dir|list>."; exit 1
fi

# --- locate tools -------------------------------------------------------------
MISSING=""
command -v minimap2 >/dev/null 2>&1 || MISSING="$MISSING minimap2"
command -v samtools >/dev/null 2>&1 || MISSING="$MISSING samtools"
command -v syri >/dev/null 2>&1 || MISSING="$MISSING syri"
if [ -n "$MISSING" ]; then
  echo "Required tool(s) not found:$MISSING"
  echo ""
  echo "Install: conda install -c bioconda minimap2 syri samtools"
  echo "Docs  : https://github.com/sschmidtu/SYRI (syri), https://github.com/lh3/minimap2"
  echo ""
  echo "Alternative (graph route): if you built a pangenome with"
  echo "  scripts/pangenome/01_pangenome_graph.sh --tool minigraph-cactus"
  echo "the emitted VCF already contains genotyped SVs — skip this script and"
  echo "go straight to scripts/pangenome/04_pan_gwas.R."
  exit 0
fi

# collect query list
QRY_LIST=()
if [ -n "$QUERY" ]; then
  QRY_LIST=("$QUERY")
elif [ -d "$QUERIES" ]; then
  for f in "$QUERIES"/*.fa "$QUERIES"/*.fasta "$QUERIES"/*.fa.gz; do
    [ -f "$f" ] && QRY_LIST+=("$f")
  done
elif [ -f "$QUERIES" ]; then
  while IFS= read -r line; do [ -n "$line" ] && [ -f "$line" ] && QRY_LIST+=("$line"); done < "$QUERIES"
fi
[ "${#QRY_LIST[@]}" -ge 1 ] || { echo "ERROR: no query assemblies found."; exit 1; }
echo "Queries: ${#QRY_LIST[@]} assembly file(s); preset=$PRESET"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/sv_log.txt"
: > "$LOG"
run() { echo "+ $*" | tee -a "$LOG"; bash -c "$*" >>"$LOG" 2>&1; }

# reference index for syri (needs .fai and optionally minimap2 index)
[ -f "${REF}.fai" ] || run "samtools faidx '$REF'"

ALL_SV="$OUTDIR/sv_all.tsv"
: > "$ALL_SV"

for q in "${QRY_LIST[@]}"; do
  s="${SAMPLE:-$(basename "$q" | sed 's/\.\(fa\|fasta\|fa\.gz\)$//')}"
  [ "${#QRY_LIST[@]}" -gt 1 ] && s=$(basename "$q" | sed 's/\.\(fa\|fasta\|fa\.gz\)$//')
  echo "--> [$s] aligning to reference (minimap2 -x $PRESET)"
  run "minimap2 -ax $PRESET -t $THREADS '$REF' '$q' \
      | samtools sort -@ $THREADS -o '$OUTDIR/$s.bam' -"
  run "samtools index '$OUTDIR/$s.bam'"

  echo "--> [$s] syri structural-variant calling"
  run "cd '$OUTDIR' && syri -c '$s.bam' -r '$REF' -q '$q' -F B --prefix '$s.'"

  # collect per-sample SV tables
  SVOUT="$OUTDIR/${s}.sv.syri.out"
  if [ -f "$SVOUT" ]; then
    awk -v smp="$s" 'BEGIN{FS=OFS="\t"} NR==1{print "sample",$0} NR>1{print smp,$0}' \
      "$SVOUT" >> "$ALL_SV"
  fi
  SUM="$OUTDIR/${s}.syri.summary"
  [ -f "$SUM" ] && { echo "    SV summary:"; sed 's/^/      /' "$SUM" | head -12; }
done

echo ""
echo "SV calling complete."
echo "  per-sample: $OUTDIR/<sample>.sv.syri.out"
echo "  merged    : $ALL_SV"
echo ""
echo "Next step:"
echo "  1. build an SV/PAV matrix (presence-absence per marker per sample), e.g."
echo "     from the merged table or from the cactus VCF;"
echo "  2. association: Rscript scripts/pangenome/04_pan_gwas.R --geno sv_matrix.csv \\"
echo "         --snp-info sv_info.csv --phenotype phenotype.csv --trait <name>"
