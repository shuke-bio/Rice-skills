#!/usr/bin/env bash
#
# 01_pangenome_graph.sh -- build a pangenome graph from multiple assemblies
#
# Description:
#   Pangenome graph construction template runner supporting three standard
#   toolchains (all fetched/installed by the user — nothing bundled):
#     --tool minigraph         fast incremental graph (graph only, no VCF)
#     --tool minigraph-cactus  graph + genotypable VCF/GFA (recommended)
#     --tool pggb              all-vs-all alignment graph (no reference-first)
#   If the selected toolchain is not installed the script prints install
#   instructions and exits.
#
# Usage:
#   # 1) assemble your panel (e.g. 10-30 high-quality rice assemblies)
#   # 2) build the graph:
#   bash scripts/pangenome/01_pangenome_graph.sh \
#        --ref IRGSP-1.0.fa --assemblies assemblies/ \
#        --tool minigraph --outdir results/pangenome/01_graph --threads 32
#
# Inputs:
#   --ref         reference genome FASTA (stays the coordinate backbone)
#   --assemblies  directory of assembly FASTAs (*.fa / *.fasta / *.fa.gz),
#                 or a file list (one path per line)
#   --tool        minigraph | minigraph-cactus | pggb (default minigraph)
#   --outdir      output directory (default results/pangenome/01_graph)
#   --threads     CPU threads (default 16)
#   --graph-name  output name (default "rice_pangenome")
# Outputs:
#   minigraph mode:         <outdir>/<name>.gfa (+ graph stats)
#   minigraph-cactus mode:  <outdir>/<name>.{gfa,gbz,vcf.gz,snarls,...}
#   pggb mode:              <outdir>/<name>.smooth.final.gfa (+ community)
#   <outdir>/graph_log.txt  full command log
#
# Notes:
#   - Install (conda):
#       minigraph        : conda install -c bioconda minigraph
#       minigraph-cactus : conda install -c bioconda cactus
#       pggb             : conda install -c bioconda pggb
#     Or run via the official Docker/Singularity images (recommended for
#     cactus — see https://github.com/ComparativeGenomicsToolkit/cactus).
#   - Resources: minigraph on 20 rice assemblies ≈ hours / tens of GB RAM;
#     minigraph-cactus is 10-50x heavier (plan a large-memory server).
#   - Choosing assemblies: use chromosome-level assemblies of diverse
#     accessions (e.g. 3K-RGP representatives + your own materials).
#   - Next step: bash scripts/pangenome/02_sv_detect.sh (graph genotyping),
#     or use the cactus --vcf output directly with the Pan-GWAS scripts.

set -u

REF=""
ASSEMBLIES=""
TOOL="minigraph"
OUTDIR="results/pangenome/01_graph"
THREADS=16
GRAPH_NAME="rice_pangenome"

while [ $# -gt 0 ]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --assemblies) ASSEMBLIES="$2"; shift 2 ;;
    --tool) TOOL="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --graph-name) GRAPH_NAME="$2"; shift 2 ;;
    -h|--help) sed -n '2,44p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

[ -n "$REF" ] && [ -n "$ASSEMBLIES" ] || {
  echo "ERROR: --ref and --assemblies are required. See --help."; exit 1; }
[ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }

# collect assembly FASTAs
ASM_LIST=()
if [ -d "$ASSEMBLIES" ]; then
  for f in "$ASSEMBLIES"/*.fa "$ASSEMBLIES"/*.fasta "$ASSEMBLIES"/*.fa.gz "$ASSEMBLIES"/*.fasta.gz; do
    [ -f "$f" ] && ASM_LIST+=("$f")
  done
elif [ -f "$ASSEMBLIES" ]; then
  while IFS= read -r line; do
    [ -n "$line" ] && [ -f "$line" ] && ASM_LIST+=("$line")
  done < "$ASSEMBLIES"
fi
if [ "${#ASM_LIST[@]}" -lt 2 ]; then
  echo "ERROR: need at least 2 assembly FASTAs (found ${#ASM_LIST[@]})."
  echo "       Provide a directory of FASTAs or a file list."
  exit 1
fi
echo "Reference : $REF"
echo "Assemblies: ${#ASM_LIST[@]} genome(s)"
echo "Toolchain : $TOOL"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/graph_log.txt"
: > "$LOG"
run() { echo "+ $*" | tee -a "$LOG"; bash -c "$*" >>"$LOG" 2>&1; }

case "$TOOL" in
  minigraph)
    if ! command -v minigraph >/dev/null 2>&1; then
      echo "minigraph not found. Install: conda install -c bioconda minigraph"
      echo "Docs: https://github.com/lh3/minigraph"
      exit 0
    fi
    echo "--> building graph with minigraph"
    run "minigraph -cxggs -t $THREADS '$REF' ${ASM_LIST[*]} > '$OUTDIR/$GRAPH_NAME.gfa'"
    echo "--> graph stats"
    if command -v gfatools >/dev/null 2>&1; then
      run "gfatools stat '$OUTDIR/$GRAPH_NAME.gfa' > '$OUTDIR/$GRAPH_NAME.stats.txt'"
    fi
    echo ""
    echo "Graph: $OUTDIR/$GRAPH_NAME.gfa"
    echo "Tip: minigraph graphs carry structural variants; convert to rGFA/BED"
    echo "     with gfatools, or genotype with vg/minigraph tools:"
    echo "     https://github.com/lh3/minigraph#calling-variants-from-a-graph"
    ;;
  minigraph-cactus)
    if ! command -v cactus-pangenome >/dev/null 2>&1 && ! command -v cactus >/dev/null 2>&1; then
      echo "minigraph-cactus (cactus) not found."
      echo ""
      echo "Install options:"
      echo "  1. conda install -c bioconda cactus"
      echo "  2. Docker/Singularity (recommended): see"
      echo "     https://github.com/ComparativeGenomicsToolkit/cactus"
      exit 0
    fi
    SEQFILE="$OUTDIR/seqfile.txt"
    {
      i=0
      for f in "${ASM_LIST[@]}"; do
        i=$((i+1))
        printf 'asm%02d\t%s\n' "$i" "$f"
      done
    } > "$SEQFILE"
    echo "--> seqfile written: $SEQFILE"
    echo "--> building pangenome with minigraph-cactus (heavy; hours)"
    run "cactus-pangenome '$OUTDIR/jobstore' '$SEQFILE' \
        --outDir '$OUTDIR/out' --outName '$GRAPH_NAME' \
        --reference '$(basename "$REF")' --vcf --gfa --gbz --giraffe \
        --maxCores $THREADS"
    echo ""
    echo "Outputs in $OUTDIR/out :"
    ls "$OUTDIR/out" 2>/dev/null | sed 's/^/  /'
    echo ""
    echo "Key files: <name>.vcf.gz (genotyped variants incl. SVs),"
    echo "           <name>.gbz + <name>.giraffe (graph for read mapping)."
    ;;
  pggb)
    if ! command -v pggb >/dev/null 2>&1; then
      echo "pggb not found. Install: conda install -c bioconda pggb"
      echo "Docs: https://github.com/pangenome/pggb"
      exit 0
    fi
    SEQFILE="$OUTDIR/input.fa"
    : > "$SEQFILE"
    for f in "$REF" "${ASM_LIST[@]}"; do
      case "$f" in
        *.gz) zcat "$f" >> "$SEQFILE" ;;
        *) cat "$f" >> "$SEQFILE" ;;
      esac
    done
    N=$(( ${#ASM_LIST[@]} + 1 ))
    echo "--> building graph with pggb ($N sequences)"
    run "pggb -i '$SEQFILE' -o '$OUTDIR/pggb_out' -n $N -t $THREADS -p 90 -s 5000"
    echo ""
    echo "Graph: $OUTDIR/pggb_out/*.smooth.final.gfa"
    ;;
  *)
    echo "ERROR: --tool must be minigraph | minigraph-cactus | pggb"; exit 1 ;;
esac

echo ""
echo "Next step:"
echo "  bash scripts/pangenome/02_sv_detect.sh --mode asm --ref $REF \\"
echo "      --queries $ASSEMBLIES --outdir results/pangenome/02_sv"
