#!/usr/bin/env bash
###############################################################################
# 08_graph_genotype.sh -- graph pangenome genotyping (vg Giraffe / PanGenie)
#
# 两条路线：
#   vg       : vg autoindex (giraffe) -> vg giraffe -> vg pack -> vg call
#   pangenie : PanGenie -i graph.index -r reads -o out
#
# 说明：
#   - 默认只生成命令计划（tool_plan.txt），加 --run 才执行；
#   - 水稻图泛基因组推荐：minigraph-cactus 产出 .gfa/.gbz，再 autoindex；
#   - 低覆盖/GBS 数据可考虑 PHG（Practical Haplotype Graph），见 docs 注释。
#
# 用法：
#   bash scripts/pangenome/08_graph_genotype.sh --graph rice.gfa --reads sample.fq \
#       --sample S01 --outdir results/pangenome/08_graph_genotype
#   bash scripts/pangenome/08_graph_genotype.sh --demo --outdir results/pangenome/08_graph_genotype_demo
###############################################################################
set -euo pipefail

GRAPH=""; READS=""; SAMPLE="sample"; OUT="results/pangenome/08_graph_genotype"
THREADS=8; PANGENIE_INDEX=""; RUN=0; DEMO=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --graph) GRAPH="$2"; shift 2;;
    --reads) READS="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --pangenie-index) PANGENIE_INDEX="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

mkdir -p "$OUT"

if [[ $DEMO -eq 1 ]]; then
  GRAPH="$OUT/demo_graph.gfa"
  READS="$OUT/demo_reads.fq"
  {
    echo "H	VN:Z:1.0"
    echo "S	1	ACGTACGTACGTACGTACGT"
    echo "S	2	TTTTACGTACGTACGTACGT"
    echo "P	ref	1,2	*"
  } > "$GRAPH"
  printf '@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n' > "$READS"
  echo "[graph] demo graph/reads -> $OUT"
fi

[[ -n "$GRAPH" && -f "$GRAPH" ]] || { echo "Need --graph or --demo" >&2; exit 1; }
if [[ -z "$READS" || ! -f "$READS" ]]; then echo "Need --reads (FASTQ/FASTA) or --demo" >&2; exit 1; fi

PREFIX="$OUT/graph"
PLAN="$OUT/tool_plan.txt"
{
  echo "# Graph genotyping plan ($(date -Iseconds))"
  echo "# graph: $GRAPH"
  echo "# reads: $READS"
  echo "# sample: $SAMPLE"
  if [[ "$GRAPH" == *.gfa ]]; then
    echo "vg autoindex --workflow giraffe -g $GRAPH -p $PREFIX -t $THREADS"
  else
    echo "# assume indexed graph: $GRAPH (gbz) or run vg autoindex manually"
  fi
  echo "vg giraffe -Z $PREFIX.gbz -f $READS -o $OUT/$SAMPLE.gam -t $THREADS"
  echo "vg pack -x $PREFIX.gbz -g $OUT/$SAMPLE.gam -o $OUT/$SAMPLE.pack"
  echo "vg call $PREFIX.gbz -k $OUT/$SAMPLE.pack -s $SAMPLE -o $OUT/$SAMPLE.vcf"
  if [[ -n "$PANGENIE_INDEX" ]]; then
    echo "# PanGenie alternative:"
    echo "PanGenie -i $PANGENIE_INDEX -r $READS -o $OUT/pangenie_$SAMPLE -t $THREADS"
  else
    echo "# PanGenie: skipped (need --pangenie-index)"
  fi
} > "$PLAN"
cat "$PLAN"

if [[ $RUN -eq 1 ]]; then
  if ! command -v vg >/dev/null 2>&1; then
    echo "[graph] vg not found; install vg (https://github.com/vgteam/vg) and re-run" >&2
    exit 1
  fi
  if [[ "$GRAPH" == *.gfa ]]; then
    vg autoindex --workflow giraffe -g "$GRAPH" -p "$PREFIX" -t "$THREADS" > "$OUT/autoindex.log" 2>&1
  fi
  vg giraffe -Z "$PREFIX.gbz" -f "$READS" -o "$OUT/$SAMPLE.gam" -t "$THREADS" > "$OUT/giraffe.log" 2>&1
  vg pack -x "$PREFIX.gbz" -g "$OUT/$SAMPLE.gam" -o "$OUT/$SAMPLE.pack" > "$OUT/pack.log" 2>&1
  vg call "$PREFIX.gbz" -k "$OUT/$SAMPLE.pack" -s "$SAMPLE" -o "$OUT/$SAMPLE.vcf" > "$OUT/call.log" 2>&1
  echo "[graph] vg pipeline finished -> $OUT/$SAMPLE.vcf"
  if [[ -n "$PANGENIE_INDEX" ]] && command -v PanGenie >/dev/null 2>&1; then
    PanGenie -i "$PANGENIE_INDEX" -r "$READS" -o "$OUT/pangenie_$SAMPLE" -t "$THREADS" \
      > "$OUT/pangenie.log" 2>&1 || echo "[graph] PanGenie failed; see pangenie.log" >&2
  fi
fi
