#!/usr/bin/env bash
###############################################################################
# 09_graph_visualize.sh -- pangenome graph visualization (odgi / vg)
#
# 输入：.gfa 或 .og
# 步骤：odgi build -> odgi sort -> odgi viz（全图）；可选 odgi extract 区域
#
# 默认只打印计划；--run 才执行。
#
# 用法：
#   bash scripts/pangenome/09_graph_visualize.sh --graph rice.gfa --outdir results/pangenome/09_viz
#   bash scripts/pangenome/09_graph_visualize.sh --demo --outdir results/pangenome/09_viz_demo
###############################################################################
set -euo pipefail

GRAPH=""; OUT="results/pangenome/09_viz"; REGION=""; RUN=0; DEMO=0; WIDTH=2000; HEIGHT=500
while [[ $# -gt 0 ]]; do
  case "$1" in
    --graph) GRAPH="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --region) REGION="$2"; shift 2;;
    --width) WIDTH="$2"; shift 2;;
    --height) HEIGHT="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then
  GRAPH="$OUT/demo.gfa"
  {
    echo "H	VN:Z:1.0"
    echo "S	1	ACGTACGTACGTACGTACGT"
    echo "S	2	TTTTACGTACGTACGTACGT"
    echo "P	ref	1,2	*"
  } > "$GRAPH"
  echo "[viz] demo graph -> $GRAPH"
fi
[[ -n "$GRAPH" && -f "$GRAPH" ]] || { echo "Need --graph or --demo" >&2; exit 1; }

PLAN="$OUT/tool_plan.txt"
{
  echo "# Graph visualization plan ($(date -Iseconds)) graph=$GRAPH"
  if [[ "$GRAPH" == *.gfa ]]; then
    echo "odgi build -g $GRAPH -o $OUT/graph.og"
    echo "odgi sort -i $OUT/graph.og -o $OUT/graph.sorted.og -Y"
    echo "odgi viz -i $OUT/graph.sorted.og -o $OUT/graph.png -x $WIDTH -y $HEIGHT"
  else
    echo "odgi viz -i $GRAPH -o $OUT/graph.png -x $WIDTH -y $HEIGHT"
  fi
  if [[ -n "$REGION" ]]; then
    echo "odgi extract -i $OUT/graph.sorted.og -o $OUT/region.og -r $REGION"
    echo "odgi viz -i $OUT/region.og -o $OUT/region.png -x $WIDTH -y $HEIGHT"
  fi
  echo "# optional: vg find -x graph.gbz -r chr1:1000-2000 | vg view -d - > region.gfa"
} > "$PLAN"
cat "$PLAN"

if [[ $RUN -eq 1 ]]; then
  command -v odgi >/dev/null 2>&1 || { echo "[viz] odgi not found" >&2; exit 1; }
  if [[ "$GRAPH" == *.gfa ]]; then
    odgi build -g "$GRAPH" -o "$OUT/graph.og" > "$OUT/build.log" 2>&1
    odgi sort -i "$OUT/graph.og" -o "$OUT/graph.sorted.og" -Y > "$OUT/sort.log" 2>&1
    odgi viz -i "$OUT/graph.sorted.og" -o "$OUT/graph.png" -x "$WIDTH" -y "$HEIGHT" > "$OUT/viz.log" 2>&1
  else
    odgi viz -i "$GRAPH" -o "$OUT/graph.png" -x "$WIDTH" -y "$HEIGHT" > "$OUT/viz.log" 2>&1
  fi
  echo "[viz] finished -> $OUT/graph.png"
fi
