#!/usr/bin/env bash
###############################################################################
# 03_hic_analysis.sh -- Hi-C / 3D genome plan-first pipeline
#
# 输入：.cool / .mcool / .hic 矩阵
# 步骤：平衡 -> zoomify -> TAD -> loops（HiCExplorer/cooler）
#
# 默认只打印计划；--run 才执行。
#
# 用法：
#   bash scripts/epigenomics/03_hic_analysis.sh --matrix sample.cool --resolution 10000 \
#       --outdir results/epigenomics/03_hic
#   bash scripts/epigenomics/03_hic_analysis.sh --demo --outdir results/epigenomics/03_hic_demo
###############################################################################
set -euo pipefail

MAT=""; RES=10000; OUT="results/epigenomics/03_hic"; RUN=0; DEMO=0; THREADS=8
while [[ $# -gt 0 ]]; do
  case "$1" in
    --matrix) MAT="$2"; shift 2;;
    --resolution) RES="$2"; shift 2;;
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
  MAT="$OUT/demo.cool"; : > "$MAT"; echo "[hic] demo placeholder -> $MAT"
fi
[[ -n "$MAT" && -f "$MAT" ]] || { echo "Need --matrix or --demo" >&2; exit 1; }

PLAN="$OUT/tool_plan.txt"
{
  echo "# Hi-C plan ($(date -Iseconds)) matrix=$MAT resolution=$RES"
  echo "cooler balance $MAT -o $OUT/balanced.cool"
  echo "cooler zoomify $OUT/balanced.cool -r 1000,5000,10000,50000,100000 -o $OUT/multi.cool"
  echo "hicConvertFormat -m $OUT/balanced.cool --inputFormat cool --outputFormat h5 -o $OUT/hic.h5 --resolutions $RES"
  echo "hicFindTADs -m $OUT/hic.h5 --outPrefix $OUT/tads --correctForMultipleTesting fdr --numberOfProcessors $THREADS"
  echo "hicDetectLoops -m $OUT/hic.h5 -o $OUT/loops.bedgraph --resolutions $RES --numberOfProcessors $THREADS"
  echo "# optional: Juicer A/B compartments, cooltools insulation"
} > "$PLAN"
cat "$PLAN"

if [[ $RUN -eq 1 ]]; then
  for t in cooler hicConvertFormat hicFindTADs hicDetectLoops; do
    command -v "$t" >/dev/null 2>&1 || echo "[hic] missing: $t" >&2
  done
  if command -v cooler >/dev/null 2>&1; then
    cooler balance "$MAT" -o "$OUT/balanced.cool" > "$OUT/balance.log" 2>&1 || echo "[hic] balance failed" >&2
    cooler zoomify "$OUT/balanced.cool" -r 1000,5000,10000,50000,100000 -o "$OUT/multi.cool" > "$OUT/zoomify.log" 2>&1 || true
  fi
  if command -v hicFindTADs >/dev/null 2>&1 && [[ -f "$OUT/hic.h5" ]]; then
    hicFindTADs -m "$OUT/hic.h5" --outPrefix "$OUT/tads" --correctForMultipleTesting fdr --numberOfProcessors "$THREADS" \
      > "$OUT/tads.log" 2>&1 || echo "[hic] TAD failed" >&2
  fi
  echo "[hic] finished; see $OUT"
fi
