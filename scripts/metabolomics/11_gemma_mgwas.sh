#!/usr/bin/env bash
###############################################################################
# 11_gemma_mgwas.sh -- GEMMA AI-REML mGWAS (plan-first)
#
# 输入：PLINK bfile 前缀 + GEMMA 表型文件（一列，无表头）
# 默认只打印计划；--run 才执行。
#
# 用法：
#   bash scripts/metabolomics/11_gemma_mgwas.sh --bfile geno --pheno met1.txt --outdir results/metabolomics/11_gemma
#   bash scripts/metabolomics/11_gemma_mgwas.sh --demo --outdir results/metabolomics/11_gemma_demo
###############################################################################
set -euo pipefail
BF=""; PH=""; OUT="results/metabolomics/11_gemma"; THREADS=8; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile) BF="$2"; shift 2;;
    --pheno) PH="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then BF="$OUT/demo"; PH="$OUT/demo.pheno"; : > "$PH"; echo "[gemma] demo placeholders -> $OUT"; fi
[[ -n "$BF" && -n "$PH" ]] || { echo "Need --bfile and --pheno, or --demo" >&2; exit 1; }
PLAN="$OUT/tool_plan.txt"
{
  echo "# GEMMA AI-REML mGWAS plan ($(date -Iseconds))"
  echo "gemma -bfile $BF -p $PH -gk 1 -o kinship"
  echo "gemma -bfile $BF -p $PH -k $OUT/kinship.cXX.txt -lmm 2 -o mgwas"
  echo "# outputs: $OUT/mgwas.assoc.txt (AI-REML LMM)"
} > "$PLAN"
cat "$PLAN"
if [[ $RUN -eq 1 ]]; then
  command -v gemma >/dev/null 2>&1 || { echo "[gemma] not found" >&2; exit 1; }
  gemma -bfile "$BF" -p "$PH" -gk 1 -o kinship > "$OUT/kinship.log" 2>&1
  gemma -bfile "$BF" -p "$PH" -k "$OUT/kinship.cXX.txt" -lmm 2 -o mgwas > "$OUT/mgwas.log" 2>&1
  echo "[gemma] finished -> $OUT/mgwas.assoc.txt"
fi
