#!/usr/bin/env bash
###############################################################################
# 08_selection_scan.sh -- iHS / CLR selection scans + Dsuite introgression (plan-first)
#
# 输入：
#   --vcf    phased VCF（iHS/CLR 需要 phased genotypes）
#   --map    genetic map（chr pos cM）
#   --sets   Dsuite sets 文件（A/B/C 三分类）
#   --tree   物种/群体树（newick）
# 输出：selscan iHS/CLR、norm、Dsuite Dtrios/Fbranch 的命令计划与结果
#
# 默认只打印计划；--run 才执行。
#
# 用法：
#   bash scripts/popgen/08_selection_scan.sh --vcf phased.vcf.gz --map chr.map \
#       --sets sets.txt --tree tree.nwk --outdir results/popgen/08_selection
#   bash scripts/popgen/08_selection_scan.sh --demo --outdir results/popgen/08_selection_demo
###############################################################################
set -euo pipefail

VCF=""; MAP=""; SETS=""; TREE=""; OUT="results/popgen/08_selection"; THREADS=8; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2;;
    --map) MAP="$2"; shift 2;;
    --sets) SETS="$2"; shift 2;;
    --tree) TREE="$2"; shift 2;;
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
  VCF="$OUT/demo.vcf"; MAP="$OUT/demo.map"; SETS="$OUT/demo_sets.txt"; TREE="$OUT/demo.nwk"
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n' > "$VCF"
  printf 'chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0|1\t0|0\n' >> "$VCF"
  printf 'chr1\t100\t0.01\nchr1\t200\t0.02\n' > "$MAP"
  printf 'S1\tA\nS2\tB\n' > "$SETS"
  printf '((A:1,B:1):1);\n' > "$TREE"
  echo "[sel] demo placeholders -> $OUT"
fi
[[ -n "$VCF" && -n "$MAP" ]] || { echo "Need --vcf and --map, or --demo" >&2; exit 1; }

PLAN="$OUT/tool_plan.txt"
{
  echo "# Selection scan plan ($(date -Iseconds))"
  echo "# vcf=$VCF map=$MAP sets=$SETS tree=$TREE"
  echo "selscan --ihs --vcf $VCF --map $MAP --out $OUT/ihs --threads $THREADS"
  echo "norm --ihs --files $OUT/ihs.ihs.out --log"
  echo "selscan --clr --vcf $VCF --map $MAP --out $OUT/clr --threads $THREADS"
  echo "norm --clr --files $OUT/clr.clr.out --log"
  if [[ -n "$SETS" && -n "$TREE" ]]; then
    echo "Dsuite Dtrios $VCF $SETS -t $TREE -o $OUT/dsuite"
    echo "Dsuite Fbranch $TREE $OUT/dsuite_tree.txt > $OUT/fbranch.txt"
  else
    echo "# Dsuite: skipped (need --sets and --tree)"
  fi
  echo "# XP-EHH already available in scripts/popgen/07_xpehh.R"
} > "$PLAN"
cat "$PLAN"

if [[ $RUN -eq 1 ]]; then
  for t in selscan norm Dsuite; do
    command -v "$t" >/dev/null 2>&1 || echo "[sel] missing: $t" >&2
  done
  if command -v selscan >/dev/null 2>&1; then
    selscan --ihs --vcf "$VCF" --map "$MAP" --out "$OUT/ihs" --threads "$THREADS" > "$OUT/ihs.log" 2>&1 || echo "[sel] iHS failed" >&2
    selscan --clr --vcf "$VCF" --map "$MAP" --out "$OUT/clr" --threads "$THREADS" > "$OUT/clr.log" 2>&1 || echo "[sel] CLR failed" >&2
    if command -v norm >/dev/null 2>&1; then
      [[ -f "$OUT/ihs.ihs.out" ]] && norm --ihs --files "$OUT/ihs.ihs.out" --log > "$OUT/ihs.norm.log" 2>&1 || true
      [[ -f "$OUT/clr.clr.out" ]] && norm --clr --files "$OUT/clr.clr.out" --log > "$OUT/clr.norm.log" 2>&1 || true
    fi
  fi
  if [[ -n "$SETS" && -n "$TREE" ]] && command -v Dsuite >/dev/null 2>&1; then
    Dsuite Dtrios "$VCF" "$SETS" -t "$TREE" -o "$OUT/dsuite" > "$OUT/dsuite.log" 2>&1 || echo "[sel] Dsuite failed" >&2
  fi
  echo "[sel] finished; see $OUT"
fi
