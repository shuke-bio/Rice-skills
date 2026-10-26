#!/usr/bin/env bash
###############################################################################
# 10_phg.sh -- Practical Haplotype Graph (PHG) plan/runner
#
# 说明：PHG 的 CLI 版本差异较大；本脚本生成标准流程命令，需按本机 PHG 版本调整。
#
# 用法：
#   bash scripts/pangenome/10_phg.sh --config phg_config.txt --outdir results/pangenome/10_phg
#   bash scripts/pangenome/10_phg.sh --demo --outdir results/pangenome/10_phg_demo
###############################################################################
set -euo pipefail
CFG=""; OUT="results/pangenome/10_phg"; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CFG="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,30p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then CFG="$OUT/demo_config.txt"; echo "# PHG config" > "$CFG"; fi
[[ -n "$CFG" ]] || { echo "Need --config or --demo" >&2; exit 1; }
PLAN="$OUT/tool_plan.txt"
{
  echo "# PHG plan ($(date -Iseconds))"
  echo "phg create -c $CFG"
  echo "phg add -c $CFG"
  echo "phg build -c $CFG"
  echo "phg map -c $CFG"
  echo "# 若使用 Docker: docker run --rm -v \$PWD:/work -w /work phg:latest phg create -c $CFG"
} > "$PLAN"
cat "$PLAN"
if [[ $RUN -eq 1 ]]; then
  command -v phg >/dev/null 2>&1 || { echo "[phg] not found" >&2; exit 1; }
  phg create -c "$CFG" > "$OUT/create.log" 2>&1
  phg add -c "$CFG" > "$OUT/add.log" 2>&1
  phg build -c "$CFG" > "$OUT/build.log" 2>&1
  phg map -c "$CFG" > "$OUT/map.log" 2>&1
  echo "[phg] finished -> $OUT"
fi
