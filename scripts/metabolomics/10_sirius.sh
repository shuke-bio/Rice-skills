#!/usr/bin/env bash
###############################################################################
# 10_sirius.sh -- SIRIUS formula/structure annotation (plan-first)
#
# 输入：mzML/mzXML 目录或单个文件
# 默认只打印计划；--run 才执行。GNPS 需网页上传，脚本只生成上传说明。
#
# 用法：
#   bash scripts/metabolomics/10_sirius.sh --mzml data/mzml --outdir results/metabolomics/10_sirius
#   bash scripts/metabolomics/10_sirius.sh --demo --outdir results/metabolomics/10_sirius_demo
###############################################################################
set -euo pipefail
MZ=""; OUT="results/metabolomics/10_sirius"; PPM=10; CAND=5; THREADS=8; RUN=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mzml) MZ="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --ppm) PPM="$2"; shift 2;;
    --candidates) CAND="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --run) RUN=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUT"
if [[ $DEMO -eq 1 ]]; then MZ="$OUT/demo.mzML"; : > "$MZ"; echo "[sirius] demo placeholder -> $MZ"; fi
[[ -n "$MZ" && -e "$MZ" ]] || { echo "Need --mzml or --demo" >&2; exit 1; }
PLAN="$OUT/tool_plan.txt"
{
  echo "# SIRIUS plan ($(date -Iseconds))"
  echo "sirius -i $MZ -o $OUT/sirius --threads $THREADS formula --ppm-max $PPM --candidates $CAND"
  echo "# optional structure: sirius ... structure --candidates $CAND"
  echo "# GNPS: convert to MGF and upload to https://gnps.ucsd.edu (FBMN); see $OUT/gnps_upload.md"
} > "$PLAN"
cat "$PLAN"
cat > "$OUT/gnps_upload.md" << 'MD'
# GNPS 上传步骤

1. 将 mzML 转为 MGF（MSConvert / msconvert）；
2. 登录 https://gnps.ucsd.edu，选择 Classical Molecular Networking 或 FBMN；
3. 上传 MGF 与 metadata（sample, group）；
4. 任务完成后下载 `network.tsv`、`library.tsv`，用 `metabolomics/09_annotate.py` 合并。
MD
if [[ $RUN -eq 1 ]]; then
  command -v sirius >/dev/null 2>&1 || { echo "[sirius] not found" >&2; exit 1; }
  sirius -i "$MZ" -o "$OUT/sirius" --threads "$THREADS" formula --ppm-max "$PPM" --candidates "$CAND" \
    > "$OUT/sirius.log" 2>&1
  echo "[sirius] finished -> $OUT/sirius"
fi
