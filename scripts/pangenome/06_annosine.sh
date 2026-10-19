#!/usr/bin/env bash
###############################################################################
# 06_annosine.sh -- SINE annotation for plant genomes (AnnoSINE)
#
# Title:    Rice-skills SINE annotation (AnnoSINE, MIT)
# Description:
#   Annotates short interspersed nuclear elements (SINEs) in plant genomes.
#   AnnoSINE builds a non-redundant SINE library from a genome assembly using
#   homology and/or structure evidence, then annotates the whole genome.
#   Upstream: github.com/yangli557/AnnoSINE (MIT), Yang et al., Plant Physiol
#   2022 "AnnoSINE: a short interspersed nuclear elements annotation tool for
#   plant genomes".
#
# Usage:
#   bash scripts/pangenome/06_annosine.sh \
#       --repo /path/to/AnnoSINE --genome assembly.fa --mode 3 \
#       --outdir results/pangenome/06_annosine
#
#   bash scripts/pangenome/06_annosine.sh --demo
#
# Inputs:
#   --repo     local clone of AnnoSINE (MIT)
#   --genome   genome assembly FASTA
#   --mode     1 homology | 2 structure | 3 hybrid (default 3)
#   --outdir   output directory
#   --demo     synthetic FASTA fixture + dependency check
#
# Outputs:
#   Seed_SINE.fa (non-redundant library), <genome>.out (annotations),
#   intermediate Step*.fa/txt files, run_info.txt
#
# Notes:
#   - Upstream NOT vendored: git clone https://github.com/yangli557/AnnoSINE
#   - Requires: python3, blastn (BLAST+), TRF, node.js (TSD_Searcher.js),
#     and the upstream Family_Seq / RNA database resources shipped in the repo.
#   - Runtime on a 400 Mb rice genome is hours; start with a chromosome.
###############################################################################
set -euo pipefail

REPO=""; GENOME=""; MODE=3; OUTDIR="results/pangenome/06_annosine"; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) REPO="$2"; shift 2;;
    --genome) GENOME="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown argument: $1" >&2; exit 1;;
  esac
done
mkdir -p "$OUTDIR"

check() { command -v "$1" >/dev/null 2>&1 && echo "$1: $(command -v "$1")" || echo "$1: MISSING"; }
{
  echo "# AnnoSINE dependency report"
  echo "# generated: $(date '+%Y-%m-%d %H:%M:%S')"
  check python3; check blastn; check trf; check node
} | tee "$OUTDIR/dependency_report.txt"

HAS_DEPS=1
for t in python3 blastn trf node; do command -v "$t" >/dev/null 2>&1 || HAS_DEPS=0; done

if [[ $DEMO -eq 1 ]]; then
  echo "== DEMO: synthetic repeat-rich FASTA =="
  FIX="$OUTDIR/_demo_fixture"; mkdir -p "$FIX"
  python3 - "$FIX/demo.fa" <<'PY'
import random, sys
random.seed(42)
unit = "GGTTCACGTACGTACCTTG"
seq = []
for i in range(60):
    seq.append("".join(random.choice("ACGT") for _ in range(400)))
    seq.append(unit * random.randint(2, 5))
open(sys.argv[1], "w").write(">demo_chr\n" + "".join(seq) + "\n")
PY
  echo "fixture: $FIX/demo.fa ($(wc -c < "$FIX/demo.fa") bytes, planted 20 bp repeats)"
  if [[ $HAS_DEPS -eq 0 ]]; then
    echo "SKIP: missing external tools (see dependency_report.txt)."
    echo "  conda install -c bioconda blast trf"
    echo "  node.js: https://nodejs.org (needed by TSD_Searcher.js)"
    exit 0
  fi
  echo "Dependencies present. Real run needs a cloned repo with Family_Seq/RNA DB:"
  echo "  bash $0 --repo <clone> --genome $FIX/demo.fa --mode 3"
  exit 0
fi

[[ -z "$REPO" || -z "$GENOME" ]] && { echo "Required: --repo and --genome (or --demo)" >&2; exit 1; }
[[ -d "$REPO" ]] || { echo "Repo not found: $REPO" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "Genome not found: $GENOME" >&2; exit 1; }
[[ $HAS_DEPS -eq 1 ]] || { echo "External tools missing; see $OUTDIR/dependency_report.txt" >&2; exit 1; }

cd "$REPO/bin"
python3 AnnoSINE.py "$MODE" "$GENOME" "$OUTDIR/" 2>&1 | tee "$OUTDIR/annosine.log"
mv -f ./*.out ./*.fa "$OUTDIR/" 2>/dev/null || true
{
  echo "script: scripts/pangenome/06_annosine.sh"
  echo "repo: $REPO"; echo "genome: $GENOME"; echo "mode: $MODE"
} > "$OUTDIR/run_info.txt"
echo "done -> $OUTDIR"
