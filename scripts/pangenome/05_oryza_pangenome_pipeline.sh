#!/usr/bin/env bash
###############################################################################
# 05_oryza_pangenome_pipeline.sh -- genus-scale Oryza pan-genome workflow runner
#
# Title:    Rice-skills Oryza pan-genome pipeline (Snakemake stage runner)
# Description:
#   Drives the genus-scale Oryza pan-genome workflow
#   (github.com/LengFeng00/oryza-pangenome-pipeline, MIT) stage by stage:
#     01_assembly -> 02_annotation -> 03_phylogeny -> 04_pan_synteny ->
#     05_sv_detection -> 06_pan_gene_te_lncrna -> 07_new_genes ->
#     08_transcriptome -> 09_popgen
#   The upstream repository is a study pipeline (Snakefile + per-stage scripts),
#   not a packaged tool; this runner validates inputs, checks dependencies and
#   invokes Snakemake with the chosen stage(s).
#
# Usage:
#   bash scripts/pangenome/05_oryza_pangenome_pipeline.sh \
#       --repo /path/to/oryza-pangenome-pipeline \
#       --config config/config.yaml --stage 04_pan_synteny \
#       --outdir results/pangenome/05_pipeline
#
#   bash scripts/pangenome/05_oryza_pangenome_pipeline.sh --demo
#
# Inputs:
#   --repo      local clone of the upstream repository (MIT)
#   --config    Snakemake config YAML (see upstream config/)
#   --stage     one of the 9 stages, or "all" (default: all)
#   --dry-run   snakemake -n only
#   --cores     CPU cores (default 4)
#   --outdir    output directory
#   --demo      fixture validation + dependency check (SKIP when snakemake absent)
#
# Outputs:
#   run_info.txt, dependency_report.txt, and the upstream stage outputs under
#   --outdir (or the repo's configured working directory).
#
# Notes:
#   - Upstream is NOT vendored: clone it yourself
#     (git clone https://github.com/LengFeng00/oryza-pangenome-pipeline).
#   - Requires: snakemake >= 7, plus the stage-specific tools declared in
#     upstream envs/ (minimap2, hifiasm, BRAKER/EDTA, OrthoFinder, SyRI, ...).
#   - LICENSE: MIT (retain LICENSE + NOTICE if you redistribute any part).
###############################################################################
set -euo pipefail

REPO=""; CONFIG=""; STAGE="all"; CORES=4; OUTDIR="results/pangenome/05_pipeline"; DRY=0; DEMO=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) REPO="$2"; shift 2;;
    --config) CONFIG="$2"; shift 2;;
    --stage) STAGE="$2"; shift 2;;
    --cores) CORES="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --dry-run) DRY=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,50p' "$0"; exit 0;;
    *) echo "Unknown argument: $1" >&2; exit 1;;
  esac
done

mkdir -p "$OUTDIR"
REPORT="$OUTDIR/dependency_report.txt"
: > "$REPORT"
say() { echo "$@" | tee -a "$REPORT"; }

say "# Oryza pan-genome pipeline dependency report"
say "# generated: $(date '+%Y-%m-%d %H:%M:%S')"

# ---- dependency checks -----------------------------------------------------
HAS_SMK=0; command -v snakemake >/dev/null 2>&1 && HAS_SMK=1
say "snakemake: $([[ $HAS_SMK -eq 1 ]] && snakemake --version || echo MISSING)"
for t in git minimap2 samtools; do
  if command -v "$t" >/dev/null 2>&1; then say "$t: $(command -v $t)"; else say "$t: MISSING"; fi
done

# ---- demo mode -------------------------------------------------------------
if [[ $DEMO -eq 1 ]]; then
  echo "== DEMO: validating workflow layout =="
  TMP="$OUTDIR/_demo_fixture"; mkdir -p "$TMP"
  # synthetic fixture: a tiny config + two fake assemblies
  cat > "$TMP/config.demo.yaml" <<'YAML'
samples:
  - Nipponbare
  - Huazhan
stages:
  - 01_assembly
  - 04_pan_synteny
YAML
  printf '>chr1\nACGTACGTACGTACGT\n' > "$TMP/Nipponbare.fa"
  printf '>chr1\nACGTACGTTCGTACGT\n' > "$TMP/Huazhan.fa"
  echo "fixture: $TMP (2 assemblies, 16 bp each)"
  if [[ $HAS_SMK -eq 0 ]]; then
    echo "SKIP: snakemake not installed."
    echo "  conda install -c bioconda snakemake>=7   # or pip install snakemake"
    echo "  Then: bash $0 --repo <clone> --config <yaml> --dry-run"
    exit 0
  fi
  echo "snakemake detected; a real run still needs upstream repo + stage tools."
  echo "  bash $0 --repo <clone> --config <yaml> --stage 04_pan_synteny --dry-run"
  exit 0
fi

# ---- real run --------------------------------------------------------------
[[ -z "$REPO" || -z "$CONFIG" ]] && { echo "Required: --repo and --config (or --demo)" >&2; exit 1; }
[[ -d "$REPO" ]] || { echo "Repo not found: $REPO" >&2; exit 1; }
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
[[ $HAS_SMK -eq 1 ]] || { echo "snakemake is required. Install: pip install snakemake" >&2; exit 1; }

cd "$REPO"
EXTRA=""; [[ $DRY -eq 1 ]] && EXTRA="-n"
if [[ "$STAGE" == "all" ]]; then
  snakemake --configfile "$CONFIG" --cores "$CORES" $EXTRA
else
  snakemake --configfile "$CONFIG" --cores "$CORES" $EXTRA "$STAGE"
fi

{
  echo "script: scripts/pangenome/05_oryza_pangenome_pipeline.sh"
  echo "repo: $REPO"
  echo "config: $CONFIG"
  echo "stage: $STAGE"
  echo "cores: $CORES"
  echo "dry_run: $DRY"
} > "$OUTDIR/run_info.txt"
echo "done -> $OUTDIR"
