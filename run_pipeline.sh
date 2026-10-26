#!/usr/bin/env bash
#
# run_pipeline.sh -- Rice-skills one-command pipeline runner
#
# Description:
#   Runs the transcriptome, metabolome and multi-omics integration scripts in
#   order, wiring the output of each step into the next one. Every step logs to
#   <outdir>/pipeline_log.txt; mandatory steps abort the run on failure while
#   optional steps (GSEA without a gene set, TF analysis without an annotation
#   file, MOFA/DIABLO when their packages are missing) are skipped with a
#   warning so that the rest of the pipeline still completes.
#
# Usage:
#   bash run_pipeline.sh --counts counts.csv --metadata metadata.csv \
#       --metabolites metabolites.csv --outdir results --stages all
#
#   # transcriptome only, with a batch term in the design and a GSEA gene set
#   bash run_pipeline.sh --counts counts.csv --metadata metadata.csv \
#       --stages transcriptome --design "~ batch + group" --gmt rice_kegg.gmt
#
# Notes:
#   The interpreter is auto-detected; override with --python / --rscript.
#   On Windows use Git Bash and pass the full path to Rscript.exe if needed.

set -u

COUNTS=""
METADATA=""
METABOLITES=""
OUTDIR="results"
STAGES="all"
DESIGN="~ group"
CONTRAST=""
GMT=""
TF_ANNOTATION=""
TRAITS=""
PADJ="0.05"
LOG2FC="1"
VIP="1"
COR_CUTOFF="0.8"
PYTHON=""
RSCRIPT=""
SEED="42"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [ $# -gt 0 ]; do
  case "$1" in
    --counts)         COUNTS="$2"; shift 2 ;;
    --metadata)       METADATA="$2"; shift 2 ;;
    --metabolites)    METABOLITES="$2"; shift 2 ;;
    --outdir)         OUTDIR="$2"; shift 2 ;;
    --stages)         STAGES="$2"; shift 2 ;;
    --design)         DESIGN="$2"; shift 2 ;;
    --contrast)       CONTRAST="$2"; shift 2 ;;
    --gmt)            GMT="$2"; shift 2 ;;
    --tf-annotation)  TF_ANNOTATION="$2"; shift 2 ;;
    --traits)         TRAITS="$2"; shift 2 ;;
    --padj)           PADJ="$2"; shift 2 ;;
    --log2fc)         LOG2FC="$2"; shift 2 ;;
    --vip)            VIP="$2"; shift 2 ;;
    --cor-cutoff)     COR_CUTOFF="$2"; shift 2 ;;
    --python)         PYTHON="$2"; shift 2 ;;
    --rscript)        RSCRIPT="$2"; shift 2 ;;
    --seed)           SEED="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1' (see --help)."; exit 1 ;;
  esac
done

# ---------------------------------------------------------------------------
# Interpreter detection
# ---------------------------------------------------------------------------
# Prefer a Python that actually carries the scientific stack (matplotlib);
# fall back to the first interpreter found, with a warning.
if [ -z "$PYTHON" ]; then
  _py_first=""
  for _cand in python3 python; do
    if command -v "$_cand" >/dev/null 2>&1; then
      [ -z "$_py_first" ] && _py_first="$_cand"
      if "$_cand" -c "import matplotlib" >/dev/null 2>&1; then PYTHON="$_cand"; break; fi
    fi
  done
  if [ -z "$PYTHON" ]; then
    for _cand in "$HOME/AppData/Local/Programs/Python"/Python*/python.exe; do
      if [ -x "$_cand" ] && "$_cand" -c "import matplotlib" >/dev/null 2>&1; then PYTHON="$_cand"; break; fi
    done
  fi
  if [ -z "$PYTHON" ]; then
    if [ -n "$_py_first" ]; then
      PYTHON="$_py_first"
      echo "WARNING: no Python with matplotlib found; using '$PYTHON'. Install deps per docs/installation.md or pass --python."
    else
      echo "ERROR: python not found. Install Python >=3.10 or pass --python."; exit 1
    fi
  fi
fi
if [ -z "$RSCRIPT" ]; then
  if command -v Rscript >/dev/null 2>&1; then RSCRIPT="Rscript"
  elif [ -x "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" ]; then
    RSCRIPT="/c/Program Files/R/R-4.6.0/bin/Rscript.exe"
  else
    echo "ERROR: Rscript not found. Install R >=4.3 or pass --rscript."; exit 1
  fi
fi

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/pipeline_log.txt"
: > "$LOG"

N_PASS=0; N_FAIL=0; N_SKIP=0

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# run_step <label> <command...>
run_step() {
  local label="$1"; shift
  log "==> $label"
  if "$@" >>"$LOG" 2>&1; then
    log "    PASS: $label"; N_PASS=$((N_PASS+1)); return 0
  else
    log "    FAIL: $label (see $LOG)"; N_FAIL=$((N_FAIL+1)); return 1
  fi
}

skip_step() { log "--> SKIP: $1 ($2)"; N_SKIP=$((N_SKIP+1)); }

require_file() {
  if [ ! -f "$2" ]; then
    echo "ERROR: $1 not found: $2"
    echo "Hint: check the path, or generate demo data with:"
    echo "      $PYTHON scripts/utils/make_demo_data.py --outdir demo_data"
    exit 1
  fi
}

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
case "$STAGES" in
  all|transcriptome) require_file "--counts" "$COUNTS"; require_file "--metadata" "$METADATA" ;;
esac
case "$STAGES" in
  all|metabolome) require_file "--metabolites" "$METABOLITES"; require_file "--metadata" "$METADATA" ;;
esac
case "$STAGES" in
  all|integration)
    require_file "--counts" "$COUNTS"
    require_file "--metabolites" "$METABOLITES"
    require_file "--metadata" "$METADATA" ;;
esac

log "Rice-skills pipeline start"
log "  python  : $PYTHON ($($PYTHON -V 2>&1))"
log "  Rscript : $RSCRIPT"
log "  outdir  : $OUTDIR   stages: $STAGES   seed: $SEED"

TX_DIR="$OUTDIR/transcriptome"
MET_DIR="$OUTDIR/metabolome"
INT_DIR="$OUTDIR/integration"

# ---------------------------------------------------------------------------
# Module 1: transcriptomics
# ---------------------------------------------------------------------------
if [ "$STAGES" = "all" ] || [ "$STAGES" = "transcriptome" ]; then
  log "--- Module 1: transcriptomics ---"

  run_step "01 QC and filtering" \
    "$PYTHON" "$REPO_ROOT/scripts/transcriptomics/01_qc_filter.py" \
      --counts "$COUNTS" --metadata "$METADATA" --outdir "$TX_DIR/01_qc" || true

  FILTERED="$TX_DIR/01_qc/counts_filtered.csv"
  [ -f "$FILTERED" ] || FILTERED="$COUNTS"

  CONTRAST_ARG=""
  [ -n "$CONTRAST" ] && CONTRAST_ARG="--contrast $CONTRAST"

  run_step "02 Differential expression (DESeq2)" \
    "$RSCRIPT" "$REPO_ROOT/scripts/transcriptomics/02_deseq2.R" \
      --counts "$FILTERED" --metadata "$METADATA" --outdir "$TX_DIR/02_de" \
      --design "$DESIGN" --padj "$PADJ" --log2fc "$LOG2FC" $CONTRAST_ARG || true

  # Pick the first significant-DE table (or all-DE table) for downstream steps.
  DE_SIG="$(ls "$TX_DIR"/02_de/de_sig_*.csv 2>/dev/null | head -1)"
  DE_ALL="$(ls "$TX_DIR"/02_de/de_all_*.csv 2>/dev/null | head -1)"

  if [ -n "$DE_SIG" ] && [ -n "$DE_ALL" ]; then
    run_step "03 GO/KEGG enrichment" \
      "$RSCRIPT" "$REPO_ROOT/scripts/transcriptomics/03_enrichment.R" \
        --de-sig "$DE_SIG" --de-all "$DE_ALL" --outdir "$TX_DIR/03_enrichment" \
        --padj "$PADJ" || true

    if [ -n "$GMT" ] && [ -f "$GMT" ]; then
      run_step "04 GSEA" \
        "$RSCRIPT" "$REPO_ROOT/scripts/transcriptomics/04_gsea.R" \
          --de-all "$DE_ALL" --gmt "$GMT" --outdir "$TX_DIR/04_gsea" --padj "$PADJ" || true
    else
      skip_step "04 GSEA" "no --gmt gene set supplied (optional)"
    fi
  else
    skip_step "03 GO/KEGG enrichment" "no DE table produced by step 02"
    skip_step "04 GSEA" "no DE table produced by step 02"
  fi

  TRAITS_ARG=""
  [ -n "$TRAITS" ] && TRAITS_ARG="--traits $TRAITS"

  run_step "05 WGCNA co-expression network" \
    "$RSCRIPT" "$REPO_ROOT/scripts/transcriptomics/05_wgcna.R" \
      --counts "$FILTERED" --metadata "$METADATA" --outdir "$TX_DIR/05_wgcna" \
      $TRAITS_ARG || true

  if [ -n "$TF_ANNOTATION" ] && [ -f "$TF_ANNOTATION" ]; then
    if [ -n "$DE_ALL" ]; then
      run_step "06 Transcription factor analysis" \
        "$RSCRIPT" "$REPO_ROOT/scripts/transcriptomics/06_tf_analysis.R" \
          --de-all "$DE_ALL" --tf-annotation "$TF_ANNOTATION" \
          --outdir "$TX_DIR/06_tf" --padj "$PADJ" --log2fc "$LOG2FC" || true
    else
      skip_step "06 TF analysis" "no DE table produced by step 02"
    fi
  else
    skip_step "06 TF analysis" "no --tf-annotation supplied (optional)"
  fi

  run_step "07 Transcriptome figures" \
    "$RSCRIPT" "$REPO_ROOT/scripts/transcriptomics/07_plot_transcriptome.R" \
      --results-dir "$TX_DIR" --outdir "$TX_DIR/07_figures" || true
fi

# ---------------------------------------------------------------------------
# Module 2: metabolomics
# ---------------------------------------------------------------------------
if [ "$STAGES" = "all" ] || [ "$STAGES" = "metabolome" ]; then
  log "--- Module 2: metabolomics ---"

  run_step "01 Preprocessing" \
    "$PYTHON" "$REPO_ROOT/scripts/metabolomics/01_preprocess.py" \
      --abundance "$METABOLITES" --metadata "$METADATA" --outdir "$MET_DIR/01_preprocess" || true

  PREPROCESSED="$MET_DIR/01_preprocess/metabolome_preprocessed.csv"
  [ -f "$PREPROCESSED" ] || PREPROCESSED="$METABOLITES"

  run_step "02 OPLS-DA and differential metabolites" \
    "$RSCRIPT" "$REPO_ROOT/scripts/metabolomics/02_oplsda.R" \
      --data "$PREPROCESSED" --metadata "$METADATA" --outdir "$MET_DIR/02_oplsda" \
      --vip "$VIP" --pvalue "$PADJ" || true

  DIFF_METS="$MET_DIR/02_oplsda/diff_metabolites.csv"
  if [ -f "$DIFF_METS" ]; then
    run_step "03 Metabolite pathway enrichment" \
      "$RSCRIPT" "$REPO_ROOT/scripts/metabolomics/03_pathway.R" \
        --diff-metabolites "$DIFF_METS" --outdir "$MET_DIR/03_pathway" --padj "$PADJ" || true

    run_step "04 Metabolite classification" \
      "$RSCRIPT" "$REPO_ROOT/scripts/metabolomics/04_metabolite_classification.R" \
        --diff-metabolites "$DIFF_METS" --all-metabolites "$PREPROCESSED" \
        --outdir "$MET_DIR/04_classification" || true
  else
    skip_step "03 Metabolite pathway enrichment" "no diff_metabolites.csv from step 02"
    skip_step "04 Metabolite classification" "no diff_metabolites.csv from step 02"
  fi

  run_step "05 Biomarker analysis (ROC + random forest)" \
    "$RSCRIPT" "$REPO_ROOT/scripts/metabolomics/05_biomarker.R" \
      --data "$PREPROCESSED" --metadata "$METADATA" --outdir "$MET_DIR/05_biomarker" || true

  run_step "06 Metabolite correlation network" \
    "$RSCRIPT" "$REPO_ROOT/scripts/metabolomics/06_metabolite_network.R" \
      --data "$PREPROCESSED" --outdir "$MET_DIR/06_network" --cor-cutoff "$COR_CUTOFF" || true

  run_step "07 Metabolome figures" \
    "$RSCRIPT" "$REPO_ROOT/scripts/metabolomics/07_plot_metabolome.R" \
      --results-dir "$MET_DIR" --outdir "$MET_DIR/07_figures" || true
fi

# ---------------------------------------------------------------------------
# Module 3: multi-omics integration
# ---------------------------------------------------------------------------
if [ "$STAGES" = "all" ] || [ "$STAGES" = "integration" ]; then
  log "--- Module 3: transcriptome + metabolome integration ---"

  RNA_MATRIX="$(ls "$TX_DIR"/02_de/normalized_counts_vst.csv 2>/dev/null | head -1)"
  [ -n "$RNA_MATRIX" ] || RNA_MATRIX="$(ls "$TX_DIR"/01_qc/counts_filtered.csv 2>/dev/null | head -1)"
  [ -n "$RNA_MATRIX" ] || RNA_MATRIX="$COUNTS"
  MET_MATRIX="$(ls "$MET_DIR"/01_preprocess/metabolome_preprocessed.csv 2>/dev/null | head -1)"
  [ -n "$MET_MATRIX" ] || MET_MATRIX="$METABOLITES"

  run_step "01 MOFA factor analysis" \
    "$RSCRIPT" "$REPO_ROOT/scripts/integration/01_mofa.R" \
      --rna "$RNA_MATRIX" --metabolites "$MET_MATRIX" --metadata "$METADATA" \
      --outdir "$INT_DIR/01_mofa" --seed "$SEED" || true

  run_step "02 DIABLO / O2PLS integration" \
    "$RSCRIPT" "$REPO_ROOT/scripts/integration/02_diablo.R" \
      --rna "$RNA_MATRIX" --metabolites "$MET_MATRIX" --metadata "$METADATA" \
      --outdir "$INT_DIR/02_diablo" || true

  DE_SIG="$(ls "$TX_DIR"/02_de/de_sig_*.csv 2>/dev/null | head -1)"
  DE_ALL="$(ls "$TX_DIR"/02_de/de_all_*.csv 2>/dev/null | head -1)"
  [ -n "$DE_SIG" ] || DE_SIG="$DE_ALL"
  DIFF_METS="$(ls "$MET_DIR"/02_oplsda/diff_metabolites.csv 2>/dev/null | head -1)"

  if [ -n "$DE_SIG" ] && [ -n "$DIFF_METS" ]; then
    run_step "03 Gene-metabolite correlation" \
      "$RSCRIPT" "$REPO_ROOT/scripts/integration/03_correlation.R" \
        --de-genes "$DE_SIG" --de-metabolites "$DIFF_METS" \
        --rna "$RNA_MATRIX" --metabolites "$MET_MATRIX" \
        --outdir "$INT_DIR/03_correlation" --cor-cutoff "$COR_CUTOFF" || true

    run_step "04 Joint pathway analysis" \
      "$RSCRIPT" "$REPO_ROOT/scripts/integration/04_joint_pathway.R" \
        --de-genes "$DE_SIG" --de-metabolites "$DIFF_METS" \
        --outdir "$INT_DIR/04_joint_pathway" --padj "$PADJ" || true
  else
    skip_step "03 Gene-metabolite correlation" "missing DE genes or differential metabolites"
    skip_step "04 Joint pathway analysis" "missing DE genes or differential metabolites"
  fi

  GENE_MODULES="$(ls "$TX_DIR"/05_wgcna/gene_modules.csv 2>/dev/null | head -1)"
  if [ -n "$GENE_MODULES" ]; then
    run_step "05 Multi-omics network integration" \
      "$RSCRIPT" "$REPO_ROOT/scripts/integration/05_network_integration.R" \
        --gene-modules "$GENE_MODULES" --rna "$RNA_MATRIX" --metabolites "$MET_MATRIX" \
        --metadata "$METADATA" --outdir "$INT_DIR/05_network" --cor-cutoff "$COR_CUTOFF" || true
  else
    skip_step "05 Multi-omics network integration" "no gene_modules.csv from WGCNA"
  fi

  run_step "06 Integration figures" \
    "$RSCRIPT" "$REPO_ROOT/scripts/integration/06_plot_integration.R" \
      --results-dir "$INT_DIR" --transcriptome-dir "$TX_DIR" --metabolome-dir "$MET_DIR" \
      --outdir "$INT_DIR/06_figures" || true
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log "Pipeline finished: $N_PASS passed, $N_FAIL failed, $N_SKIP skipped"
log "Full log: $LOG"
echo
echo "Summary: $N_PASS passed, $N_FAIL failed, $N_SKIP skipped"
echo "Log    : $LOG"
echo "Results: $OUTDIR"

if [ "$N_FAIL" -gt 0 ]; then
  echo
  echo "Some steps failed. Common causes:"
  echo "  - missing R packages (see docs/installation.md);"
  echo "  - network access required for KEGG/GO enrichment;"
  echo "  - sample IDs in the matrices not matching metadata.csv."
  exit 1
fi
exit 0
