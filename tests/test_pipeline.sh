#!/usr/bin/env bash
#
# test_pipeline.sh -- Rice-skills self-test
#
# Description:
#   Verifies that the shipped skill is structurally sound and that the
#   tooling runs on this machine, without requiring any real data:
#     1. Python syntax check for every .py file
#     2. R syntax check (parse) for every .R file
#     3. Shared-utility smoke tests (tests/test_utils.py, tests/test_utils.R)
#     4. End-to-end smoke test on synthetic demo data:
#        make_demo_data.py -> 01_qc_filter.py -> 01_preprocess.py
#   Every check reports PASS, FAIL or SKIP (SKIP when the interpreter or an
#   optional package is unavailable), and the script exits non-zero only when
#   a check actually failed.
#
# Usage:
#   bash tests/test_pipeline.sh
#   bash tests/test_pipeline.sh --keep-tmp     # keep the temporary demo data
#   bash tests/test_pipeline.sh --python /path/to/python --rscript /path/to/Rscript

set -u

KEEP_TMP=0
PYTHON_OVERRIDE=""
RSCRIPT_OVERRIDE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --keep-tmp) KEEP_TMP=1; shift ;;
    --python)   PYTHON_OVERRIDE="$2"; shift 2 ;;
    --rscript)  RSCRIPT_OVERRIDE="$2"; shift 2 ;;
    *)          echo "unknown option: $1"; exit 1 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT" || exit 1

N_PASS=0; N_FAIL=0; N_SKIP=0

pass() { printf 'PASS: %s\n' "$1"; N_PASS=$((N_PASS+1)); }
fail() { printf 'FAIL: %s\n' "$1"; N_FAIL=$((N_FAIL+1)); }
skip() { printf 'SKIP: %s (%s)\n' "$1" "$2"; N_SKIP=$((N_SKIP+1)); }

# ---------------------------------------------------------------------------
# Interpreter detection
# ---------------------------------------------------------------------------
PYTHON=""
if [ -n "$PYTHON_OVERRIDE" ]; then
  PYTHON="$PYTHON_OVERRIDE"
else
  # Prefer a Python that actually carries the scientific stack (matplotlib);
  # fall back to the first interpreter found (checks will SKIP/FAIL with hints).
  _py_first=""
  for _cand in python3 python; do
    if command -v "$_cand" >/dev/null 2>&1; then
      [ -z "$_py_first" ] && _py_first="$_cand"
      if "$_cand" -c "import matplotlib" >/dev/null 2>&1; then PYTHON="$_cand"; break; fi
    fi
  done
  if [ -z "$PYTHON" ]; then
    for _cand in "$HOME/AppData/Local/Programs/Python"/Python*/python.exe \
                 /c/Users/*/AppData/Local/Programs/Python/Python*/python.exe; do
      if [ -x "$_cand" ] && "$_cand" -c "import matplotlib" >/dev/null 2>&1; then PYTHON="$_cand"; break; fi
    done
  fi
  [ -z "$PYTHON" ] && PYTHON="$_py_first"
fi

RSCRIPT=""
if [ -n "$RSCRIPT_OVERRIDE" ]; then
  RSCRIPT="$RSCRIPT_OVERRIDE"
elif command -v Rscript >/dev/null 2>&1; then RSCRIPT="Rscript"
elif [ -x "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" ]; then
  RSCRIPT="/c/Program Files/R/R-4.6.0/bin/Rscript.exe"
fi

echo "=== Rice-skills self-test ==="
echo "Repo root: $REPO_ROOT"
echo "Python   : ${PYTHON:-not found}"
echo "Rscript  : ${RSCRIPT:-not found}"
echo

# ---------------------------------------------------------------------------
# 0. Repository layout
# ---------------------------------------------------------------------------
missing_dir=0
for d in scripts/transcriptomics scripts/metabolomics scripts/integration scripts/utils \
         scripts/gwas scripts/popgen scripts/design \
         scripts/variants scripts/pangenome scripts/qtl \
         scripts/proteomics scripts/epigenomics scripts/phenomics \
         scripts/single_cell scripts/spatial \
         config templates/figure_templates examples docs tests; do
  if [ ! -d "$d" ]; then echo "  missing directory: $d"; missing_dir=$((missing_dir+1)); fi
done
if [ "$missing_dir" -eq 0 ]; then pass "repository layout"; else fail "repository layout ($missing_dir missing)"; fi

if [ -f "config/nature_palette.json" ]; then pass "palette config present"; else fail "palette config missing"; fi
if [ -f "SKILL.md" ]; then pass "SKILL.md present"; else fail "SKILL.md missing"; fi

docs_ok=1
for f in docs/installation.md docs/usage.md docs/output_description.md \
         docs/genomics_usage.md docs/design_usage.md docs/lab_protocols.md \
         docs/protocols_molecular.md docs/research_workflow_tf.md \
         docs/variants_usage.md docs/pangenome_usage.md docs/qtl_usage.md \
         docs/proteomics_usage.md docs/epigenomics_usage.md docs/phenomics_usage.md \
         docs/single_cell_usage.md docs/spatial_usage.md docs/rice_databases.md; do
  [ -f "$f" ] || { docs_ok=0; echo "  missing doc: $f"; }
done
if [ "$docs_ok" -eq 1 ]; then pass "core documentation present (17 docs)"; else fail "documentation incomplete"; fi
if [ -f "third_party/nature-skills/LICENSE" ] && [ -f "third_party/nature-skills/NOTICE.md" ] &&
   [ -f "third_party/nature-skills/skills/nature-figure/SKILL.md" ]; then
  pass "third_party nature-skills vendored (license + notice present)"
else
  fail "third_party nature-skills vendoring incomplete (expected LICENSE, NOTICE.md, skills/nature-figure)"
fi
tp_ok=1
for f in third_party/bioskills/LICENSE third_party/bioskills/NOTICE.md \
         third_party/kdense-scientific-skills/LICENSE.md \
         third_party/kdense-scientific-skills/NOTICE.md \
         third_party/cell-skills/LICENSE-Cell \
         third_party/cell-skills/LICENSE-Cell_CNS_Figure \
         third_party/cell-skills/NOTICE.md \
         third_party/sciagent-skills/LICENSE \
         third_party/sciagent-skills/NOTICE.md; do
  [ -f "$f" ] || tp_ok=0
done
if [ "$tp_ok" -eq 1 ]; then
  pass "third_party bioskills + kdense + cell-skills + sciagent-skills vendored (licenses + notices)"
else
  fail "third_party vendoring incomplete (bioskills / kdense / cell-skills / sciagent-skills)"
fi

if [ -f "assets/annotations/rapdb_to_msu_map.csv" ] && \
   [ -f "assets/annotations/rice_TF_PlantTFDB_msu.csv" ] && \
   [ -f "assets/annotations/rice_TF_families.gmt" ]; then
  pass "bundled annotation snapshot present (RAP-DB map + PlantTFDB TF + GMT)"
else
  fail "bundled annotation snapshot incomplete (assets/annotations/)"
fi

# ---------------------------------------------------------------------------
# 1. Python syntax check
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ]; then
  if "$PYTHON" -m compileall -q scripts templates tests >/dev/null 2>&1; then
    pass "Python syntax (compileall)"
  else
    fail "Python syntax (compileall)"
  fi
else
  skip "Python syntax" "python not available"
fi

# ---------------------------------------------------------------------------
# 2. R syntax check
# ---------------------------------------------------------------------------
if [ -n "$RSCRIPT" ]; then
  r_out=$("$RSCRIPT" -e 'files <- c(list.files("scripts", pattern="[.]R$", recursive=TRUE, full.names=TRUE), list.files("tests", pattern="[.]R$", full.names=TRUE)); bad <- character(0); for (f in files) { res <- tryCatch({parse(f); "ok"}, error=function(e) conditionMessage(e)); if (res != "ok") bad <- c(bad, paste0(f, ": ", res)) }; cat(length(files), "\n"); cat(bad, sep="\n")' 2>/dev/null)
  n_r=$(printf '%s\n' "$r_out" | head -1)
  r_bad=$(printf '%s\n' "$r_out" | tail -n +2 | grep -v '^$')
  if [ -z "$r_bad" ]; then
    pass "R syntax (parse) for ${n_r} files"
  else
    fail "R syntax (parse): $r_bad"
  fi
else
  skip "R syntax" "Rscript not available"
fi

# ---------------------------------------------------------------------------
# 3. Shared utility smoke tests
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ] && [ -f tests/test_utils.py ]; then
  if "$PYTHON" tests/test_utils.py >/dev/null 2>&1; then pass "Python utility smoke test"
  else fail "Python utility smoke test (run: $PYTHON tests/test_utils.py)"; fi
else
  skip "Python utility smoke test" "python or tests/test_utils.py missing"
fi

if [ -n "$RSCRIPT" ] && [ -f tests/test_utils.R ]; then
  if "$RSCRIPT" tests/test_utils.R >/dev/null 2>&1; then pass "R utility smoke test"
  else fail "R utility smoke test (run: Rscript tests/test_utils.R)"; fi
else
  skip "R utility smoke test" "Rscript or tests/test_utils.R missing"
fi

# ---------------------------------------------------------------------------
# 4. End-to-end smoke test on synthetic demo data
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ]; then
  TMP="$(mktemp -d 2>/dev/null || echo "${TEMP:-/tmp}/rice_selftest_$$")"
  mkdir -p "$TMP"

  if "$PYTHON" scripts/utils/make_demo_data.py --outdir "$TMP/data" --n-genes 500 \
       --n-metabolites 120 >/dev/null 2>&1; then
    pass "synthetic demo data generation"
  else
    fail "synthetic demo data generation"
  fi

  if "$PYTHON" scripts/transcriptomics/01_qc_filter.py \
       --counts "$TMP/data/counts.csv" --metadata "$TMP/data/metadata.csv" \
       --outdir "$TMP/tx" >/dev/null 2>&1 &&
     [ -f "$TMP/tx/counts_filtered.csv" ] && [ -f "$TMP/tx/sample_pca.pdf" ] &&
     [ -f "$TMP/tx/sample_pca.tiff" ]; then
    pass "transcriptome QC (01_qc_filter.py) end-to-end"
  else
    fail "transcriptome QC (01_qc_filter.py) end-to-end"
  fi

  if "$PYTHON" scripts/metabolomics/01_preprocess.py \
       --abundance "$TMP/data/metabolites.csv" --metadata "$TMP/data/metadata.csv" \
       --outdir "$TMP/met" >/dev/null 2>&1 &&
     [ -f "$TMP/met/metabolome_preprocessed.csv" ] && [ -f "$TMP/met/pca_after.pdf" ] &&
     [ -f "$TMP/met/pca_after.tiff" ]; then
    pass "metabolome preprocessing (01_preprocess.py) end-to-end"
  else
    fail "metabolome preprocessing (01_preprocess.py) end-to-end"
  fi

  if [ "$KEEP_TMP" -eq 0 ]; then
    rm -f "$TMP"/*/* "$TMP"/* 2>/dev/null
    rmdir "$TMP"/tx "$TMP"/met "$TMP"/data "$TMP" 2>/dev/null
  else
    echo "      (demo data kept at: $TMP)"
  fi
else
  skip "end-to-end smoke test" "python not available"
fi

# ---------------------------------------------------------------------------
# 5. Genomics smoke test (GWAS + population genetics on synthetic data)
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ] && [ -n "$RSCRIPT" ]; then
  GTMP="$(mktemp -d 2>/dev/null || echo "${TEMP:-/tmp}/rice_geno_selftest_$$")"
  mkdir -p "$GTMP"

  if "$PYTHON" scripts/utils/make_demo_genomics.py --outdir "$GTMP" \
       --n-indica 40 --n-japonica 40 --n-admixed 5 --n-snps 800 >/dev/null 2>&1 &&
     [ -f "$GTMP/genotype.csv" ] && [ -f "$GTMP/phenotype.csv" ]; then
    pass "genomics demo data generation"
  else
    fail "genomics demo data generation"
  fi

  if "$RSCRIPT" scripts/gwas/03_association.R \
       --genotype "$GTMP/genotype.csv" --snp-info "$GTMP/snp_info.csv" \
       --phenotype "$GTMP/phenotype.csv" --trait salt_tolerance_index \
       --outdir "$GTMP/gwas" >/dev/null 2>&1 &&
     [ -f "$GTMP/gwas/association_results.csv" ]; then
    pass "GWAS association (03_association.R) end-to-end"
  else
    fail "GWAS association (03_association.R) end-to-end"
  fi

  if "$RSCRIPT" scripts/gwas/05_plot_gwas.R \
       --assoc "$GTMP/gwas/association_results.csv" --snp-info "$GTMP/snp_info.csv" \
       --outdir "$GTMP/figs" >/dev/null 2>&1 &&
     [ -f "$GTMP/figs/manhattan.pdf" ] && [ -f "$GTMP/figs/manhattan.tiff" ]; then
    pass "GWAS figures (05_plot_gwas.R)"
  else
    fail "GWAS figures (05_plot_gwas.R)"
  fi

  if "$RSCRIPT" scripts/popgen/01_pca_structure.R \
       --genotype "$GTMP/genotype.csv" --snp-info "$GTMP/snp_info.csv" \
       --phenotype "$GTMP/phenotype.csv" --group-col group \
       --outdir "$GTMP/pca" >/dev/null 2>&1 &&
     [ -f "$GTMP/pca/pca_scatter.pdf" ] && [ -f "$GTMP/pca/kinship_matrix.csv" ]; then
    pass "population structure (01_pca_structure.R)"
  else
    fail "population structure (01_pca_structure.R)"
  fi

  if "$RSCRIPT" scripts/popgen/03_selection_stats.R \
       --genotype "$GTMP/genotype.csv" --snp-info "$GTMP/snp_info.csv" \
       --phenotype "$GTMP/phenotype.csv" --group-col group \
       --group1 indica --group2 japonica \
       --outdir "$GTMP/fst" >/dev/null 2>&1 &&
     [ -f "$GTMP/fst/fst_windows.csv" ]; then
    pass "selection statistics (03_selection_stats.R)"
  else
    fail "selection statistics (03_selection_stats.R)"
  fi

  rm -f "$GTMP"/*/* "$GTMP"/*/*/* "$GTMP"/* 2>/dev/null
  rmdir "$GTMP"/gwas "$GTMP"/figs "$GTMP"/pca "$GTMP"/fst "$GTMP" 2>/dev/null
else
  skip "genomics smoke test" "python or Rscript not available"
fi

# ---------------------------------------------------------------------------
# 6. Design smoke test (sgRNA + primer design on the demo gene)
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ]; then
  DTMP="$(mktemp -d 2>/dev/null || echo "${TEMP:-/tmp}/rice_design_selftest_$$")"
  mkdir -p "$DTMP"

  if "$PYTHON" scripts/utils/make_demo_genomics.py --outdir "$DTMP/data" \
       --n-indica 5 --n-japonica 5 --n-admixed 2 --n-snps 400 >/dev/null 2>&1 &&
     [ -f "$DTMP/data/demo_gene.fa" ]; then
    pass "design demo gene generation"
  else
    fail "design demo gene generation"
  fi

  if "$PYTHON" scripts/design/01_sgrna_design.py \
       --fasta "$DTMP/data/demo_gene.fa" --gff "$DTMP/data/demo_gene.gff3" \
       --outdir "$DTMP/sgrna" >/dev/null 2>&1 &&
     [ -f "$DTMP/sgrna/sgrna_recommended_pairs.csv" ] &&
     [ -f "$DTMP/sgrna/construction_primers.csv" ]; then
    pass "sgRNA design (01_sgrna_design.py)"
  else
    fail "sgRNA design (01_sgrna_design.py)"
  fi

  if "$PYTHON" -c "import primer3" >/dev/null 2>&1; then
    if "$PYTHON" scripts/design/02_primer_design.py --mode qpcr \
         --fasta "$DTMP/data/demo_gene.fa" --gff "$DTMP/data/demo_gene.gff3" \
         --outdir "$DTMP/primers" >/dev/null 2>&1 &&
       [ -f "$DTMP/primers/primers.csv" ]; then
      pass "qPCR primer design (02_primer_design.py)"
    else
      fail "qPCR primer design (02_primer_design.py)"
    fi
  else
    skip "qPCR primer design" "primer3-py not installed (pip install primer3-py)"
  fi

  if "$PYTHON" scripts/design/03_qpcr_analysis.py --demo \
       --outdir "$DTMP/qpcr" >/dev/null 2>&1 &&
     [ -f "$DTMP/qpcr/ddct_table.csv" ] &&
     [ -f "$DTMP/qpcr/rel_expression_summary.csv" ]; then
    pass "qPCR quantification (03_qpcr_analysis.py --demo)"
  else
    fail "qPCR quantification (03_qpcr_analysis.py --demo)"
  fi

  rm -f "$DTMP"/*/*/* "$DTMP"/*/* "$DTMP"/* 2>/dev/null
  rmdir "$DTMP"/sgrna "$DTMP"/primers "$DTMP"/data "$DTMP" 2>/dev/null
else
  skip "design smoke test" "python not available"
fi

# ---------------------------------------------------------------------------
# 7. v3.0 module smoke tests (pangenome / QTL / proteomics / epigenome /
#    phenomics / meta-GWAS on small synthetic data)
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ] && [ -n "$RSCRIPT" ]; then
  VTMP="$(mktemp -d 2>/dev/null || echo "${TEMP:-/tmp}/rice_v3_selftest_$$")"
  mkdir -p "$VTMP"

  # --- pangenome: PAV analysis + Pan-GWAS ---
  if "$PYTHON" scripts/utils/make_demo_pangenome.py --outdir "$VTMP/pg" \
       --n-subgroup-a 20 --n-subgroup-b 20 --n-markers 150 >/dev/null 2>&1 &&
     "$RSCRIPT" scripts/pangenome/03_pav_analysis.R --pav "$VTMP/pg/pav_matrix.csv" \
       --meta "$VTMP/pg/marker_info.csv" --pheno "$VTMP/pg/phenotype.csv" \
       --trait grain_yield --outdir "$VTMP/pg_pav" >/dev/null 2>&1 &&
     [ -f "$VTMP/pg_pav/pav_class.csv" ] && [ -f "$VTMP/pg_pav/pav_class_bar.pdf" ]; then
    pass "pangenome PAV analysis (03_pav_analysis.R)"
  else
    fail "pangenome PAV analysis (03_pav_analysis.R)"
  fi

  if [ -f "$VTMP/pg/pav_matrix.csv" ] &&
     "$RSCRIPT" scripts/pangenome/04_pan_gwas.R --pav "$VTMP/pg/pav_matrix.csv" \
       --snp-info "$VTMP/pg/marker_info.csv" --phenotype "$VTMP/pg/phenotype.csv" \
       --trait grain_yield --outdir "$VTMP/pg_gwas" >/dev/null 2>&1 &&
     [ -f "$VTMP/pg_gwas/pan_gwas_results.csv" ]; then
    pass "pan-GWAS (04_pan_gwas.R)"
  else
    fail "pan-GWAS (04_pan_gwas.R)"
  fi

  # --- QTL mapping (eQTL cis scan) ---
  if "$PYTHON" scripts/utils/make_demo_genomics.py --outdir "$VTMP/gq" \
       --n-indica 30 --n-japonica 30 --n-admixed 10 --n-snps 600 >/dev/null 2>&1 &&
     "$PYTHON" scripts/utils/make_demo_qtl.py --genotype "$VTMP/gq/genotype.csv" \
       --snp-info "$VTMP/gq/snp_info.csv" --outdir "$VTMP/gq" \
       --n-features 60 --n-causal 25 >/dev/null 2>&1 &&
     "$RSCRIPT" scripts/qtl/01_qtl_mapping.R --genotype "$VTMP/gq/genotype.csv" \
       --snp-info "$VTMP/gq/snp_info.csv" --features "$VTMP/gq/expression.csv" \
       --feature-info "$VTMP/gq/gene_info.csv" --mode cis \
       --outdir "$VTMP/qtl" >/dev/null 2>&1 &&
     [ -f "$VTMP/qtl/qtl_results.csv" ]; then
    pass "QTL mapping (01_qtl_mapping.R, cis eQTL)"
  else
    fail "QTL mapping (01_qtl_mapping.R, cis eQTL)"
  fi

  # --- proteomics ---
  if "$PYTHON" scripts/utils/make_demo_proteomics.py --outdir "$VTMP/prot" \
       --n-proteins 400 --n-de 30 --n-ctrl 8 --n-salt 8 >/dev/null 2>&1 &&
     "$PYTHON" scripts/proteomics/01_preprocess.py --abundance "$VTMP/prot/proteins.csv" \
       --metadata "$VTMP/prot/metadata.csv" --outdir "$VTMP/prot1" >/dev/null 2>&1 &&
     "$RSCRIPT" scripts/proteomics/02_diff_protein.R \
       --matrix "$VTMP/prot1/proteins_filtered.csv" --metadata "$VTMP/prot/metadata.csv" \
       --ref control --alt salt --outdir "$VTMP/prot2" >/dev/null 2>&1 &&
     [ -f "$VTMP/prot2/protein_diff_all.csv" ]; then
    pass "proteomics pipeline (01_preprocess + 02_diff_protein)"
  else
    fail "proteomics pipeline (01_preprocess + 02_diff_protein)"
  fi

  # --- epigenome: DMR + DAR ---
  if "$PYTHON" scripts/utils/make_demo_epigenome.py --outdir "$VTMP/epi" \
       --n-meth-regions 600 --n-dmr 40 --n-peaks 600 --n-dar 40 \
       --n-ctrl 6 --n-salt 6 >/dev/null 2>&1 &&
     "$RSCRIPT" scripts/epigenomics/01_methylation_analysis.R \
       --matrix "$VTMP/epi/methylation_regions.csv" --metadata "$VTMP/epi/metadata.csv" \
       --ref control --alt salt --outdir "$VTMP/epi1" >/dev/null 2>&1 &&
     [ -f "$VTMP/epi1/dmr_all.csv" ]; then
    pass "methylation DMR analysis (01_methylation_analysis.R)"
  else
    fail "methylation DMR analysis (01_methylation_analysis.R)"
  fi

  if [ -f "$VTMP/epi/peak_counts.csv" ] &&
     "$RSCRIPT" scripts/epigenomics/02_atac_chip_analysis.R \
       --matrix "$VTMP/epi/peak_counts.csv" --metadata "$VTMP/epi/metadata.csv" \
       --ref control --alt salt --outdir "$VTMP/epi2" >/dev/null 2>&1 &&
     [ -f "$VTMP/epi2/dar_all.csv" ]; then
    pass "ATAC/ChIP differential peaks (02_atac_chip_analysis.R)"
  else
    fail "ATAC/ChIP differential peaks (02_atac_chip_analysis.R)"
  fi

  # --- phenomics: BLUP + spectral indices ---
  if "$PYTHON" scripts/utils/make_demo_phenomics.py --outdir "$VTMP/phe" \
       --n-lines 30 >/dev/null 2>&1 &&
     "$RSCRIPT" scripts/phenomics/01_phenotype_blup.R --pheno "$VTMP/phe/phenotype_long.csv" \
       --env-col env --outdir "$VTMP/phe1" >/dev/null 2>&1 &&
     [ -f "$VTMP/phe1/blups.csv" ] && [ -f "$VTMP/phe1/heritability.csv" ]; then
    pass "phenotype BLUP (01_phenotype_blup.R)"
  else
    fail "phenotype BLUP (01_phenotype_blup.R)"
  fi

  if [ -f "$VTMP/phe/spectra.csv" ] &&
     "$RSCRIPT" scripts/phenomics/02_spectral_indices.R --spectra "$VTMP/phe/spectra.csv" \
       --outdir "$VTMP/phe2" >/dev/null 2>&1 &&
     [ -f "$VTMP/phe2/spectral_indices.csv" ]; then
    pass "spectral indices (02_spectral_indices.R)"
  else
    fail "spectral indices (02_spectral_indices.R)"
  fi

  # --- meta-GWAS (tiny synthetic study tables) ---
  # NB: pass $VTMP via argv (not string interpolation): on Windows/Git-Bash
  # the temp path may contain backslashes that break Python string literals.
  "$PYTHON" - "$VTMP" <<'PYEOF' >/dev/null 2>&1
import sys
import numpy as np
import pandas as pd
tmp = sys.argv[1]
rng = np.random.default_rng(1)
n = 80
snp = [f"s{i:03d}" for i in range(n)]
pd.DataFrame({"snp_id": snp, "chr": "Chr01",
              "pos": (np.arange(n) + 1) * 1000}).to_csv(tmp + "/meta_info.csv", index=False)
for tag in ("a", "b"):
    beta = rng.normal(0, 0.2, n)
    beta[10:14] += 0.8
    pd.DataFrame({"snp_id": snp, "beta": beta, "se": 0.1}).to_csv(
        tmp + "/meta_" + tag + ".csv", index=False)
PYEOF
  if [ -f "$VTMP/meta_a.csv" ] && [ -f "$VTMP/meta_b.csv" ] &&
     "$RSCRIPT" scripts/gwas/06_meta_gwas.R \
       --studies "$VTMP/meta_a.csv,$VTMP/meta_b.csv" --study-names "studyA,studyB" \
       --snp-info "$VTMP/meta_info.csv" --outdir "$VTMP/meta" >/dev/null 2>&1 &&
     [ -f "$VTMP/meta/meta_results.csv" ]; then
    pass "meta-GWAS (06_meta_gwas.R)"
  else
    fail "meta-GWAS (06_meta_gwas.R)"
  fi

  rm -f "$VTMP"/*/*/* "$VTMP"/*/* "$VTMP"/* 2>/dev/null
  rmdir "$VTMP"/pg "$VTMP"/pg_pav "$VTMP"/pg_gwas "$VTMP"/gq "$VTMP"/qtl \
        "$VTMP"/prot "$VTMP"/prot1 "$VTMP"/prot2 "$VTMP"/epi "$VTMP"/epi1 \
        "$VTMP"/epi2 "$VTMP"/phe "$VTMP"/phe1 "$VTMP"/phe2 "$VTMP"/meta \
        "$VTMP" 2>/dev/null
else
  skip "v3.0 module smoke tests" "python or Rscript not available"
fi

# ---------------------------------------------------------------------------
# 8. v3.1 smoke tests (haplotype / colocalization / GBLUP / AMMI / XP-EHH /
#    KA-KS / promoter k-mer on small synthetic data)
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ] && [ -n "$RSCRIPT" ]; then
  WTMP="$(mktemp -d 2>/dev/null || echo "${TEMP:-/tmp}/rice_v31_selftest_$$")"
  mkdir -p "$WTMP"
  "$PYTHON" scripts/utils/make_demo_genomics.py --outdir "$WTMP/gen" \
      --n-indica 50 --n-japonica 50 --n-admixed 10 --n-snps 1200 >/dev/null 2>&1
  LEAD=$(awk -F, 'NR==2{print $1}' "$WTMP/gen/snp_info.csv")

  if "$RSCRIPT" scripts/gwas/07_haplotype.R --genotype "$WTMP/gen/genotype.csv" \
       --snp-info "$WTMP/gen/snp_info.csv" --phenotype "$WTMP/gen/phenotype.csv" \
       --trait salt_tolerance_index --region "Chr01:1-43000000" --lead-snp "$LEAD" \
       --outdir "$WTMP/hap" >/dev/null 2>&1 &&
     [ -f "$WTMP/hap/sample_haplotypes.csv" ] && [ -f "$WTMP/hap/haplotype_phenotype.pdf" ]; then
    pass "haplotype analysis (07_haplotype.R)"
  else
    fail "haplotype analysis (07_haplotype.R)"
  fi

  if [ -f "$WTMP/hap/sample_haplotypes.csv" ] &&
     "$RSCRIPT" scripts/gwas/07_haplotype.R --genotype "$WTMP/gen/genotype.csv" \
       --snp-info "$WTMP/gen/snp_info.csv" --phenotype "$WTMP/gen/phenotype.csv" \
       --trait salt_tolerance_index --region "Chr01:1-43000000" --lead-snp "$LEAD" \
       --mode snp --outdir "$WTMP/hap_snp" >/dev/null 2>&1 &&
     [ -f "$WTMP/hap_snp/group_stats.csv" ]; then
    pass "allele-dose mode (07_haplotype.R --mode snp)"
  else
    fail "allele-dose mode (07_haplotype.R --mode snp)"
  fi

  "$PYTHON" - "$WTMP" <<'PYEOF' >/dev/null 2>&1
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]
rng = np.random.default_rng(5)
n = 200
snp = [f"s{i:04d}" for i in range(1, n + 1)]
pd.DataFrame({"snp_id": snp, "chr": "Chr01", "pos": np.arange(1, n + 1) * 5000}).to_csv(
    tmp + "/c_snp_info.csv", index=False)
p_g = np.ones(n); p_g[94:100] = 10 ** (-rng.uniform(4, 9, 6))
pd.DataFrame({"snp_id": snp, "p": np.round(p_g, 8)}).to_csv(tmp + "/c_gwas.csv", index=False)
pd.DataFrame([{"feature_id": "GENE_T", "snp_id": s, "p": round(p, 8), "type": "cis"}
              for s, p in zip(snp, p_g)]).to_csv(tmp + "/c_eqtl.csv", index=False)
PYEOF
  if [ -f "$WTMP/c_gwas.csv" ] &&
     "$RSCRIPT" scripts/qtl/02_colocalization.R --gwas "$WTMP/c_gwas.csv" \
       --eqtl "$WTMP/c_eqtl.csv" --snp-info "$WTMP/c_snp_info.csv" \
       --region "Chr01:1-1000000" --outdir "$WTMP/coloc" >/dev/null 2>&1 &&
     [ -f "$WTMP/coloc/coloc_summary.csv" ]; then
    pass "GWAS x QTL colocalization (02_colocalization.R)"
  else
    fail "GWAS x QTL colocalization (02_colocalization.R)"
  fi

  if "$RSCRIPT" scripts/phenomics/03_gs_gblup.R --genotype "$WTMP/gen/genotype.csv" \
       --phenotype "$WTMP/gen/phenotype.csv" --trait salt_tolerance_index \
       --k-folds 3 --repeats 1 --outdir "$WTMP/gs" >/dev/null 2>&1 &&
     [ -f "$WTMP/gs/gebv_all.csv" ] && [ -f "$WTMP/gs/cv_lambda_summary.csv" ]; then
    pass "GBLUP genomic selection (03_gs_gblup.R)"
  else
    fail "GBLUP genomic selection (03_gs_gblup.R)"
  fi

  if "$PYTHON" scripts/utils/make_demo_phenomics.py --outdir "$WTMP/phe" \
       --n-lines 30 >/dev/null 2>&1 &&
     "$RSCRIPT" scripts/phenomics/04_ammi_gge.R --pheno "$WTMP/phe/phenotype_long.csv" \
       --env-col env --trait grain_yield --outdir "$WTMP/ammi" >/dev/null 2>&1 &&
     [ -f "$WTMP/ammi/gge_biplot.pdf" ] && [ -f "$WTMP/ammi/stability.csv" ]; then
    pass "AMMI/GGE biplot (04_ammi_gge.R)"
  else
    fail "AMMI/GGE biplot (04_ammi_gge.R)"
  fi

  if "$RSCRIPT" scripts/popgen/07_xpehh.R --genotype "$WTMP/gen/genotype.csv" \
       --snp-info "$WTMP/gen/snp_info.csv" --metadata "$WTMP/gen/phenotype.csv" \
       --group-col group --group1 indica --group2 japonica --chr Chr05 \
       --max-core 100 --max-dist 30 --outdir "$WTMP/xpehh" >/dev/null 2>&1 &&
     [ -f "$WTMP/xpehh/xpehh_sites.csv" ]; then
    pass "XP-EHH selection scan (07_xpehh.R)"
  else
    fail "XP-EHH selection scan (07_xpehh.R)"
  fi

  if "$PYTHON" scripts/utils/make_demo_evolution.py --outdir "$WTMP/evo" >/dev/null 2>&1 &&
     "$PYTHON" scripts/evolution/01_ka_ks.py --pairs "$WTMP/evo/ortholog_pairs.csv" \
       --cds "$WTMP/evo/cds.fasta" --outdir "$WTMP/ks" >/dev/null 2>&1 &&
     [ -f "$WTMP/ks/ka_ks.csv" ]; then
    pass "KA/KS evolution (01_ka_ks.py)"
  else
    fail "KA/KS evolution (01_ka_ks.py)"
  fi

  if "$PYTHON" scripts/utils/make_demo_tf.py --outdir "$WTMP/tf" >/dev/null 2>&1 &&
     "$PYTHON" scripts/tf/01_promoter_kmer.py --genome "$WTMP/tf/genome.fa" \
       --gff "$WTMP/tf/annotation.gff3" --genes "$WTMP/tf/candidate_genes.csv" \
       --outdir "$WTMP/kmer" >/dev/null 2>&1 &&
     [ -f "$WTMP/kmer/kmer_enrichment.csv" ]; then
    pass "promoter k-mer enrichment (01_promoter_kmer.py)"
  else
    fail "promoter k-mer enrichment (01_promoter_kmer.py)"
  fi

  rm -f "$WTMP"/*/*/* "$WTMP"/*/* "$WTMP"/* 2>/dev/null
  rmdir "$WTMP"/gen "$WTMP"/hap "$WTMP"/hap_snp "$WTMP"/coloc "$WTMP"/gs "$WTMP"/ammi \
        "$WTMP"/phe "$WTMP"/xpehh "$WTMP"/evo "$WTMP"/ks "$WTMP"/tf "$WTMP"/kmer \
        "$WTMP" 2>/dev/null
else
  skip "v3.1 smoke tests" "python or Rscript not available"
fi

# ---------------------------------------------------------------------------
# 9. v3.2 smoke tests (Bayesian coloc / 2D-GWAS / power / UniProt mapping /
#    InterPro domains / GRN / miRNA targets / motif-TF / QTL region SNPs /
#    genome-wide BSA scan -- all built-in --demo self-tests)
# ---------------------------------------------------------------------------
if [ -n "$PYTHON" ] && [ -n "$RSCRIPT" ]; then
  V32="$(mktemp -d 2>/dev/null || echo "${TEMP:-/tmp}/rice_v32_selftest_$$")"
  mkdir -p "$V32"

  if "$PYTHON" scripts/proteomics/04_uniprot_mapping.py --demo --outdir "$V32/map" >/dev/null 2>&1; then
    pass "UniProt mapping (04_uniprot_mapping.py --demo)"
  else
    fail "UniProt mapping (04_uniprot_mapping.py --demo)"
  fi

  if "$PYTHON" scripts/proteomics/05_interpro_domains.py --demo --outdir "$V32/ipr" >/dev/null 2>&1; then
    pass "InterPro domain enrichment (05_interpro_domains.py --demo)"
  else
    fail "InterPro domain enrichment (05_interpro_domains.py --demo)"
  fi

  if "$RSCRIPT" scripts/phenomics/05_power_analysis.R --demo --outdir "$V32/power" >/dev/null 2>&1; then
    pass "power analysis (05_power_analysis.R --demo)"
  else
    fail "power analysis (05_power_analysis.R --demo)"
  fi

  if "$RSCRIPT" scripts/qtl/03_coloc_bayesian.R --demo --outdir "$V32/coloc" >/dev/null 2>&1; then
    pass "Bayesian coloc (03_coloc_bayesian.R --demo)"
  else
    fail "Bayesian coloc (03_coloc_bayesian.R --demo)"
  fi

  if "$RSCRIPT" scripts/gwas/08_2d_association.R --demo --outdir "$V32/2d" >/dev/null 2>&1; then
    pass "2D-GWAS (08_2d_association.R --demo)"
  else
    fail "2D-GWAS (08_2d_association.R --demo)"
  fi

  if "$PYTHON" scripts/integration/07_grn_boost.py --demo --outdir "$V32/grn" >/dev/null 2>&1; then
    pass "GRN inference (07_grn_boost.py --demo)"
  else
    fail "GRN inference (07_grn_boost.py --demo)"
  fi

  if "$PYTHON" scripts/evolution/02_mirna_targets.py --demo --outdir "$V32/mirna" >/dev/null 2>&1; then
    pass "miRNA targets (02_mirna_targets.py --demo)"
  else
    fail "miRNA targets (02_mirna_targets.py --demo)"
  fi

  if "$PYTHON" scripts/tf/02_motif_tf_match.py --demo --outdir "$V32/motif" >/dev/null 2>&1; then
    pass "motif-TF matching (02_motif_tf_match.py --demo)"
  else
    fail "motif-TF matching (02_motif_tf_match.py --demo)"
  fi

  if "$PYTHON" scripts/qtl/04_qtl_region_snps.py --demo --outdir "$V32/region" >/dev/null 2>&1; then
    pass "QTL region SNP extraction (04_qtl_region_snps.py --demo)"
  else
    fail "QTL region SNP extraction (04_qtl_region_snps.py --demo)"
  fi

  if "$PYTHON" scripts/qtl/05_bsa_scan.py --demo --outdir "$V32/bsa" >/dev/null 2>&1; then
    pass "genome-wide BSA scan (05_bsa_scan.py --demo)"
  else
    fail "genome-wide BSA scan (05_bsa_scan.py --demo)"
  fi

  rm -rf "$V32" 2>/dev/null
else
  skip "v3.2 smoke tests" "python or Rscript not available"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo
echo "=== Summary: $N_PASS passed, $N_FAIL failed, $N_SKIP skipped ==="
if [ "$N_FAIL" -gt 0 ]; then exit 1; fi
exit 0
