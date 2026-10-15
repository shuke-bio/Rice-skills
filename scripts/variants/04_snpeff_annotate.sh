#!/usr/bin/env bash
#
# 04_snpeff_annotate.sh -- VCF -> functionally annotated VCF (SnpEff)
#
# Description:
#   Functional annotation of variants with SnpEff:
#     - uses a pre-built database (default Oryza_sativa, auto-downloaded
#       on first use), or
#     - builds a custom database from your own GFF3 + FASTA
#       (--build-custom, e.g. RAP-DB/Phytozome annotations).
#   Also extracts a flat annotation table (gene / effect / impact per
#   variant) and a summary report for convenient downstream filtering.
#   If SnpEff is not installed the script prints install instructions
#   and exits.
#
# Usage:
#   # standard: pre-built rice database
#   bash scripts/variants/04_snpeff_annotate.sh \
#        --vcf results/variants/03_filter/filtered.vcf.gz \
#        --outdir results/variants/04_annot
#
#   # first time: download the rice database only
#   bash scripts/variants/04_snpeff_annotate.sh --download-db --outdir data/snpeff
#
#   # custom annotation (e.g. RAP-DB IRGSP-1.0):
#   bash scripts/variants/04_snpeff_annotate.sh \
#        --vcf filtered.vcf.gz --outdir results/variants/04_annot \
#        --build-custom --gff IRGSP-1.0_annotation.gff3 --fasta IRGSP-1.0.fa \
#        --db-name rapdb_irgsp1
#
# Inputs:
#   --vcf        input VCF (.vcf.gz)
#   --outdir     output directory (default results/variants/04_annot)
#   --db         SnpEff database name (default Oryza_sativa)
#   --download-db  only download/build the database, do not annotate
#   --build-custom  build a custom database from --gff + --fasta
#   --gff, --fasta  inputs for --build-custom
#   --db-name    custom database name (with --build-custom)
#   --java-mem   Java heap (default 8g)
# Outputs:
#   <outdir>/annotated.vcf.gz (+ .tbi)       VCF with ANN/LOF annotations
#   <outdir>/ann_table.tsv                   flat table: CHROM POS REF ALT
#                                            GENE EFFECT IMPACT HGVS.c HGVS.p
#   <outdir>/snpeff_summary.html             SnpEff statistics report
#   <outdir>/annotate_log.txt                full command log
#
# Notes:
#   - Install: conda install -c bioconda snpeff   (bundles Java)
#   - dbNSFP / ClinVar-style human resources are not applicable to rice;
#     effect prediction (missense/splice/stop) is what matters here.
#   - Useful downstream filters on ann_table.tsv:
#       * HIGH impact  : grep -P "\tHIGH\t"  (stop-gained, splice donor/acceptor, frameshift)
#       * by gene      : grep "Os01g0100100" or RAP-DB gene IDs

set -u

VCF=""
OUTDIR="results/variants/04_annot"
DB="Oryza_sativa"
DOWNLOAD_ONLY=0
BUILD_CUSTOM=0
GFF=""
FASTA=""
DB_NAME=""
JAVA_MEM="8g"

while [ $# -gt 0 ]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --db) DB="$2"; shift 2 ;;
    --download-db) DOWNLOAD_ONLY=1; shift ;;
    --build-custom) BUILD_CUSTOM=1; shift ;;
    --gff) GFF="$2"; shift 2 ;;
    --fasta) FASTA="$2"; shift 2 ;;
    --db-name) DB_NAME="$2"; shift 2 ;;
    --java-mem) JAVA_MEM="$2"; shift 2 ;;
    -h|--help) sed -n '2,55p' "$0"; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

# --- locate SnpEff ------------------------------------------------------------
SNPEFF=""
if command -v snpEff >/dev/null 2>&1; then SNPEFF="snpEff"
elif [ -n "${SNPEFF_JAR:-}" ] && [ -f "$SNPEFF_JAR" ]; then SNPEFF="java -Xmx$JAVA_MEM -jar $SNPEFF_JAR"
elif [ -f "$HOME/snpEff/snpEff.jar" ]; then SNPEFF="java -Xmx$JAVA_MEM -jar $HOME/snpEff/snpEff.jar"
else
  echo "SnpEff not found."
  echo ""
  echo "Install options:"
  echo "  1. conda install -c bioconda snpeff       (recommended; bundles Java)"
  echo "  2. Download: https://pcingola.github.io/SnpEff/download/"
  echo "     then set: export SNPEFF_JAR=/path/to/snpEff.jar"
  exit 0
fi
echo "Using SnpEff: $SNPEFF"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/annotate_log.txt"
: > "$LOG"
run() { echo "+ $*" | tee -a "$LOG"; bash -c "$*" >>"$LOG" 2>&1; }

# --- database download / build ------------------------------------------------
if [ "$BUILD_CUSTOM" -eq 1 ]; then
  [ -n "$GFF" ] && [ -n "$FASTA" ] && [ -n "$DB_NAME" ] || {
    echo "ERROR: --build-custom needs --gff, --fasta and --db-name."; exit 1; }
  for f in "$GFF" "$FASTA"; do [ -f "$f" ] || { echo "ERROR: not found: $f"; exit 1; }; done
  SNPEFF_DATA_DIR="$OUTDIR/snpeff_data"
  mkdir -p "$SNPEFF_DATA_DIR/$DB_NAME"
  echo "--> copying FASTA"
  cp "$FASTA" "$SNPEFF_DATA_DIR/$DB_NAME/sequences.fa"
  echo "--> copying GFF3"
  cp "$GFF" "$SNPEFF_DATA_DIR/$DB_NAME/genes.gff"
  echo "--> building custom database: $DB_NAME"
  run "$SNPEFF build -gff3 -v -dataDir $SNPEFF_DATA_DIR $DB_NAME"
  DB="$DB_NAME"
  DATA_DIR_ARG="-dataDir $SNPEFF_DATA_DIR"
else
  # check database availability; download if missing
  if ! $SNPEFF databases 2>/dev/null | grep -q "^$DB"; then
    echo "--> downloading SnpEff database: $DB (first run only)"
    run "$SNPEFF download $DB"
    if ! $SNPEFF databases 2>/dev/null | grep -q "$DB"; then
      echo "WARNING: could not verify database '$DB'. Check spelling with:"
      echo "  snpEff databases | grep -i oryza"
      echo "Common rice databases: Oryza_sativa (Ensembl), Oryza_indica"
    fi
  else
    echo "--> database '$DB' already installed"
  fi
  DATA_DIR_ARG=""
fi

if [ "$DOWNLOAD_ONLY" -eq 1 ]; then
  echo "Database ready. Re-run with --vcf to annotate."
  echo "Tip: list available rice databases with:  snpEff databases | grep -i oryza"
  exit 0
fi

[ -n "$VCF" ] || { echo "ERROR: --vcf is required (or use --download-db). See --help."; exit 1; }
[ -f "$VCF" ] || { echo "ERROR: VCF not found: $VCF"; exit 1; }

# --- annotation ---------------------------------------------------------------
echo "--> annotating with $DB"
run "$SNPEFF -v $DATA_DIR_ARG -stats $OUTDIR/snpeff_summary.html \
    $DB '$VCF' | bgzip -c > $OUTDIR/annotated.vcf.gz"

if ! command -v bgzip >/dev/null 2>&1; then
  # bgzip unavailable: rerun without compression pipeline
  echo "    (bgzip not found — writing plain VCF and gzip-tabix-free fallback)"
  run "$SNPEFF -v $DATA_DIR_ARG -stats $OUTDIR/snpeff_summary.html $DB '$VCF' \
      > $OUTDIR/annotated.vcf"
fi

if command -v bcftools >/dev/null 2>&1 && [ -f "$OUTDIR/annotated.vcf.gz" ]; then
  run "bcftools index -t -f $OUTDIR/annotated.vcf.gz"
fi

# --- flat annotation table ----------------------------------------------------
echo "--> extracting flat ANN table"
ANNOT_VCF="$OUTDIR/annotated.vcf.gz"
[ -f "$ANNOT_VCF" ] || ANNOT_VCF="$OUTDIR/annotated.vcf"

# prefer SnpSift if available
if command -v snpSift >/dev/null 2>&1; then
  run "snpSift extractFields -s ',' -e '.' '$ANNOT_VCF' \
      CHROM POS REF ALT 'ANN[*].GENE' 'ANN[*].EFFECT' 'ANN[*].IMPACT' \
      'ANN[*].HGVS_C' 'ANN[*].HGVS_P' > $OUTDIR/ann_table.tsv"
else
  # minimal parser for the first ANN entry per record
  python3 - "$ANNOT_VCF" "$OUTDIR/ann_table.tsv" <<'PYEOF' || true
import sys, gzip
inp, out = sys.argv[1], sys.argv[2]
op = gzip.open if inp.endswith(".gz") else open
with op(inp, "rt") as fi, open(out, "w") as fo:
    fo.write("CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\tHGVS_C\tHGVS_P\n")
    for line in fi:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        info = f[7]
        ann = ""
        for kv in info.split(";"):
            if kv.startswith("ANN="):
                ann = kv[4:].split(",")[0]
                break
        if not ann:
            fo.write("\t".join(f[:5]) + "\t.\t.\t.\t.\t.\n")
            continue
        a = ann.split("|")
        gene, eff, imp = (a[3] if len(a) > 3 else ""), (a[1] if len(a) > 1 else ""), (a[2] if len(a) > 2 else "")
        hc = a[9] if len(a) > 9 else ""
        hp = a[10] if len(a) > 10 else ""
        fo.write("\t".join(f[:4] + [gene, eff, imp, hc, hp]) + "\n")
print("wrote", out)
PYEOF
fi

if [ -f "$OUTDIR/ann_table.tsv" ]; then
  N=$(($(wc -l < "$OUTDIR/ann_table.tsv") - 1))
  HIGH=$(grep -c -P "\tHIGH\t" "$OUTDIR/ann_table.tsv" 2>/dev/null || echo "?")
  echo "  annotated variants : $N"
  echo "  HIGH impact calls  : $HIGH   (stop-gained / frameshift / splice)"
fi

echo ""
echo "Annotation complete."
echo "  VCF    : $OUTDIR/annotated.vcf.gz"
echo "  Table  : $OUTDIR/ann_table.tsv"
echo "  Report : $OUTDIR/snpeff_summary.html"
echo ""
echo "Downstream: filter HIGH-impact genes for candidate lists, or join with"
echo "GWAS candidate regions (scripts/gwas/04_candidate_loci.R)."
