#!/usr/bin/env bash
#
# 02_sv_detect.sh -- structural variant (SV) detection
#
# Description:
#   Two routes:
#   [assembly]     Assembly-based SV detection with the minimap2 + SYRI standard
#                  route (default): minimap2 (asm5/asm10/asm20) -> samtools sort
#                  -> syri -c -> per-sample SV table. Also documents the
#                  graph-based alternative (genotyping via the cactus VCF).
#   [read-based]   Read-based SV detection + genotyping with HitSV (BAM input):
#                  pure LRS (ONT/PacBio-HiFi), pure SRS (Illumina/BGI-T7), or
#                  hybrid (low-depth LRS + SRS). Outputs a genotyped VCF that
#                  feeds scripts/pangenome/03_pav_analysis.R / 04_pan_gwas.R.
#
#   HitSV is GPL-3.0 and intentionally NOT bundled in this MIT-licensed repo.
#   Install it at runtime (prebuilt static binary, no compilation needed):
#     https://github.com/hitbc/HitSV/releases
#   Paper: Li G et al., Nat Biotechnol (2026), doi:10.1038/s41587-026-03332-9
#   Method notes: docs/references/hitsv_sv_calling.md
#
#   If the required tools are not installed the script prints install
#   instructions and exits (status 0).
#
# Usage:
#   # assembly route: single query assembly vs reference
#   bash scripts/pangenome/02_sv_detect.sh --ref IRGSP-1.0.fa \
#        --query sample_asm.fa --sample sampleA \
#        --outdir results/pangenome/02_sv --threads 16
#
#   # assembly route: batch (directory / list of query assemblies)
#   bash scripts/pangenome/02_sv_detect.sh --ref IRGSP-1.0.fa \
#        --queries assemblies/ --outdir results/pangenome/02_sv
#
#   # read-based route: pure short-read (Illumina/BGI-T7)
#   bash scripts/pangenome/02_sv_detect.sh --mode read-based --seq srs \
#        --ref IRGSP-1.0.fa --bam sample.sorted.bam --sample sampleA
#
#   # read-based route: pure long-read (PacBio HiFi / ONT)
#   bash scripts/pangenome/02_sv_detect.sh --mode read-based --seq lrs \
#        --ref IRGSP-1.0.fa --bam sample.hifi.bam --platform HIFI \
#        --aln-preset asm10 --sample sampleA
#
#   # read-based route: hybrid (low-depth LRS + SRS)
#   bash scripts/pangenome/02_sv_detect.sh --mode read-based --seq hybrid \
#        --ref IRGSP-1.0.fa --lrs-bam lowcov.ont.bam --srs-bam sample.srs.bam \
#        --sample sampleA
#
#   # synthetic end-to-end self-test (read-based LRS route, no input data)
#   bash scripts/pangenome/02_sv_detect.sh --demo
#
# Inputs (assembly route):
#   --ref      reference genome FASTA
#   --query    single query assembly FASTA  (or use --queries)
#   --queries  directory of FASTA files or a file list
#   --sample   sample name (single-query mode; default from filename)
#   --outdir   output directory (default results/pangenome/02_sv)
#   --threads  threads for minimap2/samtools (default 16)
#   --preset   minimap2 preset: asm5 | asm10 | asm20 (default asm10)
#
# Inputs (read-based route):
#   --mode          assembly (default) | read-based
#   --seq           lrs (default) | srs | hybrid
#   --ref           reference genome FASTA (SRS additionally needs a HitSV
#                   repeat-stat file, generated automatically)
#   --bam           SRS BAM (--seq srs) or LRS BAM (--seq lrs);
#                   must be coordinate-sorted and indexed
#   --lrs-bam       LRS BAM (--seq hybrid)
#   --srs-bam       SRS BAM (--seq hybrid)
#   --platform      LRS error-mode preset: HIFI | ONT_Q20 | ERR_PRONE | ASM
#                   (default HIFI)
#   --aln-preset    contig-to-reference alignment preset: asm5 | asm10 | asm20
#                   (default asm10; rice intraspecies is a safe default,
#                   >5% divergence use asm20)
#   --min-sv-len    minimum reported SV length in bp (default 50)
#   --sample        sample name (default from BAM filename)
#   --demo          run the synthetic self-test instead (forces read-based)
#
# Outputs (assembly route):
#   <outdir>/<sample>.sv.syri.out       SV table (per SV: type, coords, size)
#   <outdir>/<sample>.syri.summary      SV counts by type
#   <outdir>/<sample>.bam (+ .bai)      whole-genome alignment
#   <outdir>/sv_all.tsv                 merged SV table (all samples)
#   <outdir>/sv_log.txt                 full command log
#
# Outputs (read-based route):
#   <outdir>/<sample>.hitsv.vcf         SV callset with genotypes (VCF)
#   <outdir>/<sample>.sv_summary.tsv    SV counts by type
#   <outdir>/<sample>.ref.stat.txt      HitSV reference repeat-stat file
#   <outdir>/<sample>.tl.sort.bam       preprocessed SRS reads (SRS/hybrid)
#   <outdir>/sv_log.txt                 full command log
#   <outdir>/run_info.txt               provenance (params, tool paths, license note)
#
# Notes:
#   - install: conda install -c bioconda minimap2 syri samtools
#   - HitSV  : https://github.com/hitbc/HitSV (GPL-3.0; release binary or make)
#   - SYRI requires the alignment BAM sorted and indexed; both done here.
#   - asm5 = <1% divergence, asm10 = ~1-5%, asm20 = >5% (rice intraspecies
#     assemblies: asm10 is a safe default; same-species polished: asm5).
#   - Graph route (cohort with pangenome graph): genotype SVs directly from
#     the cactus VCF (already produced by 01_pangenome_graph.sh
#     --tool minigraph-cactus), then run scripts/pangenome/04_pan_gwas.R.
#   - Next step: bash scripts/pangenome/04_pan_gwas.R (SV/PAV association)
#   - read-based route: region-level runs and the population CSV workflows
#     (clustering / MEI-TR / nested variants) are available in the HitSV
#     binary directly; this wrapper covers whole-genome per-sample calling.

set -u

MODE="assembly"        # assembly | read-based
SEQ="lrs"              # read-based: lrs | srs | hybrid
BAM=""
LRS_BAM=""
SRS_BAM=""
PLATFORM="HIFI"
ALN_PRESET="asm10"
MIN_SV_LEN=50
DEMO=0
REF=""
QUERY=""
QUERIES=""
SAMPLE=""
OUTDIR="results/pangenome/02_sv"
THREADS=16
PRESET="asm10"

while [ $# -gt 0 ]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --seq) SEQ="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --lrs-bam) LRS_BAM="$2"; shift 2 ;;
    --srs-bam) SRS_BAM="$2"; shift 2 ;;
    --platform) PLATFORM="$2"; shift 2 ;;
    --aln-preset) ALN_PRESET="$2"; shift 2 ;;
    --min-sv-len) MIN_SV_LEN="$2"; shift 2 ;;
    --demo) DEMO=1; shift ;;
    --ref) REF="$2"; shift 2 ;;
    --query) QUERY="$2"; shift 2 ;;
    --queries) QUERIES="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --preset) PRESET="$2"; shift 2 ;;
    -h|--help) sed -n '2,/^set -u$/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "ERROR: unknown option '$1'"; exit 1 ;;
  esac
done

case "$MODE" in
  assembly|read-based) ;;
  *) echo "ERROR: --mode must be 'assembly' or 'read-based' (got '$MODE')"; exit 1 ;;
esac
case "$SEQ" in
  lrs|srs|hybrid) ;;
  *) echo "ERROR: --seq must be 'lrs', 'srs' or 'hybrid' (got '$SEQ')"; exit 1 ;;
esac
case "$PLATFORM" in
  HIFI|ONT_Q20|ERR_PRONE|ASM) ;;
  *) echo "ERROR: --platform must be HIFI | ONT_Q20 | ERR_PRONE | ASM (got '$PLATFORM')"; exit 1 ;;
esac
case "$ALN_PRESET" in
  asm5|asm10|asm20) ;;
  *) echo "ERROR: --aln-preset must be asm5 | asm10 | asm20 (got '$ALN_PRESET')"; exit 1 ;;
esac
case "$MIN_SV_LEN" in
  ''|*[!0-9]*) echo "ERROR: --min-sv-len must be a positive integer in bp (got '$MIN_SV_LEN')"; exit 1 ;;
esac

# ===========================================================================
# read-based route (HitSV) ---------------------------------------------------
# ===========================================================================
read_based_tools_missing() {
  # echoes the list of missing required tools; sets up install hints once
  MISSING=""
  command -v HitSV >/dev/null 2>&1    || MISSING="$MISSING HitSV"
  command -v samtools >/dev/null 2>&1 || MISSING="$MISSING samtools"
  command -v minimap2 >/dev/null 2>&1 || MISSING="$MISSING minimap2(demo)"
}

print_read_based_install() {
  echo ""
  echo "HitSV (GPL-3.0) is intentionally NOT bundled in this MIT repo."
  echo "Install at runtime:"
  echo "  prebuilt static binary: https://github.com/hitbc/HitSV/releases"
  echo "  from source           : git clone https://github.com/hitbc/HitSV \\"
  echo "                           cd HitSV && make   # needs g++/zlib"
  echo "samtools / minimap2    : conda install -c bioconda samtools minimap2"
  echo "Reference              : Li G et al., Nat Biotechnol (2026),"
  echo "                          doi:10.1038/s41587-026-03332-9"
  echo "Method notes           : docs/references/hitsv_sv_calling.md"
}

# summarize SVTYPE counts from a HitSV VCF -> TSV
sv_summary() {
  local vcf="$1" out="$2"
  {
    printf "svtype\tcount\n"
    grep -v '^#' "$vcf" 2>/dev/null | \
      awk 'BEGIN{FS=OFS="\t"}
           {type="other"; n=split($8,info,";");
            for(i=1;i<=n;i++) if(info[i] ~ /^SVTYPE=/){type=substr(info[i],8)}
            cnt[type]++}
           END{for(t in cnt) print t, cnt[t]}' | sort
  } > "$out"
}

hitsv_call() {  # args: sample ref vcf_out lrs_bam srs_bam
  local smp="$1" ref="$2" vcf="$3"
  local lrs="$4" srs="$5"
  local stat tl tl_sort

  if [ ! -f "${ref}.fai" ]; then run "samtools faidx '$ref'"; fi

  # NOTE: path elements are single-quoted inside the array so that
  # `run "HitSV call ${args[*]}"` survives paths with spaces (run passes the
  # string to bash -c). Non-path elements (presets, ints) stay unquoted.
  local args=(-r "'$ref'" -P "$ALN_PRESET" -m "$MIN_SV_LEN" -o "'$vcf'")

  if [ "$SEQ" = "lrs" ] || [ "$SEQ" = "hybrid" ]; then
    args+=(-p "$PLATFORM")
    args+=(-l "'$lrs'")
  fi
  if [ "$SEQ" = "srs" ] || [ "$SEQ" = "hybrid" ]; then
    # reference repeat-stat file (needed for SRS/hybrid calling)
    stat="$OUTDIR/${smp}.ref.stat.txt"
    run "HitSV srs_fa_stat '$ref' > '$stat'"
    # optional read preprocessing (documented SRS whole-genome flow)
    tl="$OUTDIR/${smp}.tl.bam"
    tl_sort="$OUTDIR/${smp}.tl.sort.bam"
    if [ -n "$srs" ]; then
      run "HitSV srs_trans_reads '$ref' '$srs' '$tl'"
      run "samtools sort -@ $THREADS --output-fmt=BAM -o '$tl_sort' '$tl'"
      run "samtools index '$tl_sort'"
    fi
    args+=(-n "'$srs'")
    if [ -f "$tl_sort" ]; then args+=(-L "'$tl_sort'"); fi
    args+=(-I "'$stat'")
  fi

  echo "--> [$smp] HitSV call (seq=$SEQ platform=$PLATFORM aln=$ALN_PRESET min_sv_len=$MIN_SV_LEN)"
  run "HitSV call ${args[*]}"

  if [ ! -s "$vcf" ]; then
    echo "ERROR: HitSV produced no output ($vcf). Check sv_log.txt."
    exit 1
  fi
  sv_summary "$vcf" "$OUTDIR/${smp}.sv_summary.tsv"
  echo "    SV summary:"
  sed 's/^/      /' "$OUTDIR/${smp}.sv_summary.tsv"
}

if [ "$MODE" = "read-based" ]; then
  read_based_tools_missing
  if [ -n "$MISSING" ]; then
    echo "Required tool(s) not found:$MISSING"
    print_read_based_install
    exit 0
  fi

  DEMO_DIR=""
  mkdir -p "$OUTDIR"
  LOG="$OUTDIR/sv_log.txt"
  : > "$LOG"
  run() { echo "+ $*" | tee -a "$LOG"; bash -c "$*" >>"$LOG" 2>&1; }

  if [ "$DEMO" = "1" ]; then
    mkdir -p "$OUTDIR/demo_synthetic"
    DEMO_DIR="$(cd "$OUTDIR/demo_synthetic" && pwd)"
    echo "==> synthetic self-test (read-based LRS route)"

    # -- generate a synthetic reference + a haplotype carrying a known
    #    500 bp deletion @100000 and a 1 kb novel insertion @200000,
    #    then simulate ~10x 10 kb HiFi-like reads from that haplotype
    python - "$DEMO_DIR" <<'PYEOF'
import random, sys, os
outdir = sys.argv[1]
random.seed(42)
CHR_LEN      = 300_000
DEL_START    = 100_000
DEL_LEN      = 500
INS_AT       = 200_000   # reference coordinate of the insertion point
INS_LEN      = 1000
READ_LEN     = 10_000
COV          = 10

ref = "".join(random.choice("ACGT") for _ in range(CHR_LEN))
ins = "".join(random.choice("ACGT") for _ in range(INS_LEN))

with open(os.path.join(outdir, "demo_ref.fa"), "w") as f:
    f.write(">chrSVdemo\n")
    for i in range(0, CHR_LEN, 60):
        f.write(ref[i:i+60] + "\n")

# haplotype: reference carrying the deletion, then the insertion
hap = ref[:DEL_START] + ref[DEL_START+DEL_LEN:]
ins_pos_hap = INS_AT - DEL_LEN           # coordinate shift after the deletion
hap = hap[:ins_pos_hap] + ins + hap[ins_pos_hap:]

with open(os.path.join(outdir, "demo_reads.fq"), "w") as f:
    for i in range(int(CHR_LEN * COV / READ_LEN)):
        start = random.randrange(0, len(hap) - READ_LEN)
        seq = hap[start:start+READ_LEN]
        f.write("@demo_read_%d\n%s\n+\n%s\n" % (i, seq, "I"*READ_LEN))

with open(os.path.join(outdir, "demo_truth.tsv"), "w") as f:
    f.write("DEL\t%d\t%d\n" % (DEL_START, DEL_START+DEL_LEN))
    f.write("INS\t%d\t%d\n" % (INS_AT, INS_AT))
print("[demo] reference 300 kb; haplotype carries DEL %d-%d + INS %d bp @%d"
      % (DEL_START, DEL_START+DEL_LEN, INS_LEN, INS_AT))
PYEOF

    REF="$DEMO_DIR/demo_ref.fa"
    BAM="$DEMO_DIR/demo.hifi.bam"
    SAMPLE="demo"
    SEQ="lrs"
    echo "--> [demo] aligning synthetic reads (minimap2 -ax map-hifi)"
    run "minimap2 -ax map-hifi -t $THREADS '$REF' '$DEMO_DIR/demo_reads.fq' \
        | samtools sort -@ $THREADS -o '$BAM' -"
    run "samtools index '$BAM'"
  else
    [ -n "$REF" ] || { echo "ERROR: --ref is required (read-based route). See --help."; exit 1; }
    [ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
    if [ "$SEQ" = "hybrid" ]; then
      [ -n "$LRS_BAM" ] && [ -n "$SRS_BAM" ] || {
        echo "ERROR: --seq hybrid requires --lrs-bam and --srs-bam."; exit 1; }
    else
      [ -n "$BAM" ] || { echo "ERROR: --bam is required (read-based route). See --help."; exit 1; }
      [ -f "$BAM" ] || { echo "ERROR: BAM not found: $BAM"; exit 1; }
    fi
    for b in "$BAM" "$LRS_BAM" "$SRS_BAM"; do
      [ -z "$b" ] && continue
      [ -f "$b" ] || { echo "ERROR: BAM not found: $b"; exit 1; }
      if [ ! -f "$b.bai" ] && [ ! -f "${b%.bam}.bai" ]; then
        echo "NOTE: no index for $b — creating it."
        run "samtools index '$b'"
      fi
    done
    [ -n "$SAMPLE" ] || SAMPLE="$(basename "${BAM:-$SRS_BAM}" | sed 's/\.bam$//')"
  fi

  VCF="$OUTDIR/${SAMPLE}.hitsv.vcf"

  # HitSV needs BAMs; in hybrid mode the LRS BAM is passed separately
  LRS_IN=""; SRS_IN=""
  case "$SEQ" in
    lrs)     LRS_IN="$BAM" ;;
    srs)     SRS_IN="$BAM" ;;
    hybrid)  LRS_IN="$LRS_BAM"; SRS_IN="$SRS_BAM" ;;
  esac

  hitsv_call "$SAMPLE" "$REF" "$VCF" "$LRS_IN" "$SRS_IN"

  # provenance (project contract: every script writes run_info.txt) — written
  # before the demo verification block so it exists regardless of outcome
  {
    echo "script    : 02_sv_detect.sh (mode=read-based seq=$SEQ)"
    echo "date      : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "reference : $REF"
    echo "bam_lrs   : ${LRS_IN:-NA}"
    echo "bam_srs   : ${SRS_IN:-NA}"
    echo "platform  : $PLATFORM | aln_preset: $ALN_PRESET | min_sv_len: $MIN_SV_LEN"
    echo "hitsv     : $(command -v HitSV || echo NA)"
    echo "samtools  : $(samtools --version 2>/dev/null | head -1 || echo NA)"
    echo "vcf       : $VCF"
    echo "license   : HitSV is GPL-3.0 (runtime install, not bundled)."
    echo "            Source: https://github.com/hitbc/HitSV"
    echo "paper     : Li G et al., Nat Biotechnol (2026), doi:10.1038/s41587-026-03332-9"
    if [ "$DEMO" = "1" ]; then
      echo "demo      : synthetic self-test; truth: $(cat "$DEMO_DIR/demo_truth.tsv" | tr '\n' ' ')"
    fi
  } > "$OUTDIR/run_info.txt"

  echo ""
  echo "Read-based SV calling complete."
  echo "  callset : $VCF"
  echo "  summary : $OUTDIR/${SAMPLE}.sv_summary.tsv"
  echo ""
  echo "Next step:"
  echo "  1. genotype matrix: bcftools query -f '[%GT\t]\n' '$VCF' | ... -> 0/1 matrix"
  echo "  2. PAV/association: Rscript scripts/pangenome/03_pav_analysis.R --pav <matrix.csv> ..."
  echo "     or Pan-GWAS: Rscript scripts/pangenome/04_pan_gwas.R --pav <matrix.csv> ..."

  # ---- synthetic self-test verification -----------------------------------
  if [ "$DEMO" = "1" ]; then
    DEL_START=100000; DEL_END=100500; TOL=2000
    DEL_HIT=$(awk -v ds="$DEL_START" -v de="$DEL_END" -v tol="$TOL" '
      /^#/ {next}
      {
        type=""; end=0;
        n=split($8, info, ";");
        for (i=1; i<=n; i++) {
          if (info[i] ~ /^SVTYPE=/) type = substr(info[i], 8);
          if (info[i] ~ /^END=/)    end  = substr(info[i], 5) + 0;
        }
        if (type == "DEL" && end == 0) end = $2 + 50;
        if (type == "DEL" && $2 <= de + tol && end >= ds - tol) { print; exit }
      }' "$VCF")
    SV_TOTAL=$(grep -vc '^#' "$VCF" 2>/dev/null)
    [ -n "$SV_TOTAL" ] || SV_TOTAL=0
    INS_NEAR=$(awk -v t=200000 -v tol=3000 '
      /^#/ {next}
      {
        type="";
        n=split($8, info, ";");
        for (i=1; i<=n; i++) if (info[i] ~ /^SVTYPE=/) type = substr(info[i], 8);
        if ((type == "INS" || type == "DUP") && ($2 > t - tol && $2 < t + tol)) { print; exit }
      }' "$VCF")
    echo ""
    echo "==> self-test verification (truth: DEL $DEL_START-$DEL_END, INS@200000)"
    echo "    total SV records : $SV_TOTAL"
    [ -n "$DEL_HIT" ] && echo "    DEL overlap test : PASS" \
                       || echo "    DEL overlap test : FAIL"
    if [ -n "$INS_NEAR" ]; then
      echo "    INS near-truth   : PASS (also detected)"
    else
      echo "    INS near-truth   : not called (soft check; HitSV may report it"
      echo "                       as another type or merge it — informational only)"
    fi
    if [ -n "$DEL_HIT" ]; then
      echo "==> DEMO PASS"
      exit 0
    else
      echo "==> DEMO FAIL (no DEL overlapping the synthetic deletion;"
      echo "    inspect $VCF and $OUTDIR/sv_log.txt)"
      exit 1
    fi
  fi
  exit 0
fi

# ===========================================================================
# assembly route (default) ---------------------------------------------------
# ===========================================================================

[ -n "$REF" ] || { echo "ERROR: --ref is required. See --help."; exit 1; }
[ -f "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
if [ -z "$QUERY" ] && [ -z "$QUERIES" ]; then
  echo "ERROR: provide --query <assembly.fa> or --queries <dir|list>."
  exit 1
fi

# --- locate tools -------------------------------------------------------------
MISSING=""
command -v minimap2 >/dev/null 2>&1 || MISSING="$MISSING minimap2"
command -v samtools >/dev/null 2>&1 || MISSING="$MISSING samtools"
command -v syri >/dev/null 2>&1 || MISSING="$MISSING syri"
if [ -n "$MISSING" ]; then
  echo "Required tool(s) not found:$MISSING"
  echo ""
  echo "Install: conda install -c bioconda minimap2 syri samtools"
  echo "Docs  : https://github.com/sschmidtu/SYRI (syri), https://github.com/lh3/minimap2"
  echo ""
  echo "Alternative (graph route): if you built a pangenome with"
  echo "  scripts/pangenome/01_pangenome_graph.sh --tool minigraph-cactus"
  echo "the emitted VCF already contains genotyped SVs — skip this script and"
  echo "go straight to scripts/pangenome/04_pan_gwas.R."
  exit 0
fi

# collect query list
QRY_LIST=()
if [ -n "$QUERY" ]; then
  QRY_LIST=("$QUERY")
elif [ -d "$QUERIES" ]; then
  for f in "$QUERIES"/*.fa "$QUERIES"/*.fasta "$QUERIES"/*.fa.gz; do
    [ -f "$f" ] && QRY_LIST+=("$f")
  done
elif [ -f "$QUERIES" ]; then
  while IFS= read -r line; do [ -n "$line" ] && [ -f "$line" ] && QRY_LIST+=("$line"); done < "$QUERIES"
fi
[ "${#QRY_LIST[@]}" -ge 1 ] || { echo "ERROR: no query assemblies found."; exit 1; }
echo "Queries: ${#QRY_LIST[@]} assembly file(s); preset=$PRESET"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/sv_log.txt"
: > "$LOG"
run() { echo "+ $*" | tee -a "$LOG"; bash -c "$*" >>"$LOG" 2>&1; }

# reference index for syri (needs .fai and optionally minimap2 index)
[ -f "${REF}.fai" ] || run "samtools faidx '$REF'"

ALL_SV="$OUTDIR/sv_all.tsv"
: > "$ALL_SV"

for q in "${QRY_LIST[@]}"; do
  s="${SAMPLE:-$(basename "$q" | sed 's/\.\(fa\|fasta\|fa\.gz\)$//')}"
  [ "${#QRY_LIST[@]}" -gt 1 ] && s=$(basename "$q" | sed 's/\.\(fa\|fasta\|fa\.gz\)$//')
  echo "--> [$s] aligning to reference (minimap2 -x $PRESET)"
  run "minimap2 -ax $PRESET -t $THREADS '$REF' '$q' \
      | samtools sort -@ $THREADS -o '$OUTDIR/$s.bam' -"
  run "samtools index '$OUTDIR/$s.bam'"

  echo "--> [$s] syri structural-variant calling"
  run "cd '$OUTDIR' && syri -c '$s.bam' -r '$REF' -q '$q' -F B --prefix '$s.'"

  # collect per-sample SV tables
  SVOUT="$OUTDIR/${s}.sv.syri.out"
  if [ -f "$SVOUT" ]; then
    awk -v smp="$s" 'BEGIN{FS=OFS="\t"} NR==1{print "sample",$0} NR>1{print smp,$0}' \
      "$SVOUT" >> "$ALL_SV"
  fi
  SUM="$OUTDIR/${s}.syri.summary"
  [ -f "$SUM" ] && { echo "    SV summary:"; sed 's/^/      /' "$SUM" | head -12; }
done

echo ""
echo "SV calling complete."
echo "  per-sample: $OUTDIR/<sample>.sv.syri.out"
echo "  merged    : $ALL_SV"
echo ""
echo "Next step:"
echo "  1. build an SV/PAV matrix (presence-absence per marker per sample), e.g."
echo "     from the merged table or from the cactus VCF;"
echo "  2. association: Rscript scripts/pangenome/04_pan_gwas.R --geno sv_matrix.csv \\"
echo "         --snp-info sv_info.csv --phenotype phenotype.csv --trait <name>"
