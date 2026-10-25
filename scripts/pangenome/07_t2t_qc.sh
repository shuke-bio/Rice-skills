#!/usr/bin/env bash
###############################################################################
# 07_t2t_qc.sh -- T2T assembly QC (built-in stats + optional external tools)
#
# 内置（无需外部工具）：
#   - contig 数、总长、N50、GC%、gap(N) 数、最长 contig
# 可选（需自行安装，--run-tools 才会执行）：
#   - BUSCO     : busco -i asm -l lineage -m genome -o busco -c N
#   - Merqury   : merqury.sh reads.meryl asm prefix
#   - RepeatMasker : RepeatMasker -pa N -lib rice.lib -dir out asm
#   - Liftoff   : liftoff -g anno.gff3 -o out.gff3 asm ref.fa
#
# 用法：
#   bash scripts/pangenome/07_t2t_qc.sh --assembly asm.fa --outdir results/pangenome/07_t2t_qc
#   bash scripts/pangenome/07_t2t_qc.sh --assembly asm.fa --reads reads.fq --run-tools \
#        --busco-lineage poales_odb10 --repeat-lib rice.lib --annotation anno.gff3 --reference ref.fa
#   bash scripts/pangenome/07_t2t_qc.sh --demo --outdir results/pangenome/07_t2t_qc_demo
###############################################################################
set -euo pipefail

ASM=""; OUT="results/pangenome/07_t2t_qc"; READS=""; LINEAGE="poales_odb10"
RUN_TOOLS=0; THREADS=8; REPEAT_LIB=""; ANNO=""; REF=""; DEMO=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assembly) ASM="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --reads) READS="$2"; shift 2;;
    --busco-lineage) LINEAGE="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --repeat-lib) REPEAT_LIB="$2"; shift 2;;
    --annotation) ANNO="$2"; shift 2;;
    --reference) REF="$2"; shift 2;;
    --run-tools) RUN_TOOLS=1; shift;;
    --demo) DEMO=1; shift;;
    -h|--help) sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

mkdir -p "$OUT"

if [[ $DEMO -eq 1 ]]; then
  ASM="$OUT/demo_assembly.fa"
  {
    echo ">chr1"
    printf 'ACGT%.0s' $(seq 1 20000); echo
    printf 'N%.0s' $(seq 1 200); echo
    printf 'GCGT%.0s' $(seq 1 10000); echo
    echo ">chr2"
    printf 'TTAA%.0s' $(seq 1 5000); echo
  } > "$ASM"
  echo "[t2t] demo assembly -> $ASM"
fi

[[ -n "$ASM" && -f "$ASM" ]] || { echo "Need --assembly or --demo" >&2; exit 1; }

# ---- built-in stats ----
STATS="$OUT/assembly_stats.tsv"
awk -v out="$STATS" '
BEGIN{ print "contig\tlength\tgc_pct\tgaps" > out }
/^>/ { if (name!="") print name"\t"len"\t"(len>0?100*gc/len:0)"\t"gaps > out;
       name=substr($0,2); len=0; gc=0; gaps=0; next }
{ s=$0; n=gsub(/[Nn]/,"",s); gaps+=n; L=length($0); len+=L;
  t=$0; gsub(/[^GCgc]/,"",t); gc+=length(t) }
END{ if (name!="") print name"\t"len"\t"(len>0?100*gc/len:0)"\t"gaps > out }
' "$ASM"

TOTAL=$(awk 'NR>1{s+=$2} END{print s+0}' "$STATS")
NCONTIG=$(awk 'END{print NR-1}' "$STATS")
GC=$(awk 'NR>1{g+=$3*$2; s+=$2} END{if(s>0) printf "%.2f", g/s; else print 0}' "$STATS")
GAPS=$(awk 'NR>1{g+=$4} END{print g+0}' "$STATS")
N50=$(awk 'NR>1{print $2}' "$STATS" | sort -nr | awk '{a[NR]=$1; s+=$1} END{half=s/2; c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=half){print a[i]; break}}}')
MAXLEN=$(awk 'NR>1{if($2>m)m=$2} END{print m+0}' "$STATS")

{
  echo -e "metric\tvalue"
  echo -e "assembly\t$ASM"
  echo -e "contigs\t$NCONTIG"
  echo -e "total_bp\t$TOTAL"
  echo -e "n50_bp\t$N50"
  echo -e "max_contig_bp\t$MAXLEN"
  echo -e "gc_pct\t$GC"
  echo -e "gap_bases\t$GAPS"
} > "$OUT/qc_summary.tsv"

echo "[t2t] contigs=$NCONTIG total=$TOTAL N50=$N50 GC=${GC}% gaps=$GAPS"

# ---- external tool plan ----
PLAN="$OUT/tool_plan.txt"
{
  echo "# T2T QC external tool plan (run with --run-tools)"
  echo "BUSCO: busco -i $ASM -l $LINEAGE -m genome -o busco -c $THREADS --out_path $OUT"
  if [[ -n "$READS" ]]; then
    echo "Merqury: merqury.sh <reads.meryl> $ASM $OUT/merqury"
  else
    echo "Merqury: skipped (need --reads and a meryl database)"
  fi
  if [[ -n "$REPEAT_LIB" ]]; then
    echo "RepeatMasker: RepeatMasker -pa $THREADS -lib $REPEAT_LIB -dir $OUT/repeatmasker $ASM"
  else
    echo "RepeatMasker: skipped (need --repeat-lib)"
  fi
  if [[ -n "$ANNO" && -n "$REF" ]]; then
    echo "Liftoff: liftoff -g $ANNO -o $OUT/liftoff.gff3 -u $OUT/unmapped.txt $ASM $REF"
  else
    echo "Liftoff: skipped (need --annotation and --reference)"
  fi
} > "$PLAN"
cat "$PLAN"

if [[ $RUN_TOOLS -eq 1 ]]; then
  command -v busco >/dev/null 2>&1 || { echo "[t2t] BUSCO not found; skip" >&2; }
  if command -v busco >/dev/null 2>&1; then
    busco -i "$ASM" -l "$LINEAGE" -m genome -o busco -c "$THREADS" --out_path "$OUT" \
      > "$OUT/busco.log" 2>&1 || echo "[t2t] BUSCO failed; see busco.log" >&2
  fi
  if [[ -n "$REPEAT_LIB" ]] && command -v RepeatMasker >/dev/null 2>&1; then
    mkdir -p "$OUT/repeatmasker"
    RepeatMasker -pa "$THREADS" -lib "$REPEAT_LIB" -dir "$OUT/repeatmasker" "$ASM" \
      > "$OUT/repeatmasker.log" 2>&1 || echo "[t2t] RepeatMasker failed" >&2
  fi
  if [[ -n "$ANNO" && -n "$REF" ]] && command -v liftoff >/dev/null 2>&1; then
    liftoff -g "$ANNO" -o "$OUT/liftoff.gff3" -u "$OUT/unmapped.txt" "$ASM" "$REF" \
      > "$OUT/liftoff.log" 2>&1 || echo "[t2t] Liftoff failed" >&2
  fi
  echo "[t2t] external tools finished (check logs in $OUT)"
fi
