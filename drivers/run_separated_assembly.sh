#!/usr/bin/env bash
set -uo pipefail

LABELS=""; READS=""; TYPE=""; WORK=""; HIC1=""; HIC2=""; REF=""; CHROM_SUBG=""
while [ $# -gt 0 ]; do
  case "$1" in
    --labels) LABELS=$2; shift 2;;
    --reads) READS=$2; shift 2;;
    --type) TYPE=$2; shift 2;;
    --work) WORK=$2; shift 2;;
    --hic1) HIC1=$2; shift 2;;
    --hic2) HIC2=$2; shift 2;;
    --ref) REF=$2; shift 2;;
    --chrom-subg) CHROM_SUBG=$2; shift 2;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done
for v in LABELS READS TYPE WORK; do
  [ -n "${!v}" ] || { echo "missing --$(echo "$v" | tr 'A-Z' 'a-z')"; exit 1; }
done
case "$TYPE" in
  ont)  FLYE_PRESET="--nano-hq";     MMX="map-ont";;
  hifi) FLYE_PRESET="--pacbio-hifi"; MMX="map-hifi";;
  *) echo "--type must be ont or hifi"; exit 1;;
esac

THREADS=${THREADS:-48}; PY=python3
FLYE_BIN=${FLYE_BIN:-/path/to/flye/bin}
BWA=bwa; SAMTOOLS=samtools; SEQKIT=seqkit; MINIMAP2=minimap2; YAHS=yahs
SEP=$WORK/sep_asm; mkdir -p "$SEP"
log(){ echo "### $* | $(date) ###"; }

SUBG=$(cut -f2 "$LABELS" | tail -n +2 | grep -viE 'ambiguous|unassigned|^label$' | sort -u)
log "subgenomes: $(echo $SUBG | tr '\n' ' ')"
printf "subgenome\tsize_Mb\tcontigs\tN50\tbackmap_pct\tpurity_pct\n" > "$SEP/summary.tsv"

for sg in $SUBG; do
  D=$SEP/$sg; mkdir -p "$D"; RD=$D/$sg.reads.fq.gz; ASM=$D/flye/assembly.fasta

  if [ ! -s "$RD" ]; then
    awk -v s="$sg" -F'\t' 'NR>1 && $2==s{print $1}' "$LABELS" > "$D/$sg.ids"
    "$SEQKIT" grep -f "$D/$sg.ids" "$READS" -o "$RD"
  fi

  if [ ! -s "$ASM" ]; then
    export PATH="$FLYE_BIN:$PATH"
    log "[$sg] Flye $FLYE_PRESET"
    if [ -d "$D/flye/00-assembly" ]; then
      flye $FLYE_PRESET "$RD" --keep-haplotypes --resume -o "$D/flye" -t "$THREADS"
    else
      flye $FLYE_PRESET "$RD" --keep-haplotypes -o "$D/flye" -t "$THREADS"
    fi
    [ -s "$ASM" ] || { echo "!! [$sg] Flye did not finish"; continue; }
  fi
  [ -s "$ASM.fai" ] || "$SAMTOOLS" faidx "$ASM"

  FINAL="$ASM"
  if [ -n "$HIC1" ] && [ -n "$HIC2" ]; then
    SC=$D/yahs/${sg}_scaffolds_final.fa
    if [ ! -s "$SC" ]; then
      mkdir -p "$D/yahs"; [ -s "$ASM.bwt" ] || "$BWA" index "$ASM"
      log "[$sg] Hi-C align + YaHS"
      if [ ! -s "$D/yahs/hic.bam" ]; then
        "$BWA" mem -5SP -t "$THREADS" "$ASM" "$HIC1" "$HIC2" \
          | "$SAMTOOLS" view -@ 8 -bh -F 0x904 -q 1 - \
          | "$SAMTOOLS" sort -@ 8 -m 2G -T "$D/yahs/st" -n -o "$D/yahs/hic.bam"
      fi
      ( cd "$D/yahs" && "$YAHS" "$ASM" hic.bam -o "$sg" )
    fi
    [ -s "$SC" ] && FINAL="$SC"
  fi
  [ -s "$FINAL.fai" ] || "$SAMTOOLS" faidx "$FINAL"

  read TOTBP NCTG N50 < <("$SEQKIT" stats -a -T "$FINAL" | awk 'NR==2{print $5, $4, $13}')
  SIZE_MB=$(awk -v b="$TOTBP" 'BEGIN{printf "%.0f", b/1e6}')
  BM=$("$MINIMAP2" -ax "$MMX" -t "$THREADS" "$FINAL" "$RD" 2>/dev/null \
        | "$SAMTOOLS" view -bh -F 0x900 - | "$SAMTOOLS" flagstat - \
        | awk '/ mapped \(/{gsub(/[()%]/,"",$5); print $5; exit}')
  PUR=NA
  if [ -n "$REF" ] && [ -n "$CHROM_SUBG" ]; then
    "$MINIMAP2" -cx asm10 -t "$THREADS" "$REF" "$FINAL" 2>/dev/null > "$D/asm_to_ref.paf"
    PUR=$("$PY" - "$D/asm_to_ref.paf" "$CHROM_SUBG" <<'PY'
import sys, collections
paf, csub = sys.argv[1], sys.argv[2]
sub = {}
for ln in open(csub):
    p = ln.split()
    if len(p) >= 2: sub[p[0]] = p[1]
best = {}
for ln in open(paf):
    f = ln.split('\t')
    q, nmatch, tname = f[0], int(f[9]), f[5]
    s = sub.get(tname)
    if s is None: continue
    if q not in best or nmatch > best[q][0]: best[q] = (nmatch, s)
bp = collections.Counter()
for q, (m, s) in best.items(): bp[s] += m
tot = sum(bp.values())
print(f"{100*max(bp.values())/tot:.1f}" if tot else "NA")
PY
)
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sg" "$SIZE_MB" "$NCTG" "$N50" "${BM:-NA}" "$PUR" | tee -a "$SEP/summary.tsv"
done
echo; log "DONE -> $SEP/summary.tsv"; column -t "$SEP/summary.tsv"
