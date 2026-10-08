#!/usr/bin/env bash
set -uo pipefail
export PYTHONHASHSEED=0

READS=""; HIC1=""; HIC2=""; TYPE=""; NSG=""; WORK=""; REF=""; CHROM_SUBG=""; FP=""; ASSEMBLER=flye
while [ $# -gt 0 ]; do
  case "$1" in
    --reads) READS=$2; shift 2;;
    --hic1) HIC1=$2; shift 2;;
    --hic2) HIC2=$2; shift 2;;
    --type) TYPE=$2; shift 2;;
    --nsg) NSG=$2; shift 2;;
    --work) WORK=$2; shift 2;;
    --ref) REF=$2; shift 2;;
    --chrom-subg) CHROM_SUBG=$2; shift 2;;
    --flye-preset) FP=$2; shift 2;;
    --assembler) ASSEMBLER=$2; shift 2;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done
for v in READS HIC1 HIC2 TYPE NSG WORK; do
  [ -n "${!v}" ] || { echo "missing --$(echo "$v" | tr 'A-Z' 'a-z')"; exit 1; }
done
case "$TYPE" in
  ont)  FLYE_PRESET="--nano-hq";     MMX="map-ont";;
  hifi) FLYE_PRESET="--pacbio-hifi"; MMX="map-hifi";;
  *) echo "--type must be ont or hifi"; exit 1;;
esac
if [ -n "$FP" ]; then case "$FP" in --*) FLYE_PRESET=$FP;; *) FLYE_PRESET="--$FP";; esac; fi

THREADS=${THREADS:-48}; K_PAIR=33; BETA=0.60
PY=python3; PIPE=$POLYSPLIT/pipeline
FLYE_BIN=${FLYE_BIN:-/path/to/flye/bin}
BWA=bwa; SAMTOOLS=samtools; MINIMAP2=minimap2; GSUFSORT=gsufsort; SEQKIT=seqkit
mkdir -p "$WORK"; ASM=$WORK/flye_out/assembly.fasta
if [ "$ASSEMBLER" = hifiasm ]; then
  [ "$TYPE" = hifi ] || { echo "--assembler hifiasm needs --type hifi"; exit 1; }
  ASM=$WORK/hifiasm/asm.bp.p_ctg.fa
fi
log(){ echo "### $* | $(date) ###"; }

# Stage 1: de-novo assembly
if [ "$ASSEMBLER" = hifiasm ] && [ ! -s "$ASM" ]; then
  log "Stage 1: hifiasm -l0"
  mkdir -p "$WORK/hifiasm"
  hifiasm -o "$WORK/hifiasm/asm" -t "$THREADS" -l0 "$READS" 2> "$WORK/hifiasm/hifiasm.log"
  awk '/^S/{print ">"$2"\n"$3}' "$WORK/hifiasm/asm.bp.p_ctg.gfa" > "$ASM"
  [ -s "$ASM" ] || { echo "!! hifiasm did not finish; check $WORK/hifiasm/hifiasm.log"; exit 1; }
fi
if [ ! -s "$ASM" ]; then
  export PATH="$FLYE_BIN:$PATH"
  log "Stage 1: Flye $FLYE_PRESET"
  if [ -d "$WORK/flye_out/00-assembly" ]; then
    flye $FLYE_PRESET "$READS" --keep-haplotypes --resume -o "$WORK/flye_out" -t "$THREADS"
  else
    flye $FLYE_PRESET "$READS" --keep-haplotypes -o "$WORK/flye_out" -t "$THREADS"
  fi
  [ -s "$ASM" ] || { echo "!! Flye did not finish (assembly.fasta missing); check $WORK/flye_out/flye.log"; exit 1; }
fi
[ -s "$ASM.fai" ] || "$SAMTOOLS" faidx "$ASM"

# Stage 2: Hi-C -> contig contact graph
if [ ! -s "$WORK/contacts.pkl" ]; then
  log "Stage 2: Hi-C contacts"
  [ -s "$ASM.bwt" ] || "$BWA" index "$ASM"
  cat > "$WORK/build_contacts.py" <<'PY'
import sys, pickle, collections
fai, out = sys.argv[1], sys.argv[2]
length = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(fai)}
contacts = collections.defaultdict(int); prev = None
for ln in sys.stdin:
    p = ln.split('\t'); q, c = p[0], p[2]
    if prev and prev[0] == q and prev[1] != c and c != '*' and prev[1] != '*':
        a, b = sorted((prev[1], c)); contacts[(a, b)] += 1
    prev = (q, c)
pickle.dump((dict(contacts), length), open(out, 'wb'))
PY
  "$BWA" mem -5SP -t "$THREADS" "$ASM" "$HIC1" "$HIC2" \
    | "$SAMTOOLS" view -@ 8 -bh -F 0x904 -q 1 - \
    | "$SAMTOOLS" sort -@ 8 -m 2G -T "$WORK/sorttmp" -n -o "$WORK/hic.nsort.bam"
  "$SAMTOOLS" view "$WORK/hic.nsort.bam" | "$PY" "$WORK/build_contacts.py" "$ASM.fai" "$WORK/contacts.pkl"
fi

# Stage 3: shared-33mer homoeolog edges (gsufsort GSA + LCP)
if [ ! -s "$WORK/homoeolog_edges.tsv" ]; then
  log "Stage 3: gsufsort GSA+LCP -> homoeolog edges"
  "$GSUFSORT" "$ASM" --fasta --gsa 4 8 --output "$WORK/contigs_k33"
  "$GSUFSORT" "$ASM" --fasta --lcp 1  --output "$WORK/contigs_k33"
  "$PY" "$PIPE/homoeolog_graph_from_lcp.py" --gsa "$WORK/contigs_k33.4.8.gsa" \
        --lcp "$WORK/contigs_k33.1.lcp" --fa "$ASM" --k "$K_PAIR" \
        --min-copy 2 --max-copy 6 --out "$WORK/homoeolog_edges.tsv"
fi

# Stage 4-6: Hi-C blocks -> de-chimerize -> label (K subgenomes) -> recover -> repair
if [ ! -s "$WORK/all_contig_labels_repaired.tsv" ]; then
  log "Stage 4-6: blocks/label/repair (NSG=$NSG)"
  export POLYSPLIT_NSG="$NSG" POLYSPLIT_FASTA="$ASM" \
         POLYSPLIT_CONTACTS="$WORK/contacts.pkl" POLYSPLIT_EDGES="$WORK/homoeolog_edges.tsv" \
         POLYSPLIT_TRUTH=/none
  ( cd "$WORK"
    "$PY" "$PIPE/dechimerize_structural.py"
    POLYSPLIT_OUT="$WORK/all_contig_labels_v2.tsv" \
      "$PY" "$PIPE/label_small_contigs_v2.py" "$WORK/decloud_structural_contig_labels.tsv"
    POLYSPLIT_LABELS_IN="$WORK/all_contig_labels_v2.tsv" POLYSPLIT_OUT="$WORK/all_contig_labels_repaired.tsv" \
      "$PY" "$PIPE/homoeolog_repair.py" )
  cut -f2 "$WORK/all_contig_labels_repaired.tsv" | tail -n +2 | sort | uniq -c
fi

# Stage 7: reads -> contigs
[ -s "$WORK/reads_to_contigs.paf" ] || { log "Stage 7: reads->contigs"; \
  "$MINIMAP2" -x "$MMX" -t "$THREADS" "$ASM" "$READS" > "$WORK/reads_to_contigs.paf"; }

# Stage 8: propagate contig labels -> reads
[ -s "$WORK/read_subg.tsv" ] || { log "Stage 8: propagate labels to reads"; \
  "$PY" "$PIPE/propagate_to_reads.py" --paf "$WORK/reads_to_contigs.paf" \
        --contig-labels "$WORK/all_contig_labels_repaired.tsv" --min-conf "$BETA" --weight ident \
        --out "$WORK/read_subg.tsv"; }
log "per-read subgenome labels -> $WORK/read_subg.tsv"

# Stage 9: evaluation (optional; only when a reference and chrom->subgenome map are given)
if [ -n "$REF" ] && [ -n "$CHROM_SUBG" ]; then
  log "Stage 9: evaluation vs reference"
  [ -s "$WORK/contigs_to_ref.paf" ] || "$MINIMAP2" -cx asm10 -t "$THREADS" "$REF" "$ASM" > "$WORK/contigs_to_ref.paf"
  "$PY" "$PIPE/eval_contig_labels.py" "$CHROM_SUBG" "$WORK/contigs_to_ref.paf" "$ASM.fai" \
        "$WORK/all_contig_labels_repaired.tsv" "$WORK/wg_purity.per_contig.tsv"
  [ -s "$WORK/reads_to_ref.paf" ] || "$MINIMAP2" -cx "$MMX" -t "$THREADS" "$REF" "$READS" > "$WORK/reads_to_ref.paf"
  TOTAL_READS=$("$SEQKIT" stats -T -j 8 "$READS" | awk 'NR==2{print $4}')
  "$PY" "$PIPE/allread_eval.py" --labels "$WORK/read_subg.tsv" --ref-paf "$WORK/reads_to_ref.paf" \
        --total-reads "$TOTAL_READS" --chrom-subg "$CHROM_SUBG"
fi
