#!/usr/bin/env bash
set -uo pipefail
export PYTHONHASHSEED=0

NSG=3 ; THREADS=48 ; K_PAIR=33 ; BETA=0.60
PY=python3
PKG=$POLYSPLIT
PIPE=$PKG/pipeline
DATA=$DATA/camelina/Hexaploid_data
LONG_READS=$DATA/TMP24026.ONT.fastq.gz
HIC_R1=$DATA/TMP24026.hic.R1.fastq.gz
HIC_R2=$DATA/TMP24026.hic.R2.fastq.gz
REF=$DATA/Cmicrocarpa.CN119205.fa.gz
CHROM_SUBG=$PKG/chrom_subg.hexaploid.tsv
WORK=$DATA/polysplit_run_nanohq ; mkdir -p "$WORK"
ASM=$WORK/flye_out/assembly.fasta

FLYE_BIN=${FLYE_BIN:-/path/to/flye/bin}
BWA=bwa
SAMTOOLS=samtools
SEQKIT=seqkit
MINIMAP2=minimap2
GSUFSORT=gsufsort

log(){ echo "### $* | $(date) ###"; }

if [ ! -s "$ASM" ]; then
  export PATH="$FLYE_BIN:$PATH"
  flye --version >/dev/null 2>&1 || true
  flye-samtools --version >/dev/null 2>&1 || true
  flye-modules --help >/dev/null 2>&1 || true
  if [ -d "$WORK/flye_out/00-assembly" ]; then
    log "Stage 1: Flye --nano-hq --resume (partial flye_out found)"
    flye --nano-hq "$LONG_READS" --keep-haplotypes --resume -o "$WORK/flye_out" -t "$THREADS"
  else
    log "Stage 1: Flye --nano-hq (correct preset for R10.4.1)"
    flye --nano-hq "$LONG_READS" --keep-haplotypes -t "$THREADS" -o "$WORK/flye_out"
  fi
  [ -s "$ASM" ] || { echo "!! Flye did not finish. Re-run to resume, on NUGGET."; exit 1; }
fi
[ -s "$ASM.fai" ] || "$SAMTOOLS" faidx "$ASM"
ASM_MB=$(awk '!/^>/{n+=length($0)}END{printf "%.0f", n/1e6}' "$ASM")
log "assembly: $(grep -c '^>' "$ASM") contigs, ${ASM_MB} Mb  (--nano-raw gave 241 Mb collapsed; HiFi gave 591 Mb)"

if [ -s "$WORK/total_reads.txt" ]; then TOTAL_READS=$(cat "$WORK/total_reads.txt"); else
  TOTAL_READS=$("$SEQKIT" stats -T -j 8 "$LONG_READS" | awk 'NR==2{print $4}'); echo "$TOTAL_READS" > "$WORK/total_reads.txt"
fi

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
    if prev and prev[0]==q and prev[1]!=c and c!='*' and prev[1]!='*':
        a,b = sorted((prev[1], c)); contacts[(a,b)] += 1
    prev = (q, c)
pickle.dump((dict(contacts), length), open(out, 'wb'))
print(f"[contacts] {len(contacts):,} pairs / {len(length):,} contigs", file=sys.stderr)
PY
  "$BWA" mem -5SP -t "$THREADS" "$ASM" "$HIC_R1" "$HIC_R2" \
    | "$SAMTOOLS" view -@ 8 -bh -F 0x904 -q 1 - \
    | "$SAMTOOLS" sort -@ 8 -m 2G -T "$WORK/sorttmp" -n -o "$WORK/hic.nsort.bam"
  "$SAMTOOLS" view "$WORK/hic.nsort.bam" | "$PY" "$WORK/build_contacts.py" "$ASM.fai" "$WORK/contacts.pkl"
fi

if [ ! -s "$WORK/homoeolog_edges.tsv" ]; then
  log "Stage 3: gsufsort GSA + LCP + homoeolog edges"
  "$GSUFSORT" "$ASM" --fasta --gsa 4 8 --output "$WORK/contigs_k33"
  "$GSUFSORT" "$ASM" --fasta --lcp 1  --output "$WORK/contigs_k33"
  "$PY" "$PIPE/homoeolog_graph_from_lcp.py" --gsa "$WORK/contigs_k33.4.8.gsa" \
        --lcp "$WORK/contigs_k33.1.lcp" --fa "$ASM" --k "$K_PAIR" \
        --min-copy 2 --max-copy 6 --out "$WORK/homoeolog_edges.tsv"
fi

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

[ -s "$WORK/reads_to_contigs.paf" ] || { log "Stage 7a: reads->contigs (map-ont)"; \
  "$MINIMAP2" -x map-ont -t "$THREADS" "$ASM" "$LONG_READS" > "$WORK/reads_to_contigs.paf"; }
[ -s "$WORK/reads_to_ref.paf" ] || { log "Stage 7b: reads->ref (map-ont, truth)"; \
  "$MINIMAP2" -cx map-ont -t "$THREADS" "$REF" "$LONG_READS" > "$WORK/reads_to_ref.paf"; }
[ -s "$WORK/read_subg.tsv" ] || { log "Stage 8: propagate"; \
  "$PY" "$PIPE/propagate_to_reads.py" --paf "$WORK/reads_to_contigs.paf" \
        --contig-labels "$WORK/all_contig_labels_repaired.tsv" --min-conf "$BETA" --weight ident \
        --out "$WORK/read_subg.tsv"; }

log "Stage 9a: CONTIG accuracy (K=3, --nano-hq)"
[ -s "$WORK/contigs_to_ref.paf" ] || \
  "$MINIMAP2" -cx asm10 -t "$THREADS" "$REF" "$ASM" > "$WORK/contigs_to_ref.paf"
"$PY" "$PIPE/eval_contig_labels.py" "$CHROM_SUBG" "$WORK/contigs_to_ref.paf" "$ASM.fai" \
      "$WORK/all_contig_labels_repaired.tsv" "$WORK/wg_purity.per_contig.tsv"
log "Stage 9b: READ-level accuracy"
"$PY" "$PIPE/allread_eval.py" --labels "$WORK/read_subg.tsv" --ref-paf "$WORK/reads_to_ref.paf" \
      --total-reads "$TOTAL_READS" --chrom-subg "$CHROM_SUBG"
log "DONE (hexaploid --nano-hq K=3): $WORK"
