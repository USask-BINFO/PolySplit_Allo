#!/usr/bin/env bash
set -uo pipefail
export PYTHONHASHSEED=0

# ---------------- config (edit paths only if running off nugget) ----------------
NSG=2 ; THREADS=48 ; K_PAIR=33 ; BETA=0.60
PY=python3
PKG=$POLYSPLIT
PIPE=$PKG/pipeline
DATA=$DATA/napus
LONG_READS=$DATA/ont/n99_all.fastq.gz
HIC_R1=$DATA/hic/S00E51E_r1.fastq.gz
HIC_R2=$DATA/hic/S00E51E_r2.fastq.gz
TOTAL_READS=5725172
REF=$DATA/assembly/Bnapus_N99_hifi.genome.fasta      # N99 HiFi assembly = eval truth
CHROM_SUBG=$DATA/chrom_subg.nam0.tsv                 # N1-N10 -> A, N11-N19 -> C
WORK=$DATA/polysplit_run ; mkdir -p "$WORK"
ASM=$WORK/flye_out/assembly.fasta

# tools (nugget-validated)
FLYE_BIN=${FLYE_BIN:-/path/to/flye/bin}
BWA=bwa
SAMTOOLS=samtools
MINIMAP2=minimap2
GSUFSORT=gsufsort

log(){ echo "### $* | $(date) ###"; }

# ---------------- Stage 1: de-novo assembly (Flye) ----------------
if [ ! -s "$ASM" ]; then
  export PATH="$FLYE_BIN:$PATH"
  # pre-warm the flye helper binaries (first exec over NFS can fail)
  flye --version          >/dev/null 2>&1 || true
  flye-samtools --version >/dev/null 2>&1 || true
  flye-modules --help     >/dev/null 2>&1 || true
  if [ -d "$WORK/flye_out/00-assembly" ]; then
    log "Stage 1: Flye --resume (partial flye_out found)"
    flye --nano-raw "$LONG_READS" --keep-haplotypes --resume -o "$WORK/flye_out" -t "$THREADS"
  else
    log "Stage 1: Flye assembly"
    flye --nano-raw "$LONG_READS" --keep-haplotypes -t "$THREADS" -o "$WORK/flye_out"
  fi
  [ -s "$ASM" ] || { echo "!! Flye did not finish (assembly.fasta missing). Run on NUGGET, not lode/trove; check flye_out/flye.log"; exit 1; }
fi
[ -s "$ASM.fai" ] || "$SAMTOOLS" faidx "$ASM"
log "contigs: $(grep -c '^>' "$ASM") "

# ---------------- Stage 2: Hi-C -> contig contact graph ----------------
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

# ---------------- Stage 3: shared-33mer homoeolog edges (gsufsort + LCP walk) ----------------
if [ ! -s "$WORK/homoeolog_edges.tsv" ]; then
  log "Stage 3: gsufsort GSA + LCP (two runs; combining --gsa+--lcp drops the LCP)"
  "$GSUFSORT" "$ASM" --fasta --gsa 4 8 --output "$WORK/contigs_k33"
  "$GSUFSORT" "$ASM" --fasta --lcp 1  --output "$WORK/contigs_k33"
  ls -la "$WORK/contigs_k33.4.8.gsa" "$WORK/contigs_k33.1.lcp"
  "$PY" "$PIPE/homoeolog_graph_from_lcp.py" --gsa "$WORK/contigs_k33.4.8.gsa" \
        --lcp "$WORK/contigs_k33.1.lcp" --fa "$ASM" --k "$K_PAIR" \
        --min-copy 2 --max-copy 6 --out "$WORK/homoeolog_edges.tsv"
fi

# ---------------- Stage 4-6: blocks -> de-chimerize -> label (S1/S2) -> recover -> repair ----------------
if [ ! -s "$WORK/all_contig_labels_repaired.tsv" ]; then
  log "Stage 4-6: blocks/label/repair"
  export POLYSPLIT_NSG="$NSG" POLYSPLIT_FASTA="$ASM" \
         POLYSPLIT_CONTACTS="$WORK/contacts.pkl" POLYSPLIT_EDGES="$WORK/homoeolog_edges.tsv" \
         POLYSPLIT_TRUTH=/none                                   # per-contig truth optional; eval below
  ( cd "$WORK"
    "$PY" "$PIPE/dechimerize_structural.py"
    POLYSPLIT_OUT="$WORK/all_contig_labels_v2.tsv" \
      "$PY" "$PIPE/label_small_contigs_v2.py" "$WORK/decloud_structural_contig_labels.tsv"
    POLYSPLIT_LABELS_IN="$WORK/all_contig_labels_v2.tsv" POLYSPLIT_OUT="$WORK/all_contig_labels_repaired.tsv" \
      "$PY" "$PIPE/homoeolog_repair.py" )
  cut -f2 "$WORK/all_contig_labels_repaired.tsv" | tail -n +2 | sort | uniq -c
fi

# ---------------- Stage 7: reads -> contigs PAF + reads -> ref PAF ----------------
[ -s "$WORK/reads_to_contigs.paf" ] || { log "Stage 7a: reads->contigs"; \
  "$MINIMAP2" -x map-ont -t "$THREADS" "$ASM" "$LONG_READS" > "$WORK/reads_to_contigs.paf"; }
[ -s "$WORK/reads_to_ref.paf" ] || { log "Stage 7b: reads->ref (truth)"; \
  "$MINIMAP2" -cx map-ont -t "$THREADS" "$REF" "$LONG_READS" > "$WORK/reads_to_ref.paf"; }

# ---------------- Stage 8: propagate contig labels -> reads ----------------
[ -s "$WORK/read_subg.tsv" ] || { log "Stage 8: propagate"; \
  "$PY" "$PIPE/propagate_to_reads.py" --paf "$WORK/reads_to_contigs.paf" \
        --contig-labels "$WORK/all_contig_labels_repaired.tsv" --min-conf "$BETA" --weight ident \
        --out "$WORK/read_subg.tsv"; }

# ---------------- Stage 9: evaluation (contig + read level) ----------------
log "Stage 9a: per-contig truth + CONTIG accuracy"
[ -s "$WORK/contigs_to_ref.paf" ] || \
  "$MINIMAP2" -cx asm10 -t "$THREADS" "$REF" "$ASM" > "$WORK/contigs_to_ref.paf"
"$PY" "$PIPE/eval_contig_labels.py" "$CHROM_SUBG" "$WORK/contigs_to_ref.paf" "$ASM.fai" \
      "$WORK/all_contig_labels_repaired.tsv" "$WORK/wg_purity.per_contig.tsv"

log "Stage 9b: READ-level accuracy / precision / recall / F1 (vs A/C truth)"
"$PY" "$PIPE/allread_eval.py" --labels "$WORK/read_subg.tsv" --ref-paf "$WORK/reads_to_ref.paf" \
      --total-reads "$TOTAL_READS" --chrom-subg "$CHROM_SUBG"

log "DONE: $WORK/read_subg.tsv (per-read A/C labels)"
