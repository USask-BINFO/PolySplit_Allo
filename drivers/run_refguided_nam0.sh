#!/usr/bin/env bash
set -uo pipefail
HERE=$POLYSPLIT
A=$DATA          # napus ref-guided scripts live here
N=$A/NAM0
SK=seqkit
SAMTOOLS=samtools
PY=python3
GSUFSORT=gsufsort-64  # 64-bit pos: plain gsufsort leaves the high 4 bytes of pos uninitialized

WORK=$N/polysplit_run
RGW=$WORK/refguided
CHROM=$N/n99_chr.fa                                      # N99 chromosome-level reference (N1..N19)
CHROM_SUBG=$N/chrom_subg.nam0.tsv                        # N<i> <TAB> A|C
READS=$N/ont/n99_all.fastq.gz                            # 96 Gb ONT reads
R2R=$WORK/reads_to_ref.paf                               # reads -> N99 ref (truth, already built)
GSA=$RGW/chrom64_k33.4.8.gsa
LCP=$RGW/chrom64_k33.1.lcp
log(){ echo "### $* | $(date) ###"; }
mkdir -p "$RGW"

for f in "$CHROM" "$CHROM_SUBG" "$READS" "$R2R"; do
  [ -s "$f" ] || { echo "!! missing input: $f"; exit 1; }
done

# ---------------- Stage 0: subgenome fastas + fai + gsufsort index ----------------
log "Stage 0: subgenome fastas + reference index"
[ -s "$CHROM.fai" ] || "$SAMTOOLS" faidx "$CHROM"
if [ ! -s "$RGW/A.fa" ]; then
  "$SAMTOOLS" faidx "$CHROM" $(awk '$2=="A"{print $1}' "$CHROM_SUBG") > "$RGW/A.fa"
fi
if [ ! -s "$RGW/C.fa" ]; then
  "$SAMTOOLS" faidx "$CHROM" $(awk '$2=="C"{print $1}' "$CHROM_SUBG") > "$RGW/C.fa"
fi
echo "A chroms: $(grep -c '^>' "$RGW/A.fa")   C chroms: $(grep -c '^>' "$RGW/C.fa")"
# gsufsort-64: GSA and LCP capped at k in one call
if [ ! -s "$GSA" ] || [ ! -s "$LCP" ]; then
  "$GSUFSORT" "$CHROM" --fasta --gsa 4 8 --lcp 1 --trlcp 33 --output "$RGW/chrom64_k33" --time
fi
[ -s "$GSA" ] && [ -s "$LCP" ] || { echo "!! GSA/LCP build failed"; exit 2; }

# total read count (only affects the T3 line; our headline is correct/n_tru). Cached.
TOTAL_CACHE=$N/n99_read_count.txt
if [ ! -s "$TOTAL_CACHE" ]; then
  log "counting reads (one pass; cached -> $TOTAL_CACHE)"
  "$SK" stats -T "$READS" | awk 'NR==2{print $4}' > "$TOTAL_CACHE"
fi
TOTAL=$(cat "$TOTAL_CACHE"); echo "total reads = $TOTAL"

# ---------------- Stage 1: Select_signatures ----------------
log "Stage 1: Select_signatures (gsufsort-64 GSA, k=33, max-copy 3)"
if [ ! -s "$RGW/sig_k33.A.bed" ] || [ ! -s "$RGW/sig_k33.C.bed" ]; then
  "$PY" "$A/Select_signatures.py" --k 33 --gsa "$GSA" --lcp "$LCP" \
    --ac-fa "$CHROM" --a-fa "$RGW/A.fa" --c-fa "$RGW/C.fa" --fai "$CHROM.fai" \
    --max-copy 3 --target 200000000000 --seed 1 --out-prefix "$RGW/sig_k33"
fi
wc -l "$RGW/sig_k33.A.bed" "$RGW/sig_k33.C.bed"
[ -s "$RGW/sig_k33.A.bed" ] && [ -s "$RGW/sig_k33.C.bed" ] || { echo "!! BEDs empty -- abort"; exit 3; }

# ---------------- Stage 2: extract signature k-mer sequences ----------------
log "Stage 2: seqkit subseq --bed"
[ -s "$RGW/A_sig_raw.fa" ] || "$SK" subseq --bed "$RGW/sig_k33.A.bed" "$CHROM" > "$RGW/A_sig_raw.fa" 2> "$RGW/subseq_A.err"
[ -s "$RGW/C_sig_raw.fa" ] || "$SK" subseq --bed "$RGW/sig_k33.C.bed" "$CHROM" > "$RGW/C_sig_raw.fa" 2> "$RGW/subseq_C.err"
echo "raw sig seqs: A=$(grep -c '^>' "$RGW/A_sig_raw.fa")  C=$(grep -c '^>' "$RGW/C_sig_raw.fa")"

# ---------------- Stage 3: clean signature sets ----------------
log "Stage 3: make_clean_signature_sets"
if [ ! -s "$RGW/A_sig_k33.txt" ] || [ ! -s "$RGW/C_sig_k33.txt" ]; then
  "$PY" "$A/make_clean_signature_sets.py" "$RGW/A_sig_raw.fa" "$RGW/C_sig_raw.fa" \
        "$RGW/A_sig_k33.txt" "$RGW/C_sig_k33.txt"
fi
n1=$(wc -l < "$RGW/A_sig_k33.txt"); n2=$(wc -l < "$RGW/C_sig_k33.txt")
echo "clean sig 33-mers: A=$n1  C=$n2"
{ [ "$n1" -gt 1000 ] && [ "$n2" -gt 1000 ]; } || { echo "!! sig sets too small -- abort"; exit 4; }

# ---------------- Stage 4: classify the ONT reads ----------------
log "Stage 4: classify_reads_by_kmers (k=33, min-hits 3, ratio 3.0)"
if [ ! -s "$RGW/read_labels_refguided.tsv" ]; then
  "$PY" "$A/classify_reads_by_kmers.py" "$RGW/A_sig_k33.txt" "$RGW/C_sig_k33.txt" \
        "$READS" 33 "$RGW/read_labels_refguided.tsv" --min-hits 3 --ratio 3.0 --min-len 2000 --stride 1
fi

# ---------------- Stage 5: evaluate vs A/C chromosome truth ----------------
log "Stage 5: EVALUATE (upper bound) vs A/C truth"
cut -f1,4 "$RGW/read_labels_refguided.tsv" > "$RGW/read_labels_refguided.2col.tsv"
"$PY" "$HERE/pipeline/allread_eval.py" --labels "$RGW/read_labels_refguided.2col.tsv" \
      --ref-paf "$R2R" --total-reads "$TOTAL" --chrom-subg "$CHROM_SUBG"
echo
echo "TABLE read acc = correct / chromosome-truth  (NOT the printed lenient T1;"
echo "ambiguous + unassigned reads that have truth count as errors -- same convention as the other rows)."
log "DONE refguided NAM0 (upper bound)"
