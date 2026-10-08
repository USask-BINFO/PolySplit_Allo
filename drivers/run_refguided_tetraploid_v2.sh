#!/usr/bin/env bash
set -uo pipefail
HERE=$POLYSPLIT
A=$DATA          # napus ref-guided scripts
SK=seqkit
PY=python3
RGW=$DATA/camelina/Tetraploid_data/polysplit_run/refguided
WORK=$DATA/camelina/Tetraploid_data/polysplit_run
CHROM_SUBG=$HERE/chrom_subg.tetraploid.tsv
READS=$DATA/camelina/Tetraploid_data/CN119243.ONT.fastq.gz
TOTAL=1647060
log(){ echo "### $* | $(date) ###"; }

log "Select_signatures (gsufsort-64 GSA, max-copy 3)"
"$PY" "$A/Select_signatures.py" --k 33 --gsa "$RGW/chrom64_k33.4.8.gsa" --lcp "$RGW/chrom64_k33.1.lcp" \
  --ac-fa "$RGW/chrom.fa" --a-fa "$RGW/S1.fa" --c-fa "$RGW/S2.fa" --fai "$RGW/chrom.fa.fai" \
  --max-copy 3 --target 200000000000 --seed 1 --out-prefix "$RGW/sig_k33"
wc -l "$RGW/sig_k33.A.bed" "$RGW/sig_k33.C.bed"
[ -s "$RGW/sig_k33.A.bed" ] && [ -s "$RGW/sig_k33.C.bed" ] || { echo "!! BEDs empty -- abort"; exit 2; }

log "extract signature k-mer sequences (seqkit subseq --bed)"
"$SK" subseq --bed "$RGW/sig_k33.A.bed" "$RGW/chrom.fa" > "$RGW/S1_sig_raw.fa" 2> "$RGW/subseq_A.err"
"$SK" subseq --bed "$RGW/sig_k33.C.bed" "$RGW/chrom.fa" > "$RGW/S2_sig_raw.fa" 2> "$RGW/subseq_C.err"
echo "raw sig seqs: S1=$(grep -c '^>' "$RGW/S1_sig_raw.fa")  S2=$(grep -c '^>' "$RGW/S2_sig_raw.fa")"

log "clean signature sets"
"$PY" "$A/make_clean_signature_sets.py" "$RGW/S1_sig_raw.fa" "$RGW/S2_sig_raw.fa" \
      "$RGW/S1_sig_k33.txt" "$RGW/S2_sig_k33.txt"
n1=$(wc -l < "$RGW/S1_sig_k33.txt"); n2=$(wc -l < "$RGW/S2_sig_k33.txt")
echo "clean sig 33-mers: S1=$n1  S2=$n2"
{ [ "$n1" -gt 1000 ] && [ "$n2" -gt 1000 ]; } || { echo "!! sig sets too small -- abort"; exit 3; }

log "classify reads (k=33, min-hits 3, ratio 3.0)"
"$PY" "$A/classify_reads_by_kmers.py" "$RGW/S1_sig_k33.txt" "$RGW/S2_sig_k33.txt" \
      "$READS" 33 "$RGW/read_labels_refguided.tsv" --min-hits 3 --ratio 3.0 --min-len 2000 --stride 1

log "EVALUATE: UPPER BOUND vs S1/S2 truth"
cut -f1,4 "$RGW/read_labels_refguided.tsv" > "$RGW/read_labels_refguided.2col.tsv"
"$PY" "$HERE/pipeline/allread_eval.py" --labels "$RGW/read_labels_refguided.2col.tsv" \
      --ref-paf "$WORK/reads_to_ref.paf" --total-reads "$TOTAL" --chrom-subg "$CHROM_SUBG"
log "DONE refguided (upper bound)"
