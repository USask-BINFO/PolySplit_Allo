#!/usr/bin/env bash
set -uo pipefail
HERE=$POLYSPLIT
A=$DATA
PY=python3
SK=seqkit
SAMTOOLS=samtools
GSUFSORT=gsufsort-64
DATA=$DATA/camelina/Hexaploid_data
WORK=$DATA/polysplit_run
RGW=$WORK/refguided
CHROM_SUBG=$HERE/chrom_subg.hexaploid.tsv          # Chr01-06->S1, Chr07-13->S2, Chr14-20->S3
REF=$DATA/Cmicrocarpa.CN119205.fa.gz               # eval reference = truth
READS=$DATA/TMP24026.ONT.fastq.gz
TOTAL=1903171                                       # cached total ONT reads ($WORK/total_reads.txt)
K=33
mkdir -p "$RGW"
log(){ echo "### $* | $(date) ###"; }

# ---- Stage 0: per-subgenome chromosome FASTAs (chrom.fa = doc order for the index) ----
if [ ! -s "$RGW/chrom.fa.fai" ]; then
  log "build chrom.fa + S1/S2/S3.fa from reference (chromosomes with subgenome truth only)"
  awk '$2=="S1"{print $1}' "$CHROM_SUBG" > "$RGW/s1.ids"
  awk '$2=="S2"{print $1}' "$CHROM_SUBG" > "$RGW/s2.ids"
  awk '$2=="S3"{print $1}' "$CHROM_SUBG" > "$RGW/s3.ids"
  cat "$RGW"/s1.ids "$RGW"/s2.ids "$RGW"/s3.ids > "$RGW/chrom.ids"
  "$SK" grep -f "$RGW/chrom.ids" "$REF" > "$RGW/chrom.fa" 2>/dev/null
  "$SK" grep -f "$RGW/s1.ids"    "$REF" > "$RGW/S1.fa"    2>/dev/null
  "$SK" grep -f "$RGW/s2.ids"    "$REF" > "$RGW/S2.fa"    2>/dev/null
  "$SK" grep -f "$RGW/s3.ids"    "$REF" > "$RGW/S3.fa"    2>/dev/null
  "$SAMTOOLS" faidx "$RGW/chrom.fa"
  echo "chrom.fa seqs: $(grep -c '^>' "$RGW/chrom.fa")  (S1 $(grep -c '^>' "$RGW/S1.fa") / S2 $(grep -c '^>' "$RGW/S2.fa") / S3 $(grep -c '^>' "$RGW/S3.fa"))"
fi

# ---- Stage 1: gsufsort GSA + LCP (k=33) ----
GSA="$RGW/chrom_k33.4.8.gsa"; LCP="$RGW/chrom_k33.1.lcp"
if [ ! -s "$GSA" ] || [ ! -s "$LCP" ]; then
  log "gsufsort-64 index (GSA 4,8 + 1-byte trLCP 33)"
  "$GSUFSORT" "$RGW/chrom.fa" --fasta --gsa 4 8 --lcp 1 --trlcp "$K" --output "$RGW/chrom_k33" --time
  ls -la "$GSA" "$LCP" 2>/dev/null
  [ -s "$GSA" ] && [ -s "$LCP" ] || { echo "!! gsufsort did not produce GSA/LCP -- abort (run on NUGGET)"; exit 2; }
fi

# ---- Stage 2: select subgenome-specific k=33 signatures (K=3) ----
if [ ! -s "$RGW/sig_k33.S1.bed" ]; then
  log "Select_signatures_k3 (max-copy 3, keep all)"
  "$PY" "$A/Select_signatures_k3.py" --k "$K" --gsa "$GSA" --lcp "$LCP" \
    --ac-fa "$RGW/chrom.fa" --s1-fa "$RGW/S1.fa" --s2-fa "$RGW/S2.fa" --s3-fa "$RGW/S3.fa" \
    --fai "$RGW/chrom.fa.fai" --max-copy 3 --target 200000000000 --seed 1 --out-prefix "$RGW/sig_k33"
fi
wc -l "$RGW/sig_k33.S1.bed" "$RGW/sig_k33.S2.bed" "$RGW/sig_k33.S3.bed"
for g in S1 S2 S3; do [ -s "$RGW/sig_k33.$g.bed" ] || { echo "!! $g BED empty -- abort"; exit 3; }; done

# ---- Stage 3: signature sequences -> clean (drop k-mers shared across subgenomes) ----
if [ ! -s "$RGW/S1_sig_k33.txt" ]; then
  log "extract signature sequences (seqkit subseq --bed) + clean"
  for g in S1 S2 S3; do
    "$SK" subseq --bed "$RGW/sig_k33.$g.bed" "$RGW/chrom.fa" > "$RGW/${g}_sig_raw.fa" 2> "$RGW/subseq_$g.err"
  done
  echo "raw sig seqs: S1=$(grep -c '^>' "$RGW/S1_sig_raw.fa") S2=$(grep -c '^>' "$RGW/S2_sig_raw.fa") S3=$(grep -c '^>' "$RGW/S3_sig_raw.fa")"
  "$PY" "$A/make_clean_signature_sets_k3.py" "$RGW/S1_sig_raw.fa" "$RGW/S2_sig_raw.fa" "$RGW/S3_sig_raw.fa" \
        "$RGW/S1_sig_k33.txt" "$RGW/S2_sig_k33.txt" "$RGW/S3_sig_k33.txt"
fi
n1=$(wc -l < "$RGW/S1_sig_k33.txt"); n2=$(wc -l < "$RGW/S2_sig_k33.txt"); n3=$(wc -l < "$RGW/S3_sig_k33.txt")
echo "clean sig 33-mers: S1=$n1 S2=$n2 S3=$n3"
{ [ "$n1" -gt 1000 ] && [ "$n2" -gt 1000 ] && [ "$n3" -gt 1000 ]; } || { echo "!! a sig set too small -- abort (S1/S2 may be too similar for ref-specific kmers)"; exit 4; }

# ---- Stage 4: classify reads by signature hits (k=33, min-hits 3, ratio 3.0) ----
if [ ! -s "$RGW/read_labels_refguided.tsv" ]; then
  log "classify reads (K=3)"
  "$PY" "$A/classify_reads_k3.py" "$RGW/S1_sig_k33.txt" "$RGW/S2_sig_k33.txt" "$RGW/S3_sig_k33.txt" \
        "$READS" "$K" "$RGW/read_labels_refguided.tsv" --min-hits 3 --ratio 3.0 --min-len 2000 --stride 1
fi
cut -f2 "$RGW/read_labels_refguided.tsv" | tail -n +2 | sort | uniq -c

# ---- Stage 5: evaluate vs S1/S2/S3 truth (Acc/Prec/Rec/F1) ----
log "EVALUATE reference-guided vs S1/S2/S3 truth"
cut -f1,5 "$RGW/read_labels_refguided.tsv" > "$RGW/read_labels_refguided.2col.tsv"
"$PY" "$HERE/pipeline/allread_eval.py" --labels "$RGW/read_labels_refguided.2col.tsv" \
      --ref-paf "$WORK/reads_to_ref.paf" --total-reads "$TOTAL" --chrom-subg "$CHROM_SUBG"
log "DONE refguided hexaploid (precision ceiling)"
