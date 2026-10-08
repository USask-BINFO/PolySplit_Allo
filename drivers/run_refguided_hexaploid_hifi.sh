#!/usr/bin/env bash
set -uo pipefail
HERE=$POLYSPLIT
A=$DATA
PY=python3
SK=seqkit
DATA=$DATA/camelina/Hexaploid_data
SIG=$DATA/polysplit_run/refguided                  # read-independent signatures (built on ONT run)
WORK=$DATA/polysplit_run_hifi                      # HiFi run: truth paf lives here
RGW=$WORK/refguided                                # HiFi reference-guided outputs (separate dir)
CHROM_SUBG=$HERE/chrom_subg.hexaploid.tsv
READS=$DATA/TMP24026.hifi.fastq.gz                 # HiFi reads (8.6 GB)
R2R=$WORK/reads_to_ref.paf                         # HiFi reads -> reference (truth)
K=33
log(){ echo "### $* | $(date) ###"; }
mkdir -p "$RGW"

for f in "$SIG/S1_sig_k33.txt" "$SIG/S2_sig_k33.txt" "$SIG/S3_sig_k33.txt" "$READS" "$R2R"; do
  [ -s "$f" ] || { echo "!! missing input: $f"; exit 1; }
done
echo "reusing signatures: S1=$(wc -l < "$SIG/S1_sig_k33.txt")  S2=$(wc -l < "$SIG/S2_sig_k33.txt")  S3=$(wc -l < "$SIG/S3_sig_k33.txt")"

# total HiFi read count (T3 line only; headline is correct/n_tru). Cached.
TOTAL_CACHE=$WORK/total_reads.txt
if [ ! -s "$TOTAL_CACHE" ]; then
  log "counting HiFi reads (cached -> $TOTAL_CACHE)"
  "$SK" stats -T "$READS" | awk 'NR==2{print $4}' > "$TOTAL_CACHE"
fi
TOTAL=$(cat "$TOTAL_CACHE"); echo "total HiFi reads = $TOTAL"

# ---- classify HiFi reads by signature hits (K=3, same params as the ONT run) ----
log "classify HiFi reads (K=3, min-hits 3, ratio 3.0)"
if [ ! -s "$RGW/read_labels_refguided.tsv" ]; then
  "$PY" "$A/classify_reads_k3.py" "$SIG/S1_sig_k33.txt" "$SIG/S2_sig_k33.txt" "$SIG/S3_sig_k33.txt" \
        "$READS" "$K" "$RGW/read_labels_refguided.tsv" --min-hits 3 --ratio 3.0 --min-len 2000 --stride 1
fi
cut -f5 "$RGW/read_labels_refguided.tsv" | tail -n +2 | sort | uniq -c

# ---- evaluate vs S1/S2/S3 truth (Acc/Prec/Rec/F1) ----
log "EVALUATE reference-guided (HiFi) vs S1/S2/S3 truth"
cut -f1,5 "$RGW/read_labels_refguided.tsv" > "$RGW/read_labels_refguided.2col.tsv"
"$PY" "$HERE/pipeline/allread_eval.py" --labels "$RGW/read_labels_refguided.2col.tsv" \
      --ref-paf "$R2R" --total-reads "$TOTAL" --chrom-subg "$CHROM_SUBG"
echo
echo "TABLE read acc = correct / chromosome-truth (same convention as the other rows)."
log "DONE refguided hexaploid HiFi"
