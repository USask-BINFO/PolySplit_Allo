#!/usr/bin/env bash
set -uo pipefail
HERE=$POLYSPLIT
A=$DATA
PY=python3
SK=seqkit
DATA=$DATA/camelina/Tetraploid_data
SIG=$DATA/polysplit_run/refguided                  # read-independent S1/S2 signatures (built on ONT run)
WORK=$DATA/polysplit_run_hifi                      # HiFi run: truth paf + total_reads live here
RGW=$WORK/refguided                                # HiFi reference-guided outputs
CHROM_SUBG=$HERE/chrom_subg.tetraploid.tsv
READS=$DATA/CN119243.hifi.fastq.gz                 # HiFi reads
R2R=$WORK/reads_to_ref.paf                          # HiFi reads -> reference (truth)
K=33
log(){ echo "### $* | $(date) ###"; }
mkdir -p "$RGW"

for f in "$SIG/S1_sig_k33.txt" "$SIG/S2_sig_k33.txt" "$READS" "$R2R"; do
  [ -s "$f" ] || { echo "!! missing input: $f"; exit 1; }
done
echo "reusing signatures: S1=$(wc -l < "$SIG/S1_sig_k33.txt")  S2=$(wc -l < "$SIG/S2_sig_k33.txt")"

TOTAL_CACHE=$WORK/total_reads.txt                   # the HiFi run already cached this
if [ ! -s "$TOTAL_CACHE" ]; then "$SK" stats -T "$READS" | awk 'NR==2{print $4}' > "$TOTAL_CACHE"; fi
TOTAL=$(cat "$TOTAL_CACHE"); echo "total HiFi reads = $TOTAL"

# ---- classify HiFi reads by signature hits (K=2, same params as ONT) ----
log "classify HiFi reads (K=2, min-hits 3, ratio 3.0)"
if [ ! -s "$RGW/read_labels_refguided.tsv" ]; then
  "$PY" "$A/classify_reads_by_kmers.py" "$SIG/S1_sig_k33.txt" "$SIG/S2_sig_k33.txt" \
        "$READS" "$K" "$RGW/read_labels_refguided.tsv" --min-hits 3 --ratio 3.0 --min-len 2000 --stride 1
fi
cut -f4 "$RGW/read_labels_refguided.tsv" | tail -n +2 | sort | uniq -c

# ---- evaluate vs S1/S2 truth ----
log "EVALUATE reference-guided (HiFi) vs S1/S2 truth"
cut -f1,4 "$RGW/read_labels_refguided.tsv" > "$RGW/read_labels_refguided.2col.tsv"
"$PY" "$HERE/pipeline/allread_eval.py" --labels "$RGW/read_labels_refguided.2col.tsv" \
      --ref-paf "$R2R" --total-reads "$TOTAL" --chrom-subg "$CHROM_SUBG"
echo
echo "TABLE read acc = correct / chromosome-truth (same convention as the other rows)."
log "DONE refguided tetraploid HiFi"
