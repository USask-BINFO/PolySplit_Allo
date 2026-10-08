#!/usr/bin/env bash
set -uo pipefail
HERE=$POLYSPLIT
A=$DATA
N=$A/NAM0
PY=python3
SK=seqkit
MINIMAP2=minimap2

WORK=$N/polysplit_run
SIG=$WORK/refguided                                # A/C signatures built by the ONT run
RGW=$WORK/refguided_hifi                           # HiFi reference-guided outputs (separate dir)
CHROM=$N/n99_chr.fa
CHROM_SUBG=$N/chrom_subg.nam0.tsv
READS=$N/pacbio/hifi_reads.fastq.gz                # HiFi reads (19.7 Gb)
R2R=$WORK/reads_to_ref.hifi.paf                    # HiFi reads -> reference (truth) -- built here
K=33
THREADS=${THREADS:-32}
log(){ echo "### $* | $(date) ###"; }
mkdir -p "$RGW"

for f in "$SIG/A_sig_k33.txt" "$SIG/C_sig_k33.txt" "$CHROM" "$CHROM_SUBG" "$READS"; do
  [ -s "$f" ] || { echo "!! missing input: $f  (did the ONT NAM0 reference-guided run finish?)"; exit 1; }
done
echo "reusing signatures: A=$(wc -l < "$SIG/A_sig_k33.txt")  C=$(wc -l < "$SIG/C_sig_k33.txt")"

# ---- build HiFi reads -> reference truth alignment (the only piece NAM0 HiFi lacks) ----
if [ ! -s "$R2R" ]; then
  log "minimap2 HiFi reads -> N99 reference (truth PAF)"
  "$MINIMAP2" -x map-hifi -t "$THREADS" "$CHROM" "$READS" > "$R2R" 2> "$RGW/minimap2_hifi.err"
fi
[ -s "$R2R" ] || { echo "!! truth PAF build failed -- see $RGW/minimap2_hifi.err"; exit 2; }

TOTAL_CACHE=$N/n99_hifi_read_count.txt
if [ ! -s "$TOTAL_CACHE" ]; then "$SK" stats -T "$READS" | awk 'NR==2{print $4}' > "$TOTAL_CACHE"; fi
TOTAL=$(cat "$TOTAL_CACHE"); echo "total HiFi reads = $TOTAL"

# ---- classify HiFi reads by signature hits (K=2) ----
log "classify HiFi reads (K=2, min-hits 3, ratio 3.0)"
if [ ! -s "$RGW/read_labels_refguided.tsv" ]; then
  "$PY" "$A/classify_reads_by_kmers.py" "$SIG/A_sig_k33.txt" "$SIG/C_sig_k33.txt" \
        "$READS" "$K" "$RGW/read_labels_refguided.tsv" --min-hits 3 --ratio 3.0 --min-len 2000 --stride 1
fi
cut -f4 "$RGW/read_labels_refguided.tsv" | tail -n +2 | sort | uniq -c

# ---- evaluate vs A/C truth ----
log "EVALUATE reference-guided (HiFi) vs A/C truth"
cut -f1,4 "$RGW/read_labels_refguided.tsv" > "$RGW/read_labels_refguided.2col.tsv"
"$PY" "$HERE/pipeline/allread_eval.py" --labels "$RGW/read_labels_refguided.2col.tsv" \
      --ref-paf "$R2R" --total-reads "$TOTAL" --chrom-subg "$CHROM_SUBG"
echo
echo "TABLE read acc = correct / chromosome-truth (same convention as the other rows)."
log "DONE refguided NAM0 HiFi"
