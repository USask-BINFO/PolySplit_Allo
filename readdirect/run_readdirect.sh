#!/usr/bin/env bash
set -uo pipefail
READS=$1; PAF=$2; CHROM=$3; NSG=$4; W=$5
RD=$POLYSPLIT/readdirect
PY=python3
KMC=$KMC_BIN/kmc; KMCT=$KMC_BIN/kmc_tools; KMCD=$KMC_BIN/kmc_dump
K=21; CI=3; CS=100000; THREADS=16; MEM=64
mkdir -p "$W/kmc_tmp"

if [ ! -s "$W/reads_k$K.kmc_pre" ]; then
  echo "$READS" > "$W/reads.lst"
  "$KMC" -k$K -ci$CI -cs$CS -t$THREADS -m$MEM -fq @"$W/reads.lst" "$W/reads_k$K" "$W/kmc_tmp"
fi
[ -s "$W/reads_k$K.histo.txt" ] || "$KMCT" transform "$W/reads_k$K" histogram "$W/reads_k$K.histo.txt" -cx$CS
read LO HI < <("$PY" "$RD/pick_band.py" "$W/reads_k$K.histo.txt" "$CI" 300)
echo "band [$LO,$HI]"
[ -s "$W/band_k$K.txt" ] || "$KMCD" -ci$LO -cx$HI "$W/reads_k$K" "$W/band_k$K.txt"

"$PY" "$RD/read_cluster_eval.py" --reads "$READS" --band "$W/band_k$K.txt" \
      --paf "$PAF" --chrom-subg "$CHROM" --k $K --nsg $NSG \
      --max-feats 3000000 --sample-reads 80000 --stride 1 --svd 40
