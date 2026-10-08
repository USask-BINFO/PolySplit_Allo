#!/usr/bin/env bash
set -euo pipefail
BASE=$DATA/camelina/Tetraploid_data/polycracker_baseline
DRIVER=$DATA/polycracker_compare/run_polycracker_napus.sh
FLYE=$DATA/camelina/Tetraploid_data/polysplit_run/flye_out/assembly.fasta
GENOME=Cmi4x_flye.fa
WORK="$BASE/pc_work"
mkdir -p "$WORK"

# stage the SAME assembly PolySplit used (real file, not symlink -- udocker mounts $WORK)
[ -s "$WORK/$GENOME" ] || cp "$FLYE" "$WORK/$GENOME"
echo "staged $GENOME ($(du -h "$WORK/$GENOME" | cut -f1)) in $WORK"
md5sum "$FLYE" "$WORK/$GENOME"   # confirm identical to PolySplit's input

# run polyCRACKER: genome, n_subgenomes=2, chunk=100kb  (identical config to napus)
WORK="$WORK" bash "$DRIVER" "$GENOME" 2 100000

echo "DONE -> $WORK/analysisOutputs_${GENOME%.fa}_n2/"
echo "next: bash run_eval.sh"
