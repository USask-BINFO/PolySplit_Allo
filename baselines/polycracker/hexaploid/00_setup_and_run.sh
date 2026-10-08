#!/usr/bin/env bash
set -euo pipefail
BASE=$DATA/camelina/Hexaploid_data/polycracker_baseline
DRIVER=$DATA/polycracker_compare/run_polycracker_napus.sh
ASM=$DATA/camelina/Hexaploid_data/hifi_asm_flye/assembly.fasta
GENOME=CmiT1_hifi.fa
WORK="$BASE/pc_work"
mkdir -p "$WORK"

# stage the SAME assembly PolySplit used (real file -- udocker mounts $WORK, no symlinks)
[ -s "$WORK/$GENOME" ] || cp "$ASM" "$WORK/$GENOME"
echo "staged $GENOME ($(du -h "$WORK/$GENOME" | cut -f1)) in $WORK"
md5sum "$ASM" "$WORK/$GENOME"     # confirm identical to PolySplit's HiFi input

# run polyCRACKER: genome, n_subgenomes=3 (hexaploid), chunk=100kb (same config as napus/tetraploid)
WORK="$WORK" bash "$DRIVER" "$GENOME" 3 100000

echo "DONE -> $WORK/analysisOutputs_${GENOME%.fa}_n3/"
echo "next: bash run_eval.sh"
