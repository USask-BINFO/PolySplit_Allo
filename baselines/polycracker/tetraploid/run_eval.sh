#!/usr/bin/env bash
set -uo pipefail
BASE=$DATA/camelina/Tetraploid_data/polycracker_baseline
PY=python3
# auto-discover the clusterResults dir (the SpectralClustering..._n3 name can vary)
CLUSTER="$(find "$BASE/pc_work" -type d -name clusterResults 2>/dev/null | head -1)"

if [ -z "$CLUSTER" ] || [ -z "$(ls "$CLUSTER"/subgenome_*.txt 2>/dev/null)" ]; then
  echo "!! no clusterResults/subgenome_*.txt found yet under $BASE/pc_work"
  echo "   -> polyCRACKER (00_setup_and_run.sh) has not finished. Outputs so far:"
  ls -d "$BASE"/pc_work/analysisOutputs_* 2>/dev/null || echo "   (none yet)"
  exit 2
fi
echo "using clusterResults: $CLUSTER"

echo "############ CONTIG-level oracle accuracy ############"
"$PY" "$BASE/eval_polycracker_tetraploid.py" "$CLUSTER"

echo
echo "############ per-contig S1/S2 labels (oracle) ############"
"$PY" "$BASE/make_polycracker_contig_labels.py" "$CLUSTER" "$BASE/polycracker_contig_labels.tsv"

echo
echo "############ READ-level accuracy ############"
bash "$BASE/eval_reads.sh"
