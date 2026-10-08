#!/usr/bin/env bash
set -uo pipefail
W=$1; CHROM=$2
PIPE=$POLYSPLIT/pipeline; RG=$POLYSPLIT/refguided; PY=python3
PAF=$W/reads_to_ref.paf
R2C=$W/reads_to_contigs.paf
LAB=$W/all_contig_labels_repaired.tsv
BLOCKS=$W/decloud_structural_contig_labels.tsv
FASTA=$W/flye_out/assembly.fasta
TMP=${TMPDIR:-/tmp}
export POLYSPLIT_TRUTH=/nonexistent

echo "== t / rho (PolySplit-Ref, re-threshold cached hits) =="
"$PY" "$RG/rethreshold_hits.py" --hits "$W/refguided/read_labels_refguided.tsv" \
      --paf "$PAF" --chrom-subg "$CHROM"

echo "== beta (PolySplit read propagation) =="
for B in 0.50 0.55 0.60 0.70 0.80 0.90; do
  "$PY" "$PIPE/propagate_to_reads.py" --paf "$R2C" --contig-labels "$LAB" \
        --min-conf "$B" --weight ident --out "$TMP/sw_b.tsv" >/dev/null 2>&1
  echo "  beta=$B  acc=$("$PY" "$PIPE/score_reads_strict.py" --labels "$TMP/sw_b.tsv" --paf "$PAF" --chrom-subg "$CHROM")"
  rm -f "$TMP/sw_b.tsv"
done

echo "== alpha (small-contig recovery floor; re-run recovery -> repair -> propagation, beta=0.60) =="
for AL in 0.40 0.50 0.55 0.60 0.70 0.80; do
  POLYSPLIT_HIC_FRAC="$AL" POLYSPLIT_CONTACTS="$W/contacts.pkl" POLYSPLIT_FASTA="$FASTA" \
    POLYSPLIT_OUT="$TMP/sw_v2.tsv" "$PY" "$PIPE/label_small_contigs_v2.py" "$BLOCKS" >/dev/null 2>&1
  POLYSPLIT_CONTACTS="$W/contacts.pkl" POLYSPLIT_EDGES="$W/homoeolog_edges.tsv" \
    POLYSPLIT_LABELS_IN="$TMP/sw_v2.tsv" POLYSPLIT_OUT="$TMP/sw_rep.tsv" \
    "$PY" "$PIPE/homoeolog_repair.py" >/dev/null 2>&1
  "$PY" "$PIPE/propagate_to_reads.py" --paf "$R2C" --contig-labels "$TMP/sw_rep.tsv" \
        --min-conf 0.60 --weight ident --out "$TMP/sw_a.tsv" >/dev/null 2>&1
  echo "  alpha=$AL  acc=$("$PY" "$PIPE/score_reads_strict.py" --labels "$TMP/sw_a.tsv" --paf "$PAF" --chrom-subg "$CHROM")"
  rm -f "$TMP/sw_v2.tsv" "$TMP/sw_rep.tsv" "$TMP/sw_a.tsv"
done
