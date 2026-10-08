#!/usr/bin/env bash
set -uo pipefail
RD=$1; ASM=$2; NSG=$3; OUT=$4
PY=python3; PIPE=$POLYSPLIT/pipeline
for COND in default contacts2 contacts5 contacts10 tau5e4 tau2e5; do
  W=$OUT/$COND; mkdir -p "$W"
  [ -s "$W/read_subg.tsv" ] && continue
  CONTACTS=$RD/contacts.pkl; PAIR=100000
  case $COND in
    contacts*) "$PY" "$PIPE/filter_contacts.py" "$RD/contacts.pkl" "$W/contacts.pkl" "${COND#contacts}"; CONTACTS=$W/contacts.pkl;;
    tau5e4) PAIR=50000;;
    tau2e5) PAIR=200000;;
  esac
  ( cd "$W"
    export PYTHONPATH=$PIPE POLYSPLIT_NSG=$NSG POLYSPLIT_FASTA=$ASM POLYSPLIT_CONTACTS=$CONTACTS \
           POLYSPLIT_EDGES=$RD/homoeolog_edges.tsv POLYSPLIT_TRUTH=/none POLYSPLIT_PAIR_MIN=$PAIR
    "$PY" "$PIPE/dechimerize_structural.py" > dechim.log 2>&1
    POLYSPLIT_OUT=$W/all_contig_labels_v2.tsv "$PY" "$PIPE/label_small_contigs_v2.py" decloud_structural_contig_labels.tsv > recover.log 2>&1
    POLYSPLIT_LABELS_IN=$W/all_contig_labels_v2.tsv POLYSPLIT_OUT=$W/all_contig_labels_repaired.tsv \
      "$PY" "$PIPE/homoeolog_repair.py" > repair.log 2>&1
    "$PY" "$PIPE/propagate_to_reads.py" --paf "$RD/reads_to_contigs.paf" --contig-labels all_contig_labels_repaired.tsv \
          --min-conf 0.60 --weight ident --out read_subg.tsv > propagate.log 2>&1 )
  echo "$COND done"
done
