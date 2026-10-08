#!/usr/bin/env bash
set -uo pipefail
REF=$1; CS=$2; ORG=$3; OUT=$4; T=${5:-16}
PY=python3; PIPE=$POLYSPLIT/pipeline
MINIMAP2=minimap2; SAMTOOLS=samtools; WINDOWMASKER=windowmasker
mkdir -p "$OUT/windows"; cd "$OUT"

if [[ "$REF" == *.gz ]]; then [ -s ref.fa ] || zcat "$REF" > ref.fa; else ln -sf "$REF" ref.fa; fi
[ -s ref.fa.fai ] || "$SAMTOOLS" faidx ref.fa

# homoeolog identity of each window
[ -s windows/win.done ] || { "$PY" "$PIPE/ref_windows.py" ref.fa "$CS" windows 10000 && touch windows/win.done; }
for W in windows/win.*.fa; do
  S=$(basename "$W" .fa); S=${S#win.}
  [ -s "windows/cross.$S.paf" ] && continue
  cat $(ls windows/chroms.*.fa | grep -v "chroms.$S.fa") > "windows/other.$S.fa"
  "$MINIMAP2" -cx asm20 --secondary=no -t "$T" "windows/other.$S.fa" "$W" > "windows/cross.$S.paf"
  rm -f "windows/other.$S.fa"
done
[ -s windows/windows.tsv ] || "$PY" "$PIPE/window_identity.py" windows

# repeats
[ -s wm.counts ] || "$WINDOWMASKER" -mk_counts -in ref.fa -out wm.counts -infmt fasta -sformat obinary
[ -s wm.intervals ] || "$WINDOWMASKER" -ustat wm.counts -in ref.fa -outfmt interval -out wm.intervals

# unplaced scaffolds vs organelle genomes
cut -f1 ref.fa.fai | grep -vxFf <(cut -f1 "$CS") > unplaced.names
[ -s unplaced_vs_organelles.paf ] || { "$SAMTOOLS" faidx ref.fa -r unplaced.names > unplaced.fa
  "$MINIMAP2" -cx asm10 -t "$T" "$ORG" unplaced.fa > unplaced_vs_organelles.paf; }
