#!/usr/bin/env bash
set -uo pipefail
WORK=$DATA/camelina/Tetraploid_data/subphaser_baseline
RUN=$DATA/camelina/Tetraploid_data/polysplit_run
PIPE=$POLYSPLIT/pipeline
PY=python3
CONDA=$CONDA_PREFIX/etc/profile.d/conda.sh
SPENV=subphaser_core
SP=$SUBPHASER

ASM=$RUN/flye_out/assembly.fasta
HICBAM=$RUN/hic.nsort.bam
TRUTH=$RUN/wg_purity.per_contig.tsv
R2C=$RUN/reads_to_contigs.paf
R2R=$RUN/reads_to_ref.paf
CHROM=$POLYSPLIT/chrom_subg.tetraploid.tsv
TOTAL=1647060
MINLEN=${MINLEN:-5000000}
THREADS=${THREADS:-16}

YAHS=yahs
SAMTOOLS=samtools
MINIMAP2=minimap2

SCAF=$WORK/cmi4x_scaffolds_final.fa
AGP=$WORK/cmi4x_scaffolds_final.agp
log(){ echo "### $* | $(date) ###"; }

if [ ! -s "$SCAF" ]; then
  log "Stage 1: YaHS"
  [ -s "$ASM.fai" ] || "$SAMTOOLS" faidx "$ASM"
  "$YAHS" "$ASM" "$HICBAM" -o "$WORK/cmi4x"
fi
[ -s "$SCAF.fai" ] || "$SAMTOOLS" faidx "$SCAF"
log "scaffolds: $(grep -c '^>' "$SCAF")   >= ${MINLEN} bp: $(awk -v m=$MINLEN '$2>=m' "$SCAF.fai" | wc -l)"

if [ ! -s "$WORK/scaf_pairs.tsv" ]; then
  log "Stage 2: self-align chromosome-scale scaffolds"
  awk -v m=$MINLEN '$2>=m{print $1}' "$SCAF.fai" > "$WORK/big.txt"
  "$SAMTOOLS" faidx "$SCAF" $(cat "$WORK/big.txt") > "$WORK/big.fa"
  "$MINIMAP2" -x asm10 -t "$THREADS" "$WORK/big.fa" "$WORK/big.fa" > "$WORK/scaf_self.paf"
  "$PY" - "$WORK/scaf_self.paf" "$WORK/scaf_pairs.tsv" <<'PY'
import sys
from collections import defaultdict
paf, out = sys.argv[1], sys.argv[2]
agg = defaultdict(int)
for ln in open(paf):
    p = ln.split("\t")
    if len(p) < 10: continue
    q, t, m = p[0], p[5], int(p[9])
    if q == t: continue
    a, b = sorted((q, t))
    agg[(a, b)] += m
with open(out, "w") as o:
    for (a, b), s in sorted(agg.items(), key=lambda kv: -kv[1]):
        o.write(f"{s}\t{a}\t{b}\n")
print(f"[scaf_pairs] {len(agg)} scaffold pairs -> {out}", file=sys.stderr)
PY
fi

if [ ! -s "$WORK/pairs.config" ]; then
  log "Stage 3: build SubPhaser config"
  "$PY" "$WORK/build_config.py" "$SCAF.fai" "$WORK/scaf_pairs.tsv" "$WORK/pairs.config" "$MINLEN" 2
fi
[ -s "$WORK/pairs.config" ] || { echo "!! empty config -- lower MINLEN (few big scaffolds?)"; exit 2; }

SGOUT=$(find "$WORK" -name '*chrom-subgenome.tsv' 2>/dev/null | head -1)
if [ -z "$SGOUT" ]; then
  log "Stage 4: SubPhaser -just_core"
  set +u; source "$CONDA"; conda activate "$SPENV"; set -u
  command -v subphaser >/dev/null 2>&1 || pip install -e "$SP"
  ( cd "$WORK"
    subphaser -i "$SCAF" -c "$WORK/pairs.config" -nsg 2 -just_core -pre cmi4x_ \
      -o phase-results -tmpdir tmp_subphaser -p "$THREADS" )
  set +u; conda deactivate; set -u
  SGOUT=$(find "$WORK" -name '*chrom-subgenome.tsv' 2>/dev/null | head -1)
fi
[ -n "$SGOUT" ] || { echo "!! SubPhaser produced no *chrom-subgenome.tsv -- check tmp_subphaser/log"; exit 3; }
log "SubPhaser output: $SGOUT"

log "Stage 5: CONTIG accuracy"
"$PY" "$WORK/score_subphaser_tetraploid.py" "$SGOUT" "$AGP" "$WORK/subphaser_contig_labels.tsv" "$TRUTH"

log "Stage 6: READ accuracy"
"$PY" "$PIPE/propagate_to_reads.py" --paf "$R2C" \
      --contig-labels "$WORK/subphaser_contig_labels.tsv" --min-conf 0.6 --weight ident \
      --out "$WORK/read_subphaser.tsv"
"$PY" "$PIPE/allread_eval.py" --labels "$WORK/read_subphaser.tsv" \
      --ref-paf "$R2R" --total-reads "$TOTAL" --chrom-subg "$CHROM"
echo
echo "TABLE: contig acc = the Stage-5 number; read acc = correct/chromosome-truth"
echo "(use correct/n_tru, NOT the printed T2, since SubPhaser leaves reads unlabelled -- same as polyCRACKER)."
log "DONE subphaser baseline"
