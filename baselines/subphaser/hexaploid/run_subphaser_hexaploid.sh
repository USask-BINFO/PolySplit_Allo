#!/usr/bin/env bash
set -uo pipefail
HX=$DATA/camelina/Hexaploid_data
WORK=$HX/subphaser_baseline
RUN=$HX/polysplit_run_hifi                       # HiFi PolySplit run (Hi-C BAM, read alignments, truth)
PIPE=$POLYSPLIT/pipeline
PY=python3
CONDA=$CONDA_PREFIX/etc/profile.d/conda.sh
SPENV=subphaser_core
SP=$SUBPHASER

ASM=$HX/hifi_asm_flye/assembly.fasta             # the HiFi assembly PolySplit used
HICBAM=$RUN/hic.nsort.bam                        # name-sorted bwa -5SP, aligned to the HiFi assembly
TRUTH=$RUN/wg_purity.per_contig.tsv
R2C=$RUN/reads_to_contigs.paf
R2R=$RUN/reads_to_ref.paf
CHROM=$POLYSPLIT/chrom_subg.hexaploid.tsv
TOTAL=925873
NSG=3
MINLEN=${MINLEN:-5000000}
THREADS=${THREADS:-16}

YAHS=yahs
SAMTOOLS=samtools
MINIMAP2=minimap2

SCAF=$WORK/CmiT1_scaffolds_final.fa
AGP=$WORK/CmiT1_scaffolds_final.agp
log(){ echo "### $* | $(date) ###"; }
mkdir -p "$WORK"

# ---------------- Stage 1: YaHS scaffolding ----------------
if [ ! -s "$SCAF" ]; then
  log "Stage 1: YaHS"
  [ -s "$ASM.fai" ] || "$SAMTOOLS" faidx "$ASM"
  "$YAHS" "$ASM" "$HICBAM" -o "$WORK/CmiT1"
fi
[ -s "$SCAF.fai" ] || "$SAMTOOLS" faidx "$SCAF"
log "scaffolds: $(grep -c '^>' "$SCAF")   >= ${MINLEN} bp: $(awk -v m=$MINLEN '$2>=m' "$SCAF.fai" | wc -l)"

# ---------------- Stage 2: self-align big scaffolds -> scaf_pairs.tsv ----------------
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
    a, b = sorted((q, t)); agg[(a, b)] += m
with open(out, "w") as o:
    for (a, b), s in sorted(agg.items(), key=lambda kv: -kv[1]):
        o.write(f"{s}\t{a}\t{b}\n")
print(f"[scaf_pairs] {len(agg)} scaffold pairs -> {out}", file=sys.stderr)
PY
fi

# ---------------- Stage 3: SubPhaser homoeolog-triple config (NSG=3) ----------------
if [ ! -s "$WORK/pairs.config" ]; then
  log "Stage 3: build SubPhaser config (triples)"
  "$PY" "$WORK/build_config.py" "$SCAF.fai" "$WORK/scaf_pairs.tsv" "$WORK/pairs.config" "$MINLEN" "$NSG"
fi
[ -s "$WORK/pairs.config" ] || { echo "!! empty config -- lower MINLEN (few big scaffolds?)"; exit 2; }

# ---------------- Stage 4: SubPhaser -just_core -nsg 3 ----------------
# the R heatmap step may fail; the tsv is written before it
SGOUT=$(find "$WORK" -name '*chrom-subgenome.tsv' 2>/dev/null | head -1)
if [ -z "$SGOUT" ]; then
  log "Stage 4: SubPhaser -just_core -nsg $NSG"
  set +u; source "$CONDA"; conda activate "$SPENV"; set -u
  command -v subphaser >/dev/null 2>&1 || pip install -e "$SP"
  ( cd "$WORK"
    subphaser -i "$SCAF" -c "$WORK/pairs.config" -nsg "$NSG" -just_core -pre CmiT1_ \
      -o phase-results -tmpdir tmp_subphaser -p "$THREADS" )
  set +u; conda deactivate; set -u
  SGOUT=$(find "$WORK" -name '*chrom-subgenome.tsv' 2>/dev/null | head -1)
fi
[ -n "$SGOUT" ] || { echo "!! SubPhaser produced no *chrom-subgenome.tsv -- check tmp_subphaser/log"; exit 3; }
log "SubPhaser output: $SGOUT"

# ---------------- Stage 5: contig-level score ----------------
log "Stage 5: CONTIG accuracy"
"$PY" "$WORK/score_subphaser_hexaploid.py" "$SGOUT" "$AGP" "$WORK/subphaser_contig_labels.tsv" "$TRUTH"

# ---------------- Stage 6: read-level (identical Step-7 path) ----------------
log "Stage 6: READ accuracy"
"$PY" "$PIPE/propagate_to_reads.py" --paf "$R2C" \
      --contig-labels "$WORK/subphaser_contig_labels.tsv" --min-conf 0.6 --weight ident \
      --out "$WORK/read_subphaser.tsv"
"$PY" "$PIPE/allread_eval.py" --labels "$WORK/read_subphaser.tsv" \
      --ref-paf "$R2R" --total-reads "$TOTAL" --chrom-subg "$CHROM"
echo
echo "TABLE: contig acc = Stage-5 number; read acc = correct/n_tru (NOT the printed T2,"
echo "since SubPhaser leaves reads unlabelled -- same convention as polyCRACKER)."
log "DONE subphaser hexaploid baseline"
