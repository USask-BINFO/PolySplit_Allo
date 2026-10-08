#!/usr/bin/env bash
set -uo pipefail
ASM=""; BAM=""; NSG=""; RUN=""; SCORE=""; WORK=""; GROUPS=""; PRE=scaf
while [ $# -gt 0 ]; do
  case "$1" in
    --asm) ASM=$2; shift 2;; --hic-bam) BAM=$2; shift 2;; --nsg) NSG=$2; shift 2;;
    --run-dir) RUN=$2; shift 2;; --score) SCORE=$2; shift 2;; --work) WORK=$2; shift 2;;
    --groups) GROUPS=$2; shift 2;; --prefix) PRE=$2; shift 2;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done
PY=python3; PIPE=$POLYSPLIT/pipeline; THREADS=${THREADS:-16}; MINLEN=${MINLEN:-5000000}
CONDA=$CONDA_PREFIX/etc/profile.d/conda.sh; SPENV=subphaser_core
mkdir -p "$WORK"; cd "$WORK"
SCAF=${PRE}_scaffolds_final.fa

[ -s "$SCAF" ] || yahs "$ASM" "$BAM" -o "$PRE"
[ -s "$SCAF.fai" ] || samtools faidx "$SCAF"

if [ -z "$GROUPS" ]; then
  if [ ! -s groups.config ]; then
    awk -v m=$MINLEN '$2>=m{print $1}' "$SCAF.fai" > big.txt
    samtools faidx "$SCAF" $(cat big.txt) > big.fa
    minimap2 -x asm10 -t "$THREADS" big.fa big.fa > scaf_self.paf
    awk -F'\t' '$1!=$6{a=($1<$6)?$1:$6; b=($1<$6)?$6:$1; s[a"\t"b]+=$10} END{for(k in s) print s[k]"\t"k}' scaf_self.paf \
      | sort -k1,1nr > scaf_pairs.tsv
    "$PY" "$POLYSPLIT/baselines/subphaser/hexaploid/build_config.py" "$SCAF.fai" scaf_pairs.tsv groups.config "$MINLEN" "$NSG"
  fi
  GROUPS=$WORK/groups.config
fi

SGOUT=$(find "$WORK" -name '*chrom-subgenome.tsv' | head -1)
if [ -z "$SGOUT" ]; then
  set +u; source "$CONDA"; conda activate "$SPENV"; set -u
  subphaser -i "$SCAF" -c "$GROUPS" -nsg "$NSG" -just_core -pre ${PRE}_ -o phase-results -tmpdir tmp_subphaser -p "$THREADS"
  set +u; conda deactivate; set -u
  SGOUT=$(find "$WORK" -name '*chrom-subgenome.tsv' | head -1)
fi
[ -n "$SGOUT" ] || { echo "!! SubPhaser produced no *chrom-subgenome.tsv"; exit 3; }

"$PY" "$SCORE" "$SGOUT" ${PRE}_scaffolds_final.agp subphaser_contig_labels.tsv "$RUN/wg_purity.per_contig.tsv"
"$PY" "$PIPE/propagate_to_reads.py" --paf "$RUN/reads_to_contigs.paf" --contig-labels subphaser_contig_labels.tsv \
      --min-conf 0.6 --weight ident --out read_subphaser.tsv
echo "read labels -> $WORK/read_subphaser.tsv"
