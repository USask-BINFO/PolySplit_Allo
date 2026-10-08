#!/usr/bin/env bash
set -euo pipefail
SP=$SUBPHASER
CONDA=$CONDA_PREFIX/etc/profile.d/conda.sh
ENV=subphaser_core    # MINIMAL core env (the full SubPhaser.yaml has stale, unsolvable pins)
HERE=$(cd "$(dirname "$0")" && pwd)
SCAF=$DATA/Napus_nam_denovo/sa_read_correct/hic_ac_cluster/yahs/napus_scaffolds_final.fa
CFG="$HERE/denovo_pairs_sg.config"
THREADS=${THREADS:-16}

source "$CONDA"
if ! conda env list | grep -qiE "(^|/)$ENV([ /]|$)"; then
  echo "[env] creating MINIMAL core env '$ENV' (jellyfish + sci stack; skips the unsolvable R/circos/perl deps) ..."
  conda create -y -n "$ENV" -c conda-forge -c bioconda \
    python=3.9 jellyfish scikit-learn numpy scipy pandas biopython matplotlib networkx xopen
fi
conda activate "$ENV"
command -v subphaser >/dev/null 2>&1 || { echo "[install] pip install $SP"; pip install -e "$SP"; }
# NOTE: -just_core writes *.chrom-subgenome.tsv BEFORE the (R) heatmap step, so even if the heatmap
# errors for lack of R in this minimal env, the per-scaffold assignment we need is already on disk.

cd "$HERE"
DT=$(date +"%y%m%d%H%M")
echo "[subphaser] -i $(basename "$SCAF")  -c $(basename "$CFG")  -nsg 2 -just_core"
subphaser -i "$SCAF" -c "$CFG" -nsg 2 -just_core -pre napusscaf_ \
    -o phase-results_denovo -tmpdir tmp_subphaser -p "$THREADS" \
    2>&1 | tee "subphaser_denovo.log.$DT"

echo
echo "DONE. Now score (base env):"
echo "  conda deactivate"
echo "  SG=\$(ls phase-results_denovo/*chrom-subgenome.tsv | head -1)"
echo "  python3 score_subphaser_denovo.py \"\$SG\" subphaser_contig_labels.tsv"
echo "  # read-level (same Step-7 path as polyCRACKER):"
echo "  cd .. && python3 propagate_to_reads.py --paf reads_to_contigs.paf \\"
echo "      --contig-labels baselines/subphaser_contig_labels.tsv --min-conf 0.6 --weight ident \\"
echo "      --out read_AC_subphaser.tsv && \\"
echo "  python3 allread_eval.py --labels read_AC_subphaser.tsv \\"
echo "      --ref-paf $DATA/Napus_nam_denovo/hybrid_method_code/long_reads_to_ref.primary.paf \\"
echo "      --total-reads 2301470"
