#!/usr/bin/env bash
set -uo pipefail
# ultra-long ONT + Hi-C, SRR27463913 and SRR27463924
D=$DATA/millet
bash $POLYSPLIT/drivers/run_polysplit.sh --reads $D/AJ8_ont.fastq --hic1 $D/AJ8_hic_R1.fastq --hic2 $D/AJ8_hic_R2.fastq \
     --type ont --nsg 2 --flye-preset nano-hq --work $D/polysplit_run_ont \
     --ref $D/GCA_038442765.1.fna --chrom-subg $POLYSPLIT/tables/chrom_subg.millet.tsv
