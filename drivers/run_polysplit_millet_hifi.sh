#!/usr/bin/env bash
set -uo pipefail
# HiFi + Hi-C, SRR27463934 and SRR27463924
D=$DATA/millet
bash $POLYSPLIT/drivers/run_polysplit.sh --reads $D/AJ8_hifi.fastq --hic1 $D/AJ8_hic_R1.fastq --hic2 $D/AJ8_hic_R2.fastq \
     --type hifi --nsg 2 --work $D/polysplit_run_hifi \
     --ref $D/GCA_038442765.1.fna --chrom-subg $POLYSPLIT/tables/chrom_subg.millet.tsv
