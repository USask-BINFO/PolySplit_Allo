#!/usr/bin/env bash
set -uo pipefail
# HiFi + Omni-C, SRR32133635 and SRR32133636
D=$DATA/coffee
bash $POLYSPLIT/drivers/run_polysplit.sh --reads $D/Typica_hifi.fastq --hic1 $D/Typica_omnic_R1.fastq --hic2 $D/Typica_omnic_R2.fastq \
     --type hifi --nsg 2 --work $D/polysplit_run_hifi \
     --ref $D/GCA_049114775.1.fna --chrom-subg $POLYSPLIT/tables/chrom_subg.coffee.tsv
