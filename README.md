# PolySplit-Allo

Reference-free subgenome separation of allopolyploid long reads by assembly-first
homoeolog phasing. This repository accompanies the paper accepted at IEEE BIBM 2026.

The paper evaluates the method on five allopolyploid genomes of four species: *Brassica napus* NAM0 (4x),
*Camelina microcarpa* (4x and 6x), *Coffea arabica* 'Kona Typica' (4x) and proso millet, *Panicum miliaceum*
AJ8 (4x).

<p align="center"><img src="figures/Denovo_assembly_polysplit_V2.png" alt="Full PolySplit pipeline schematic" width="100%"></p>

## Layout
- `pipeline/`  core PolySplit stages
  - `homoeolog_graph_from_lcp.py`  shared-33mer homoeolog edges from a gsufsort GSA+LCP
  - `dechimerize_structural.py`    Hi-C blocks + structural de-chimerization (imports `decloud_blocks_v2.py`)
  - `label_small_contigs_v2.py`    repeat-composition contrast labelling + small-contig recovery
  - `homoeolog_repair.py`          homoeolog-consistency repair
  - `propagate_to_reads.py`        identity-weighted contig-label -> read propagation
  - `allread_eval.py`, `eval_contig_labels.py`  read- and contig-level evaluation
  - `read_accounting.py`           read accounting and read-level scoring by MAPQ and homoeolog identity (Tables S8, S9)
  - `ref_windows.py`, `window_identity.py`  homoeolog identity of 10-kb reference windows
  - `unlabelled_reads.py`          what the reads without a truth label are (Table S8)
  - `validate_edges.py`            homoeolog edges and Hi-C normalisation checked against the reference (Tables S10, S11)
  - `filter_contacts.py`, `score_sensitivity.py`  Hi-C contact-count and tau_H sensitivity (Table S11)
- `refguided/`  reference-guided signature baseline (K=2 and K=3)
- `drivers/`    `run_polysplit.sh` (generalized: `--type ont|hifi`, `--nsg K`, `--assembler flye|hifiasm`) and `run_separated_assembly.sh`; plus per-dataset scripts that reproduce the paper runs (PolySplit, reference-guided, sensitivity, timing)
  - `run_reference_annotation.sh`  window homoeolog identity, WindowMasker repeats and organelle check of unplaced scaffolds (evaluation only)
  - `sensitivity_blocks.sh`        re-run from the Hi-C blocks onward with a minimum contact count or another tau_H
- `baselines/`  wrappers + cluster-to-subgenome mapping + read-eval for the comparison methods
  - `polycracker/{napus,tetraploid,hexaploid}/`  run polyCRACKER, map clusters to subgenomes (best 1:1), score
  - `subphaser/{napus,tetraploid,hexaploid}/`     YaHS scaffold + SubPhaser (scaffold-first), score
  - `subphaser/run_scaffold_first.sh`             the same on any assembly, with automatic or given homoeolog groups
  - `subphaser/mixed_scaffolds.py`, `check_groups.py`, `true_groups_config.py`  mixed scaffolds, correctness of the
    automatic homoeolog groups, and the true groups from the reference (Table S12)
  These call the external tools polyCRACKER and SubPhaser, which must be installed separately.
- `readdirect/`  assembly-free control: cluster raw reads by single-copy k-mer incidence and score vs truth (shows the assembly step is necessary)
- `tables/`     `datasets.template.tsv` (run directories and label files per dataset), `make_table_s7.py`, `make_tables_s8_s12.py`,
  and the chromosome -> subgenome maps
- `figures/`    figures of the paper and supplement; `make_kmer_survival.py`, `make_edge_weights.py`, `make_confusion.py` draw Figs S3-S6

## Dependencies
Python 3 (numpy, scipy, scikit-learn, networkx, matplotlib) and the external tools Flye, bwa, samtools,
minimap2, gsufsort (32- and 64-bit builds), KMC, YaHS, SubPhaser (baseline only), hifiasm (optional) and
WindowMasker from BLAST+ (evaluation only).

## Usage
Inputs: pooled long reads (one FASTQ), a Hi-C/Omni-C read pair, and the number of subgenomes
`K`. No reference genome is needed to run the tool; a reference is only used to score accuracy.

**1. Configure.** Edit `config.sh` to point `DATA`, `POLYSPLIT`, and `FLYE_BIN` at your
environment (`SUBPHASER` and `KMC_BIN` are only needed for the baselines and the read-direct
control), then source it:
```
source config.sh
```

**2. Separate the reads into subgenomes.**
```
bash drivers/run_polysplit.sh \
     --reads reads.fq.gz --hic1 hic_R1.fq.gz --hic2 hic_R2.fq.gz \
     --type ont --nsg 2 --work out/
```
- `--type`  `ont` (Flye `--nano-hq`, minimap `map-ont`) or `hifi` (`--pacbio-hifi`, `map-hifi`)
- `--nsg`   number of subgenomes to separate (2 for a tetraploid, 3 for a hexaploid)
- `--work`  output directory
- `--flye-preset nano-raw`  optional; overrides the default assembly preset
- `--assembler hifiasm`  optional; pre-assemble with hifiasm (`-l0`) instead of Flye (HiFi reads only)

This runs the full pipeline (assemble -> Hi-C blocks -> homoeolog pairing -> subgenome
labelling -> propagate) and writes per-read subgenome labels to `out/read_subg.tsv`.

**3. Assemble each subgenome (optional).**
```
bash drivers/run_separated_assembly.sh \
     --labels out/read_subg.tsv --reads reads.fq.gz \
     --type ont --work out/ --hic1 hic_R1.fq.gz --hic2 hic_R2.fq.gz
```
Splits the reads by label, assembles each subgenome with Flye, scaffolds with YaHS, and writes
`out/sep_asm/summary.tsv` (size, contigs, N50, read back-mapping rate, subgenome purity).

**Scoring accuracy (benchmarking only).** Add `--ref REF.fa --chrom-subg map.tsv` to either
command, where `map.tsv` is two columns, chromosome and its subgenome (`chrom<TAB>subgenome`).
`run_polysplit.sh` then also reports read-level accuracy, precision, recall, and F1.

**Reproducing the paper.** The per-dataset `drivers/run_polysplit_<name>.sh` scripts run the
exact commands (paths and presets) used for each dataset in the manuscript (`nam0`, `tetraploid`, `hexaploid`,
`coffee_hifi`, `millet_hifi`, `millet_ont`). The tetraploid and
NAM0 ONT runs used `--nano-raw`, reproducible with `run_polysplit.sh ... --flye-preset nano-raw`.

**Supplementary tables S7-S12.** Copy `tables/datasets.template.tsv`, fill in the paths, and run (all outputs go to
`results/`; `organelles.fa` holds the chloroplast and mitochondrial genomes listed in the script header):
```
bash drivers/run_reference_annotation.sh REF.fa tables/chrom_subg.<x>.tsv organelles.fa results/reference/<ref_id>
python3 pipeline/read_accounting.py datasets.tsv results
python3 pipeline/unlabelled_reads.py datasets.tsv results
python3 pipeline/validate_edges.py datasets.tsv results
bash drivers/sensitivity_blocks.sh RUN_DIR ASSEMBLY K results/sensitivity/<ref_id>_<ont|hifi>
python3 pipeline/score_sensitivity.py datasets.tsv results
python3 tables/make_table_s7.py datasets.tsv
python3 tables/make_tables_s8_s12.py datasets.tsv results
```
For Table S12, write `mixed_scaffolds.py` and `check_groups.py` output to `results/scaffold_first/<ref_id>_hifi/`
(`mixed.json`, `groups.json`), run `run_scaffold_first.sh --groups` with the config from `true_groups_config.py`, and
score its read labels with `read_accounting.py --label <ref_id>_hifi "SubPhaser (true groups)" read_subphaser.tsv`.

Every stage is idempotent: re-running skips any stage whose output already exists.

## Data
Sequencing inputs and the evaluation references are not redistributed here.

| Genome | Reads | Reference |
|---|---|---|
| *B. napus* NAM0 (4x) | see the manuscript's data-availability statement | NAM0 reference of the SKBnNAM panel |
| *C. microcarpa* (4x, 6x) | EBI-ENA PRJEB96055 | EBI-ENA PRJEB96055 |
| *C. arabica* 'Kona Typica' (4x) | NCBI PRJNA1215012: SRR32133635 (HiFi), SRR32133636 (Omni-C) | GCA_049114775.1 |
| *P. miliaceum* AJ8 (4x) | NCBI PRJNA1059665: SRR27463934 (HiFi), SRR27463913 (ONT), SRR27463924 (Hi-C) | GCA_038442765.1 |
