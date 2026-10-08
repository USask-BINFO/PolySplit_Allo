#!/usr/bin/env bash
set -euo pipefail

GENOME=${1:?usage: run_polycracker_napus.sh <genome.fa> <n_subgenomes> [chunk_bp]}
NSUBG=${2:?need n_subgenomes (2 for A/C, 3 for triplication)}
CHUNK=${3:-100000}
IMG=sgordon/polycracker-miniconda:1.0.3
WORK=${WORK:-$(pwd)/pc_work}
mkdir -p "$WORK"
[ -s "$WORK/$GENOME" ] || { echo "ERROR: put $GENOME in $WORK first"; exit 1; }

cat > "$WORK/config_polyCRACKER.txt" <<EOF
blastPath = ./blast_files/
kmercountPath = ./kmercount_files/
fastaPath = ./fasta_files/
bedPath = ./bed_files/
genome = $GENOME
local = 1
BB = 1
n_subgenomes = $NSUBG
n_dimensions = $((NSUBG+2))
splitFasta = 1
preFilter = 0
splitFastaLineLength = $CHUNK
writeKmer = 1
kmerLength = 26
kmer2Fasta = 1
kmer_low_count = 30
use_high_count = 0
kmer_high_count = 2000000
sampling_sensitivity = 1
original = 0
writeBlast = 1
k_search_length = 13
runBlastParallel = 0
blastMemory = 40
threads = 16
blast2bed = 1
generateClusteringMatrix = 1
lowMemory = 0
minChunkSize = $CHUNK
removeNonChunk = 1
minChunkThreshold = 0
tfidf = 1
perfect_mode = 0
transformData = 1
reduction_techniques = tsne
transformMetric = linear
ClusterAll = 1
clusterMethods = SpectralClustering
grabAllClusters = 1
n_neighbors = 20
metric = cosine
weighted_nn = 0
mst = 0
extract = 1
diff_kmer_threshold = 20
default_kmercount_value = 3
diff_sample_rate = 1
unionbed_threshold = 10,2
bootstrap = 0
EOF

OUTDIR="analysisOutputs_${GENOME%.fa}_n${NSUBG}"
echo "[polycracker] $GENOME  n_subgenomes=$NSUBG  chunk=$CHUNK  -> $WORK/$OUTDIR"
udocker pull "$IMG"
RUNDIR="pcrun_${GENOME%.fa}_n${NSUBG}"
udocker run --rm -w / -v "$WORK":/work "$IMG" bash -lc "
  set -e ; export PATH=/opt/conda/bin:\$PATH
  export _JAVA_OPTIONS='-Xms3G -Xmx${MEM:-40}G'
  export TMPDIR=/work/$RUNDIR/tmp
  if [ ! -d /work/$RUNDIR ] || [ -n \"${FRESH:-}\" ]; then
    rm -rf /work/$RUNDIR ; cp -R /workdir/polycracker /work/$RUNDIR ;
  fi
  cd /work/$RUNDIR
  mkdir -p ./fasta_files ./kmercount_files ./blast_files ./bed_files ./analysisOutputs ./tmp
  cp /work/$GENOME ./fasta_files/ ;
  cp /work/config_polyCRACKER.txt ./config_polyCRACKER.txt ;
  nextflow run polycracker.nf -process.echo true -resume ;
  rm -rf /work/$OUTDIR ; cp -R analysisOutputs /work/$OUTDIR
"
echo "DONE -> $WORK/$OUTDIR (cluster results + extractedSubgenomes/*.fa)"
