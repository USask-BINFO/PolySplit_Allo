#!/usr/bin/env python3
import sys, os, glob
from collections import defaultdict

CLUSTER_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "pc_work/analysisOutputs_NAM29_AC_n2/SpectralClusteringmain_tsne_2_n3/clusterResults"
OUT = sys.argv[2] if len(sys.argv) > 2 else "polycracker_contig_labels.tsv"
TRUTH = "$DATA/Napus_nam_denovo/sa_read_correct/purity_compare/wg_purity.per_contig.tsv"

# oracle cluster -> subgenome (from eval_polycracker_denovo.py: sg0 majority A, sg1 majority C)
ORACLE = {"subgenome_0": "A", "subgenome_1": "C"}

truth = {}
with open(TRUTH) as f:
    next(f)
    for ln in f:
        c, L, a, cc, mf, kls = ln.rstrip("\n").split("\t")
        truth[c] = {"pure_A": "A", "pure_C": "C"}.get(kls, "chimeric")

def chunk_to_contig(chunk):
    contig, start, end = chunk.rsplit("_", 2)
    return contig, int(end) - int(start)

contig_cluster_bp = defaultdict(lambda: defaultdict(int))
for fp in sorted(glob.glob(os.path.join(CLUSTER_DIR, "subgenome_*.txt"))):
    name = os.path.basename(fp).replace(".txt", "")
    for ln in open(fp):
        ch = ln.strip()
        if not ch: continue
        contig, bp = chunk_to_contig(ch)
        contig_cluster_bp[contig][name] += bp

n = defaultdict(int)
with open(OUT, "w") as o:
    o.write("contig\tlabel\ttruth\n")
    for contig, cb in contig_cluster_bp.items():
        best = max(cb.items(), key=lambda kv: kv[1])[0]
        label = ORACLE[best]
        n[label] += 1
        o.write(f"{contig}\t{label}\t{truth.get(contig, '?')}\n")

print(f"[polycracker labels] wrote {OUT}: {sum(n.values())} contigs  (A={n['A']}, C={n['C']})", file=sys.stderr)
