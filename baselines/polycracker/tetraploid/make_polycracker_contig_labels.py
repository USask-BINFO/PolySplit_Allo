#!/usr/bin/env python3
import sys, os, glob
from collections import defaultdict

CLUSTER_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "pc_work/analysisOutputs_Cmi4x_flye_n2/SpectralClusteringmain_tsne_2_n3/clusterResults"
OUT = sys.argv[2] if len(sys.argv) > 2 else "polycracker_contig_labels.tsv"
TRUTH = "$DATA/camelina/Tetraploid_data/polysplit_run/wg_purity.per_contig.tsv"

truth = {}
with open(TRUTH) as f:
    next(f)
    for ln in f:
        c, L, a, cc, mf, kls = ln.rstrip("\n").split("\t")
        truth[c] = {"pure_S1": "S1", "pure_S2": "S2"}.get(kls, "chimeric")

def chunk_to_contig(chunk):
    contig, start, end = chunk.rsplit("_", 2)
    return contig, int(end) - int(start)

# read clusters
contig_cluster_bp = defaultdict(lambda: defaultdict(int))
cluster_true_bp = defaultdict(lambda: defaultdict(int))   # cluster -> {S1:bp, S2:bp}
for fp in sorted(glob.glob(os.path.join(CLUSTER_DIR, "subgenome_*.txt"))):
    name = os.path.basename(fp).replace(".txt", "")
    for ln in open(fp):
        ch = ln.strip()
        if not ch: continue
        contig, bp = chunk_to_contig(ch)
        contig_cluster_bp[contig][name] += bp
        t = truth.get(contig)
        if t in ("S1", "S2"):
            cluster_true_bp[name][t] += bp

# dynamic oracle: each cluster -> its true-majority subgenome
ORACLE = {}
for name, tb in cluster_true_bp.items():
    ORACLE[name] = "S1" if tb.get("S1", 0) >= tb.get("S2", 0) else "S2"
for name in list(contig_cluster_bp):
    for cl in contig_cluster_bp[name]:
        ORACLE.setdefault(cl, "S1")   # cluster with no pure-contig bp: arbitrary
print(f"[oracle] {ORACLE}", file=sys.stderr)

n = defaultdict(int)
with open(OUT, "w") as o:
    o.write("contig\tlabel\ttruth\n")
    for contig, cb in contig_cluster_bp.items():
        best = max(cb.items(), key=lambda kv: kv[1])[0]
        label = ORACLE[best]
        n[label] += 1
        o.write(f"{contig}\t{label}\t{truth.get(contig, '?')}\n")

print(f"[polycracker labels] wrote {OUT}: {sum(n.values())} contigs  (S1={n['S1']}, S2={n['S2']})",
      file=sys.stderr)
