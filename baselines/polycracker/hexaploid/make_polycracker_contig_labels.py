#!/usr/bin/env python3
import sys, os, glob
from collections import defaultdict

CLUSTER_DIR = sys.argv[1]
OUT = sys.argv[2]
TRUTH = sys.argv[3] if len(sys.argv) > 3 else \
    "$DATA/camelina/Hexaploid_data/polysplit_run_hifi/wg_purity.per_contig.tsv"

SUBS = ("S1", "S2", "S3")
truth = {}
with open(TRUTH) as f:
    next(f)
    for ln in f:
        c, L, a, cc, mf, kls = ln.rstrip("\n").split("\t")
        truth[c] = {f"pure_{s}": s for s in SUBS}.get(kls, "chimeric")


def chunk_to_contig(chunk):
    contig, start, end = chunk.rsplit("_", 2)
    return contig, int(end) - int(start)


contig_cluster_bp = defaultdict(lambda: defaultdict(int))
cluster_true_bp = defaultdict(lambda: defaultdict(int))
for fp in sorted(glob.glob(os.path.join(CLUSTER_DIR, "subgenome_*.txt"))):
    name = os.path.basename(fp).replace(".txt", "")
    for ln in open(fp):
        ch = ln.strip()
        if not ch:
            continue
        contig, bp = chunk_to_contig(ch)
        contig_cluster_bp[contig][name] += bp
        t = truth.get(contig)
        if t in SUBS:
            cluster_true_bp[name][t] += bp

for name, tb in cluster_true_bp.items():
    print(f"[cluster {name}] truth bp: {dict(tb)}", file=sys.stderr)

n = defaultdict(int)
with open(OUT, "w") as o:
    o.write("contig\tlabel\ttruth\n")
    for contig, cb in contig_cluster_bp.items():
        best = max(cb.items(), key=lambda kv: kv[1])[0]
        n[best] += 1
        o.write(f"{contig}\t{best}\t{truth.get(contig, '?')}\n")

print(f"[polycracker labels] wrote {OUT}: {sum(n.values())} contigs  {dict(n)}", file=sys.stderr)
