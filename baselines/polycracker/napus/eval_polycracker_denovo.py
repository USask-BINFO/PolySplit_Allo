#!/usr/bin/env python3
import sys, os, glob
from collections import defaultdict

CLUSTER_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "pc_work/analysisOutputs_NAM29_AC_n2/SpectralClusteringmain_tsne_2_n3/clusterResults"
TRUTH = "$DATA/Napus_nam_denovo/sa_read_correct/purity_compare/wg_purity.per_contig.tsv"

# ---- truth: contig -> 'A'/'C' (pure only) ; clen for all contigs ----
truth = {}
clen = {}
total_assembly_bp = 0
with open(TRUTH) as f:
    next(f)
    for ln in f:
        c, L, a, cc, mf, kls = ln.rstrip("\n").split("\t")
        L = int(L); clen[c] = L; total_assembly_bp += L
        if kls == "pure_A": truth[c] = "A"
        elif kls == "pure_C": truth[c] = "C"

def chunk_to_contig(chunk):
    contig, start, end = chunk.rsplit("_", 2)
    return contig, int(end) - int(start)

files = sorted(glob.glob(os.path.join(CLUSTER_DIR, "subgenome_*.txt")))
if not files:
    sys.exit("no subgenome_*.txt in %s" % CLUSTER_DIR)

cluster_chunks = {}
contig_cluster_bp = defaultdict(lambda: defaultdict(int))
for fp in files:
    name = os.path.basename(fp).replace(".txt", "")
    rows = []
    for ln in open(fp):
        ch = ln.strip()
        if not ch: continue
        contig, bp = chunk_to_contig(ch)
        rows.append((contig, bp))
        contig_cluster_bp[contig][name] += bp
    cluster_chunks[name] = rows

# ---- (1) cluster-level oracle on clustered chunks ----
print("== clusters ==")
clustered_bp = oracle_correct = scored_bp = 0
cluster_label = {}
for name, rows in cluster_chunks.items():
    a = sum(bp for c, bp in rows if truth.get(c) == "A")
    cc = sum(bp for c, bp in rows if truth.get(c) == "C")
    allbp = sum(bp for c, bp in rows)
    clustered_bp += allbp
    maj = "A" if a >= cc else "C"
    cluster_label[name] = maj
    oracle_correct += max(a, cc); scored_bp += (a + cc)
    print(f"  {name:<14} chunks={len(rows):>5}  {allbp/1e6:6.1f} Mb   "
          f"trueA={a/1e6:6.1f}  trueC={cc/1e6:6.1f}  -> oracle label {maj}")

print()
print(f"(1) cluster-oracle accuracy (clustered, pure-contig bp): "
      f"{100*oracle_correct/max(1,scored_bp):.2f}%   ({scored_bp/1e6:.0f} Mb scored)")
print(f"(2) genome coverage clustered: {100*clustered_bp/total_assembly_bp:.1f}%   "
      f"({clustered_bp/1e6:.0f} of {total_assembly_bp/1e6:.0f} Mb)")

# ---- (3) contig-level oracle over the WHOLE assembly (unchunked = wrong) ----
contig_correct = contig_total = n_unassigned = unassigned_bp = 0
for c, t in truth.items():
    L = clen[c]; contig_total += L
    cb = contig_cluster_bp.get(c)
    if not cb:
        n_unassigned += 1; unassigned_bp += L; continue
    best = max(cb.items(), key=lambda kv: kv[1])[0]
    if cluster_label[best] == t:
        contig_correct += L
print()
print(f"(3) contig-oracle accuracy over ALL pure contigs (unchunked=wrong): "
      f"{100*contig_correct/max(1,contig_total):.2f}%   "
      f"({contig_total/1e6:.0f} Mb; {n_unassigned} contigs / {unassigned_bp/1e6:.0f} Mb unchunked)")
