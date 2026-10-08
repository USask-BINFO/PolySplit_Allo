#!/usr/bin/env python3
import sys, os, glob, bisect
from collections import defaultdict

CLUSTER_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "pc_work/analysisOutputs_NAM29_AC_scaf_n2/SpectralClusteringmain_tsne_2_n3/clusterResults"
TRUTH = "$DATA/Napus_nam_denovo/sa_read_correct/purity_compare/wg_purity.per_contig.tsv"
AGP = "$DATA/Napus_nam_denovo/sa_read_correct/hic_ac_cluster/yahs/napus_scaffolds_final.agp"

# ---- truth: contig -> 'A'/'C' (pure only); total assembly bp ----
truth = {}
total_assembly_bp = 0
with open(TRUTH) as f:
    next(f)
    for ln in f:
        c, L, a, cc, mf, kls = ln.rstrip("\n").split("\t")
        total_assembly_bp += int(L)
        if kls == "pure_A": truth[c] = "A"
        elif kls == "pure_C": truth[c] = "C"

# ---- AGP: scaffold -> sorted [(scaf_start, scaf_end, contig)] for W (contig) rows ----
scaf_iv = defaultdict(list)
with open(AGP) as f:
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 9 or p[4] != "W":
            continue
        scaf_iv[p[0]].append((int(p[1]), int(p[2]), p[5]))
for s in scaf_iv:
    scaf_iv[s].sort()
scaf_starts = {s: [iv[0] for iv in ivs] for s, ivs in scaf_iv.items()}

def chunk_to_contig(chunk):
    scaf, start, end = chunk.rsplit("_", 2)
    cs, ce = int(start), int(end)          # 0-based-ish; 1bp offset vs AGP is negligible
    bp = ce - cs
    ivs = scaf_iv.get(scaf)
    if not ivs:
        return None, bp
    # scan intervals overlapping [cs, ce]; start a little left of cs
    i = max(0, bisect.bisect_right(scaf_starts[scaf], cs) - 2)
    ov = defaultdict(int)
    for j in range(i, len(ivs)):
        a, b, contig = ivs[j]
        if a > ce:
            break
        o = min(ce, b) - max(cs, a)
        if o > 0:
            ov[contig] += o
    if not ov:
        return None, bp
    best = max(ov.items(), key=lambda kv: kv[1])[0]
    return best, bp

files = sorted(glob.glob(os.path.join(CLUSTER_DIR, "subgenome_*.txt")))
if not files:
    sys.exit("no subgenome_*.txt in %s" % CLUSTER_DIR)

cluster_chunks = {}
for fp in files:
    name = os.path.basename(fp).replace(".txt", "")
    rows = []
    for ln in open(fp):
        ch = ln.strip()
        if not ch: continue
        contig, bp = chunk_to_contig(ch)
        rows.append((contig, bp))
    cluster_chunks[name] = rows

print("== clusters (scaffold run) ==")
clustered_bp = oracle_correct = scored_bp = 0
for name, rows in cluster_chunks.items():
    a = sum(bp for c, bp in rows if truth.get(c) == "A")
    cc = sum(bp for c, bp in rows if truth.get(c) == "C")
    allbp = sum(bp for c, bp in rows)
    clustered_bp += allbp
    maj = "A" if a >= cc else "C"
    oracle_correct += max(a, cc); scored_bp += (a + cc)
    print(f"  {name:<14} chunks={len(rows):>5}  {allbp/1e6:6.1f} Mb   "
          f"trueA={a/1e6:6.1f}  trueC={cc/1e6:6.1f}  -> oracle label {maj}")

print()
print(f"(1) cluster-oracle accuracy (clustered, pure-contig bp): "
      f"{100*oracle_correct/max(1,scored_bp):.2f}%   ({scored_bp/1e6:.0f} Mb scored)")
print(f"(2) genome coverage clustered: {100*clustered_bp/total_assembly_bp:.1f}%   "
      f"({clustered_bp/1e6:.0f} of {total_assembly_bp/1e6:.0f} Mb)")
