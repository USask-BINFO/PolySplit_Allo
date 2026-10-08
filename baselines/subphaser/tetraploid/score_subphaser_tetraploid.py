#!/usr/bin/env python3
import sys, re
from collections import defaultdict
from itertools import permutations

SGOUT  = sys.argv[1]
AGP    = sys.argv[2]
OUTLAB = sys.argv[3] if len(sys.argv) > 3 else "subphaser_contig_labels.tsv"
TRUTH  = sys.argv[4] if len(sys.argv) > 4 else \
    "$DATA/camelina/Tetraploid_data/polysplit_run/wg_purity.per_contig.tsv"

# truth: contig -> S1/S2 (pure only)
ctruth = {}
for ln in open(TRUTH):
    if ln.startswith("contig\t"):
        continue
    c, L, a, cc, mf, kls = ln.rstrip("\n").split("\t")
    if kls.startswith("pure_"):
        ctruth[c] = kls.split("_", 1)[1]

# AGP: scaffold -> [(contig, bp)]
scaf_contigs = defaultdict(list)
for ln in open(AGP):
    p = ln.rstrip("\n").split("\t")
    if len(p) < 9 or p[4] != "W":
        continue
    scaf_contigs[p[0]].append((p[5], int(p[2]) - int(p[1]) + 1))

def norm(chrom):
    m = re.search(r"scaffold_\d+", chrom)
    return m.group(0) if m else chrom

# SubPhaser scaffold -> SG
scaf_sg = {}
for ln in open(SGOUT):
    if ln.startswith("#") or not ln.strip():
        continue
    p = ln.split()
    if len(p) < 2:
        continue
    scaf_sg[norm(p[0])] = p[1]

# oracle: best 1:1 SG-cluster -> subgenome
sg_ac = defaultdict(lambda: defaultdict(int))
for s, sg in scaf_sg.items():
    for contig, bp in scaf_contigs.get(s, []):
        t = ctruth.get(contig)
        if t:
            sg_ac[sg][t] += bp
sgs = list(sg_ac.keys())
truth_labels = sorted({t for d in sg_ac.values() for t in d})
sg_label = max(
    (dict(zip(sgs, perm)) for perm in permutations(truth_labels, min(len(sgs), len(truth_labels)))),
    key=lambda m: sum(sg_ac[s].get(m.get(s), 0) for s in sgs))
print("== SubPhaser clusters (best 1:1 assignment) ==")
for sg, d in sg_ac.items():
    comp = "  ".join(f"{t}={d[t]/1e6:.1f}Mb" for t in sorted(d))
    print(f"  {sg}: {comp}  -> {sg_label.get(sg,'?')}  "
          f"(scaffolds: {sorted(s for s,g in scaf_sg.items() if g==sg)})")

correct = total = 0
with open(OUTLAB, "w") as o:
    o.write("contig\tlabel\ttruth\n")
    for s, sg in scaf_sg.items():
        lab = sg_label[sg]
        for contig, bp in scaf_contigs.get(s, []):
            o.write(f"{contig}\t{lab}\t{ctruth.get(contig,'?')}\n")
            t = ctruth.get(contig)
            if t:
                total += bp
                if t == lab:
                    correct += bp
print()
print(f"(1) SubPhaser CONTIG accuracy (scaffold SG->oracle S1/S2 applied to contigs): "
      f"{100*correct/max(1,total):.2f}%   ({total/1e6:.0f} Mb over {len(scaf_sg)} scaffolds)")
print(f"    wrote {OUTLAB} -> propagate_to_reads.py for the read-level number.")
