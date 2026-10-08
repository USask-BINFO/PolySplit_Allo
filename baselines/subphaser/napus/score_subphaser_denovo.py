#!/usr/bin/env python3
import sys, re
from collections import defaultdict

if len(sys.argv) < 2:
    sys.exit("usage: score_subphaser_denovo.py <*.chrom-subgenome.tsv> [out_labels.tsv]")
SGOUT = sys.argv[1]
OUTLAB = sys.argv[2] if len(sys.argv) > 2 else "subphaser_contig_labels.tsv"
TRUTH = "$DATA/Napus_nam_denovo/sa_read_correct/purity_compare/wg_purity.per_contig.tsv"
AGP = "$DATA/Napus_nam_denovo/sa_read_correct/hic_ac_cluster/yahs/napus_scaffolds_final.agp"

# truth: contig -> A/C
ctruth = {}
for ln in open(TRUTH):
    if ln.startswith("contig\t"): continue
    c, L, a, cc, mf, kls = ln.rstrip("\n").split("\t")
    if kls == "pure_A": ctruth[c] = "A"
    elif kls == "pure_C": ctruth[c] = "C"

# AGP: scaffold -> [(contig, bp)] ; and contig -> scaffold
scaf_contigs = defaultdict(list)
for ln in open(AGP):
    p = ln.rstrip("\n").split("\t")
    if len(p) < 9 or p[4] != "W": continue
    scaf_contigs[p[0]].append((p[5], int(p[2]) - int(p[1]) + 1))

def norm(chrom):
    m = re.search(r"scaffold_\d+", chrom)
    return m.group(0) if m else chrom

# SubPhaser scaffold -> SG
scaf_sg = {}
for ln in open(SGOUT):
    if ln.startswith("#") or not ln.strip(): continue
    p = ln.split()
    if len(p) < 2: continue
    scaf_sg[norm(p[0])] = p[1]

# oracle: SG -> A/C by true-majority bp of member scaffolds' contigs
sg_ac = defaultdict(lambda: [0, 0])
for s, sg in scaf_sg.items():
    for contig, bp in scaf_contigs.get(s, []):
        t = ctruth.get(contig)
        if t == "A": sg_ac[sg][0] += bp
        elif t == "C": sg_ac[sg][1] += bp
# best 1:1 assignment of the two clusters to {A,C} (standard clustering-accuracy oracle;
# avoids both clusters collapsing onto the same subgenome) -- same rule used for polyCRACKER
from itertools import permutations
sgs = list(sg_ac.keys())
sg_label = max(
    (dict(zip(sgs, perm)) for perm in permutations(["A", "C"], len(sgs))),
    key=lambda m: sum(sg_ac[s][0] if m[s] == "A" else sg_ac[s][1] for s in sgs))
print("== SubPhaser clusters (best 1:1 assignment) ==")
for sg, (a, c) in sg_ac.items():
    print(f"  {sg}: trueA={a/1e6:.1f} Mb  trueC={c/1e6:.1f} Mb  -> {sg_label[sg]}  "
          f"(scaffolds: {sorted(s for s,g in scaf_sg.items() if g==sg)})")

# (1) contig accuracy + (2) emit contig labels
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
                if t == lab: correct += bp
print()
print(f"(1) SubPhaser CONTIG accuracy (scaffold SG->oracle A/C, applied to contigs): "
      f"{100*correct/max(1,total):.2f}%   ({total/1e6:.0f} Mb scored over {len(scaf_sg)} scaffolds)")
print(f"    wrote {OUTLAB} -> feed to propagate_to_reads.py for the matched READ-level number.")
