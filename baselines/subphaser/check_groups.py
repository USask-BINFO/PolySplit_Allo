#!/usr/bin/env python3
import sys, json
from collections import defaultdict, Counter

cfg, agp, pur, c2r = sys.argv[1:5]
sgout = sys.argv[5] if len(sys.argv) > 5 else None
truth = {}
for ln in open(pur):
    x = ln.rstrip("\n").split("\t")
    if x[0] != "contig" and x[5].startswith("pure_"):
        truth[x[0]] = x[5][5:]
chrom = defaultdict(Counter)
for ln in open(c2r):
    p = ln.split("\t"); chrom[p[0]][p[5]] += int(p[9])
sub, chr_ = defaultdict(Counter), defaultdict(Counter)
for ln in open(agp):
    x = ln.rstrip("\n").split("\t")
    if ln.startswith("#") or x[4] != "W":
        continue
    L = int(x[7]) - int(x[6]) + 1
    if x[5] in truth:
        sub[x[0]][truth[x[5]]] += L
    if chrom[x[5]]:
        chr_[x[0]][chrom[x[5]].most_common(1)[0][0]] += L
called = {}
if sgout:
    for ln in open(sgout):
        if not ln.startswith("#"):
            x = ln.split(); called[x[0]] = x[1]
groups = []
for ln in open(cfg):
    g = ln.split()
    if not g:
        continue
    subs = [sub[s].most_common(1)[0][0] if sub[s] else "?" for s in g]
    groups.append(dict(scaffolds=g, subgenomes=subs, chromosomes=[chr_[s].most_common(1)[0][0] if chr_[s] else "?" for s in g],
                       subphaser=[called.get(s, "-") for s in g], one_per_subgenome=len(set(subs)) == len(subs)))
print(json.dumps(dict(correct=sum(g["one_per_subgenome"] for g in groups), groups=len(groups), detail=groups), indent=1))
