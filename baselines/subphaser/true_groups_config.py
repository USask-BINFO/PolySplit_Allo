#!/usr/bin/env python3
import sys, glob
from collections import defaultdict, Counter
import numpy as np
from scipy.optimize import linear_sum_assignment

pre, c2r, cs_path, wdir, out = sys.argv[1:6]
MINLEN = int(sys.argv[6]) if len(sys.argv) > 6 else 5_000_000
cs = {p[0]: p[1] for p in (ln.split() for ln in open(cs_path)) if len(p) >= 2}
subs = sorted(set(cs.values())); chroms = {s: sorted(c for c in cs if cs[c] == s) for s in subs}
cnt = defaultdict(Counter)
for paf in glob.glob(f"{wdir}/cross.*.paf"):
    seen = set()
    for ln in open(paf):
        p = ln.split("\t")
        if (p[0], p[5]) not in seen:
            seen.add((p[0], p[5])); cnt[p[0].rsplit(":", 1)[0]][p[5]] += 1


def match(sa, sb):
    A, B = chroms[sa], chroms[sb]
    M = np.array([[cnt[a][b] + cnt[b][a] for b in B] for a in A], float)
    r, c = linear_sum_assignment(-M)
    return {A[i]: B[j] for i, j in zip(r, c) if M[i, j] > 0}


maps = {s: match(subs[0], s) for s in subs[1:]}
groups = [g for g in ([a] + [maps[s].get(a) for s in subs[1:]] for a in chroms[subs[0]]) if all(g)]
cc = defaultdict(Counter)
for ln in open(c2r):
    p = ln.split("\t"); cc[p[0]][p[5]] += int(p[9])
ctop = {c: d.most_common(1)[0][0] for c, d in cc.items() if d}
slen = {l.split()[0]: int(l.split()[1]) for l in open(f"{pre}.fa.fai")}
schr = defaultdict(Counter)
for ln in open(f"{pre}.agp"):
    x = ln.rstrip("\n").split("\t")
    if not ln.startswith("#") and x[4] == "W" and x[5] in ctop:
        schr[x[0]][ctop[x[5]]] += int(x[7]) - int(x[6]) + 1
best = {}
for s, d in schr.items():
    if slen.get(s, 0) >= MINLEN and d:
        ch = d.most_common(1)[0][0]
        if ch not in best or slen[s] > slen[best[ch]]:
            best[ch] = s
with open(out, "w") as o:
    for g in groups:
        r = [best.get(c) for c in g]
        if all(r) and len(set(r)) == len(r):
            o.write("\t".join(r) + "\n")
print(f"homoeologous chromosome groups: {groups}", file=sys.stderr)
