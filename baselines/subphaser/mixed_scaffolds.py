#!/usr/bin/env python3
import sys, json
from collections import defaultdict

agp, pur = sys.argv[1], sys.argv[2]
truth = {}
for ln in open(pur):
    x = ln.rstrip("\n").split("\t")
    if x[0] != "contig" and x[5].startswith("pure_"):
        truth[x[0]] = x[5][5:]
comp = defaultdict(lambda: defaultdict(int)); size = defaultdict(int)
for ln in open(agp):
    if ln.startswith("#"):
        continue
    x = ln.rstrip("\n").split("\t")
    if x[4] == "W":
        L = int(x[7]) - int(x[6]) + 1; size[x[0]] += L
        if x[5] in truth:
            comp[x[0]][truth[x[5]]] += L
    else:
        size[x[0]] += int(x[2]) - int(x[1]) + 1
rows = []
for s in sorted((s for s in size if size[s] >= 10_000_000), key=lambda s: -size[s]):
    d = comp[s]; tot = sum(d.values())
    rows.append((s, size[s], 1 - max(d.values()) / tot if tot else 0.0, dict(d)))
mixed = [r for r in rows if r[2] >= 0.05]
anch = sum(sum(r[3].values()) for r in rows)
print(json.dumps(dict(scaffolds_ge10Mb=len(rows), mixed_ge5pct=len(mixed),
                      mixed_bp_Mb=round(sum(r[1] for r in mixed) / 1e6, 1), big_bp_Mb=round(sum(r[1] for r in rows) / 1e6, 1),
                      minority_bp_frac_in_big=round(sum(r[2] * sum(r[3].values()) for r in rows) / anch, 3) if anch else None,
                      mixed=[(r[0], round(r[1] / 1e6, 1), round(100 * r[2], 1)) for r in mixed]), indent=1))
