#!/usr/bin/env python3
import sys
from collections import defaultdict

FAI = sys.argv[1]; PAIRS = sys.argv[2]; OUT = sys.argv[3]
MINLEN = int(sys.argv[4]) if len(sys.argv) > 4 else 5_000_000
NSG = int(sys.argv[5]) if len(sys.argv) > 5 else 2

length = {}
for ln in open(FAI):
    p = ln.split("\t"); length[p[0]] = int(p[1])
big = sorted([s for s, L in length.items() if L >= MINLEN], key=lambda s: -length[s])
bigset = set(big)
sys.stderr.write(f"[config] {len(big)} scaffolds >= {MINLEN/1e6:.0f} Mb; NSG={NSG}\n")

score = defaultdict(float)
for ln in open(PAIRS):
    p = ln.rstrip("\n").split("\t")
    if len(p) < 3:
        continue
    s = float(p[0]); a, b = p[1], p[2]
    if a in bigset and b in bigset and a != b:
        score[(a, b)] = score[(b, a)] = max(score[(a, b)], s)

used = set(); groups = []
for a in big:                                            # largest scaffolds anchor groups first
    if a in used:
        continue
    partners = sorted([x for x in big if x not in used and x != a and score.get((a, x), 0) > 0],
                      key=lambda x: -score[(a, x)])
    if len(partners) >= NSG - 1:
        grp = [a] + partners[:NSG - 1]
        groups.append(grp); used.update(grp)

# leftover scaffolds with no scored partners: pack into NSG-rows so SubPhaser still sees them
leftover = [s for s in big if s not in used]
for i in range(0, len(leftover) - (len(leftover) % NSG), NSG):
    groups.append(leftover[i:i + NSG]); used.update(leftover[i:i + NSG])
dropped = [s for s in big if s not in used]

with open(OUT, "w") as o:
    for g in groups:
        o.write("\t".join(g) + "\n")
sys.stderr.write(f"[config] wrote {OUT}: {len(groups)} groups of {NSG}\n")
for g in groups:
    sys.stderr.write("  " + "  ".join(g) + "\n")
if dropped:
    sys.stderr.write(f"[config] WARNING leftover dropped (< NSG): {dropped}\n")
