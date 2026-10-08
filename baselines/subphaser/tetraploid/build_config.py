#!/usr/bin/env python3
import sys
from collections import defaultdict

FAI   = sys.argv[1]
PAIRS = sys.argv[2]
OUT   = sys.argv[3]
MINLEN = int(sys.argv[4]) if len(sys.argv) > 4 else 5_000_000
NSG    = int(sys.argv[5]) if len(sys.argv) > 5 else 2

length = {}
for ln in open(FAI):
    p = ln.split("\t"); length[p[0]] = int(p[1])
big = sorted([s for s, L in length.items() if L >= MINLEN], key=lambda s: -length[s])
bigset = set(big)
sys.stderr.write(f"[config] {len(big)} scaffolds >= {MINLEN/1e6:.0f} Mb\n")
if NSG != 2:
    sys.exit("this builder writes 2-col (nsg=2) configs only; for nsg=3 group into triples.")

cand = []
for ln in open(PAIRS):
    p = ln.rstrip("\n").split("\t")
    if len(p) < 3:
        continue
    score = float(p[0]); a, b = p[1], p[2]
    if a in bigset and b in bigset and a != b:
        cand.append((score, a, b))
cand.sort(reverse=True)

used = set(); pairs = []
for score, a, b in cand:
    if a in used or b in used:
        continue
    pairs.append((a, b)); used.add(a); used.add(b)

# pair any leftover big scaffolds with each other so all are included (even-count rows only)
unmatched = [s for s in big if s not in used]
for i in range(0, len(unmatched) - 1, 2):
    pairs.append((unmatched[i], unmatched[i + 1])); used.update(unmatched[i:i+2])
leftover = [s for s in big if s not in used]

with open(OUT, "w") as o:
    for a, b in pairs:
        o.write(f"{a}\t{b}\n")
sys.stderr.write(f"[config] wrote {OUT}: {len(pairs)} pairs\n")
for a, b in pairs:
    sys.stderr.write(f"  {a:<14} {b}\n")
if leftover:
    sys.stderr.write(f"[config] WARNING leftover (odd count, dropped): {leftover}\n")
