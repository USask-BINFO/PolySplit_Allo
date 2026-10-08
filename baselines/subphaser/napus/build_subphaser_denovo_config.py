#!/usr/bin/env python3
import sys
from collections import defaultdict

FAI = "$DATA/Napus_nam_denovo/sa_read_correct/hic_ac_cluster/yahs/napus_scaffolds_final.fa.fai"
PAIRS = "$DATA/Napus_nam_denovo/sa_read_correct/hic_ac_cluster/baselines/scaf_pairs.tsv"
OUT = "$DATA/Napus_nam_denovo/sa_read_correct/hic_ac_cluster/baselines/denovo_pairs_sg.config"
MINLEN = 10_000_000

length = {}
for ln in open(FAI):
    p = ln.split("\t"); length[p[0]] = int(p[1])
big = sorted([s for s, L in length.items() if L >= MINLEN], key=lambda s: -length[s])
bigset = set(big)
print(f"[config] {len(big)} scaffolds >= {MINLEN/1e6:.0f} Mb", file=sys.stderr)

# greedy max-score matching restricted to big scaffolds
cand = []
for ln in open(PAIRS):
    p = ln.rstrip("\n").split("\t")
    if len(p) < 3: continue
    score = float(p[0]); a, b = p[1], p[2]
    if a in bigset and b in bigset and a != b:
        cand.append((score, a, b))
cand.sort(reverse=True)

used = set(); pairs = []
for score, a, b in cand:
    if a in used or b in used: continue
    pairs.append((a, b)); used.add(a); used.add(b)

unmatched = [s for s in big if s not in used]
# pair leftover scaffolds with each other (no self-aln link) so all are included
for i in range(0, len(unmatched) - 1, 2):
    pairs.append((unmatched[i], unmatched[i+1]))
    used.add(unmatched[i]); used.add(unmatched[i+1])
leftover = [s for s in big if s not in used]

with open(OUT, "w") as o:
    for a, b in pairs:
        o.write(f"{a}\t{b}\n")
print(f"[config] wrote {OUT}: {len(pairs)} pairs", file=sys.stderr)
if leftover:
    print(f"[config] WARNING leftover (odd count, dropped from 2-col config): {leftover}", file=sys.stderr)
for a, b in pairs:
    print(f"  {a:<12} {b}", file=sys.stderr)
