#!/usr/bin/env python3
import sys, gzip, os

ref, cs, out = sys.argv[1], sys.argv[2], sys.argv[3]
W = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
subg = {p[0]: p[1] for p in (ln.split() for ln in open(cs)) if len(p) >= 2}
os.makedirs(out, exist_ok=True)
fh_c, fh_w = {}, {}


def flush(name, seq):
    s = subg.get(name)
    if s is None:
        return
    seq = "".join(seq)
    if s not in fh_c:
        fh_c[s] = open(f"{out}/chroms.{s}.fa", "w"); fh_w[s] = open(f"{out}/win.{s}.fa", "w")
    fh_c[s].write(f">{name}\n")
    for i in range(0, len(seq), 100):
        fh_c[s].write(seq[i:i + 100] + "\n")
    for st in range(0, len(seq), W):
        w = seq[st:st + W]
        if len(w) >= W // 2 and w.upper().count("N") <= 0.1 * len(w):
            fh_w[s].write(f">{name}:{st}-{st + len(w)}\n{w}\n")


name, seq = None, []
with (gzip.open if ref.endswith(".gz") else open)(ref, "rt") as f:
    for ln in f:
        if ln.startswith(">"):
            if name is not None:
                flush(name, seq)
            name, seq = ln[1:].split()[0], []
        else:
            seq.append(ln.strip())
if name is not None:
    flush(name, seq)
for h in list(fh_c.values()) + list(fh_w.values()):
    h.close()
