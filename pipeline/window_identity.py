#!/usr/bin/env python3
import sys, glob
from collections import defaultdict

D = sys.argv[1]
hits = defaultdict(list)
for paf in glob.glob(f"{D}/cross.*.paf"):
    for ln in open(paf):
        p = ln.split("\t")
        de = next((float(x[5:]) for x in p[12:] if x.startswith("de:f:")), None)
        if de is not None:
            hits[p[0]].append((int(p[2]), int(p[3]), 1 - de))
with open(f"{D}/windows.tsv", "w") as out:
    out.write("window\tsubg\tlen\tcov\tident\n")
    for fa in sorted(glob.glob(f"{D}/win.*.fa")):
        s = fa.split("win.")[-1][:-3]
        for ln in open(fa):
            if not ln.startswith(">"):
                continue
            w = ln[1:].strip()
            a, b = map(int, w.rsplit(":", 1)[1].split("-")); L = b - a
            cov = 0; cs = ce = -1; wsum = isum = 0.0
            for qs, qe, idt in sorted(hits.get(w, [])):
                wsum += qe - qs; isum += (qe - qs) * idt
                if qs > ce:
                    cov += max(0, ce - cs); cs, ce = qs, qe
                else:
                    ce = max(ce, qe)
            cov += max(0, ce - cs)
            out.write(f"{w}\t{s}\t{L}\t{cov / L:.3f}\t{isum / wsum if wsum else float('nan'):.4f}\n")
