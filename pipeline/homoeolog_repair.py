#!/usr/bin/env python3
import sys, os, pickle
from collections import defaultdict

NONLAB = {"", "ambiguous", "unassigned", "?", "chimeric", "label", "NA"}
CONTACTS = os.environ.get("POLYSPLIT_CONTACTS", "contacts.pkl")
TRUTH = os.environ.get("POLYSPLIT_TRUTH", "../purity_compare/wg_purity.per_contig.tsv")
EDGES = os.environ.get("POLYSPLIT_EDGES", "homoeolog_edges.tsv")
LABELS_IN = os.environ.get("POLYSPLIT_LABELS_IN", "all_contig_labels_v2.tsv")
OUT = os.environ.get("POLYSPLIT_OUT", "all_contig_labels_repaired.tsv")
EDGE_MIN = int(os.environ.get("POLYSPLIT_PAIR_MIN", "100000"))   # strong homoeolog twins only (= tau_H)
MARGIN = 1.0        # flip iff same-label twin weight > other-label twin weight

with open(CONTACTS, "rb") as f:
    _, lengths = pickle.load(f)
ctruth = {}
if os.path.exists(TRUTH):                                   # truth is eval-only; tolerate absence
    with open(TRUTH) as f:
        next(f)
        for ln in f:
            c, L, a, cc_, mf, kls = ln.rstrip("\n").split("\t")
            if kls.startswith("pure_"): ctruth[c] = kls.split("_", 1)[1]


def load_labels():
    lab = {}
    with open(LABELS_IN) as f:
        next(f)
        for ln in f:
            p = ln.split("\t")
            if len(p) >= 2 and p[1] not in NONLAB:
                lab[p[0]] = p[1]
    return lab


def acc(lab, tag):
    okbp = badbp = 0
    for c, l in lab.items():
        if c in ctruth:
            if l == ctruth[c]: okbp += lengths.get(c, 0)
            else: badbp += lengths.get(c, 0)
    if okbp + badbp:
        print(f"  [{tag}] contig bp acc {100*okbp/(okbp+badbp):.2f}%")


def repair(edge_min, margin=1.0):
    lab = load_labels()
    labels_present = sorted(set(lab.values()))
    # same strongest-partner filter as in the labelling step
    FRAC = float(os.environ.get("POLYSPLIT_EDGE_FRAC", "0.5"))
    raw = []; mx = {}
    with open(EDGES) as f:
        next(f)
        for ln in f:
            a, b, w = ln.rstrip("\n").split("\t"); w = int(w)
            raw.append((a, b, w))
            if w > mx.get(a, 0): mx[a] = w
            if w > mx.get(b, 0): mx[b] = w
    # pair_mode.txt = mutual: use only edges among the K-1 strongest of both contigs
    nk = max(1, int(os.environ.get("POLYSPLIT_NSG", str(len(labels_present)))) - 1)
    ws = defaultdict(list)
    for a, b, w in raw:
        ws[a].append(w); ws[b].append(w)
    kth = {c: sorted(v, reverse=True)[:nk][-1] for c, v in ws.items()}
    mode = os.environ.get("POLYSPLIT_MUTUAL_PAIRS", "auto")
    pm = os.environ.get("POLYSPLIT_PAIR_MODE", "pair_mode.txt")
    MUTUAL = mode == "1" or (mode == "auto" and os.path.exists(pm) and open(pm).read().strip() == "mutual")
    nbrs = defaultdict(list)
    for a, b, w in raw:
        if MUTUAL and not (w >= kth[a] and w >= kth[b]):
            continue
        if w >= edge_min and a in lab and b in lab and w >= FRAC * mx[a] and w >= FRAC * mx[b]:
            nbrs[a].append((b, w)); nbrs[b].append((a, w))
    nflip = 0
    for _ in range(20):
        changed = False
        for c in list(nbrs):
            tw = defaultdict(float)                          # strong-twin weight per subgenome
            for p, w in nbrs[c]:
                tw[lab[p]] += w
            cur = lab[c]
            same = tw.get(cur, 0.0)
            other = sum(v for l, v in tw.items() if l != cur)
            if same > other * margin and same > 0 and len(labels_present) > 1:
                target = min(labels_present, key=lambda l: tw.get(l, 0.0))   # twins occupy this least
                if target != cur:
                    lab[c] = target; nflip += 1; changed = True
        if not changed:
            break
    return lab, nflip


base = load_labels()
lab, nf = repair(EDGE_MIN, MARGIN)
with open(OUT, "w") as o:
    o.write("contig\tlabel\tsource\n")
    for c, l in lab.items():
        o.write(f"{c}\t{l}\trepair\n")
print(f"repair: FIXED edge_min={EDGE_MIN//1000}k margin={MARGIN} (no truth tuning) -> {nf} flips",
      file=sys.stderr)
print(f"-> {OUT}", file=sys.stderr)
print("[eval only] contig accuracy (bp) vs truth:")
acc(base, "before repair")
acc(lab, "after repair")
