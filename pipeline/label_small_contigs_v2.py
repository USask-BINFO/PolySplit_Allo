#!/usr/bin/env python3
import pickle, sys, os
from collections import defaultdict
import numpy as np

K, MC = 15, 2
HIC_FRAC = float(os.environ.get("POLYSPLIT_HIC_FRAC", "0.55"))
SAMPLE_PER_SG, CAP_BP = 30, 3_000_000
CONTACTS = os.environ.get("POLYSPLIT_CONTACTS", "contacts.pkl")
FASTA = os.environ.get("POLYSPLIT_FASTA", "../flye_out/assembly.fasta")
TRUTH = os.environ.get("POLYSPLIT_TRUTH", "../purity_compare/wg_purity.per_contig.tsv")
OUT = os.environ.get("POLYSPLIT_OUT", "all_contig_labels_v2.tsv")

MAP = np.full(256, -1, np.int64)
for ch, v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    MAP[ch] = v


def repeat_counts(seq):
    a = np.frombuffer(seq.encode(), np.uint8); code = MAP[a]
    if code.size < K:
        return {}
    sw = np.lib.stride_tricks.sliding_window_view(code, K)
    valid = (sw >= 0).all(axis=1)
    w = (np.int64(4) ** np.arange(K - 1, -1, -1)).astype(np.int64)
    fwd = (sw.astype(np.int64) @ w)[valid]
    if fwd.size == 0:
        return {}
    rc = np.zeros_like(fwd)
    for j in range(K):
        sym = (fwd >> np.int64(2 * j)) & np.int64(3)
        rc |= (np.int64(3) - sym) << np.int64(2 * (K - 1 - j))
    canon = np.minimum(fwd, rc)
    u, c = np.unique(canon, return_counts=True)
    keep = c >= MC
    return dict(zip(u[keep].tolist(), c[keep].tolist()))


# ---------- inputs ----------
with open(CONTACTS, "rb") as f:
    contacts, lengths = pickle.load(f)
lab = {}
LABELS_IN = sys.argv[1] if len(sys.argv) > 1 else "decloud_v2_contig_labels.tsv"
with open(LABELS_IN) as f:
    next(f)
    for ln in f:
        p = ln.rstrip("\n").split("\t"); lab[p[0]] = p[1]
truth = {}
if os.path.exists(TRUTH):                                  # truth is eval-only; tolerate absence
    with open(TRUTH) as f:
        next(f)
        for ln in f:
            c, L, a, cc_, mf, kls = ln.rstrip("\n").split("\t")
            if kls.startswith("pure_"): truth[c] = kls.split("_", 1)[1]

labels_present = sorted(set(lab.values()))

# ---------- (1) Hi-C linkage (argmax over K subgenomes) ----------
score = defaultdict(lambda: defaultdict(float))
for (a, b), w in contacts.items():
    la, lb = lengths.get(a, 0), lengths.get(b, 0)
    if la == 0 or lb == 0:
        continue
    nw = w / ((la * lb) ** 0.5)
    if a not in lab and b in lab:
        score[a][lab[b]] += nw
    elif b not in lab and a in lab:
        score[b][lab[a]] += nw

new, diffuse = {}, []
for c, s in score.items():
    tot = sum(s.values())
    if tot <= 0:
        diffuse.append(c); continue
    best, bw = max(s.items(), key=lambda kv: kv[1])
    if bw / tot >= HIC_FRAC:
        new[c] = best
    else:
        diffuse.append(c)
print(f"[hic] called {len(new):,}  diffuse->composition {len(diffuse):,}", file=sys.stderr)

# ---------- (2) composition fallback: nearest of K per-subgenome centroids ----------
big = {sg: sorted([c for c in lab if lab[c] == sg], key=lambda c: -lengths.get(c, 0))[:SAMPLE_PER_SG]
       for sg in labels_present}
sampled = set().union(*big.values()) if big else set()
want = set(diffuse) | sampled
seqs = {}; cn = None; buf = []
with open(FASTA) as f:
    for line in f:
        if line[0] == ">":
            if cn in want:
                seqs[cn] = "".join(buf)[:CAP_BP] if cn in sampled else "".join(buf)
            cn = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if cn in want:
        seqs[cn] = "".join(buf)

rc = {c: repeat_counts(seqs[c]) for c in want if c in seqs}
occ = defaultdict(int)
for c in sampled:
    for km in rc.get(c, {}):
        occ[km] += 1
feats = [km for km, o in occ.items() if o >= 8]
fidx = {km: i for i, km in enumerate(feats)}


def vec(c):
    v = np.zeros(len(feats)); d = rc.get(c, {}); tot = sum(d.values()) or 1
    for km, n in d.items():
        if km in fidx:
            v[fidx[km]] = n / tot
    return v


centroids = {sg: np.mean([vec(c) for c in big[sg] if c in rc] or [np.zeros(len(feats))], axis=0)
             for sg in labels_present}
comp = {}
if feats:
    for c in diffuse:
        if c not in rc or not rc[c]:
            continue
        v = vec(c)
        comp[c] = min(labels_present, key=lambda sg: np.linalg.norm(v - centroids[sg]))
print(f"[comp] called {len(comp):,} of {len(diffuse):,} diffuse", file=sys.stderr)


def acc(d, tag):
    ok = bad = okbp = badbp = 0
    for c, l in d.items():
        if c in truth:
            if l == truth[c]: ok += 1; okbp += lengths.get(c, 0)
            else: bad += 1; badbp += lengths.get(c, 0)
    if ok + bad:
        print(f"  [{tag}] contigs {ok+bad:,}  acc {100*ok/(ok+bad):.1f}%  "
              f"bp {100*okbp/max(1,okbp+badbp):.1f}%", file=sys.stderr)


acc(new, "hic@%.2f" % HIC_FRAC)
acc(comp, "composition")

with open(OUT, "w") as o:
    o.write("contig\tlabel\tsource\n")
    for c, l in lab.items(): o.write(f"{c}\t{l}\tblock\n")
    for c, l in new.items(): o.write(f"{c}\t{l}\thic\n")
    for c, l in comp.items():
        if c not in new: o.write(f"{c}\t{l}\tcomp\n")
print(f"-> {OUT} (block {len(lab)} + hic {len(new)} + comp {len(comp)})", file=sys.stderr)
