#!/usr/bin/env python3
import sys, os, pickle
from collections import defaultdict, Counter
import numpy as np
import networkx as nx
from networkx.algorithms.community import louvain_communities

MAP = np.full(256, -1, np.int64)
for ch, v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    MAP[ch] = v
MIN_LEN, MIN_BLOCK, MIN_SUB = 50000, 3000000, 1000000
K, MC_CONTIG, MC_BLOCK = 15, 2, 10
RES, SUBRES = 2.5, 5.0
PAIR_MIN = int(os.environ.get("POLYSPLIT_PAIR_MIN", "100000"))
MAX_ITERS = 4
SWEEP = [(100000, 1), (50000, 2), (30000, 3), (20000, 3), (20000, 5)]


def contig_repeat_counts(seq, k, mc):
    a = np.frombuffer(seq.encode(), np.uint8); code = MAP[a]
    if code.size < k:
        return {}
    sw = np.lib.stride_tricks.sliding_window_view(code, k)
    valid = (sw >= 0).all(axis=1)
    w = (np.int64(4) ** np.arange(k - 1, -1, -1)).astype(np.int64)
    fwd = (sw.astype(np.int64) @ w)[valid]
    if fwd.size == 0:
        return {}
    rc = np.zeros_like(fwd)
    for j in range(k):
        sym = (fwd >> np.int64(2 * j)) & np.int64(3)
        rc |= (np.int64(3) - sym) << np.int64(2 * (k - 1 - j))
    canon = np.minimum(fwd, rc)
    u, c = np.unique(canon, return_counts=True)
    keep = c >= mc
    return dict(zip(u[keep].tolist(), c[keep].tolist()))


def enforce_pairs(score, pairs, nb, iters=500):
    s = score - np.median(score); abss = np.abs(s)
    lab = (s >= 0).astype(int)
    w_c = np.median(abss) if abss.size else 1.0
    nbrs = defaultdict(list)
    for i, j in pairs:
        nbrs[i].append(j); nbrs[j].append(i)
    for _ in range(iters):
        ch = False
        for i in range(nb):
            aligned = (lab[i] == 1) == (s[i] >= 0)
            d = -abss[i] if aligned else abss[i]
            sat = sum(1 for j in nbrs[i] if lab[i] != lab[j])
            sat_f = sum(1 for j in nbrs[i] if (1 - lab[i]) != lab[j])
            if d + w_c * (sat_f - sat) > 1e-9:
                lab[i] = 1 - lab[i]; ch = True
        if not ch:
            break
    return lab


def k_colourable(edges, k, max_steps=2000000):
    adj = defaultdict(set)
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    order = sorted(adj, key=lambda n: -len(adj[n])); col = {}; steps = [0]

    def go(i):
        if i == len(order): return True
        steps[0] += 1
        if steps[0] > max_steps: return True
        n = order[i]; used = {col[m] for m in adj[n] if m in col}
        for c in range(k):
            if c in used: continue
            col[n] = c
            if go(i + 1): return True
            del col[n]
            if c == len(set(col.values())): break
        return False
    sys.setrecursionlimit(max(10000, len(order) + 100))
    return go(0)


def mutual_top_pairs(bpw, k):
    part = defaultdict(list)
    for (i, j), w in bpw.items():
        part[i].append((w, j)); part[j].append((w, i))
    topk = {i: {j for w, j in sorted(v, reverse=True)[:k - 1]} for i, v in part.items()}
    return [(i, j) for (i, j) in bpw if j in topk[i] and i in topk[j]]


def dechimerize(blocks, G, hedges, lengths, edge_min, min_edges):
    nsplit = 0
    for _ in range(MAX_ITERS):
        blk_of = {c: i for i, com in enumerate(blocks) for c in com}
        cnt = Counter()
        for (a, b), w in hedges.items():
            if w >= edge_min and a in blk_of and b in blk_of and blk_of[a] == blk_of[b]:
                cnt[blk_of[a]] += 1
        chim = {i for i, c in cnt.items() if c >= min_edges}
        if not chim:
            break
        new = []
        for i, com in enumerate(blocks):
            if i in chim:
                nsplit += 1
                sub = G.subgraph(com)
                subs = [list(c) for c in louvain_communities(sub, weight="weight",
                                                             resolution=SUBRES, seed=1)]
                kept = [s for s in subs if sum(lengths[c] for c in s) >= MIN_SUB]
                if len(kept) >= 2:
                    new.extend(kept)
                    rem = [c for s in subs if s not in kept for c in s]
                    if rem:
                        new[-1].extend(rem)
                else:
                    new.append(com)
            else:
                new.append(com)
        if len(new) == len(blocks):
            break
        blocks = new
    return blocks, nsplit


def label_eval(blocks, cc, lengths, top, truth, tag, nsg=2):
    blk_of = {c: i for i, com in enumerate(blocks) for c in com}
    nb = len(blocks)
    bcnt = [defaultdict(int) for _ in range(nb)]
    for c in blk_of:
        for code, n in cc.get(c, {}).items():
            bcnt[blk_of[c]][code] += n
    occ = Counter()
    for d in bcnt:
        for code, n in d.items():
            if n >= MC_BLOCK: occ[code] += 1
    feats = [code for code, o in occ.items() if 2 <= o <= nb - 1]
    fidx = {code: j for j, code in enumerate(feats)}
    F = np.zeros((nb, len(feats)))
    for i, d in enumerate(bcnt):
        tot = sum(d.values()) or 1
        for code, n in d.items():
            if code in fidx: F[i, fidx[code]] = n / tot
    Fz = (F - F.mean(0)) / (F.std(0) + 1e-9)
    bpw = defaultdict(float)
    for a, b, w in top:
        if a in blk_of and b in blk_of and blk_of[a] != blk_of[b]:
            i, j = blk_of[a], blk_of[b]; bpw[(min(i, j), max(i, j))] += w
    pairs = list(bpw)
    if len(pairs) < 2:
        return None, 0, 0
    D = np.array([Fz[i] - Fz[j] for i, j in pairs])
    val, vec = np.linalg.eigh(D @ D.T)
    axis = D.T @ vec[:, -1]; axis /= (np.linalg.norm(axis) + 1e-9)
    mode = os.environ.get("POLYSPLIT_MUTUAL_PAIRS", "auto")
    colourable = k_colourable(pairs, nsg)
    mutual = mode == "1" or (mode == "auto" and not colourable)
    cpairs = mutual_top_pairs(bpw, nsg) if mutual else pairs
    with open("pair_mode.txt", "w") as o:
        o.write("mutual\n" if mutual else "all\n")
    print(f"  pairs: {len(pairs)} block pairs, {nsg}-colourable={colourable} -> using {len(cpairs)} ({'mutual top-(K-1)' if mutual else 'all'})")
    if nsg == 2:
        lab = enforce_pairs(Fz @ axis, cpairs, nb)
    else:
        from sklearn.cluster import KMeans
        axes = D.T @ vec[:, -(nsg - 1):]
        axes /= (np.linalg.norm(axes, axis=0, keepdims=True) + 1e-9)
        emb = Fz @ axes
        lab = np.zeros(nb, dtype=int)
        while len(set(lab.tolist())) < nsg:
            sizes = {c: int(np.sum(lab == c)) for c in set(lab.tolist())}
            tgt = max(sizes, key=sizes.get)
            idx = np.where(lab == tgt)[0]
            if len(idx) < 2:
                break
            sub = KMeans(n_clusters=2, n_init=10, random_state=1).fit_predict(emb[idx])
            lab[idx[sub == 1]] = max(lab.tolist()) + 1
    cl = {c: f"S{int(lab[blk_of[c]]) + 1}" for c in blk_of}
    from itertools import permutations
    bw = np.array([sum(lengths[c] for c in com) for com in blocks], float)
    yb = []
    for com in blocks:
        bp = defaultdict(int)
        for c in com:
            if c in truth: bp[truth[c]] += lengths[c]
        yb.append(max(bp, key=bp.get) if bp else None)
    clusters = sorted({int(x) for x in lab}); tls = sorted({t for t in yb if t})
    bestm, bacc = {}, 0.0
    if tls:
        bestm = max((dict(zip(clusters, perm))
                     for perm in permutations(tls, min(len(clusters), len(tls)))),
                    key=lambda m: sum(bw[i] for i in range(nb) if m.get(int(lab[i])) == yb[i]))
        bacc = sum(bw[i] for i in range(nb) if bestm.get(int(lab[i])) == yb[i]) / bw.sum()
    cw = sum(lengths[c] for c in cl if c in truth)
    cacc = (sum(lengths[c] for c in cl if c in truth
                and bestm.get(int(lab[blk_of[c]])) == truth[c]) / cw) if cw else 0.0
    print(f"  [{tag:>22}] {nb:>3} blocks  block {100*bacc:5.1f}%  contig(bp) {100*cacc:5.1f}%")
    return cl, bacc, cacc


def main():
    with open("contacts.pkl", "rb") as f:
        contacts, lengths = pickle.load(f)
    truth = {}
    with open("../purity_compare/wg_purity.per_contig.tsv") as f:
        next(f)
        for line in f:
            c, L, a, cc_, mf, kls = line.rstrip("\n").split("\t")
            if kls.startswith("pure_"): truth[c] = kls.split("_", 1)[1]
    top = []; hedges = {}
    with open("homoeolog_edges.tsv") as f:
        next(f)
        for line in f:
            a, b, w = line.rstrip("\n").split("\t"); w = int(w)
            if w >= PAIR_MIN: top.append((a, b, w))
            if w >= 20000: hedges[(a, b) if a < b else (b, a)] = w

    nodes = [c for c in lengths if lengths[c] >= MIN_LEN]
    G = nx.Graph(); G.add_nodes_from(nodes)
    for (a, b), w in contacts.items():
        if lengths.get(a, 0) >= MIN_LEN and lengths.get(b, 0) >= MIN_LEN:
            G.add_edge(a, b, weight=w / ((lengths[a] * lengths[b]) ** 0.5))
    base_blocks = [com for com in (list(c) for c in louvain_communities(G, weight="weight",
                                   resolution=RES, seed=1)) if sum(lengths[x] for x in com) >= MIN_BLOCK]

    want = {c for com in base_blocks for c in com}
    seqs = {}; cur, cn = [], None
    with open("../flye_out/assembly.fasta") as f:
        for line in f:
            if line[0] == ">":
                if cn in want: seqs[cn] = "".join(cur)
                cn = line[1:].split()[0]; cur = []
            else: cur.append(line.strip())
        if cn in want: seqs[cn] = "".join(cur)
    cc = {c: contig_repeat_counts(seqs[c], K, MC_CONTIG) for c in want}

    print("config (edge_min, min_edges) -> result:")
    label_eval(base_blocks, cc, lengths, top, truth, "baseline (no decloud)")
    best = (0, None, None)
    for emin, mine in SWEEP:
        bl, ns = dechimerize([list(b) for b in base_blocks], G, hedges, lengths, emin, mine)
        cl, ba, ca = label_eval(bl, cc, lengths, top, truth, f"emin={emin//1000}k,n>={mine} ({ns} split)")
        if cl and ca > best[0]:
            best = (ca, cl, (emin, mine))
    if best[1]:
        with open("decloud_v2_contig_labels.tsv", "w") as o:
            o.write("contig\tlabel\ttruth\n")
            for c in best[1]:
                o.write(f"{c}\t{best[1][c]}\t{truth.get(c,'NA')}\n")
        print(f"\nBEST config {best[2]} contig {100*best[0]:.1f}% -> decloud_v2_contig_labels.tsv",
              file=sys.stderr)


if __name__ == "__main__":
    main()
