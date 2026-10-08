#!/usr/bin/env python3
import sys, pickle, os
from collections import defaultdict
import networkx as nx
from networkx.algorithms.community import louvain_communities
from decloud_blocks_v2 import (contig_repeat_counts, label_eval,
                               MIN_LEN, MIN_BLOCK, MIN_SUB, K, MC_CONTIG,
                               RES, SUBRES, PAIR_MIN, MAX_ITERS)

NSG = int(os.environ.get("POLYSPLIT_NSG", "2"))
FASTA = os.environ.get("POLYSPLIT_FASTA", "../flye_out/assembly.fasta")
TRUTH = os.environ.get("POLYSPLIT_TRUTH", "../purity_compare/wg_purity.per_contig.tsv")
CONTACTS = os.environ.get("POLYSPLIT_CONTACTS", "contacts.pkl")
EDGES = os.environ.get("POLYSPLIT_EDGES", "homoeolog_edges.tsv")


def structural_dechimerize(blocks, G, strong_edges, lengths, verbose=True):
    nsplit = 0
    for it in range(MAX_ITERS):
        changed = False
        new = []
        for com in blocks:
            sub = G.subgraph(com)
            if sub.number_of_edges() == 0 or sub.size(weight="weight") == 0:
                parts = [list(com)]
            else:
                parts = [list(c) for c in louvain_communities(sub, weight="weight",
                                                              resolution=SUBRES, seed=1)]
            big = [p for p in parts if sum(lengths[c] for c in p) >= MIN_SUB]
            seam = None
            if len(big) >= 2:
                part_of = {c: i for i, p in enumerate(big) for c in p}
                for a, b, w in strong_edges:
                    if a in part_of and b in part_of and part_of[a] != part_of[b]:
                        seam = (a, b, w)
                        break
            if seam is not None:
                nsplit += 1
                changed = True
                assigned = {c for p in big for c in p}
                result = [list(p) for p in big]
                for c in com:
                    if c in assigned:
                        continue
                    best, bw = 0, -1.0
                    for i, p in enumerate(big):
                        wsum = sum(G[c][d]["weight"] for d in p if G.has_edge(c, d))
                        if wsum > bw:
                            bw, best = wsum, i
                    result[best].append(c)
                new.extend(result)
                if verbose:
                    sizes = " + ".join(f"{sum(lengths[c] for c in p)/1e6:.1f}Mb" for p in result)
                    a, b, w = seam
                    print(f"  [split it{it}] block {sum(lengths[c] for c in com)/1e6:.1f}Mb "
                          f"-> {len(result)} pure pieces [{sizes}]   seam: {a}–{b} ({w:,} 33-mers)")
            else:
                new.append(com)
        blocks = new
        if not changed:
            break
    return blocks, nsplit


def main():
    with open(CONTACTS, "rb") as f:
        contacts, lengths = pickle.load(f)
    truth = {}
    if os.path.exists(TRUTH):
        with open(TRUTH) as f:
            next(f)
            for line in f:
                c, L, a, cc_, mf, kls = line.rstrip("\n").split("\t")
                if kls.startswith("pure_"): truth[c] = kls.split("_", 1)[1]
    FRAC = float(os.environ.get("POLYSPLIT_EDGE_FRAC", "0.5"))
    raw = []; mx = {}
    with open(EDGES) as f:
        next(f)
        for line in f:
            a, b, w = line.rstrip("\n").split("\t"); w = int(w)
            raw.append((a, b, w))
            if w > mx.get(a, 0): mx[a] = w
            if w > mx.get(b, 0): mx[b] = w
    top = [(a, b, w) for (a, b, w) in raw
           if w >= PAIR_MIN and w >= FRAC * mx[a] and w >= FRAC * mx[b]]

    nodes = [c for c in lengths if lengths[c] >= MIN_LEN]
    G = nx.Graph(); G.add_nodes_from(nodes)
    for (a, b), w in contacts.items():
        if lengths.get(a, 0) >= MIN_LEN and lengths.get(b, 0) >= MIN_LEN:
            G.add_edge(a, b, weight=w / ((lengths[a] * lengths[b]) ** 0.5))
    base_blocks = [com for com in (list(c) for c in louvain_communities(G, weight="weight",
                                   resolution=RES, seed=1)) if sum(lengths[x] for x in com) >= MIN_BLOCK]

    want = {c for com in base_blocks for c in com}
    seqs = {}; cur, cn = [], None
    with open(FASTA) as f:
        for line in f:
            if line[0] == ">":
                if cn in want: seqs[cn] = "".join(cur)
                cn = line[1:].split()[0]; cur = []
            else: cur.append(line.strip())
        if cn in want: seqs[cn] = "".join(cur)
    cc = {c: contig_repeat_counts(seqs[c], K, MC_CONTIG) for c in want}

    print(f"base blocks: {len(base_blocks)}  ({len(top)} strong homoeolog edges, >= {PAIR_MIN//1000}k)\n")
    print("=== STRUCTURAL de-chimerize (no tuned count, no truth) ===")
    bl, ns = structural_dechimerize([list(b) for b in base_blocks], G, top, lengths)
    print(f"\n{ns} splits  ->  {len(base_blocks)} -> {len(bl)} blocks\n")
    print("config                          blocks  block%  contig(bp)%")
    label_eval(base_blocks, cc, lengths, top, truth, "baseline (no decloud)", nsg=NSG)
    cl, _, _ = label_eval(bl, cc, lengths, top, truth, "structural de-chimerize", nsg=NSG)

    out = "decloud_structural_contig_labels.tsv"
    with open(out, "w") as o:
        o.write("contig\tlabel\ttruth\n")
        for c in cl:
            o.write(f"{c}\t{cl[c]}\t{truth.get(c, 'NA')}\n")
    print(f"\n-> {out}  ({len(cl):,} block contigs labelled; truth-free de-chimerize)")


if __name__ == "__main__":
    main()
