#!/usr/bin/env python3
import argparse, sys
from collections import Counter
import numpy as np


def fasta_names(fa):
    names = []
    with open(fa) as f:
        for line in f:
            if line[0] == ">":
                names.append(line[1:].split()[0])
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gsa", required=True)
    ap.add_argument("--lcp", required=True)
    ap.add_argument("--fa", required=True, help="contigs fasta (defines doc order)")
    ap.add_argument("--k", type=int, default=33)
    ap.add_argument("--min-copy", type=int, default=2)
    ap.add_argument("--max-copy", type=int, default=6, help="skip k-mers in > this many places (repeats)")
    ap.add_argument("--sample", type=int, default=1, help="process every Nth qualifying block")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    names = fasta_names(args.fa)
    gsa = np.memmap(args.gsa, dtype=[("doc", "<u4"), ("pos", "<u8")], mode="r")
    lcp = np.memmap(args.lcp, dtype=np.uint8, mode="r")
    n = len(gsa); assert len(lcp) == n
    docs = np.asarray(gsa["doc"])
    print(f"[sa] {n:,} suffixes; {len(names):,} contigs", file=sys.stderr)

    bnd = np.flatnonzero(np.asarray(lcp) < args.k)        # block starts (new k-mer)
    sizes = np.diff(np.append(bnd, n))
    nc = len(names)
    edges = Counter()

    # --- FAST vectorized path for copy==2 k-mers (a 33-mer shared by exactly 2 contigs = a clean
    #     homoeolog/paralog pair); this is the bulk of the discriminating signal ---
    s2 = bnd[sizes == 2]
    da, db = docs[s2].astype(np.int64), docs[s2 + 1].astype(np.int64)
    diff = da != db
    da, db = da[diff], db[diff]
    lo = np.minimum(da, db); hi = np.maximum(da, db)
    key = lo * nc + hi
    uk, cnt = np.unique(key, return_counts=True)
    for kk, c in zip(uk.tolist(), cnt.tolist()):
        edges[(kk // nc, kk % nc)] += c
    print(f"[copy2] {len(s2):,} copy-2 k-mers -> {len(edges):,} pairs", file=sys.stderr)

    # --- loop the rarer copy 3..max_copy blocks (few) ---
    if args.max_copy >= 3:
        sel = np.flatnonzero((sizes >= 3) & (sizes <= args.max_copy))
        if args.sample > 1:
            sel = sel[::args.sample]
        for bi in sel:
            s = int(bnd[bi]); e = s + int(sizes[bi])
            d = np.unique(docs[s:e])
            for x in range(d.size):
                for y in range(x + 1, d.size):
                    a, b = int(d[x]), int(d[y])
                    edges[(a, b) if a < b else (b, a)] += 1
        print(f"[copy3-{args.max_copy}] processed {len(sel):,} blocks; total {len(edges):,} pairs",
              file=sys.stderr)

    with open(args.out, "w") as o:
        o.write("contigA\tcontigB\tshared_kmers\n")
        for (a, b), w in edges.items():
            if a < nc and b < nc:
                o.write(f"{names[a]}\t{names[b]}\t{w}\n")
    print(f"[done] {len(edges):,} homoeolog/paralog edges -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
