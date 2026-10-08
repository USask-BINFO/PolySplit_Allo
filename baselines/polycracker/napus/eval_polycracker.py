#!/usr/bin/env python3
import argparse, glob, os, re, sys
from collections import defaultdict, Counter
from itertools import permutations
import bisect


def parse_header(h):
    h = h[1:].split()[0]
    m = re.match(r"^(.*)_(\d+)_(\d+)$", h)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    # fallback: scaffold:start-end
    m = re.match(r"^(.*):(\d+)-(\d+)$", h)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return (h, None, None)


def load_chunks(subg_fastas):
    chunks = []
    for ci, fa in enumerate(sorted(subg_fastas)):
        with open(fa) as f:
            for line in f:
                if line.startswith(">"):
                    contig, s, e = parse_header(line.rstrip("\n"))
                    chunks.append((ci, contig, s, e))
    return chunks


def load_geneblk(paths):
    t = defaultdict(list)
    for p in paths:
        with open(p) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            ci = {n: i for i, n in enumerate(hdr)}
            cc = ci.get("nam_contig", 3); cp = ci.get("start", 4); cs = ci.get("subgenome", 1)
            for line in f:
                x = line.rstrip("\n").split("\t")
                t[x[cc]].append((int(x[cp]), int(x[cs])))
    for c in t:
        t[c].sort()
    return t


def chunk_truth_AC(contig):
    return "A" if ".A" in contig else ("C" if ".C" in contig else None)


def chunk_truth_trip(contig, s, e, truth):
    arr = truth.get(contig)
    if not arr or s is None:
        return None
    starts = [p for p, _ in arr]
    lo = bisect.bisect_left(starts, s); hi = bisect.bisect_right(starts, e)
    votes = Counter(arr[i][1] for i in range(lo, hi))
    return votes.most_common(1)[0][0] if votes else None


def best_perm_accuracy(pairs, classes):
    preds = sorted({p for p, _ in pairs})
    best = 0.0; bestmap = None
    # try all assignments of cluster->class (clusters may be <= len(classes))
    for perm in permutations(classes, len(preds)):
        mp = dict(zip(preds, perm))
        acc = sum(mp[p] == t for p, t in pairs) / len(pairs)
        if acc > best:
            best, bestmap = acc, mp
    return best, bestmap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subg-fastas", nargs="+", required=True,
                    help="polyCRACKER extractedSubgenomes/*.fa (one per predicted subgenome)")
    ap.add_argument("--mode", choices=["AC", "triplication"], required=True)
    ap.add_argument("--geneblk", nargs="*", default=[],
                    help="geneblk_A.gene_loci.tsv geneblk_C.gene_loci.tsv (triplication mode)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    chunks = load_chunks(args.subg_fastas)
    print(f"[chunks] {len(chunks):,} chunks across {len(args.subg_fastas)} predicted subgenomes",
          file=sys.stderr)

    if args.mode == "AC":
        classes = ["A", "C"]
        pairs = [(ci, chunk_truth_AC(c)) for (ci, c, s, e) in chunks if chunk_truth_AC(c)]
    else:
        truth = load_geneblk(args.geneblk)
        classes = [1, 2, 3]
        pairs = []
        for (ci, c, s, e) in chunks:
            tt = chunk_truth_trip(c, s, e, truth)
            if tt:
                pairs.append((ci, tt))

    if not pairs:
        print("[error] no chunks could be matched to truth (check header format / geneblk paths)",
              file=sys.stderr); sys.exit(1)

    acc, mp = best_perm_accuracy(pairs, classes)
    base = Counter(t for _, t in pairs).most_common(1)[0][1] / len(pairs)
    # bin composition
    comp = defaultdict(Counter)
    for ci, t in pairs:
        comp[ci][t] += 1

    lines = [f"polyCRACKER {args.mode}: chunks evaluated = {len(pairs):,}",
             f"  best-permutation accuracy = {acc:.3f}   (majority baseline {base:.3f}, "
             f"chance {1/len(classes):.3f})",
             f"  cluster->class map = {mp}",
             "  cluster composition (cluster: trueclass=count):"]
    for ci in sorted(comp):
        lines.append(f"    cluster {ci}: " + ", ".join(f"{k}={v}" for k, v in comp[ci].most_common()))
    out = "\n".join(lines)
    print(out)
    if args.out:
        open(args.out, "w").write(out + "\n")


if __name__ == "__main__":
    main()
