#!/usr/bin/env python3
import argparse, sys
from collections import defaultdict, Counter

NONLAB = {"", "ambiguous", "unassigned", "?", "chimeric", "label"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True, help="minimap2 reads -> assembly PAF")
    ap.add_argument("--contig-labels", default="contig_labels.tsv")
    ap.add_argument("--compare", default=None, help="optional reference-guided read labels to cross-check")
    ap.add_argument("--min-conf", type=float, default=0.6, help="min vote fraction to keep a read")
    ap.add_argument("--weight", choices=["ident", "sum"], default="ident",
                    help="'ident' = matched bases x identity (down-weights weaker homoeolog alignments); "
                         "'sum' = matched bases only")
    ap.add_argument("--out", default="read_subg.tsv")
    args = ap.parse_args()

    clab = {}
    with open(args.contig_labels) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1] not in NONLAB:
                clab[p[0]] = p[1]
    labels = sorted(set(clab.values()))
    print(f"[labels] {len(clab)} labelled contigs across {len(labels)} subgenomes: {labels}", file=sys.stderr)

    score = defaultdict(lambda: defaultdict(float))   # read -> {subgenome: weight}
    n = 0
    with open(args.paf) as f:
        for line in f:
            p = line.split("\t")
            if len(p) < 11:
                continue
            n += 1
            lab = clab.get(p[5])
            if lab is None:
                continue
            try:
                nm = int(p[9]); bl = int(p[10])
            except ValueError:
                continue
            w = nm * (nm / bl) if (args.weight == "ident" and bl) else nm
            score[p[0]][lab] += w
    print(f"[paf] {n:,} alignment lines; {len(score):,} reads hit a labelled contig", file=sys.stderr)

    calls = {}; counts = Counter()
    with open(args.out, "w") as o:
        o.write("read\tlabel\tconf\tweights\n")
        for r, sc in score.items():
            tot = sum(sc.values())
            if tot == 0:
                continue
            lab, w = max(sc.items(), key=lambda kv: kv[1])
            conf = w / tot
            wstr = ";".join(f"{k}:{v:.0f}" for k, v in sorted(sc.items()))
            if conf < args.min_conf:
                lab = "ambiguous"
            else:
                calls[r] = lab
            counts[lab] += 1
            o.write(f"{r}\t{lab}\t{conf:.3f}\t{wstr}\n")
    print(f"[reads] {dict(counts)}  (min-conf {args.min_conf})", file=sys.stderr)

    if args.compare:                                   # cross-check vs another label set (any names)
        ref = {}
        with open(args.compare) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[1] not in NONLAB:
                    ref[p[0]] = p[1]
        shared = [(calls[r], ref[r]) for r in calls if r in ref]
        if shared:
            from itertools import permutations
            cl = sorted({x for x, _ in shared}); tl = sorted({y for _, y in shared})
            best = max((dict(zip(cl, perm)) for perm in permutations(tl, min(len(cl), len(tl)))),
                       key=lambda m: sum(m.get(x) == y for x, y in shared), default={})
            agree = sum(best.get(x) == y for x, y in shared) / len(shared)
            print(f"[validate] {len(shared):,} reads in both; best-map agreement {100*agree:.1f}%",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
