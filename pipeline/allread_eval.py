#!/usr/bin/env python3
import argparse, sys
from collections import defaultdict, Counter
from itertools import permutations

NONLAB = {"", "ambiguous", "unassigned", "?", "chimeric", "label"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True, help="predicted read labels TSV (read<TAB>label...)")
    ap.add_argument("--ref-paf", required=True, help="reads -> reference PAF (for ground truth)")
    ap.add_argument("--total-reads", type=int, required=True)
    ap.add_argument("--chrom-subg", default=None,
                    help="optional TSV (reference_chrom <TAB> subgenome). If omitted, the subgenome is "
                         "the first char after the last '.' of the chrom name (B. napus A/C style).")
    a = ap.parse_args()

    cmap = {}
    if a.chrom_subg:
        for ln in open(a.chrom_subg):
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                cmap[p[0]] = p[1]

    def subg(t):
        if cmap:                       # explicit map given: only mapped chroms carry a truth
            return cmap.get(t)         # None for unmapped targets (e.g. unplaced scaffolds)
        tok = t.split(".")[-1]         # no map -> derive subgenome from chrom name (B. napus A/C)
        return tok[:1] if tok[:1].isalpha() else None

    lab = {}
    with open(a.labels) as f:
        next(f)
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                lab[p[0]] = p[1]

    tscore = defaultdict(lambda: defaultdict(int))     # read -> {subgenome: matched bases}
    seen = set()                                       # reads aligning to ANY reference sequence
    with open(a.ref_paf) as f:
        for ln in f:
            p = ln.split("\t")
            if len(p) < 11:
                continue
            seen.add(p[0])
            sg = subg(p[5])
            if sg is None:                             # mapped only to unplaced/un-mapped seq -> no truth
                continue
            try:
                nm = int(p[9])
            except ValueError:
                continue
            tscore[p[0]][sg] += nm
    truth = {r: max(d.items(), key=lambda kv: kv[1])[0] for r, d in tscore.items() if d}
    n_map, n_tru = len(seen), len(truth)

    pred_labels = sorted({lab[r] for r in lab if lab[r] not in NONLAB})
    truth_labels = sorted(set(truth.values()))
    cb = [(lab[r], truth[r]) for r in lab if lab[r] in pred_labels and r in truth]
    best = max((dict(zip(pred_labels, perm))
                for perm in permutations(truth_labels, min(len(pred_labels), len(truth_labels)))),
               key=lambda m: sum(m.get(x) == y for x, y in cb), default={})
    mp = lambda l: best.get(l, l)

    called = [r for r in lab if lab[r] in pred_labels and r in truth]
    correct = sum(1 for r in called if mp(lab[r]) == truth[r])
    labeled_truth = sum(1 for r in lab if r in truth)

    print(f"labels: {a.labels.split('/')[-1]}   K(pred)={len(pred_labels)}  truth={truth_labels}")
    print(f"  predicted->truth map : {best}")
    print(f"  called+truth={len(called):,}  correct={correct:,}  total reads={a.total_reads:,}")
    print(f"  reads: chromosome-truth={n_tru:,}  scaffold-only/no-truth={n_map-n_tru:,}  unmapped={a.total_reads-n_map:,}")
    print(f"  T1 (called)                 : {100*correct/max(1,len(called)):.2f}%")
    print(f"  T2 (over chromosome-truth)  : {100*correct/max(1,labeled_truth):.2f}%")
    print(f"  T3 (over all reads)         : {100*correct/a.total_reads:.2f}%  (incl. scaffold-only/unmapped as errors)")

    cols = truth_labels + ["ambiguous", "unassigned"]
    conf = {t: Counter() for t in truth_labels}
    for r, t in truth.items():
        if r in lab:
            l = lab[r]
            conf[t][mp(l) if l in pred_labels else "ambiguous"] += 1
        else:
            conf[t]["unassigned"] += 1
    print("  confusion (rows=truth, cols=predicted):")
    print("    " + "\t".join([""] + cols))
    for t in truth_labels:
        print("    " + "\t".join([t] + [str(conf[t][c]) for c in cols]))


if __name__ == "__main__":
    main()
