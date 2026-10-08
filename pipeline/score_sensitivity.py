#!/usr/bin/env python3
import os, re, sys, glob, json, pickle
from collections import Counter
from itertools import permutations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_accounting import read_config, load_labels, best_map, score

NONLAB = {"", "ambiguous", "unassigned", "?", "chimeric", "label", "NA"}


def contig_labels(path):
    d = {}
    with open(path) as f:
        next(f)
        for ln in f:
            x = ln.rstrip("\n").split("\t")
            if len(x) > 1 and x[1] not in NONLAB:
                d[x[0]] = x[1]
    return d


def matched(a, b, weight=None):
    """best one-to-one agreement between two labelings (weighted by contig length if given)"""
    common = [c for c in a if c in b]
    la, lb = sorted(set(a.values())), sorted(set(b.values()))
    best = 0.0
    for perm in permutations(lb, len(la)) if len(la) <= len(lb) else [lb]:
        m = dict(zip(la, perm))
        best = max(best, sum((weight or {}).get(c, 1) for c in common if m.get(a[c]) == b[c]))
    tot = sum((weight or {}).get(c, 1) for c in common)
    return 100 * best / tot if tot else float("nan")


cfg, outdir = sys.argv[1], sys.argv[2]
rows = []
for d in read_config(cfg):
    runs = sorted(glob.glob(f"{outdir}/sensitivity/{d['key']}/*/read_subg.tsv"))
    if not runs:
        continue
    truth_c, L = {}, {}
    for ln in open(f"{d['run_dir']}/wg_purity.per_contig.tsv"):
        x = ln.rstrip("\n").split("\t")
        if x[0] != "contig":
            L[x[0]] = int(x[1])
            if x[5].startswith("pure_"):
                truth_c[x[0]] = x[5][5:]
    reads = pickle.load(open(f"{outdir}/reads/cache/{d['key']}.pkl", "rb"))
    truth = {r: v[0] for r, v in reads.items() if v[0] is not None}
    default = contig_labels(f"{d['run_dir']}/all_contig_labels_repaired.tsv")
    for rs in runs:
        W = os.path.dirname(rs)
        pred = contig_labels(f"{W}/all_contig_labels_repaired.tsv")
        m = re.search(r"(\d+) splits\s+->\s+(\d+) -> (\d+) blocks", open(f"{W}/dechim.log").read())
        fl = re.search(r"-> (\d+) flips", open(f"{W}/repair.log").read())
        lab = load_labels(rs)
        rows.append(dict(key=d["key"], setting=os.path.basename(W), blocks=m and int(m.group(3)), splits=m and int(m.group(1)),
                         flips=fl and int(fl.group(1)),
                         contig_acc=round(matched({c: l for c, l in pred.items() if c in truth_c}, truth_c, L), 2),
                         read_acc=round(score(lab, truth, best_map(lab, truth))["acc"], 2),
                         same_labels_as_default=round(matched(pred, default), 2)))
        print("\t".join(str(v) for v in rows[-1].values()), flush=True)
json.dump(rows, open(f"{outdir}/sensitivity/summary.json", "w"), indent=1)
