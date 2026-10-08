#!/usr/bin/env python3
import sys, os, json, pickle, bisect, argparse
from collections import defaultdict, Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_accounting import read_config


def repeat_index(path):
    iv = defaultdict(list); cur = None
    for ln in open(path):
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
        else:
            a, b = ln.split(" - "); iv[cur].append((int(a), int(b) + 1))
    idx = {}
    for s, L in iv.items():
        L.sort(); st = np.array([x[0] for x in L]); en = np.array([x[1] for x in L]); idx[s] = (st, en)
    return idx


def masked(idx, s, a, b):
    if s not in idx:
        return 0.0
    st, en = idx[s]
    tot = sum(max(0, min(en[k], b) - max(st[k], a)) for k in range(max(bisect.bisect_right(st, a) - 1, 0), bisect.bisect_left(st, b)))
    return tot / max(1, b - a)


def run(d, outdir):
    R = f"{outdir}/reference/{d['ref_id']}"
    idx = repeat_index(f"{R}/wm.intervals")
    seqlen = {l.split()[0]: int(l.split()[1]) for l in open(f"{R}/ref.fa.fai")}
    unpl = {l.strip() for l in open(f"{R}/unplaced.names")}
    hits = defaultdict(list)
    for ln in open(f"{R}/unplaced_vs_organelles.paf"):
        p = ln.split("\t")
        if int(p[10]) and int(p[9]) / int(p[10]) >= 0.9:
            hits[p[0]].append((int(p[2]), int(p[3])))
    organellar = set()
    for s, h in hits.items():
        cov = 0; ce = -1
        for a, b in sorted(h):
            if b > ce:
                cov += b - max(a, ce); ce = b
        if cov / seqlen[s] >= 0.5:
            organellar.add(s)
    reads = pickle.load(open(f"{outdir}/reads/cache/{d['key']}.pkl", "rb"))
    bp = Counter()
    for v in reads.values():
        if v[1] is not None and v[4]:
            bp[v[1]] += v[4]
    dep = {s: bp[s] / seqlen[s] for s in bp if seqlen.get(s)}
    chrom_dep = float(np.median([dep[s] for s in dep if s not in unpl]))
    highcopy = {s for s in unpl if dep.get(s, 0) >= 3 * chrom_dep} - organellar
    cat = Counter(); rep = defaultdict(list)
    for truth, tgt, mid, mq, ql in reads.values():
        if tgt is None:
            c = "mapped, no primary record"
        elif truth is not None:
            c = "truth, primary on unplaced" if tgt in unpl else ("truth, MAPQ>=20" if mq >= 20 else "truth, MAPQ<20")
        else:
            c = ("no truth: organellar scaffold" if tgt in organellar else
                 "no truth: high-copy scaffold" if tgt in highcopy else "no truth: other unplaced")
        cat[c] += 1
        if tgt is not None and mid is not None:
            rep[c].append(masked(idx, tgt, max(0, mid - 5000), mid + 5000))
    T = d["total"]
    out = dict(key=d["key"], total=T, unmapped=T - len(reads), chromosome_median_depth=round(chrom_dep, 1),
               n_organellar_scaffolds=len(organellar), organellar_Mb=round(sum(seqlen[s] for s in organellar) / 1e6, 2),
               n_highcopy_scaffolds=len(highcopy), highcopy_Mb=round(sum(seqlen[s] for s in highcopy) / 1e6, 2),
               categories={c: dict(reads=n, pct_of_all=round(100 * n / T, 2),
                                   mean_repeat_frac=round(float(np.mean(rep[c])), 3) if rep[c] else None)
                           for c, n in cat.most_common()},
               repeat_frac_chromosomes=round(float(np.mean([masked(idx, s, 0, seqlen[s]) for s in seqlen if s not in unpl])), 3),
               repeat_frac_unplaced=round(sum(masked(idx, s, 0, seqlen[s]) * seqlen[s] for s in unpl) /
                                          max(1, sum(seqlen[s] for s in unpl)), 3))
    os.makedirs(f"{outdir}/unlabelled", exist_ok=True)
    json.dump(out, open(f"{outdir}/unlabelled/{d['key']}.json", "w"), indent=1)
    print(d["key"], json.dumps(out["categories"]), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("config"); ap.add_argument("outdir"); ap.add_argument("--only")
    a = ap.parse_args()
    for d in read_config(a.config):
        if not a.only or d["key"] == a.only:
            run(d, a.outdir)
