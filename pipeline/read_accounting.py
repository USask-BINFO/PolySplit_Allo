#!/usr/bin/env python3
import sys, os, pickle, subprocess, json, argparse
from collections import defaultdict, Counter
from itertools import permutations

NONLAB = {"", "ambiguous", "unassigned", "?", "chimeric", "label"}
WIN = 10000
BINS = [("no homoeolog detected", None), ("<92%", 0.92), ("92-95%", 0.95), ("95-97%", 0.97),
        ("97-99%", 0.99), (">=99%", 1.01)]
UNPLACED = "primary not on anchored chromosome"


def read_config(path):
    rows = []
    for ln in open(path):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 11 or p[0].startswith("#") or p[0].lower() == "name":
            continue
        name, reads, run_dir, cmap, ref_id, ref, asm, total, rg, pc, sp = p[:11]
        rows.append(dict(key=f"{ref_id}_{reads.lower()}", name=name, reads=reads, run_dir=run_dir, chrom_map=cmap,
                         ref_id=ref_id, reference=ref, assembly=asm, total=int(total),
                         labels={m: f for m, f in [("PolySplit-Allo", f"{run_dir}/read_subg.tsv"), ("PolySplit-Ref", rg),
                                                   ("polyCRACKER", pc), ("SubPhaser", sp)] if f != "-"}))
    return rows


def load_map(path):
    return {p[0]: p[1] for p in (ln.split() for ln in open(path)) if len(p) >= 2}


def parse_paf(paf, cmap, cache):
    """read -> (truth, primary target, primary midpoint, primary MAPQ, read length)"""
    if os.path.exists(cache):
        return pickle.load(open(cache, "rb"))
    sums = defaultdict(lambda: defaultdict(int)); prim = {}; mapped = set()
    proc = subprocess.Popen(["cut", "-f1-12,17", paf], stdout=subprocess.PIPE, text=True, bufsize=1 << 24)
    for ln in proc.stdout:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 12:
            continue
        try:
            qlen, ts, te, nm, mq = int(p[1]), int(p[7]), int(p[8]), int(p[9]), int(p[11])
        except ValueError:
            continue
        r, tgt = p[0], p[5]
        mapped.add(r)
        if tgt in cmap:
            sums[r][cmap[tgt]] += nm
        if len(p) > 12 and p[12] == "tp:A:P" and (r not in prim or nm > prim[r][0]):
            prim[r] = (nm, tgt, (ts + te) // 2, mq, qlen)
    proc.wait()
    out = {}
    for r in mapped:
        sm = sums.get(r)
        truth = max(sm.items(), key=lambda kv: kv[1])[0] if sm else None
        pr = prim.get(r)
        out[r] = (truth,) + (pr[1:] if pr else (None, None, None, None))
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    pickle.dump(out, open(cache, "wb"))
    return out


def load_labels(path):
    """read -> subgenome label, 'ambiguous' or 'unassigned'. Files with *_hits columns (PolySplit-Ref):
    a no-call with zero hits is unassigned. Reads absent from a file are unassigned."""
    lab = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        hits = [i for i, h in enumerate(hdr) if h.endswith("_hits")]
        li = hdr.index("label") if "label" in hdr else 1
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= li:
                continue
            l = p[li]
            if l in NONLAB:
                l = "unassigned" if hits and sum(int(p[i]) for i in hits) == 0 else "ambiguous"
            lab[p[0]] = l
    return lab


def best_map(lab, truth):
    pl = sorted({v for v in lab.values() if v not in ("ambiguous", "unassigned")})
    tl = sorted(set(truth.values()))
    cnt = Counter((lab[r], t) for r, t in truth.items() if r in lab and lab[r] in pl)
    best, bs = {}, -1
    for perm in permutations(tl, min(len(pl), len(tl))):
        m = dict(zip(pl, perm)); s = sum(c for (x, y), c in cnt.items() if m.get(x) == y)
        if s > bs:
            best, bs = m, s
    return best


def score(lab, truth, mp, subset=None):
    tl = sorted(set(truth.values()))
    conf = {t: Counter() for t in tl}
    n = c = 0
    for r, t in truth.items():
        if subset is not None and r not in subset:
            continue
        l = lab.get(r, "unassigned")
        col = l if l in ("ambiguous", "unassigned") else mp.get(l, "ambiguous")
        conf[t][col] += 1; n += 1; c += (col == t)
    P, R = [], []
    for s in tl:
        pred = sum(conf[t][s] for t in tl); tot = sum(conf[s].values())
        P.append(conf[s][s] / pred if pred else 0.0); R.append(conf[s][s] / tot if tot else 0.0)
    mP, mR = sum(P) / len(P), sum(R) / len(R)
    return dict(n=n, correct=c, acc=100 * c / n if n else float("nan"), prec=100 * mP, rec=100 * mR,
                f1=100 * (2 * mP * mR / (mP + mR) if mP + mR else 0.0), conf={t: dict(conf[t]) for t in tl})


def load_windows(path):
    wins, allw = defaultdict(dict), defaultdict(set)
    with open(path) as f:
        next(f)
        for ln in f:
            w, s, L, cov, ident = ln.rstrip("\n").split("\t")
            chrom, rng = w.rsplit(":", 1); wi = int(rng.split("-")[0]) // WIN
            allw[chrom].add(wi); wins[chrom][wi] = (float(ident), float(cov))
    return wins, allw


def ident_bin(wins, allw, chrom, mid):
    if chrom not in allw or mid // WIN not in allw[chrom]:
        return None
    h = wins[chrom].get(mid // WIN)
    if h is None or not h[1] >= 0.5 or h[0] != h[0]:
        return BINS[0][0]
    return next(name for name, ub in BINS[1:] if h[0] < ub)


def run(d, outdir, extra):
    cmap = load_map(d["chrom_map"])
    reads = parse_paf(f"{d['run_dir']}/reads_to_ref.paf", cmap, f"{outdir}/reads/cache/{d['key']}.pkl")
    truth = {r: v[0] for r, v in reads.items() if v[0] is not None}
    tprim = {r: cmap[v[1]] for r, v in reads.items() if v[1] in cmap}
    tp20 = {r: t for r, t in tprim.items() if reads[r][3] >= 20}
    hq = {r for r in truth if r in tp20}; lq = {r for r in truth if r in tprim and r not in tp20}
    res = dict(key=d["key"], total=d["total"], mapped=len(reads), unmapped=d["total"] - len(reads),
               unplaced_only=sum(1 for v in reads.values() if v[0] is None), truth_reads=len(truth),
               truth_by_subg=dict(Counter(truth.values())),
               truth_primary_on_unplaced=sum(1 for r in truth if r not in tprim),
               truth_sum_vs_prim_agree=sum(truth[r] == tprim[r] for r in truth if r in tprim),
               truth_sum_vs_prim_n=sum(1 for r in truth if r in tprim),
               mapq_ge20_by_subg=dict(Counter(truth[r] for r in hq)), mapq_lt20_by_subg=dict(Counter(truth[r] for r in lq)),
               truth_p20_n=len(tp20), truth_sum_vs_p20_agree=sum(truth.get(r) == t for r, t in tp20.items()))
    wpath = f"{outdir}/reference/{d['ref_id']}/windows/windows.tsv"
    rb = None
    if os.path.exists(wpath):
        wins, allw = load_windows(wpath)
        rb = {r: (ident_bin(wins, allw, reads[r][1], reads[r][2]) if r in tprim else None) or UNPLACED for r in truth}
        res["ident_bin_counts"] = dict(Counter(rb.values()))
        res["ident_bin_counts_mapq_ge20"] = dict(Counter(rb[r] for r in hq))
        res["ident_bin_counts_mapq_lt20"] = dict(Counter(rb[r] for r in lq))
    res["methods"] = {}
    for m, path in list(d["labels"].items()) + extra:
        if not os.path.exists(path):
            continue
        lab = load_labels(path); mp = best_map(lab, truth)
        mr = dict(map=mp, T_sum=score(lab, truth, mp), T_prim20=score(lab, tp20, best_map(lab, tp20)),
                  T_sum_mapq_ge20=score(lab, truth, mp, hq), T_sum_mapq_lt20=score(lab, truth, mp, lq))
        fate = Counter()
        for r, v in reads.items():
            if v[0] is None:
                l = lab.get(r, "unassigned"); fate["assigned" if l not in ("ambiguous", "unassigned") else l] += 1
        mr["no_truth_fate"] = dict(fate)
        if rb is not None:
            mr["by_ident_bin"] = {b: score(lab, truth, mp, {r for r, x in rb.items() if x == b})
                                  for b in [b for b, _ in BINS] + [UNPLACED] if b in res["ident_bin_counts"]}
        res["methods"][m] = mr
        s = mr["T_sum"]
        print(f"{d['key']:16s} {m:28s} acc={s['acc']:.2f} P={s['prec']:.1f} R={s['rec']:.1f} F1={s['f1']:.1f}", flush=True)
    os.makedirs(f"{outdir}/reads", exist_ok=True)
    json.dump(res, open(f"{outdir}/reads/{d['key']}.json", "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config"); ap.add_argument("outdir")
    ap.add_argument("--only", default=None, help="process one dataset KEY")
    ap.add_argument("--label", nargs=3, action="append", default=[], metavar=("KEY", "NAME", "FILE"),
                    help="score an extra read-label file for dataset KEY")
    a = ap.parse_args()
    for d in read_config(a.config):
        if a.only and d["key"] != a.only:
            continue
        run(d, a.outdir, [(n, f) for k, n, f in a.label if k == d["key"]])


if __name__ == "__main__":
    main()
