#!/usr/bin/env python3
import sys, os, json, pickle, glob, argparse
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_accounting import read_config, load_map

WIN, SLACK, MIN_LEN, FRAC = 10000, 20000, 50000, 0.5
THRESH = [100, 1000, 10000, 30000, 100000, 300000, 1000000]


def run(d, outdir):
    rd, cs = d["run_dir"], load_map(d["chrom_map"])
    pure = {}
    for ln in open(f"{rd}/wg_purity.per_contig.tsv"):
        x = ln.rstrip("\n").split("\t")
        if x[0] != "contig" and x[5].startswith("pure_"):
            pure[x[0]] = x[5][5:]
    mb = defaultdict(Counter); blocks = defaultdict(list)
    for ln in open(f"{rd}/contigs_to_ref.paf"):
        p = ln.split("\t")
        if p[5] in cs:
            mb[p[0]][p[5]] += int(p[9]); blocks[(p[0], p[5])].append((int(p[7]), int(p[8])))
    anch = {c: mb[c].most_common(1)[0][0] for c in pure if c in mb}
    whits = defaultdict(list); cnt = defaultdict(Counter)
    for paf in glob.glob(f"{outdir}/reference/{d['ref_id']}/windows/cross.*.paf"):
        for ln in open(paf):
            p = ln.split("\t"); ch, rng = p[0].rsplit(":", 1)
            whits[(ch, int(rng.split("-")[0]) // WIN)].append((p[5], int(p[7]), int(p[8]))); cnt[ch][p[5]] += 1
    partner = {ch: {t for t, n in c.items() if n >= 0.10 * sum(c.values())} for ch, c in cnt.items()}

    def syntenic(a, b):
        for x, y in ((a, b), (b, a)):
            X, Y = anch[x], anch[y]; yb = blocks[(y, Y)]
            for s, e in blocks[(x, X)]:
                for wi in range(s // WIN, e // WIN + 1):
                    for t, hs, he in whits.get((X, wi), ()):
                        if t == Y and any(hs < ye + SLACK and he > ys - SLACK for ys, ye in yb):
                            return True
        return False

    cache = {}

    def klass(a, b):
        if (a, b) not in cache:
            cache[(a, b)] = ("unanchored" if a not in anch or b not in anch else
                             "same subgenome" if pure[a] == pure[b] else
                             "syntenic homoeolog" if syntenic(a, b) else
                             "homoeologous chromosomes, outside mapped windows"
                             if anch[b] in partner.get(anch[a], ()) or anch[a] in partner.get(anch[b], ())
                             else "other cross-subgenome")
        return cache[(a, b)]

    raw, mx = [], defaultdict(int)
    with open(f"{rd}/homoeolog_edges.tsv") as f:
        next(f)
        for ln in f:
            a, b, w = ln.rstrip("\n").split("\t"); w = int(w)
            raw.append((a, b, w)); mx[a] = max(mx[a], w); mx[b] = max(mx[b], w)
    out = dict(key=d["key"], n_edges=len(raw), sweep=[], weights_by_class={})
    for t in THRESH:
        for rule in (False, True):
            sel = [(a, b, w) for a, b, w in raw if w >= t and (not rule or (w >= FRAC * mx[a] and w >= FRAC * mx[b]))]
            ws = [w for _, _, w in sel]
            out["sweep"].append(dict(thresh=t, strongest_partner_rule=rule, n=len(sel), classes=dict(Counter(klass(a, b) for a, b, _ in sel)),
                                     w_median=int(np.median(ws)) if ws else None, w_min=min(ws) if ws else None, w_max=max(ws) if ws else None))
    for a, b, w in raw:
        if w >= 100:
            out["weights_by_class"].setdefault(klass(a, b), []).append(w)

    contacts, L = pickle.load(open(f"{rd}/contacts.pkl", "rb"))
    keep = {c for c in anch if L.get(c, 0) >= MIN_LEN}
    n_, g_, cis_ = [], [], []; top_raw, top_norm = {}, {}
    for (a, b), n in contacts.items():
        if a in keep and b in keep:
            g = (L[a] * L[b]) ** 0.5
            n_.append(n); g_.append(g); cis_.append(anch[a] == anch[b])
            for x, y in ((a, b), (b, a)):
                if n > top_raw.get(x, (0, None))[0]: top_raw[x] = (n, y)
                if n / g > top_norm.get(x, (0, None))[0]: top_norm[x] = (n / g, y)
    n_, g_, cis_ = np.array(n_, float), np.array(g_), np.array(cis_)
    w_ = n_ / g_
    hic = dict(pairs=int(len(n_)), cis_pairs=int(cis_.sum()), contigs=len(keep))
    for nm, m in (("same-chromosome", cis_), ("different-chromosome", ~cis_)):
        hic[f"spearman_raw_vs_len[{nm}]"] = round(float(spearmanr(n_[m], g_[m])[0]), 3)
        hic[f"spearman_norm_vs_len[{nm}]"] = round(float(spearmanr(w_[m], g_[m])[0]), 3)
    hic["auc_same_vs_diff_chrom_raw"] = round(float(roc_auc_score(cis_, n_)), 3)
    hic["auc_same_vs_diff_chrom_norm"] = round(float(roc_auc_score(cis_, w_)), 3)
    hic["top_partner_same_chrom_raw"] = round(100 * np.mean([anch[x] == anch[y] for x, (_, y) in top_raw.items()]), 1)
    hic["top_partner_same_chrom_norm"] = round(100 * np.mean([anch[x] == anch[y] for x, (_, y) in top_norm.items()]), 1)
    out["hic"] = hic
    os.makedirs(f"{outdir}/edges", exist_ok=True)
    json.dump(out, open(f"{outdir}/edges/{d['key']}.json", "w"), indent=1)
    s = next(x for x in out["sweep"] if x["thresh"] == 100000 and x["strongest_partner_rule"])
    print(d["key"], "retained edges:", s["n"], s["classes"], flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("config"); ap.add_argument("outdir"); ap.add_argument("--only")
    a = ap.parse_args()
    for d in read_config(a.config):
        if not a.only or d["key"] == a.only:
            run(d, a.outdir)
