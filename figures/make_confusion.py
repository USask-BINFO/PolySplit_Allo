#!/usr/bin/env python3
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
from read_accounting import read_config

NAME = {"A": "S1 (A)", "C": "S2 (C)"}
rows = read_config(sys.argv[1]); outdir, figdir = sys.argv[2], sys.argv[3]
genomes = list(dict.fromkeys(d["ref_id"] for d in rows))


def label(d):
    w = d["name"].split()
    return r"$\it{" + r"\ ".join(w[:2]) + "}$ " + " ".join(w[2:]).replace("x)", r"$\times$)")


def matrix(key, method):
    s = json.load(open(f"{outdir}/reads/{key}.json"))["methods"][method]["T_sum"]
    tl = sorted(s["conf"])
    M = np.array([[s["conf"][t].get(c, 0) for c in tl + ["ambiguous", "unassigned"]] for t in tl], float)
    return M, [NAME.get(t, t) for t in tl], s["acc"]


def panel(ax, M, tl, title):
    P = 100 * M / M.sum(1, keepdims=True)
    ax.imshow(P, cmap="Blues", vmin=0, vmax=100, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = int(M[i, j])
            ax.text(j, i, f"{v / 1000:.1f}k" if v >= 1000 else str(v), ha="center", va="center", fontsize=8,
                    color="white" if P[i, j] > 55 else "#222")
    ax.set_yticks(range(len(tl))); ax.set_yticklabels(tl, fontsize=8.5)
    ax.set_xticks(range(len(tl) + 2)); ax.set_xticklabels(tl + ["amb.", "unas."], fontsize=8, rotation=30, ha="right")
    ax.set_title(title, fontsize=10)


methods = ["PolySplit-Ref", "PolySplit-Allo", "polyCRACKER", "SubPhaser"]
fig, axes = plt.subplots(len(genomes), 4, figsize=(15, 3.35 * len(genomes)), squeeze=False)
for r, g in enumerate(genomes):
    d = next(x for x in rows if x["ref_id"] == g and x["reads"].lower() == "hifi")
    for c, m in enumerate(methods):
        M, tl, acc = matrix(d["key"], m)
        panel(axes[r, c], M, tl, f"{acc:.1f}% (HiFi)")
    axes[r, 0].annotate(label(d), xy=(-0.42, 0.5), xycoords="axes fraction", ha="center", va="center", rotation=90, fontsize=11)
for c, m in enumerate(methods):
    axes[0, c].annotate(m + (" (scaffold-first)" if m == "SubPhaser" else ""), xy=(0.5, 1.22), xycoords="axes fraction",
                        ha="center", fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0.03, 0, 1, 0.96]); fig.savefig(f"{figdir}/FigS4_confusion_all.png", dpi=150, bbox_inches="tight"); plt.close(fig)

pairs = [d for g in genomes for rt in ("hifi", "ont") for d in rows if d["ref_id"] == g and d["reads"].lower() == rt]
fig, axes = plt.subplots(len(pairs), 2, figsize=(8.5, 2.85 * len(pairs)), squeeze=False)
for r, d in enumerate(pairs):
    for c, m in enumerate(["PolySplit-Ref", "PolySplit-Allo"]):
        M, tl, acc = matrix(d["key"], m)
        panel(axes[r, c], M, tl, f"{m}: {acc:.1f}% ({d['reads']})")
    axes[r, 0].annotate(label(d), xy=(-0.45, 0.5), xycoords="axes fraction", ha="center", va="center", rotation=90, fontsize=10)
fig.tight_layout(rect=[0.04, 0, 1, 1]); fig.savefig(f"{figdir}/FigS5_confusion_chem.png", dpi=150, bbox_inches="tight")
