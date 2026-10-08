#!/usr/bin/env python3
import sys, os, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
from read_accounting import read_config


def identity(paf, nlines=2_000_000):
    out = subprocess.run(f"head -n {nlines} '{paf}' | cut -f13-22", shell=True, capture_output=True, text=True).stdout
    vals = []
    for ln in out.splitlines():
        tags = ln.split("\t")
        de = next((t[5:] for t in tags if t.startswith("de:f:")), None)
        if "tp:A:P" in tags and de is not None:
            vals.append(1 - float(de))
    return float(np.mean(vals))


def title(name):
    w = name.split()
    return r"$\it{" + r"\ ".join(w[:2]) + "}$ " + " ".join(w[2:]).replace("x)", r"$\times$)")


rows = read_config(sys.argv[1])
genomes = list(dict.fromkeys(d["ref_id"] for d in rows))
k = np.arange(5, 51)
fig, axes = plt.subplots(1, len(genomes), figsize=(4.7 * len(genomes), 4.2), sharey=True)
for ax, g in zip(np.atleast_1d(axes), genomes):
    for d in [d for d in rows if d["ref_id"] == g]:
        hifi = d["reads"].lower() == "hifi"; col = "#2b7bba" if hifi else "#c0392b"
        p = identity(f"{d['run_dir']}/reads_to_ref.paf"); s33 = 100 * p ** 33
        ax.plot(k, 100 * p ** k, color=col, lw=2.2, label=f"{d['reads']} ({100 * p:.1f}%)")
        ax.plot([33], [s33], "o", color=col)
        ax.annotate(f"{s33:.1f}%" if s33 < 10 else f"{s33:.0f}%", (34, s33), xytext=(34, s33 + (-10 if hifi else 8)), color=col, fontsize=10)
        ax.set_title(title(d["name"]), fontsize=12)
    ax.axvline(33, color="grey", ls=":", lw=1)
    ax.set_xlabel("$k$-mer length $k$"); ax.set_xlim(5, 50); ax.set_ylim(0, 100); ax.grid(alpha=0.3)
    ax.legend(loc="center right", fontsize=9)
np.atleast_1d(axes)[0].set_ylabel("error-free $k$-mers per read (%, = identity$^k$)")
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=150, bbox_inches="tight")
