#!/usr/bin/env python3
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
from read_accounting import read_config

CL = [("syntenic homoeolog", "#2166ac", "syntenic homoeologs (different subgenomes)"),
      ("homoeologous chromosomes, outside mapped windows", "#67a9cf", "homoeologous chromosomes, other region"),
      ("other cross-subgenome", "#d1e5f0", "other cross-subgenome"),
      ("same subgenome", "#d6604d", "same subgenome (paralog / repeat)"),
      ("unanchored", "#bababa", "unanchored contig")]
rows = read_config(sys.argv[1]); outdir = sys.argv[2]
genomes = list(dict.fromkeys(d["ref_id"] for d in rows))
bins = np.logspace(2, 6.6, 47)
fig, axes = plt.subplots(len(genomes), 2, figsize=(11, 3.7 * len(genomes)), sharex=True, squeeze=False)
for r, g in enumerate(genomes):
    for c, rt in enumerate(["ont", "hifi"]):
        ax = axes[r, c]
        d = next((x for x in rows if x["ref_id"] == g and x["reads"].lower() == rt), None)
        if d is None:
            ax.axis("off"); continue
        w = json.load(open(f"{outdir}/edges/{d['key']}.json"))["weights_by_class"]
        ax.hist([np.array(w.get(k, []), float) for k, _, _ in CL], bins=bins, stacked=True, color=[x[1] for x in CL],
                label=[x[2] for x in CL], edgecolor="white", linewidth=0.3)
        ax.axvline(1e5, color="k", ls="--", lw=1.3, label=r"$\tau_H=10^5$")
        ax.set_xscale("log"); ax.set_yscale("log")
        nm = d["name"].split()
        ax.set_title(r"$\it{" + r"\ ".join(nm[:2]) + "}$ " + " ".join(nm[2:]).replace("x)", r"$\times$)") + f", {d['reads']}", fontsize=11)
        if c == 0:
            ax.set_ylabel("contig pairs")
        if r == len(genomes) - 1:
            ax.set_xlabel(r"shared 33-mers $w_{\mathrm{kmer}}$ (33-mers occurring 2 to 6 times in the assembly)")
h, l = axes[0, 1].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.045))
fig.tight_layout(); fig.savefig(sys.argv[3], dpi=150, bbox_inches="tight")
