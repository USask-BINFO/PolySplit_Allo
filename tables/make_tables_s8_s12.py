#!/usr/bin/env python3
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
from read_accounting import read_config

rows = read_config(sys.argv[1]); O = sys.argv[2]
SG = {"A": "S1 (A)", "C": "S2 (C)"}
f = lambda x: f"{x:,}"
pc = lambda a, b: f"{100 * a / b:.2f}" if 100 * a / b < 1 else f"{100 * a / b:.1f}"
gn = lambda d: "*" + " ".join(d["name"].split()[:2]) + "* " + " ".join(d["name"].split()[2:]).replace("x)", "×)")
J = lambda sub, d: json.load(open(f"{O}/{sub}/{d['key']}.json"))
out = []; P = out.append

P("### Table S8. Read accounting, truth labels and accuracy by mapping quality\n")
P("Every input read falls into exactly one class. *Truth-labelled* reads align to the anchored chromosomes; their truth "
  "subgenome is the one receiving the most matched bases over all their alignments (primary and secondary). No MAPQ "
  "filter is applied, so these are exactly the reads used to compute accuracy. Reads that align only to unplaced "
  "scaffolds, or not at all, have no truth label and are not scored. Percentages are of all input reads.\n")
P("| Genome | reads | input reads | truth-labelled (used for accuracy) | by subgenome | only unplaced scaffolds | unaligned |")
P("|---|---|---|---|---|---|---|")
for d in rows:
    x = J("reads", d); T = x["total"]
    bys = "; ".join(f"{SG.get(k, k)} {f(v)}" for k, v in sorted(x["truth_by_subg"].items()))
    P(f"| {gn(d)} | {d['reads']} | {f(T)} | {f(x['truth_reads'])} ({pc(x['truth_reads'], T)}%) | {bys} | "
      f"{f(x['unplaced_only'])} ({pc(x['unplaced_only'], T)}%) | {f(x['unmapped'])} ({pc(x['unmapped'], T)}%) |")
P("\n**Reads kept despite low mapping quality, and accuracy with and without them.** *MAPQ<20* counts truth-labelled "
  "reads whose primary alignment has MAPQ below 20; they are included in every reported accuracy. The last columns give "
  "accuracy on the truth-labelled reads with primary MAPQ ≥ 20 and < 20, and accuracy under the stricter rule "
  "'subgenome of the primary alignment, MAPQ ≥ 20 only'.\n")
P("| Genome | reads | MAPQ<20 by subgenome | method | accuracy (reported) | MAPQ ≥ 20 | MAPQ < 20 | primary MAPQ ≥ 20 rule |")
P("|---|---|---|---|---|---|---|---|")
for d in rows:
    x = J("reads", d)
    lq = "; ".join(f"{SG.get(k, k)} {f(v)}" for k, v in sorted(x["mapq_lt20_by_subg"].items()))
    for m in ["PolySplit-Allo", "PolySplit-Ref"]:
        s = x["methods"][m]
        P(f"| {gn(d)} | {d['reads']} | {lq} | {m} | {s['T_sum']['acc']:.1f} | {s['T_sum_mapq_ge20']['acc']:.1f} | "
          f"{s['T_sum_mapq_lt20']['acc']:.1f} | {s['T_prim20']['acc']:.1f} |")
P("\n**What the reads without truth are.** Unplaced scaffolds were aligned to the *B. napus* and *C. sativa* chloroplast "
  "and mitochondrial genomes (NC_016734.1, NC_008285.1, NC_029337.1, PQ165104.1; minimap2 asm10, ≥90% identity); a "
  "scaffold is *organellar* if ≥50% of it is covered. *High-copy* scaffolds have read depth ≥3× the chromosome median. "
  "In brackets: WindowMasker-masked (repeat) fraction of the 10-kb reference window around the read's primary alignment.\n")
P("| Genome | reads | truth, MAPQ ≥ 20 | truth, MAPQ < 20 | truth via secondary alignment (primary on unplaced) | no truth: organellar | no truth: high-copy | no truth: other unplaced | unaligned |")
P("|---|---|---|---|---|---|---|---|---|")
for d in rows:
    u = J("unlabelled", d); e = u["categories"]
    cell = lambda k: f"{e[k]['pct_of_all']:.1f}% ({e[k]['mean_repeat_frac']:.2f})" if k in e else "0"
    P(f"| {gn(d)} | {d['reads']} | {cell('truth, MAPQ>=20')} | {cell('truth, MAPQ<20')} | {cell('truth, primary on unplaced')} | "
      f"{cell('no truth: organellar scaffold')} | {cell('no truth: high-copy scaffold')} | {cell('no truth: other unplaced')} | "
      f"{pc(u['unmapped'], u['total'])}% |")

P("\n### Table S9. Read accuracy in conserved and divergent regions\n")
P("Each anchored chromosome was cut into 10-kb windows and each window aligned to the chromosomes of the other "
  "subgenome(s) (minimap2 asm20). Homoeolog identity = length-weighted identity of these alignments; windows with less "
  "than half of their bases aligned have no detected homoeolog (divergent or subgenome-specific sequence). A read takes "
  "the value of the window under the midpoint of its primary alignment. Cells give accuracy (%) of the truth-labelled "
  "reads in each bin; *% reads* is the share of truth-labelled reads, *% of MAPQ<20* the share of low-MAPQ reads in the bin.\n")
BINS = ["no homoeolog detected", "<92%", "92-95%", "95-97%", "97-99%", ">=99%", "primary not on anchored chromosome"]
P("| Genome | reads | homoeolog identity | % reads | % of MAPQ<20 | PolySplit-Allo | PolySplit-Ref | polyCRACKER | SubPhaser |")
P("|---|---|---|---|---|---|---|---|---|")
for d in rows:
    x = J("reads", d); tot = sum(x["ident_bin_counts"].values())
    lo = x["ident_bin_counts_mapq_lt20"]; lt = sum(lo.values())
    for b in BINS:
        cells = [f"{s['acc']:.1f}" if (s := x["methods"].get(m, {}).get("by_ident_bin", {}).get(b)) else "–"
                 for m in ["PolySplit-Allo", "PolySplit-Ref", "polyCRACKER", "SubPhaser"]]
        lab = b.replace(">=", "≥").replace("primary not on anchored chromosome", "primary on unplaced scaffold")
        P(f"| {gn(d)} | {d['reads']} | {lab} | {pc(x['ident_bin_counts'].get(b, 0), tot)} | {pc(lo.get(b, 0), lt)} | " + " | ".join(cells) + " |")

P("\n### Table S10. Validation of shared-33-mer homoeolog edges against the reference\n")
P("Each contig pair is classified with the reference: *syntenic homoeolog* (different subgenomes, and a 10-kb window of "
  "one contig has its best cross-subgenome alignment inside the other contig's reference interval), *same subgenome* "
  "(paralogs, shared repeats, duplicated haplotigs) or *unanchored* (a contig that is not a pure, anchored contig). "
  "Percentages are of pairs whose two contigs are anchored. The last block is the edge set PolySplit-Allo uses: "
  "w ≥ τ_H = 10^5 and the strongest-partner rule (an edge must be at least half as strong as the strongest edge of each "
  "of its two contigs).\n")
P("| Genome | reads | w ≥ 10^4: edges | same subgenome | w ≥ 10^5: edges | syntenic homoeolog | same subgenome | retained (τ_H + rule): edges | syntenic homoeolog | median w (min–max) |")
P("|---|---|---|---|---|---|---|---|---|---|")
for d in rows:
    sw = {(s["thresh"], s["strongest_partner_rule"]): s for s in J("edges", d)["sweep"]}
    def fr(s, k):
        a = sum(v for kk, v in s["classes"].items() if kk != "unanchored")
        return f"{100 * s['classes'].get(k, 0) / a:.1f}%" if a else "–"
    s4, s5, sr = sw[(10000, False)], sw[(100000, False)], sw[(100000, True)]
    P(f"| {gn(d)} | {d['reads']} | {f(s4['n'])} | {fr(s4, 'same subgenome')} | {f(s5['n'])} | {fr(s5, 'syntenic homoeolog')} | "
      f"{fr(s5, 'same subgenome')} | {sr['n']} | {fr(sr, 'syntenic homoeolog')} | {f(sr['w_median'])} ({f(sr['w_min'])}–{f(sr['w_max'])}) |")

P("\n### Table S11. Hi-C normalisation and sensitivity of the block and pairing thresholds\n")
P("**Length normalisation.** Contig pairs with both contigs ≥50 kb, pure and anchored. Spearman ρ between the contact "
  "measure and √(ℓ_a ℓ_b); AUC = area under the receiver operating characteristic curve for telling same-chromosome "
  "from different-chromosome pairs; *top partner* = % of contigs whose strongest Hi-C partner lies on the same chromosome.\n")
P("| Genome | reads | pairs (same chrom.) | ρ raw / norm., same chrom. | ρ raw / norm., different chrom. | AUC raw / norm. | top partner raw / norm. |")
P("|---|---|---|---|---|---|---|")
for d in rows:
    h = J("edges", d)["hic"]
    P(f"| {gn(d)} | {d['reads']} | {f(h['pairs'])} ({f(h['cis_pairs'])}) | {h['spearman_raw_vs_len[same-chromosome]']:.2f} / "
      f"{h['spearman_norm_vs_len[same-chromosome]']:.2f} | {h['spearman_raw_vs_len[different-chromosome]']:.2f} / "
      f"{h['spearman_norm_vs_len[different-chromosome]']:.2f} | {h['auc_same_vs_diff_chrom_raw']:.2f} / "
      f"{h['auc_same_vs_diff_chrom_norm']:.2f} | {h['top_partner_same_chrom_raw']:.1f} / {h['top_partner_same_chrom_norm']:.1f} |")
if os.path.exists(f"{O}/sensitivity/summary.json"):
    by = {(r["key"], r["setting"]): r for r in json.load(open(f"{O}/sensitivity/summary.json"))}
    P("\n**Threshold sensitivity.** PolySplit-Allo re-run from the Hi-C blocks onward on the same assemblies, Hi-C contacts "
      "and shared-33-mer edges, changing one setting at a time: a minimum raw contact count per contig pair (default: none) "
      "or τ_H (default 10^5). Cells: contig accuracy / read accuracy (%).\n")
    P("| Genome | reads | default | ≥2 contacts | ≥5 contacts | ≥10 contacts | τ_H = 5×10^4 | τ_H = 2×10^5 |")
    P("|---|---|---|---|---|---|---|---|")
    for d in rows:
        cells = [f"{r['contig_acc']:.1f} / {r['read_acc']:.1f}" if (r := by.get((d["key"], s))) else "–"
                 for s in ["default", "contacts2", "contacts5", "contacts10", "tau5e4", "tau2e5"]]
        P(f"| {gn(d)} | {d['reads']} | " + " | ".join(cells) + " |")

P("\n### Table S12. Why scaffold-first separation failed (HiFi assemblies)\n")
P("YaHS scaffolds of the same Flye assembly used by PolySplit-Allo. *Mixed* = scaffold ≥10 Mb with ≥5% of its anchored "
  "length from another subgenome. *Correct groups* = homoeolog groups built without a reference (scaffold self-alignment) "
  "that hold one scaffold per subgenome. The last columns give read accuracy of SubPhaser with these automatic groups "
  "(Table II) and with the true homoeologous groups taken from the reference, next to PolySplit-Allo.\n")
P("| Genome | scaffolds ≥10 Mb | mixed (Mb) | correct automatic groups | SubPhaser, automatic groups | SubPhaser, true groups | PolySplit-Allo |")
P("|---|---|---|---|---|---|---|")
for d in [d for d in rows if d["reads"].lower() == "hifi"]:
    sf = f"{O}/scaffold_first/{d['key']}"
    if not os.path.exists(f"{sf}/mixed.json"):
        continue
    mx, gr, m = json.load(open(f"{sf}/mixed.json")), json.load(open(f"{sf}/groups.json")), J("reads", d)["methods"]
    tg = m.get("SubPhaser (true groups)", {}).get("T_sum", {}).get("acc")
    P(f"| {gn(d)} | {mx['scaffolds_ge10Mb']} | {mx['mixed_ge5pct']} ({mx['mixed_bp_Mb']}) | {gr['correct']} / {gr['groups']} | "
      f"{m['SubPhaser']['T_sum']['acc']:.1f} | {'–' if tg is None else f'{tg:.1f}'} | {m['PolySplit-Allo']['T_sum']['acc']:.1f} |")
print("\n".join(out))
