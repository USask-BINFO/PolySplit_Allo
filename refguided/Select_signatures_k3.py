#!/usr/bin/env python3
import argparse
import numpy as np
import random


def fasta_names_in_order(fa_path):
    names = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].strip().split()[0])
    return names


def read_fai_lengths(fai_path):
    lens = {}
    with open(fai_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            lens[fields[0]] = int(fields[1])
    return lens


def reservoir_add(res, seen_count, item, target, rng):
    seen_count += 1
    if len(res) < target:
        res.append(item)
    else:
        j = rng.randrange(seen_count)
        if j < target:
            res[j] = item
    return seen_count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--gsa", required=True)
    ap.add_argument("--lcp", required=True)
    ap.add_argument("--ac-fa", required=True, help="all reference chromosomes, FASTA order = doc order")
    ap.add_argument("--s1-fa", required=True)
    ap.add_argument("--s2-fa", required=True)
    ap.add_argument("--s3-fa", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--max-copy", type=int, default=3)
    ap.add_argument("--target", type=int, default=200000000000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    K = args.k
    rng = random.Random(args.seed)

    names_all = fasta_names_in_order(args.ac_fa)
    grp_names = [set(fasta_names_in_order(args.s1_fa)),
                 set(fasta_names_in_order(args.s2_fa)),
                 set(fasta_names_in_order(args.s3_fa))]
    lengths = read_fai_lengths(args.fai)

    gsa = np.memmap(args.gsa, dtype=[("doc", "<u4"), ("pos", "<u8")], mode="r")
    lcp = np.memmap(args.lcp, dtype=np.uint8, mode="r")
    n = len(gsa)
    assert len(lcp) == n

    max_doc = int(gsa["doc"].max())
    doc_group = np.full(max_doc + 1, -1, dtype=np.int8)
    for doc_id, name in enumerate(names_all):
        if doc_id > max_doc:
            break
        for g in range(3):
            if name in grp_names[g]:
                doc_group[doc_id] = g
                break

    res = [[], [], []]                 # signature reservoirs per subgenome
    seen = [0, 0, 0]
    cnt = [0, 0, 0]                     # copies of current block per subgenome
    rep = [None, None, None]           # representative (contig, pos0) per subgenome

    def flush_block():
        for g in range(3):
            others = sum(cnt[h] for h in range(3) if h != g)
            if cnt[g] > 0 and others == 0 and cnt[g] <= args.max_copy and rep[g] is not None:
                contig, pos0 = rep[g]
                seen[g] = reservoir_add(res[g], seen[g],
                                        (contig, pos0, pos0 + K, cnt[g]), args.target, rng)
        for g in range(3):
            cnt[g] = 0
            rep[g] = None

    for i in range(n):
        if i > 0 and lcp[i] < K:
            flush_block()
        doc = int(gsa[i]["doc"])
        grp = doc_group[doc] if doc <= max_doc else -1
        if grp < 0:
            continue
        if doc >= len(names_all):
            continue
        contig = names_all[doc]
        pos0 = int(gsa[i]["pos"])
        clen = lengths.get(contig, None)
        if clen is None or pos0 > clen - K:
            continue
        cnt[grp] += 1
        if rep[grp] is None:
            rep[grp] = (contig, pos0)
    flush_block()

    print("DONE  k =", K)
    for g, tag in enumerate(["S1", "S2", "S3"]):
        bed = f"{args.out_prefix}.{tag}.bed"
        with open(bed, "w") as out:
            for j, (ctg, s, e, copy) in enumerate(res[g]):
                out.write(f"{ctg}\t{s}\t{e}\t{tag}_sig_{j:07d}|copy={copy}\n")
        print(f"{tag} candidates seen: {seen[g]}  selected: {len(res[g])} -> {bed}")


if __name__ == "__main__":
    main()
