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
    ap.add_argument("--ac-fa", required=True)
    ap.add_argument("--a-fa", required=True)
    ap.add_argument("--c-fa", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--max-copy", type=int, default=5)
    ap.add_argument("--target", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    K = args.k
    rng = random.Random(args.seed)

    names_ac = fasta_names_in_order(args.ac_fa)
    names_A  = set(fasta_names_in_order(args.a_fa))
    names_C  = set(fasta_names_in_order(args.c_fa))
    lengths  = read_fai_lengths(args.fai)

    gsa = np.memmap(args.gsa, dtype=[("doc","<u4"), ("pos","<u8")], mode="r")
    lcp = np.memmap(args.lcp, dtype=np.uint8, mode="r")
    n = len(gsa)
    assert len(lcp) == n

    max_doc = int(gsa["doc"].max())
    doc_group = np.full(max_doc + 1, -1, dtype=np.int8)

    # Map doc_id -> contig name using FASTA order (0-based). Leave extras as -1.
    for doc_id, name in enumerate(names_ac):
        if doc_id > max_doc:
            break
        if name in names_A:
            doc_group[doc_id] = 0
        elif name in names_C:
            doc_group[doc_id] = 1

    # Reservoirs store tuples: (contig, start0, end0, copy)
    A_res, C_res = [], []
    seenA = 0
    seenC = 0

    # Current block state
    cntA = 0
    cntC = 0
    repA = None  # (contig, pos0)
    repC = None  # (contig, pos0)

    def flush_block():
        nonlocal seenA, seenC, cntA, cntC, repA, repC
        # A-only
        if cntA > 0 and cntC == 0 and cntA <= args.max_copy and repA is not None:
            contig, pos0 = repA
            seenA = reservoir_add(A_res, seenA, (contig, pos0, pos0 + K, cntA), args.target, rng)
        # C-only
        if cntC > 0 and cntA == 0 and cntC <= args.max_copy and repC is not None:
            contig, pos0 = repC
            seenC = reservoir_add(C_res, seenC, (contig, pos0, pos0 + K, cntC), args.target, rng)

        cntA = 0
        cntC = 0
        repA = None
        repC = None

    for i in range(n):
        if i > 0 and lcp[i] < K:
            flush_block()

        doc = int(gsa[i]["doc"])
        grp = doc_group[doc] if doc <= max_doc else -1
        if grp < 0:
            continue

        # doc -> contig name
        if doc < len(names_ac):
            contig = names_ac[doc]
        else:
            continue

        pos0 = int(gsa[i]["pos"])
        clen = lengths.get(contig, None)
        if clen is None:
            continue

        # only consider positions where we can take a full k-mer
        if pos0 > clen - K:
            continue

        if grp == 0:
            cntA += 1
            if repA is None:
                repA = (contig, pos0)
        else:
            cntC += 1
            if repC is None:
                repC = (contig, pos0)

    flush_block()

    # Write BED
    A_bed = args.out_prefix + ".A.bed"
    C_bed = args.out_prefix + ".C.bed"

    with open(A_bed, "w") as out:
        for j, (ctg, s, e, copy) in enumerate(A_res):
            out.write(f"{ctg}\t{s}\t{e}\tA_sig_{j:07d}|copy={copy}\n")

    with open(C_bed, "w") as out:
        for j, (ctg, s, e, copy) in enumerate(C_res):
            out.write(f"{ctg}\t{s}\t{e}\tC_sig_{j:07d}|copy={copy}\n")

    print("DONE")
    print("k =", K)
    print("A candidates seen:", seenA, "A selected:", len(A_res), "->", A_bed)
    print("C candidates seen:", seenC, "C selected:", len(C_res), "->", C_bed)

if __name__ == "__main__":
    main()

