#!/usr/bin/env python3
import argparse
import sys
import gzip


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def iter_reads_auto(path):
    with open_text(path) as f:
        first = f.read(1)
        f.seek(0)
        if first == ">":
            rid = None
            buf = []
            for line in f:
                if line.startswith(">"):
                    if rid is not None:
                        yield rid, "".join(buf)
                    rid = line[1:].strip().split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
            if rid is not None:
                yield rid, "".join(buf)
        else:
            while True:
                h = f.readline()
                if not h:
                    break
                seq = f.readline().strip()
                f.readline()
                f.readline()
                yield h[1:].strip().split()[0], seq


def build_base_table():
    tbl = [-1] * 256
    for ch, v in [(b"A", 0), (b"C", 1), (b"G", 2), (b"T", 3), (b"a", 0), (b"c", 1), (b"g", 2), (b"t", 3)]:
        tbl[ch[0]] = v
    return tbl


def kmer_canon_code(seq, k, tbl):
    if len(seq) != k:
        return None
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    fwd = rev = 0
    for ch in seq.encode("ascii", "ignore"):
        b = tbl[ch]
        if b < 0:
            return None
        fwd = ((fwd << 2) | b) & mask
        rev = (rev >> 2) | ((3 - b) << shift)
    return fwd if fwd < rev else rev


def load_codes(path, k, tbl):
    S = set()
    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            if len(line) < k:
                continue
            if len(line) != k:
                line = line[:k]
            c = kmer_canon_code(line, k, tbl)
            if c is not None:
                S.add(c)
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("S1_sig"); ap.add_argument("S2_sig"); ap.add_argument("S3_sig")
    ap.add_argument("reads"); ap.add_argument("k", type=int); ap.add_argument("out")
    ap.add_argument("--min-hits", type=int, default=3)
    ap.add_argument("--ratio", type=float, default=3.0)
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--progress", type=int, default=200000)
    args = ap.parse_args()

    k = args.k
    tbl = build_base_table()
    sets = [load_codes(args.S1_sig, k, tbl), load_codes(args.S2_sig, k, tbl), load_codes(args.S3_sig, k, tbl)]
    # drop any code present in >1 set (defensive; clean step should already ensure disjoint)
    from collections import Counter
    seen = Counter()
    for S in sets:
        for x in S:
            seen[x] += 1
    lab = {}
    for g, S in enumerate(sets):
        for x in S:
            if seen[x] == 1:
                lab[x] = g
    print(f"signature codes: S1={len(sets[0])} S2={len(sets[1])} S3={len(sets[2])}  unique-used={len(lab)}",
          file=sys.stderr)

    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    names = ["S1", "S2", "S3"]
    with open(args.out, "w", buffering=1 << 20) as o:
        o.write("read_id\tS1_hits\tS2_hits\tS3_hits\tlabel\n")
        n_reads = 0
        for rid, seq in iter_reads_auto(args.reads):
            n_reads += 1
            L = len(seq)
            if L < k or (args.min_len and L < args.min_len):
                o.write(f"{rid}\t0\t0\t0\tambiguous\n")
            else:
                fwd = rev = valid = 0
                h = [0, 0, 0]
                for i, ch in enumerate(seq.encode("ascii", "ignore")):
                    b = tbl[ch]
                    if b < 0:
                        fwd = rev = valid = 0
                        continue
                    fwd = ((fwd << 2) | b) & mask
                    rev = (rev >> 2) | ((3 - b) << shift)
                    valid += 1
                    if valid >= k:
                        start = i - k + 1
                        if args.stride > 1 and (start % args.stride != 0):
                            continue
                        code = fwd if fwd < rev else rev
                        g = lab.get(code)
                        if g is not None:
                            h[g] += 1
                order = sorted(range(3), key=lambda g: h[g], reverse=True)
                top, second = order[0], order[1]
                label = "ambiguous"
                if h[top] >= args.min_hits and h[top] >= args.ratio * (h[second] + 1):
                    label = names[top]
                o.write(f"{rid}\t{h[0]}\t{h[1]}\t{h[2]}\t{label}\n")
            if args.progress and (n_reads % args.progress == 0):
                print("processed", n_reads, file=sys.stderr)


if __name__ == "__main__":
    main()
