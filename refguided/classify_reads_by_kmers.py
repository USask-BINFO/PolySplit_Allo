#!/usr/bin/env python3
import argparse
import gzip
import sys

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def iter_reads_auto(path):
    with open_text(path) as f:
        first = f.readline()
        if not first:
            return

        # FASTA
        if first.startswith(">"):
            rid = first[1:].strip().split()[0]
            chunks = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    yield rid, "".join(chunks)
                    rid = line[1:].strip().split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if rid is not None:
                yield rid, "".join(chunks)

        # FASTQ (assumes 4-line records; standard SRA FASTQ is like this)
        elif first.startswith("@"):
            h = first
            while True:
                seq = f.readline()
                plus = f.readline()
                qual = f.readline()
                if not qual:
                    return
                rid = h.strip().split()[0][1:]
                yield rid, seq.strip()
                h = f.readline()
                if not h:
                    return
        else:
            raise ValueError("Unknown format: expected FASTA '>' or FASTQ '@' as first character.")

def build_base_table():
    tbl = [-1] * 256
    for ch, v in [(b"A",0),(b"C",1),(b"G",2),(b"T",3),(b"a",0),(b"c",1),(b"g",2),(b"t",3)]:
        tbl[ch[0]] = v
    return tbl

def kmer_canon_code(seq, k, tbl):
    if len(seq) != k:
        return None
    bseq = seq.encode("ascii", "ignore")
    mask = (1 << (2*k)) - 1
    shift = 2 * (k - 1)
    fwd = 0
    rev = 0
    for ch in bseq:
        b = tbl[ch]
        if b < 0:
            return None
        fwd = ((fwd << 2) | b) & mask
        rev = (rev >> 2) | ((3 - b) << shift)
    return fwd if fwd < rev else rev

def load_signature_codes(path, k, tbl):
    S = set()
    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            # allow either: pure k-mer lines, or fasta sequences lines
            if len(line) < k:
                continue
            if len(line) != k:
                # if someone gave longer lines, ignore (safer) or take first k
                line = line[:k]
            code = kmer_canon_code(line, k, tbl)
            if code is not None:
                S.add(code)
    return S

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("A_sig", help="A signatures (txt or fasta; one kmer per line)")
    ap.add_argument("C_sig", help="C signatures (txt or fasta; one kmer per line)")
    ap.add_argument("reads", help="Reads FASTA/FASTQ (optionally .gz)")
    ap.add_argument("k", type=int, help="k-mer length (e.g. 21)")
    ap.add_argument("out", help="output TSV")
    ap.add_argument("--min-hits", type=int, default=3)
    ap.add_argument("--ratio", type=float, default=3.0)
    ap.add_argument("--min-len", type=int, default=0, help="Skip scanning reads shorter than this; output 0/0 ambiguous")
    ap.add_argument("--stride", type=int, default=1, help="Check only every stride-th k-mer start (1 = all). Faster if >1.")
    ap.add_argument("--progress", type=int, default=100000, help="Print progress every N reads (0 disables)")
    args = ap.parse_args()

    k = args.k
    tbl = build_base_table()

    # Load signatures as canonical integer codes
    A = load_signature_codes(args.A_sig, k, tbl)
    C = load_signature_codes(args.C_sig, k, tbl)

    # Remove overlaps (non-discriminative)
    both = A & C
    if both:
        A -= both
        C -= both

    # Build one lookup dict (one hash lookup per k-mer instead of two)
    # 0 => A, 1 => C
    lab = {}
    for x in A:
        lab[x] = 0
    for x in C:
        lab[x] = 1

    mask = (1 << (2*k)) - 1
    shift = 2 * (k - 1)

    with open(args.out, "w", buffering=1<<20) as o:
        o.write("read_id\tA_hits\tC_hits\tlabel\n")

        n_reads = 0
        for rid, seq in iter_reads_auto(args.reads):
            n_reads += 1
            L = len(seq)

            # quick skip
            if L < k or (args.min_len and L < args.min_len):
                o.write(f"{rid}\t0\t0\tambiguous\n")
                if args.progress and (n_reads % args.progress == 0):
                    print("processed", n_reads, file=sys.stderr)
                continue

            bseq = seq.encode("ascii", "ignore")

            fwd = 0
            rev = 0
            valid = 0
            Ah = 0
            Ch = 0

            # scan with rolling codes
            for i, ch in enumerate(bseq):
                b = tbl[ch]
                if b < 0:
                    # reset on N or non-ACGT
                    fwd = 0
                    rev = 0
                    valid = 0
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
                    if g == 0:
                        Ah += 1
                    elif g == 1:
                        Ch += 1

            # decision rule
            label = "ambiguous"
            if Ah >= args.min_hits or Ch >= args.min_hits:
                if Ah >= args.ratio * (Ch + 1):
                    label = "A"
                elif Ch >= args.ratio * (Ah + 1):
                    label = "C"

            o.write(f"{rid}\t{Ah}\t{Ch}\t{label}\n")

            if args.progress and (n_reads % args.progress == 0):
                print("processed", n_reads, file=sys.stderr)

if __name__ == "__main__":
    main()
