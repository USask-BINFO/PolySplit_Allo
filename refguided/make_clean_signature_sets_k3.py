#!/usr/bin/env python3
import sys
from collections import Counter

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def canon(s):
    rc = revcomp(s)
    return rc if rc < s else s


def read_fasta_kmers(path):
    out = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            s = line.upper()
            if "N" in s:
                continue
            out.add(canon(s))
    return out


sets = [read_fasta_kmers(sys.argv[1]), read_fasta_kmers(sys.argv[2]), read_fasta_kmers(sys.argv[3])]
# count how many subgenome sets each canonical k-mer appears in; keep only count==1
seen = Counter()
for S in sets:
    for kmer in S:
        seen[kmer] += 1
removed = 0
clean = []
for g, S in enumerate(sets):
    keep = {kmer for kmer in S if seen[kmer] == 1}
    removed += len(S) - len(keep)
    clean.append(keep)

for g, out_path in enumerate(sys.argv[4:7]):
    with open(out_path, "w") as o:
        for s in sorted(clean[g]):
            o.write(s + "\n")
    print(f"S{g+1} clean: {len(clean[g])}")
print("removed non-unique (in >=2 subgenomes):", removed)
