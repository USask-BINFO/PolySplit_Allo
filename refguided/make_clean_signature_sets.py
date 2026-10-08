import sys

COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]

def canon(s: str) -> str:
    rc = revcomp(s)
    return rc if rc < s else s

def read_fasta_kmers(path):
    kmers = []
    with open(path) as f:
        seq = None
        for line in f:
            line = line.strip()
            if not line: 
                continue
            if line.startswith(">"):
                continue
            s = line.upper()
            kmers.append(s)
    return kmers

A_raw = read_fasta_kmers(sys.argv[1])
C_raw = read_fasta_kmers(sys.argv[2])

A = set()
for s in A_raw:
    if "N" in s: 
        continue
    A.add(canon(s))

C = set()
for s in C_raw:
    if "N" in s:
        continue
    C.add(canon(s))

# remove overlaps (anything appearing in both is not discriminative)
both = A & C
A -= both
C -= both

with open(sys.argv[3], "w") as out:
    for s in sorted(A):
        out.write(s + "\n")

with open(sys.argv[4], "w") as out:
    for s in sorted(C):
        out.write(s + "\n")

print("A clean:", len(A))
print("C clean:", len(C))
print("removed overlaps:", len(both))

