#!/usr/bin/env python3
import sys, pickle

contacts, lengths = pickle.load(open(sys.argv[1], "rb"))
f = int(sys.argv[3])
kept = {p: n for p, n in contacts.items() if n >= f}
pickle.dump((kept, lengths), open(sys.argv[2], "wb"))
print(f"[contacts >= {f}] kept {len(kept):,} of {len(contacts):,} contig pairs", file=sys.stderr)
