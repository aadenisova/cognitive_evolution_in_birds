#!/usr/bin/env python3
import sys
from Bio import SeqIO

STOP_CODONS = {"TAA", "TAG", "TGA"}

def internal_stop_codon(seq):
    seq = seq.replace("-", "")
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            continue
        if codon in STOP_CODONS:
            return True
    return False

def bad_gaps(seq):
    return seq.count("-") / len(seq) > 0.5

def has_macse_error(seq):
    return "!" in seq or "#" in seq

def main():
    if len(sys.argv) != 3:
        print("Usage: filter_macse.py input.fa output.fa")
        sys.exit(1)

    inp = sys.argv[1]
    out = sys.argv[2]

    records = list(SeqIO.parse(inp, "fasta"))
    cleaned = []

    for r in records:
        s = str(r.seq)

        if has_macse_error(s):
            print(f"[DROP] {r.id}: MACSE error (! or #)")
            continue

        if len(s) % 3 != 0:
            print(f"[DROP] {r.id}: length %3 != 0")
            continue

        if internal_stop_codon(s):
            print(f"[DROP] {r.id}: internal stop codon")
            continue

        if bad_gaps(s):
            print(f"[DROP] {r.id}: >50% gaps")
            continue

        cleaned.append(r)

    if len(cleaned) < 4:
        print("[WARN] Too few sequences after filtering")

    SeqIO.write(cleaned, out, "fasta")
    print(f"[OK] Written cleaned alignment to {out}")

if __name__ == "__main__":
    main()
