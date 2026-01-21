#!/usr/bin/env python3
import sys
from Bio import SeqIO

# UAA, UAG, and UGA
STOP_CODONS = {"TAA", "TAG", "TGA"}

def is_valid_seq(seq):
    seq = seq.upper().replace("-", "N")
    return all(c in "ACGTN" for c in seq)

def has_internal_stop(seq):
    seq = seq.upper().replace("-", "N")
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            continue
        if "N" not in codon and codon in STOP_CODONS:
            return True
    return False

def gap_fraction(seq):
    return seq.count("-") / len(seq)

def main():
    if len(sys.argv) != 2:
        print("Usage: check_cds_qc.py alignment.fa", file=sys.stderr)
        sys.exit(2)

    filename = sys.argv[1]
    bad = False

    for record in SeqIO.parse(filename, "fasta"):
        name = record.id
        seq = str(record.seq)

        # 1) длина должна быть кратна 3
        if len(seq) % 3 != 0:
            print(f"[ERROR] {name}: length {len(seq)} not divisible by 3")
            bad = True

        # 2) странные символы
        if not is_valid_seq(seq):
            print(f"[ERROR] {name}: invalid characters in sequence")
            bad = True

        # 3) внутренние стоп-кодоны
        if has_internal_stop(seq):
            print(f"[ERROR] {name}: internal stop codon detected")
            bad = True

        # 4) слишком много гэпов (если нужно)
        if gap_fraction(seq) > 0.5:
            print(f"[WARN] {name}: >50% gaps — likely bad sequence")

    if bad:
        sys.exit(1)
    else:
        print("[OK] QC passed")

if __name__ == "__main__":
    main()
