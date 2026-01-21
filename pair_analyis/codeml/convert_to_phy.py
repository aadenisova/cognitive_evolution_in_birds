#!/usr/bin/env python3
from Bio import SeqIO
import argparse
from pathlib import Path

ALIGN_DIR = "alignments_5"

parser = argparse.ArgumentParser()
parser.add_argument(
    "--strict", 
    action="store_true",
    help="Require presence of all species listed in SPECIES_NAMES"
)
args = parser.parse_args()

# Словарь для замены названий
SPECIES_NAMES = {
    "GCA_014839755.1": "Sylvia_borin",
    "GCF_015227805.2": "Hirundo_rustica", 
    # "GCA_965140915.1": "Aythya_marila",
    # "GCA_964211825.1": "Aythya_ferina",
    # "GCA_020800305.1": "Porphyrio_hochstetteri",
    # "GCA_964237585.1": "Gallinula_chloropus",
    "GCA_964417175.1": "Larus_argentatus",
    "GCA_963932325.2": "Larus_fuscus",
    "GCF_023634155.1": "Falco_peregrinus",
    "GCF_015220075.1": "Falco_rusticolus",
    "GCF_009650955.1": "Corvus_moneduloides",
    "GCA_014706295.1": "Lycocorax_pyrrhopterus",
    "GCF_036417665.1": "Passer_domesticus",
    "GCF_047830755.1": "Zonotrichia_albicollis",
    "GCF_947461875.1": "Haliaeetus_albicilla",
    "GCF_964188355.1": "Buteo_buteo",
}

def fasta_to_phylip(fasta_path, out_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    n = len(records)
    L = len(records[0].seq)

    # check codon alignment validity
    if any(len(r.seq) != L for r in records):
        raise ValueError("All sequences must be same length")
    if L % 3 != 0:
        raise ValueError("Alignment length not divisible by 3")
    
    # write phylip
    with open(out_path, "w") as out:
        out.write(f"{n} {L}\n")
        for r in records:
            # Извлекаем ID и заменяем на видовое название
            seq_id = r.id

            for gca_id, species_name in SPECIES_NAMES.items():
                if gca_id in seq_id:
                    seq_id = species_name
                    name = seq_id[:30].ljust(30)     # PAML требует ≤30 символов
                    seq = str(r.seq).replace("\n", "")
                    out.write(f"{name} {seq[:-3]}\n")      # Пробел между названием и последовательностью

    

# convert all *.fa files
ogs = sorted([f for f in Path(ALIGN_DIR).glob("*.macse.clean.fa") if f.stat().st_size > 0])

for fa in ogs:
    records = list(SeqIO.parse(fa, "fasta"))
    
    # Если включён strict-режим — проверяем число видов
    if args.strict:
        expected = len(SPECIES_NAMES)
        present = 0
        for r in records:
            if any(gca in r.id for gca in SPECIES_NAMES):
                present += 1

        if present != expected:
            print(f"[SKIP] {fa.name}: only {present}/{expected} species present")
            continue

    out = fa.with_suffix(".phy")
    fasta_to_phylip(fa, out)
    print("converted:", out)