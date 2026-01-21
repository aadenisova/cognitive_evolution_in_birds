#!/usr/bin/env python3
from Bio import SeqIO
import argparse
from pathlib import Path

ALIGN_DIR = "alignments"

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

def rename_fasta(fasta_path, out_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    # write phylip
    with open(out_path, "w") as out:
        for r in records:
            r.id = SPECIES_NAMES['_'.join(r.id.split("_")[:2])]
            r.description = ""
        SeqIO.write(records, out_path, 'fasta')

# convert all *.fa files
ogs = sorted([f for f in Path(ALIGN_DIR).glob("*.prot.aln.fa") if f.stat().st_size > 0])

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

    out = fa.with_suffix(".renamed.fa")
    rename_fasta(fa, out)