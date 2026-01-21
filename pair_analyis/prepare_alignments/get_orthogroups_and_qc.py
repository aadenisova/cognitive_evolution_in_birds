#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import SeqIO
import pandas as pd

DATA=sys.argv[1]  # Results_Nov25
PATH_TO_OG = f"longest_isoforms_proteomes_clean_for_5/OrthoFinder/Results_{DATA}"

OG_FILE = Path(f"{PATH_TO_OG}/Orthogroups/Orthogroups.txt")
OG_SC_FILE = Path(f"{PATH_TO_OG}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt")
PROT_ALL = Path("proteins_all.fa")
CDS_ALL = Path("cds_all.fa")

OUT_DIR = Path("per_OG_cleaned")
OUT_DIR.mkdir(exist_ok=True)

# ---- helpers ----

def species_from_id(seq_id: str) -> str:
    """
    Превращает что-то типа:
    Phalacrocorax_carbo_GENE12345
    XM_00123_Passer_domesticus
    в нормальный species name.

    Переделай под свои ID, если нужно!
    """
    parts = seq_id.split("_")
    return "_".join(parts[0:2])   # Phalacrocorax_carbo


# ---- load sequences ----
prot = {r.id: r for r in SeqIO.parse(PROT_ALL, "fasta")}
cds  = {r.id: r for r in SeqIO.parse(CDS_ALL,  "fasta")}

# ---- load single-copy OGs ----
sc_ogs = set()
for line in open(OG_SC_FILE):
    og = line.strip()
    sc_ogs.add(og)


qc_rows = []

# ---- parse OGs ----
for line in open(OG_FILE):
    if ":" not in line:
        continue

    og, rest = line.split(":", 1)
    og = og.strip()
    
    if og not in sc_ogs:
        continue
    
    ids = rest.strip().split()
    
    print(ids)

    # --- species → protein/CDS (guarantee unique species)
    species_prot = {}
    species_cds  = {}

    for id_ in ids:
        print(id_)
        sp = species_from_id(id_)
        print(sp)

        if id_ in prot and sp not in species_prot:
            species_prot[sp] = prot[id_]

        if id_ in cds and sp not in species_cds:
            species_cds[sp] = cds[id_]

    print(f'cds: {cds[id_]}')
    print(f'cds: {prot[id_]}')

    # QC
    miss_p = len(ids) - len(species_prot)
    miss_c = len(ids) - len(species_cds)

    prot_out = list(species_prot.values())
    cds_out  = list(species_cds.values())

    qc_rows.append({
        "OG": og,
        "species_in_OG": len(ids),
        "unique_species": len(species_prot),
        "missing_proteins": miss_p,
        "missing_cds": miss_c,
        "avg_protein_len": sum(len(r.seq) for r in prot_out)/len(prot_out) if prot_out else 0,
        "avg_cds_len": sum(len(r.seq) for r in cds_out)/len(cds_out) if cds_out else 0,
    })

    # Write
    if prot_out and cds_out:
        SeqIO.write(prot_out, OUT_DIR / f"{og}.prot.fa", "fasta")
        SeqIO.write(cds_out,  OUT_DIR / f"{og}.cds.fa",  "fasta")

# ---- write QC ----
pd.DataFrame(qc_rows).to_csv("OG_QC_stats.tsv", sep="\t", index=False)
print("Done.")
