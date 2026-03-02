from pathlib import Path
from Bio import SeqIO
import pandas as pd
from sys import argv

PATH_TO_PROJECT="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

Prefix = argv[1]
cds_dir = argv[2]
alignment_file=argv[3]

SPECIES_TABLE = f"{PATH_TO_PROJECT}/src/whole_tree_analysis/cds/data/inno_aligned_with_TOGA.csv"

CDS_DIR = Path(f"{PATH_TO_PROJECT}/initial_data/{cds_dir}")
OUT_FASTA_DIR = Path(f"{PATH_TO_PROJECT}/tmp/{Prefix}_clean_fasta")
OUT_FASTA_DIR.mkdir(exist_ok=True)


gene = alignment_file.split(".fa")[0].split("/")[-1]
species_df = pd.read_csv(SPECIES_TABLE, sep=",")
# species_df = species_df.sample(frac=0.1, random_state=21)
# print(species_df.shape)

species_df.index = species_df["Assembly name"]
allowed_species = species_df["sci_name_2025"].to_dict()

def extract_species(header: str) -> str:
    if header.startswith("vs_"):
        return header.split()[0].replace("vs_", "")
    raise ValueError(f"Unexpected header: {header}")

def is_cds_valid(seq: str):
    seq = seq.replace("-", "").upper()
    length_ok = (len(seq) % 3 == 0)
    stop_codons = ["TAA", "TAG", "TGA"]
    stops = sum(seq[i:i+3] in stop_codons for i in range(0, len(seq)-2, 3))
    return length_ok, stops

def has_error(seq):
    return "!" in seq

records = list(SeqIO.parse(alignment_file, "fasta"))
species_map = dict()

for rec in records:
    if rec.id == "REFERENCE":
        continue
    try:
        sp = extract_species(rec.id)
    except ValueError:
        continue

    if sp in allowed_species:
        sp_name = allowed_species[sp]
        if sp_name in species_map:
            rec = min([species_map[sp_name], rec], key=lambda r: r.seq.count("-"))

        if not has_error(rec) and is_cds_valid(rec.seq)[0]:
            rec.id = sp_name
            rec.description = ""

            species_map[sp_name] = rec

kept_records = species_map.values()

fasta_out = OUT_FASTA_DIR / (gene + ".fa")
SeqIO.write(kept_records, fasta_out, "fasta")
