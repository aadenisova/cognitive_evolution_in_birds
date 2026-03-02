from Bio import SeqIO
from pathlib import Path
from sys import argv

PREFIX=argv[1] #"TOGA_ALL"
GENE_NAME=argv[2] #"ATP11A_rna-XM_015277887.2"

PATH = "/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp"

out_dir = f"{PATH}/{PREFIX}_polished_phy"
out_dir_stat = f"{PATH}/{PREFIX}_polished_phy_stat"
Path(out_dir).mkdir(exist_ok = True)
Path(out_dir_stat).mkdir(exist_ok = True)

# /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/tmp/TOGA_ALL_prank_phylip/ATP11A_rna-XM_015277887.2.phy
alignment_file = f"{PATH}/{PREFIX}_prank_fa/{GENE_NAME}.fa"
records = list(SeqIO.parse(alignment_file, "fasta"))

gap_threshold = 0.5

kept_records = []
sp_names = ""
for rec in records:
    seq = rec.seq
    gap_proportion = len([i for i in seq if i == "-"])/len(seq)
    print(f'{rec.id}: {gap_proportion}')

    if gap_proportion < gap_threshold:
        kept_records.append(rec)
        sp_names+=f"{rec.id},"
sp_names=sp_names.strip(",")

out_name = f"{out_dir}/{GENE_NAME}.phy"
SeqIO.write(kept_records, out_name, "phylip-relaxed")

print(f"Polished alignment saved to {out_name}")

out_name_stats = f"{out_dir_stat}/{GENE_NAME}.txt"
with open(out_name_stats, "w") as f:
    f.write(f"{GENE_NAME}\t{len(rec.seq)}\t{len(kept_records)}\t{sp_names}\n")

print(f"Stats saved to {out_name_stats}")