#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=longest_isoforms_parallel
#SBATCH --output=longest_isoforms_parallel_%j.out
#SBATCH --error=longest_isoforms_parallel_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --mail-user=savouriess2112@gmail.com

data=$1

# #change_to_cfg
PATH_TO_FILE=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025

# source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

# # conda activate base
# # . ${PATH_TO_FILE}/src/pair_analyis/filter_protein_seqs.sh

# conda activate orthofinder

# # # предполагается, что в ./proteomes по одному файлу на вид
# orthofinder -f longest_isoforms_proteomes_clean_for_5 -t 64 -a 10

# # conda activate seqkit
# # # GCA_964417175.1_rna-XM_074588430.1
# # cat longest_isoforms_cds/*_cds.fa > cds_all.fa
# # cat longest_isoforms_proteomes_clean/*_proteins.fa > proteins_all.fa

# # # cat longest_isoforms_cds/*_cds.fa | seqkit rmdup -s > cds_all.fa
# # # cat longest_isoforms_proteomes_clean/*_proteins.fa | seqkit rmdup -s > proteins_all.fa

# # seqkit seq -i proteins_all.fa > proteins_all.indexed.fa
# # seqkit seq -i cds_all.fa > cds_all.indexed.fa

# conda activate base

# rm -r per_OG_cleaned
python $PATH_TO_FILE/src/pair_analyis/get_orthogroups_and_qc.py $data

# some_file.txt
# for file in per_OG_cleaned/*.prot.fa
# do
#     num=$(grep '>' $file | wc -l)
#     if [ "$num" -eq "12" ]; then
#         echo $file
#         echo $file $num >> some_file.txt
#     fi
# done

