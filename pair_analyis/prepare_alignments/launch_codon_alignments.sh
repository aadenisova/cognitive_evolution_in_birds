#!/bin/bash
#SBATCH --job-name=alignOG
#SBATCH --output=logs/align_%A_%a.out
#SBATCH --error=logs/align_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=100
#SBATCH --mem=8G
#SBATCH --array=1-100

# OG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" og_ids.txt)

# cd alignments

# source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

# conda activate mafft
# mafft --auto ../per_OG_cleaned/${OG}.prot.fa > ${OG}.prot.aln.fa

# pal2nal.pl \
#     ${OG}.prot.aln.fa \
#     ../per_OG_cleaned/${OG}.cds.fa \
#     -output fasta \
#     > ${OG}.codon_aln.fa

PATH_TO_PAIRS=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/

# Папка с логами
mkdir -p logs

SLURM_ARRAY_TASK_ID=1
# список всех CDS файлов
ls ${PATH_TO_PAIRS}/genome_data/ncbi_dataset/data/per_OG/*.cds.fa > filelist_cds.txt

# выбираем файл, соответствующий номеру задачи
cds=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filelist_cds.txt)
rm filelist_cds.txt

echo "[INFO] Task $SLURM_ARRAY_TASK_ID running for $cds"

# Запускаем обработку
. ${PATH_TO_PAIRS}/src/pair_analyis/get_codon_alignmnets.sh "$cds"