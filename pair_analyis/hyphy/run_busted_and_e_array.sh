#!/bin/bash
#SBATCH --job-name=absrel_array
#SBATCH --output=logs/absrel_%A_%a.out
#SBATCH --error=logs/absrel_%A_%a.err
#SBATCH --partition=hpc_l40s
#SBATCH --array=60-99        
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# export PATH=$PATH:/ru-auth/local/home/adenisova/programs/hyphy

PATH_TO_PAIRS="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

TREE_FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/pruned_12_species_analysis_inno.nwk"

RESULT_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/busted_results_inno"
RESULT_DIR_E="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/busted_e_results_inno"

BATCH_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/batches"

mkdir -p "$RESULT_DIR"
mkdir -p "$RESULT_DIR_E"
mkdir -p logs

echo $TREE_FILE
# Определяем текущий batch-файл
# BATCH_FILE=$(printf "$BATCH_DIR/batch_%03d.txt" 1)
BATCH_FILE=$(printf "$BATCH_DIR/batch_%03d.txt" "$SLURM_ARRAY_TASK_ID")

echo ">>> Processing batch: $BATCH_FILE"

# Запускаем aBSREL для всех alignment'ов в batch-файле
while read aln; do
    og=$(basename "$aln" .macse.clean.phy)

    out="$RESULT_DIR/${og}.json"
    echo ">>> Running Busted on $og"
    hyphy busted \
        --alignment "$aln" \
        --tree "$TREE_FILE" \
        --branches TEST \
        --rates 3 \
        --grid-size 250 \
        --starting-points 5 \
        --output "$out" \
        --kill-zero-lengths Yes \
        --code Universal

    out="$RESULT_DIR_E/${og}.json"
    echo ">>> Running Busted E on $og"
    hyphy busted \
        --error-sink Yes \
        --alignment "$aln" \
        --tree "$TREE_FILE" \
        --branches TEST \
        --rates 3 \
        --grid-size 250 \
        --starting-points 5 \
        --output "$out" \
        --kill-zero-lengths Yes \
        --code Universal

done < "$BATCH_FILE"
