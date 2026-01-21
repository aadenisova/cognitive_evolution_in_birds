#!/bin/bash
#SBATCH --job-name=absrel_array
#SBATCH --output=logs/absrel_%A_%a.out
#SBATCH --error=logs/absrel_%A_%a.err
#SBATCH --array=0-99        
#SBATCH --time=24:00:00
#SBATCH --partition=hpc_l40s
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
conda activate hyphy

PATH_TO_PAIRS="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

# ALIGN_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/alignments_5"
TREE_FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/pruned_12_species_analysis_inno.nwk"
RESULT_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/absrel_results_inno_12"
BATCH_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/batches_absrel"

mkdir -p "$RESULT_DIR"
mkdir -p logs

echo $TREE_FILE
# Определяем текущий batch-файл
# BATCH_FILE=$(printf "$BATCH_DIR/batch_%03d.txt" 1)
BATCH_FILE=$(printf "$BATCH_DIR/batch_%d.txt" "$SLURM_ARRAY_TASK_ID")

echo ">>> Processing batch: $BATCH_FILE"

# Запускаем aBSREL для всех alignment'ов в batch-файле
while read aln; do
    og=$(basename "$aln" .macse.clean.phy)
    out="$RESULT_DIR/${og}.json"

    echo ">>> Running aBSREL on $og"

    hyphy absrel \
        --alignment "$aln" \
        --tree "$TREE_FILE" \
        --branches TEST \
        --output "$out" \
        --code Universal
        --multiple-hits TRUE

done < "$BATCH_FILE"
