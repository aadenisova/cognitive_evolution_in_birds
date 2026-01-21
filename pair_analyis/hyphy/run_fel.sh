#!/bin/bash
#SBATCH --job-name=absrel_array
#SBATCH --output=logs/absrel_%A_%a.out
#SBATCH --error=logs/absrel_%A_%a.err
#SBATCH --array=0-99        
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
conda activate hyphy

PATH_TO_PAIRS="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

ALIGN_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/alignments"
TREE_FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/pruned_12_species_analysis_noninno.nwk"
RESULT_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/fel_results_correct_noninno"
BATCH_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/batches"

mkdir -p "$RESULT_DIR"
mkdir -p logs

echo $TREE_FILE
# Определяем текущий batch-файл
# BATCH_FILE=$(printf "$BATCH_DIR/batch_%03d.txt" 1)
BATCH_FILE=$(printf "$BATCH_DIR/batch_%03d.txt" "$SLURM_ARRAY_TASK_ID")

echo ">>> Processing batch: $BATCH_FILE"

# Запускаем aBSREL для всех alignment'ов в batch-файле
while read aln; do

RESULT_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/fel_results_correct_noninno"

PATH_TO_PAIRS="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"
TREE_FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/pruned_12_species_analysis_inno.nwk"

"((((((Zonotrichia_albicollis:18.0594,Passer_domesticus:18.0594):9.83095,(Sylvia_borin:20.0422,Hirundo_rustica{TEST}:20.0422):7.84812):2.81251,(Lycocorax_pyrrhopterus:18.5451,Corvus_moneduloides{TEST}:18.5451):12.1577):28.0818,(Falco_rusticolus:4.37196,Falco_peregrinus{TEST}:4.37196):54.4126):4.12531,(Buteo_buteo{TEST}:9.95686,Haliaeetus_albicilla:9.95686):52.9531):2.18132,(Larus_fuscus:0.420387,Larus_argentatus:0.420387):64.6709);" > tree.nwk

aln=alignments/OG0001890.macse.clean.phy
og=$(basename "$aln" .macse.clean.phy)
out="${og}.json"

echo ">>> Running aBSREL on $og"

hyphy meme \
    --alignment "$aln" \
    --tree "$TREE_FILE" \
    --branches TEST \
    --output "$out" > "$out".meme

hyphy fel \
    --alignment "$aln" \
    --tree "$TREE_FILE" \
    --branches TEST \
    --output "$out" \
    --code Universal \
    --srv yes \
    --multiple-hits TRUE

done < "$BATCH_FILE"
