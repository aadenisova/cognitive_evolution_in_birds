#!/bin/bash
#SBATCH --job-name=absrel
#SBATCH --output=logs/absrel_%A_%a.out
#SBATCH --error=logs/absrel_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-user=savouriess2112@gmail.com

# source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
# conda activate hyphy

export PATH=$PATH:/ru-auth/local/home/adenisova/programs/hyphy

PATH_TO_PAIRS="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

ALIGN_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/alignments"
TREE_FILE="$PATH_TO_PAIRS/hyphy_tutorial/pruned_12_species_analysis_inno.nwk"
RESULT_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/busted_e_results_whole"

mkdir -p "$RESULT_DIR"

for aln in "$ALIGN_DIR"/*.macse.clean.phy; do
    og=$(basename "$aln" .macse.clean.phy)
    out="$RESULT_DIR/${og}.json"

    echo ">>> Running aBSREL on $og"

    # 
    hyphy busted \
        --error-sink Yes \
        --alignment "$aln" \
        --tree "$TREE_FILE" \
        --branches TEST \
        --save-fit "${out%.json}.fit" \
        --rates 3 \
        --output "$out" \
        --code Universal

    # hyphy absrel \
    #     --alignment "$aln" \
    #     --tree "$TREE_FILE" \
    #     --branches TEST \
    #     --output "$out" \
    #     --code Universal
done

hyphy remove-duplicates \
    --msa alignments/OG0001351.macse.clean.phy \
    --tree /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/genome_data/ncbi_dataset/data/pruned_12_species_analysis_inno_no_TEST.nwk \
    --output OG0001351.macse.clean.no_dups.phy \
    ENV="DATA_FILE_PRINT_FORMAT=2"

hyphy remove-duplicates \
    /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/genome_data/ncbi_dataset/data/alignments/OG0001351.macse.clean.phy 

# source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
# conda activate hyphy

# # Параметры
# PATH_TO_PAIRS=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/

# # /lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/hyphy_tutorial/pruned_12_species_analysis_inno.nwk
# TREE_FILE="${PATH_TO_PAIRS}/hyphy_tutorial/pruned_12_species_analysis_inno.nwk"          # ваше дерево с пометками
# ALIGNMENT_FILE="${PATH_TO_PAIRS}/genome_data/ncbi_dataset/data/alignments/OG0001342.macse.clean.phy"    # кодонное выравнивание
# OUTPUT_FILE="absrel_results.json"  # выходной файл

# # https://github.com/veg/hyphy-analyses/tree/master/remove-duplicates
# # Запуск aBSREL
# hyphy absrel \
#     --alignment $ALIGNMENT_FILE \
#     --tree $TREE_FILE \
#     --branches TEST \
#     --output $OUTPUT_FILE \
#     --code Universal

# Альтернативный вариант с более подробными настройками
# hyphy absrel \
#     --alignment $ALIGNMENT_FILE \
#     --tree $TREE_FILE \
#     --branches All \
#     --output $OUTPUT_FILE \
#     --code Vertebrate-mtDNA \
#     --p-value 0.05 \
#     --rates 3 \
#     --lrt Yes