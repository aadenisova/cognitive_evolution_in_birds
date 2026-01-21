#!/bin/bash
#SBATCH --job-name=meme
#SBATCH --output=logs/meme_%A_%a.out
#SBATCH --error=logs/meme_%A_%a.err
#SBATCH --partition=hpc_l40s
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-user=savouriess2112@gmail.com
#SBATCH --array=1-20 

# source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
# conda activate hyphy

# export PATH=$PATH:/ru-auth/local/home/adenisova/programs/hyphy

PATH_TO_PAIRS="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

ALIGN_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/alignments"
TREE_FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/pruned_12_species_analysis_inno.nwk"
RESULT_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/meme_results"
FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/OG_selected_absrel.txt"

# Создаем директорию для результатов
mkdir -p "$RESULT_DIR"

# Получаем общее количество строк в файле
TOTAL_LINES=$(wc -l < "$FILE")

# Проверяем, что номер задачи не превышает количество строк
if [ $SLURM_ARRAY_TASK_ID -gt $TOTAL_LINES ]; then
    echo "Номер задачи $SLURM_ARRAY_TASK_ID превышает количество OG в файле ($TOTAL_LINES)"
    exit 0  # Завершаем без ошибки
fi

# Получаем конкретный OG для этой задачи
og=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE")

# Формируем пути к файлам
aln="${ALIGN_DIR}/${og}.macse.clean.phy"
out="$RESULT_DIR/${og}.json"

echo "=== Обработка задачи $SLURM_ARRAY_TASK_ID: $og ==="
echo "Файл выравнивания: $aln"
echo "Выходной файл: $out"

# Проверяем существование файла
if [[ ! -f "$aln" ]]; then
    echo "ОШИБКА: Файл $aln не найден"
    exit 1
fi

# Запускаем MEME
echo ">>> Running MEME on $og"

hyphy meme \
    --alignment "$aln" \
    --tree "$TREE_FILE" \
    --code Universal \
    --branches TEST \
    --output "$out"

echo ">>> Завершено для $og"

# #!/bin/bash
# #SBATCH --job-name=meme
# #SBATCH --output=logs/meme_%A_%a.out
# #SBATCH --error=logs/meme_%A_%a.err
# #SBATCH --time=24:00:00
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=2G
# #SBATCH --mail-user=savouriess2112@gmail.com

# # source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
# # conda activate hyphy

# # export PATH=$PATH:/ru-auth/local/home/adenisova/programs/hyphy

# PATH_TO_PAIRS="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025"

# ALIGN_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/alignments"
# TREE_FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/pruned_12_species_analysis_inno.nwk"
# RESULT_DIR="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/meme_results"
# FILE="$PATH_TO_PAIRS/genome_data/ncbi_dataset/data/OG_selected_absrel.txt"

# mkdir -p "$RESULT_DIR"

# while read og; do
#     aln="${ALIGN_DIR}/${og}.macse.clean.phy"
#     if [[ -f "$aln" ]]; then
#         echo "Обрабатываю $og"
#         out="$RESULT_DIR/${og}.json"

#         echo ">>> Running MEME on $og"

#         hyphy meme \
#             --alignment "$aln" \
#             --tree "$TREE_FILE" \
#             --code Universal \
#             --branches TEST \
#             --output "$out"
#     else
#         echo "Файл $aln не найден"
#     fi
# done < $FILE