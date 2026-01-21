#!/bin/bash

#SBATCH --job-name=codeml_array
#SBATCH --array=0-99
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/codeml_%A_%a.out
#SBATCH --mail-user=savouriess2112@gmail.com

mkdir -p logs

# module load python/3.10
# module load paml

export PATH=$PATH:/ru-auth/local/home/adenisova/programs/paml/bin
PATH_TO_SRC=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/pair_analyis/codeml/run_by_batches/

BATCH=$(printf "batches/batch_%03d.txt" ${SLURM_ARRAY_TASK_ID})


SPS_INNO=(
Sylvia_borin
Passer_domesticus
Corvus_moneduloides
Falco_peregrinus
Larus_argentatus
Buteo_buteo
)


# Цикл по всем видам
for sp in "${SPS_INNO[@]}"; do
    # Создаем имя файла для каждого вида
    echo "Обработка вида: $sp"
    
    if [[ -f "$BATCH" ]]; then
        python3 $PATH_TO_SRC/run_batch.py "$BATCH" "$sp"
    else
        echo "Batch file not found: $BATCH"
    fi

done


