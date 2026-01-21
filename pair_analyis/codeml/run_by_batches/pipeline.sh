
PATH_TO_SRC=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/pair_analyis/codeml/run_by_batches/

python3 $PATH_TO_SRC/prep_batches.py
sbatch $PATH_TO_SRC/run_array.sh