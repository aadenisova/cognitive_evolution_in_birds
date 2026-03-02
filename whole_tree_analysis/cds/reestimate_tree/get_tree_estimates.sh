#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=tree_estimate
#SBATCH --output=logs/tree_estimate_%j.out
#SBATCH --error=logs/tree_estimate_%j.err
#SBATCH --nodes=1
#SBATCH --partition=hpc_l40_b 
#SBATCH --ntasks=1
#SBATCH --mem=4G  
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=savouriess2112@gmail.com
#SBATCH --array=1-100  

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh
conda activate iqtree

Prefix="TOGA_ALL"

PATH_TO_DATA=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/
tree=$PATH_TO_DATA/src/whole_tree_analysis/cds/data/${Prefix}_tree_aligned_with_TOGA.newick
alignment=$PATH_TO_DATA/tmp/${Prefix}_polished_phy/
stats=$PATH_TO_DATA/tmp/${Prefix}_polished_phy_stat
FOR_RER=$PATH_TO_DATA/tmp/${Prefix}_polished_for_RER

NEW_TREES=$PATH_TO_DATA/tmp/${Prefix}_new_tree
NEW_TREES_REEST=$PATH_TO_DATA/tmp/${Prefix}_new_tree_reest
mkdir -p $NEW_TREES $NEW_TREES_REEST $FOR_RER

GENES_TO_FILTER=$PATH_TO_DATA/tmp/genes_to_filter.txt

# ALL_FILES=($alignment/*.phy)
# TOTAL_FILES=${#ALL_FILES[@]}

TOTAL_FILES=$(wc -l $GENES_TO_FILTER | awk '{print $1}')
echo $TOTAL_FILES

FILES_PER_JOB=160
# SLURM_ARRAY_TASK_ID=1
START_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) * FILES_PER_JOB )) 
END_IDX=$(( START_IDX + FILES_PER_JOB - 1 ))

for (( i=START_IDX; i<=END_IDX && i<TOTAL_FILES; i++ )); do
  # ALIGN_FILE="${ALL_FILES[$i]}"
  # RELATIVE_PATH="${ALIGN_FILE#$alignment/}"
  # GENENAME="${RELATIVE_PATH%.phy}"
  i=$(($i + 1 ))
  # sed "${NUM}q;d" file
  GENENAME=$(sed $i'q;d' $GENES_TO_FILTER)

  # echo $GENENAME

  sp_names=$(awk '{print $4}' $stats/$GENENAME.txt)
  # echo $sp_names

  # echo $ALIGN_FILE

  conda activate base 
  new_tree_path=$NEW_TREES/$GENENAME.newick
  python $PATH_TO_DATA/src/whole_tree_analysis/cds/reestimate_tree/prune_tree.py \
    $tree \
    $sp_names \
    $new_tree_path

  conda activate iqtree

  cd $NEW_TREES_REEST
  iqtree2 \
    -s $alignment/$GENENAME.phy \
    --seqtype CODON \
    -te $new_tree_path \
    -m GY \
    -pre $GENENAME

  conda activate base 
  python $PATH_TO_DATA/src/whole_tree_analysis/cds/reestimate_tree/compare_topology.py $GENENAME.treefile $new_tree_path

  echo -e "$GENENAME\t$(cat $NEW_TREES_REEST/$GENENAME.treefile)" > $FOR_RER/$GENENAME.txt

done

echo "Задача $SLURM_ARRAY_TASK_ID завершена"