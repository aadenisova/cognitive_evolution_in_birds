#!/bin/bash

## FILTER ALIGNMENTS 
prefix="TOGA_ALL_selected"
# to get filtered alignments from TOGA
sbatch prepare_alignment/alignment_filtration.sh
# to concat all stats in single file
cat ${prefix}_polished_phy_stat/* > ${prefix}_polished_phy_stat.tsv
# to check filtered alignments + get genes_to_filter.txt for tree estimation
prepare_alignment/analyse_polished_align.ipynb

## REESTIMATE TREES
# to get new trees for each gene (genes_to_filter.txt is used to run for some genes)
sbatch reestimate_tree/get_tree_estimates.sh

## RERConverge
# Adapt tree format
RERConverge/last_prep.ipynb
# Find correlation between RER and traits
RERConverge/find_correlation.R

## RERConverge analysis
RERConverge/analyse_RER.ipynb