#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%j.log

OG="$1"
ALIGN="$2"
TREE="$3"
TEMPLATE_ALT="$4"
TEMPLATE_NULL="$5"

mkdir -p results/$OG/alt
mkdir -p results/$OG/null

python run_one_og.py "$OG" "$ALIGN" "$TREE" "$TEMPLATE_ALT" "$TEMPLATE_NULL"