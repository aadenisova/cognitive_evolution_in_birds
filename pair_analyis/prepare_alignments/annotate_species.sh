#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=liftoff_batch
#SBATCH --output=liftoff_batch_%j.out
#SBATCH --error=liftoff_batch_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G  
#SBATCH --cpus-per-task=16 
#SBATCH --mail-user=savouriess2112@gmail.com

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

# Функция для выполнения liftoff
run_liftoff() {
    local ref_gff=$1
    local ref_fasta=$2
    local target_fasta=$3
    local output_dir=$4
    local prefix=$5
    
    echo "=== Processing $prefix ==="
    
    # Активируем окружение liftoff
    conda activate liftOff
    
    # Определяем митохондриальную хромосому и исключаем её
    mito_chrom=$(grep '>' "$ref_fasta" | grep -i mito | sed 's/>//; s/,.*//' | awk '{print $1}')
    if [ -n "$mito_chrom" ]; then
        echo "Excluding mitochondrial chromosome: $mito_chrom"
        grep -v "^${mito_chrom}" "$ref_gff" > "${ref_gff%.*}_no_mito.gff"
        ref_gff_used="${ref_gff%.*}_no_mito.gff"
    else
        echo "No mitochondrial chromosome found, using original GFF"
        ref_gff_used="$ref_gff"
    fi
    
    # # Создаем выходную директорию
    # mkdir -p "$output_dir"
    
    # Запускаем liftoff
    liftoff -g "$ref_gff_used" \
        -p 16 \
        -o "$output_dir/lifted_annotations_genomic.gff" \
        -u "$output_dir/unmapped_features_$prefix.txt" \
        -infer_genes \
        -cds \
        "$target_fasta" "$ref_fasta"
    
    # Извлекаем features
    conda activate gffread
    
    gffread "$output_dir/lifted_annotations_genomic.gff" \
        -x "$output_dir/cds_sequences.fa" \
        -g "$target_fasta"
        
    gffread "$output_dir/lifted_annotations_genomic.gff" \
        -w "$output_dir/transcripts.fa" \
        -g "$target_fasta"
        
    gffread "$output_dir/lifted_annotations_genomic.gff" \
        -y "$output_dir/proteins.fa" \
        -g "$target_fasta"
    
    echo "=== Completed $prefix ==="
    echo ""
}

# Основной код
echo "Starting batch liftoff..."

# # 1. GCF_009819655.1 (Sylvia atricapilla) -> GCA_014839755.1 (Sylvia borin)
# run_liftoff \
#     "GCF_009819655.1/genomic.gff" \
#     "GCF_009819655.1/GCF_009819655.1_bSylAtr1.pri_genomic.fna" \
#     "GCA_014839755.1/GCA_014839755.1_bSylBor1.pri_genomic.fna" \
#     "GCA_014839755.1" \
#     "bSylBor1"

# # 2. GCF_009819795.1 (Aythya fuligula) -> GCA_965140915.1 (Aythya marila)
# run_liftoff \
#     "GCF_009819795.1/genomic.gff" \
#     "GCF_009819795.1/GCF_009819795.1_bAytFul2.pri_genomic.fna" \
#     "GCA_965140915.1/GCA_965140915.1_bAytMar2.hap1.1_genomic.fna" \
#     "GCA_965140915.1" \
#     "bAytMar2"

# # 3. GCF_009819795.1 (Aythya fuligula) -> GCA_964211825.1 (Aythya ferina)
# run_liftoff \
#     "GCF_009819795.1/genomic.gff" \
#     "GCF_009819795.1/GCF_009819795.1_bAytFul2.pri_genomic.fna" \
#     "GCA_964211825.1/GCA_964211825.1_bAytFer1.hap1.1_genomic.fna" \
#     "GCA_964211825.1" \
#     "bAytFer1"

# 4. GCF_009650955.1 (Corvus_moneduloides) -> GCA_014706295.1 (Lycocorax_pyrrhopterus)
# run_liftoff \
#     "GCF_009650955.1/genomic.gff" \
#     "GCF_009650955.1/GCF_009650955.1_bCorMon1.pri_genomic.fna" \
#     "GCA_014706295.1/GCA_014706295.1_ASM1470629v1_genomic.fna" \
#     "GCA_014706295.1" \
#     "ASM1470629v1"

# # 5. GCF_964199755.1 (Larus michahellis) -> GCA_964417175.1 (Larus_argentatus)
# run_liftoff \
#     "GCF_964199755.1/genomic.gff" \
#     "GCF_964199755.1/GCF_964199755.1_bLarMic1.1_genomic.fna" \
#     "GCA_964417175.1/GCA_964417175.1_bLarArg3.hap1.1_genomic.fna" \
#     "GCA_964417175.1" \
#     "bLarArg3"

# 6. GCF_964199755.1 (Larus michahellis) -> GCA_964417175.1 (Larus_fuscus)
run_liftoff \
    "GCF_964199755.1/genomic.gff" \
    "GCF_964199755.1/GCF_964199755.1_bLarMic1.1_genomic.fna" \
    "GCA_963932325.2/GCA_963932325.2_bLarFus1.hap2.2_genomic.fna" \
    "GCA_963932325.2" \
    "bLarFus1"


echo "All liftoff operations completed!"