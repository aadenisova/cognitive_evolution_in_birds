#!/bin/bash
#SBATCH --job-name=ortho_align
#SBATCH --output=ortho_align_%j.log
#SBATCH --error=ortho_align_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=24:00:00

set -e
PATH_TO_SRC=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/src/pair_analyis/
PATH_TO_MAC=/ru-auth/local/home/adenisova/programs/

# Директории

OGdir="per_OG_cleaned"
ALIGNdir="alignments_5"
ALIGNdir_final="alignments_final_5"

mkdir -p "$ALIGNdir"
mkdir -p "$ALIGNdir_final"

# Функция для обработки одного OG
process_og() {
    local cds="$1"
    local og=$(basename "$cds" .cds.fa)
    local prot="$OGdir/$og.prot.fa"

    # Если нет protein файла — пропускаем
    if [ ! -f "$prot" ]; then
        echo "[WARN] No protein file for $og" >&2
        return 0
    fi

    source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

    conda activate mafft
    # Выравнивание аминокислот
    mafft --quiet --auto "$prot" > "$ALIGNdir/$og.prot.aln.fa"

    # # # Кодонное выравнивание
    # pal2nal.pl \
    #     "$ALIGNdir/$og.prot.aln.fa" \
    #     "$cds" \
    #     -output fasta \
    #     > "$ALIGNdir/$og.codon_aln.fa"

    # java -jar ${PATH_TO_MAC}/macse_v2.07.jar -prog alignSequences \
    #     -seq_ref "$prot" \
    #     -seq "$cds" \
    #     -out_NT "$ALIGNdir/$og.macse.fa"

    java -jar ${PATH_TO_MAC}/macse_v2.07.jar -prog reportGapsAA2NT \
        -align_AA "$ALIGNdir/$og.prot.aln.fa" \
        -seq "$cds" \
        -out_NT "$ALIGNdir/$og.macse.fa"

    # mafft --quiet --auto  per_OG_cleaned/OG0001704.prot.fa >  per_OG_cleaned/OG0001704.prot.aln.fa  
    # java -jar ${PATH_TO_MAC}/macse_v2.07.jar -prog reportGapsAA2NT \
    #     -align_AA per_OG_cleaned/OG0001704.prot.aln.fa \
    #     -seq per_OG_cleaned/OG0001704.cds.fa \
    #     -out_NT OG0001704.macse.fa

    # per_OG_cleaned/OG0001704.prot.fa

    # QC: проверить, что длина кратна 3 и нет стопов
    conda activate base
    # python ${PATH_TO_SRC}/check_cds_qc.py "$ALIGNdir/$og.codon_aln.fa" || {
    #     echo "[WARN] QC failed for $og"
    #     return 0
    # }
    # 2. Фильтрация MACSE ошибок (!, #, >50% gaps)
    python ${PATH_TO_SRC}/filter_macse.py \
        "$ALIGNdir/$og.macse.fa" \
        "$ALIGNdir/$og.macse.clean.fa"

    # conda activate mafft
    # # Trim плохих позиций
    # trimal \
    #     -in "$ALIGNdir/$og.codon_aln.fa" \
    #     -out "$ALIGNdir/$og.codon_aln.trim.fa" \
    #     -automated1

    # # Удалить слишком плохие последовательности
    # trimal \
    #     -in "$ALIGNdir/$og.codon_aln.trim.fa" \
    #     -out "$ALIGNdir_final/$og.codon_aln.final.fa" \
    #     -resoverlap 0.5 \
    #     -seqoverlap 60

    # trimal -in "$ALIGNdir/$og.codon_aln.fa" \
    #     -out "$ALIGNdir/$og.codon_aln.trim.fa" -gt 0.8

    # echo "[OK] trimmed $og"

    echo "[OK] $og"
}

export -f process_og
export PATH_TO_SRC PATH_TO_MAC
export OGdir ALIGNdir ALIGNdir_final 

# Запускаем параллельно для всех OG
find "$OGdir" -name "*.cds.fa" | parallel -j "$SLURM_CPUS_PER_TASK" process_og

echo "[DONE] All orthogroups processed"



# #!/bin/bash
# cds="$1"

# PATH_TO_PAIRS=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/

# OGdir="${PATH_TO_PAIRS}/genome_data/ncbi_dataset/data/per_OG"
# ALIGNdir="${PATH_TO_PAIRS}/codeml/alignments"

# mkdir -p "$ALIGNdir"

# og=$(basename "$cds" .cds.fa)
# prot="$OGdir/$og.prot.fa"

# echo "[INFO] Processing $og"

# if [ ! -f "$prot" ]; then
#     echo "[WARN] No protein file for $og" >&2
#     exit 0
# fi

# source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

# conda activate mafft
# mafft --quiet --auto "$prot" > "$ALIGNdir/$og.prot.aln.fa"

# # Кодонное выравнивание
# pal2nal.pl \
#     "$ALIGNdir/$og.prot.aln.fa" \
#     "$cds" \
#     -output fasta \
#     > "$ALIGNdir/$og.codon_aln.fa"

# echo "[OK] $og"