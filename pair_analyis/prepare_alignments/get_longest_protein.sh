#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=longest_isoforms_parallel
#SBATCH --output=longest_isoforms_parallel_%j.out
#SBATCH --error=longest_isoforms_parallel_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --mail-user=savouriess2112@gmail.com

module load parallel   # если нужно — иначе убрать

mkdir -p longest_isoforms_gff
mkdir -p longest_isoforms_cds
mkdir -p longest_isoforms_proteomes

DIRS=(
GCA_014839755.1
GCF_015227805.2
GCF_036417665.1
GCF_047830755.1
GCF_009650955.1
GCA_014706295.1
GCF_023634155.1
GCF_015220075.1
# GCA_964417175.1
# GCA_963932325.2
GCF_947461875.1
GCF_964188355.1
)

source /vggpfs/fs3/vgl/store/adenisova/anaconda3/etc/profile.d/conda.sh

process_dir() {
    DIR="$1"
    echo "=== Processing $DIR ==="

    cd $DIR #|| exit 1

    # определяем GFF + genome
    if [[ -f "lifted_annotations_genomic.gff" ]]; then
        GFF="lifted_annotations_genomic.gff"
        GENOME=$(ls *genomic.fna)
    elif [[ -f "genomic.gff" ]]; then
        GFF="genomic.gff"
        GENOME=$(ls $DIR* | grep _genomic.fna)
    else
        echo "No valid GFF in $DIR"
        cd ..
        return
    fi

    echo $DIR $GENOME

    # run AGAT
    source ~/.bashrc
    conda activate agat
    agat_sp_keep_longest_isoform.pl \
        -g "$GFF" \
        -o longest_isoforms.gff

    sed -E "s/ID=([^;]+)/ID=${DIR}_\1/; s/Parent=([^;]+)/Parent=${DIR}_\1/" \
    longest_isoforms.gff > longest_isoforms_renamed.gff

    ls

    # run gffread
    conda activate gffread
    gffread longest_isoforms_renamed.gff \
        -g $GENOME \
        -x cds_longest.fa \
        -y proteins_longest.fa

    # copy results
    cp longest_isoforms_renamed.gff ../longest_isoforms_gff/"${DIR}_longest.gff"
    cp cds_longest.fa ../longest_isoforms_cds/"${DIR}_cds.fa"
    cp proteins_longest.fa ../longest_isoforms_proteomes/"${DIR}_proteins.fa"

    cd ..
}

export -f process_dir

parallel -j $SLURM_CPUS_PER_TASK process_dir ::: "${DIRS[@]}"
