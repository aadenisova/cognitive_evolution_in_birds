PATH_TO_PAIRS=/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/

# Inno
SPS_INNO=(
Hirundo_rustica
Passer_domesticus
Corvus_moneduloides
Falco_peregrinus
Larus_argentatus
Buteo_buteo
)

TREE_FILE="${PATH_TO_PAIRS}/initial_data/trees/pruned_12_species_analysis.newick"

# Цикл по всем видам
for sp in "${SPS_INNO[@]}"; do
    # Создаем имя файла для каждого вида
    OUTPUT_FILE="${PATH_TO_PAIRS}/initial_data/trees_labeled/tree_labelled_fg_inno_${sp}.nwk"
    
    echo "Обработка вида: $sp"
    echo "Сохранение в: $OUTPUT_FILE"
    
    python \
        ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
        $TREE_FILE \
        $sp \
        $OUTPUT_FILE
done

python \
    ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
    ${PATH_TO_PAIRS}/initial_data/trees/pruned_10_species_analysis.newick \
    Hirundo_rustica,Passer_domesticus,Corvus_moneduloides,Falco_peregrinus,Buteo_buteo \
    pruned_10_species_analysis_inno.nwk

python \
    ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
    ${PATH_TO_PAIRS}/initial_data/trees/pruned_10_species_analysis.newick \
    Sylvia_borin,Zonotrichia_albicollis,Lycocorax_pyrrhopterus,Falco_rusticolus,Haliaeetus_albicilla \
    pruned_10_species_analysis_noninno.nwk

# python \
#     ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
#     ${PATH_TO_PAIRS}/initial_data/trees/pruned_12_species_analysis.newick \
#     Hirundo_rustica,Passer_domesticus,Corvus_moneduloides,Falco_peregrinus,Larus_argentatus,Buteo_buteo \
#     pruned_12_species_analysis_inno.nwk

# python \
#     ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
#     ${PATH_TO_PAIRS}/initial_data/trees/pruned_12_species_analysis.newick \
#     Sylvia_borin,Zonotrichia_albicollis,Lycocorax_pyrrhopterus,Falco_rusticolus,Larus_fuscus,Haliaeetus_albicilla \
#     pruned_12_species_analysis_noninno.nwk

# python \
#     ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
#     ${PATH_TO_PAIRS}/initial_data/trees/pruned_8_species_analysis_01.newick \
#     $SPS_INNO \
#     tree_labelled_fg_inno.nwk

# python \
#     ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
#     ${PATH_TO_PAIRS}/initial_data/trees/pruned_8_species_analysis.newick \
#     Zonotrichia_albicollis,Sylvia_borin,Falco_rusticolus,Aythya_ferina \
#     tree_labelled_fg_non_inno.nwk
