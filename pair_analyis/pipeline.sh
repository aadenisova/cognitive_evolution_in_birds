annotate_species.sh

get_longest_protein.sh

# need to paste date to get OrthoFinder results
get_orthogroups.sh Dec02
get_codon_alignmnets.sh

codeml/convert_to_phy.py

PATH_TO_PAIRS=/vggpfs/fs3/vgl/store/adenisova/Inno_2025/pair_analyis

python \
    ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
    ${PATH_TO_PAIRS}/initial_data/trees/pruned_8_species_analysis.newick \
    Zonotrichia_albicollis,Sylvia_borin,Falco_rusticolus,Aythya_ferina \
    tree_labelled_fg_non_inno.nwk

python \
    ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
    ${PATH_TO_PAIRS}/initial_data/trees/pruned_8_species_analysis.newick \
    Zonotrichia_albicollis,Sylvia_borin,Falco_rusticolus,Aythya_ferina \
    tree_labelled_fg_non_inno.nwk

python \
    ${PATH_TO_PAIRS}/src/pair_analyis/codeml/label_tree.py \
    ${PATH_TO_PAIRS}/initial_data/trees/pruned_8_species_analysis.newick \
    Zonotrichia_albicollis,Sylvia_borin,Falco_rusticolus,Aythya_ferina \
    tree_labelled_fg_non_inno.nwk

codeml/label_tree.py

codeml/convert_to_phy.py --strict 