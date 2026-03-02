# conda activate r-environment
# cat TOGA_ALL_selected_polished_for_RER/* > TOGA_ALL_new_tree_reest.txt

library(RERconverge)
library(tidyverse)
library(ape)

matchNodesInject=function (tr1, tr2){
  if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
    #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
    return(FALSE)
    stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
  }
  commontiplabels <- intersect(tr1$tip,tr2$tip)
  if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
    return(FALSE)
    stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
  }
  #if(RF.dist(tr1,tr2)>0){
  #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
  #}
  
  toRm=setdiff(tr2$tip.label, tr1$tip.label)
  desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                              x)$tip.label)
  names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
  desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                              x)$tip.label)
  names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
  Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                           c("tr1", "tr2")))
  for (i in 1:length(desc.tr1)) {
    Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
    for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in%
                                          desc.tr2[[j]]))
      Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
  }
  
  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  if(any(table(Nodes[,2])>1)){
    return(FALSE)
    stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
  }
  
  return(TRUE)
}


path="/lustre/fs5/jarv_lab/scratch/adenisova/Inno_2025/"
treeFile = paste(path, "tmp/checked_TOGA_ALL_new_tree_reest.txt", sep="/")

prefix <- "TOGA_ALL_selected"

inno_rates <- read_csv(paste(path, "/src/whole_tree_analysis/cds/data/inno_aligned_with_TOGA.csv", sep = "/"))

tree_df = read_tsv(
    treeFile,
    col_names = c("gene", "tree")
)

ourTrees=readTrees(
    treeFile,
    max.read = 500   
)

masterTree <- ourTrees$masterTree
tree_df$good <- logical(nrow(tree_df))

for (i in seq_len(nrow(tree_df))) {
  tr <- read.tree(text = tree_df$tree[i])

  tree_df$good[i] <- matchNodesInject(
    tr,
    masterTree
  )
}
table(tree_df$good)
filtered_df <- tree_df[tree_df$good, c("gene", "tree")]

write.table(
  filtered_df,
  file = "filtered_trees.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

ourTrees=readTrees(
    "filtered_trees.tsv"
)
ourTrees$numTrees

mamRERw=RERconverge::getAllResiduals(
    ourTrees,
    useSpecies=names(sp_vect),
    transform = "sqrt", 
    weighted = T,
    scale = T,
)

inno_rates %>% colnames()
cols_of_interest <- c("FOODINNO2025_ResEff", "TECHINNO2025_ResEff", "Relative_brain_size", "Wing_length")

for (col in cols_of_interest) {
    print(col)

    sp_vect <- setNames(pull(inno_rates[col]), inno_rates$sci_name_2025)
    charpaths=char2Paths(sp_vect, ourTrees)

    res=RERconverge::getAllCor(
        mamRERw, 
        charpaths, 
        method = "p", 
        min.pos = 0, 
        winsorizeRER = 3,
    ) 

    ordered_results = res[order(res$P),]
    head(ordered_results)
    ordered_results$gene <- row.names(ordered_results)

    ordered_results %>% head()

    write_tsv(
        ordered_results,
        paste(paste(path,"tmp", prefix,sep="/"), col, "RERConverge_result.tsv", sep="_"),
    )
}
