# conda activate r-environment

require("rphast")

ref <- "GCA_048771995.1.CM109824.1"

# upload data
align <- read.msa("maf_by_chr/1.CM109824.1.maf")
birdsTree <- read.newick.tree("roadies_v1.1.16b_renamed.nwk")
feats <- read.feat("tmp/gff_by_chromosome/CM109824.1.cds.gff")

# get whole chromosome region
wholeChrom <- feat(
    seq=ref, 
    src=".", 
    feature="all", 
    start=align$offset, 
    end=align$offset+ncol.msa(align, ref)
)

# get 4d sites and estimate neutral model
align4d <- get4d.msa(align, feats)

neutralMod <- phyloFit(
    align4d, 
    tree=birdsTree, 
    subst.mod="REV",
    EM=TRUE
)
summary(neutralMod)

# phastCons 
pc <- phastCons(align, neutralMod, expected.length=6,target.coverage=0.125, viterbi=TRUE)
consElements <- pc$most.conserved

# phyloP
pp <- phyloP(neutralMod, align, method="LRT", mode="CONACC")

# plotResults
codingFeats <- feats[feats$feature=="CDS",]
geneTrack <- as.track.feat(codingFeats, "genes", is.gene=TRUE)
consElTrack <- as.track.feat(consElements, "phastCons most conserved", col="red")
phastConsScoreTrack <- as.track.wig(wig=pc$post.prob.wig, name="phastCons post prob", col="red", ylim=c(0, 1))
phyloPTrack <- as.track.wig(coord=pp$coord, score=pp$score, name="phyloP score",col="blue", smooth=TRUE, horiz.line=0)


align$offset+ncol.msa(align, "GCA_048771995.1")
pdf("tracks_region_CM109824.1.pdf", 
    width = 10,  
    height = 6)

plot.track(
    list(
        geneTrack, 
        consElTrack, 
        phastConsScoreTrack, 
        phyloPTrack
    ),
    # xlim=c(3000000, 3080000), 
    cex.labels=1.25, 
    cex.axis=1.25, 
    cex.lab=1.5
)

dev.off()
