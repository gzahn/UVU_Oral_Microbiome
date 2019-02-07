# Oral 16S Analyses
# add phylogenetic tree info and explore with phyloseq

# packages ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(phangorn)
library(msa)
library(microbiome)
library(RColorBrewer)

# Read original ps object ####
ps = readRDS("16S/output/phyloseq_object_16S_noncontam.RDS")

ps@sam_data$Project_Name

# subset to oral microbiome and taxa with more than 5 hits ####
ps = subset_samples(ps, Project_Name == "Oral_Mycobiome")
ps = subset_taxa(ps,(colSums(otu_table(ps)) > 5))


# Build Phylogenetic Tree ####
seqs <- rownames(tax_table(ps))
names(seqs) <- seqs # This propagates to the tip labels of the tree
# alignment
alignment <- msa(seqs,method = "Muscle", type = "dna")
# saveRDS(alignment,"~/Desktop/oral_dna_alignment_muscle.RDS")
phang.align = as.phyDat(alignment, type = "DNA")

# distance max likelihood
dm <- dist.ml(phang.align)
# neighbor-joining tree
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeNJ$tip.label <- seqs

fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!

(unique(names(seqs)))


fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    control = phangorn::pml.control(trace = 1L))


saveRDS(fitGTR, "16S/output/oral_fitGTR2.RDS") # This is the new tree using optim.pml
fitGTR = readRDS("16S/output/oral_fitGTR2.RDS") # loading new tree. does it work??
# write.tree(fitGTR$tree, file = "~/Desktop/oral_bact_tree.nwk")

# bs = bootstrap.pml(fitGTR, bs=100, optNni=TRUE, multicore=TRUE)

detach("package:phangorn", unload=TRUE)

# add tree to phyloseq object ####
ps2 = phyloseq(tax_table(tax_table(ps)),
               otu_table(otu_table(ps)),
               sample_data(sample_data(ps)),
               phy_tree(fitGTR$tree))
saveRDS(ps2, "~/Desktop/oral_phyloseq_object_w_tree.RDS")
