# Oral 16S Analyses

# packages ####
library(tidyverse)
library(phyloseq)
library(vegan)
library(microbiome)
library(RColorBrewer)

#Load phyloseq object with tree
ps2 = readRDS("./16S/output/oral_phyloseq_object_w_tree.RDS")

# Explore and prepare metadata ####
meta = as.data.frame(sample_data(ps2))

glimpse(meta)
# condense metadata
OF_Variables = c("Diet","ActivityLevel","Age","Gender","Ethnicity",
                 "WeeklyChewTobacco","WeeklyAlcohol","WeeklyCigarettes","WeeklyVape","WeeklyDrugs")

# subset to new ps object
OF_meta = select(meta,OF_Variables)
ps2 = phyloseq(otu_table(otu_table(ps2)),
         tax_table(tax_table(ps2)),
         sample_data(OF_meta),
         phy_tree(phy_tree(ps2)))

# Remove eukaryotes left over
ps2 = subset_taxa(ps2, Kingdom == "Bacteria")

# check for missing metadata
complete.cases(ps2@sam_data)

glimpse(ps2@sam_data)

# convert to appropriate classes
ps2@sam_data$Diet = factor(ps2@sam_data$Diet)
ps2@sam_data$ActivityLevel[ps2@sam_data$ActivityLevel == "NA"] <- "0"
ps2@sam_data$ActivityLevel = factor(ps2@sam_data$ActivityLevel, ordered = TRUE, levels = c("0","1-2","3-4","5-7"))
ps2@sam_data$Age = as.numeric(ps2@sam_data$Age)
ps2@sam_data$Gender = as.factor(ps2@sam_data$Gender)
ps2@sam_data$Ethnicity = as.factor(ps2@sam_data$Ethnicity)
ps2@sam_data$WeeklyChewTobacco = as.numeric(ps2@sam_data$WeeklyChewTobacco)
ps2@sam_data$WeeklyAlcohol = as.numeric(ps2@sam_data$WeeklyAlcohol)
ps2@sam_data$WeeklyCigarettes = as.numeric(ps2@sam_data$WeeklyCigarettes)
ps2@sam_data$WeeklyVape = as.numeric(ps2@sam_data$WeeklyVape)
ps2@sam_data$WeeklyDrugs = as.numeric(ps2@sam_data$WeeklyDrugs)

glimpse(ps2@sam_data)

# save summary info about metadata
sink("./16S/output/OF_metadata_summaries.txt")
for(i in names(ps2@sam_data)){
  print(summary(ps2@sam_data[,i], main = i))
}
sink(NULL)

# Convert data to compositional
ps_ra = microbiome::transform(ps2, "compositional")

# Save taxonomy to SeqID and Seq dictionary
tax = as.data.frame(ps_ra@tax_table@.Data)
taxIDs = paste0("k__",tax[,1],";p__",tax[,2],";c__",tax[,3],";o__",tax[,4],
       ";f__",tax[,5],";g__",tax[,6],";s__",tax[,7])

seq_to_tax = data.frame(SeqID = paste0("ASV_", seq(1:length(tax[,1]))),
                       Tax = taxIDs, Seq = colnames(ps_ra@otu_table))

# rename phyloseq "otu" table
taxa_names(ps_ra) <- seq_to_tax$SeqID


# find core members of oral microbiome ####
# Core with compositionals:
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(ps_ra, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")
ggsave("./16S/output/CoreSize_vs_RelativeAbundance.png", dpi = 300)

prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette
gray <- gray(seq(0,1,length=10))
p <- plot_core(ps_ra, plot.type = "heatmap", colours = gray,
               prevalences = prevalences, detections = detections) +
  xlab("Detection Threshold (Relative Abundance (%))")
print(p)    


# Core heatmap
prevalences <- seq(.05, .5, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

p <- plot_core(ps_ra, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE)
print(p)
ggsave("./16S/output/Core_Heatmap.png", dpi=300)

# prevalence
taxa_prevalence = (prevalence(ps_ra, detection = 0.001, sort = TRUE))

# core members at >= 0.2 sample prevalence 0.01% detection threshold
core.taxa.standard <- core_members(ps_ra, detection = 0.001, prevalence = .1)

# Total core abundance in each sample (sum of abundances of the core members):
ps_core <- core(ps_ra, detection = 0.001, prevalence = .1)
# subset to only samples containing core microbiome
ps_core_samples = subset_samples(ps_core,rowSums(otu_table(ps_core)) > 0)

nmds_core = ordinate(ps_core_samples, "NMDS")
plot_ordination(ps_core_samples, nmds_core, color = "Diet") +
  stat_ellipse()
ggsave("./16S/output/NMDS_core_taxa-and-samples.png", dpi=300)

tiff("./16S/output/NMDS_stressplot.tiff")
stressplot(nmds_core)
dev.off()

permanova_core = adonis(otu_table(ps_core_samples) ~ ps_core_samples@sam_data$Diet + ps_core_samples@sam_data$WeeklyAlcohol +
         ps_core_samples@sam_data$ActivityLevel + ps_core_samples@sam_data$Age + ps_core_samples@sam_data$Gender + ps_core_samples@sam_data$Ethnicity)

sink("./16S/output/permanova_table.txt")
print("Permanova test on core taxa/samples:", quote = FALSE)
permanova_core
sink(NULL)

# names of core taxa
core_IDs = which(seq_to_tax$SeqID %in% taxa(ps_core) )
core_taxa_names = as.character(unique(seq_to_tax[core_IDs,"Tax"]))

# Core taxa broken up by gender
ps_male = subset_samples(ps_ra,Gender == "Male")
ps_female = subset_samples(ps_ra,Gender == "Female")

# Core plots
plot_core(ps_male, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")
plot_core(ps_female, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")

# Total core abundance in each sample (sum of abundances of the core members):
ps_core_m <- core(ps_male, detection = 0.001, prevalence = .1)
ps_core_f <- core(ps_female, detection = 0.001, prevalence = .1)

p.malecore <- plot_core(ps_male, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE) + ggtitle("Core - Male subjects")
print(p.malecore)
ggsave("./16S/output/Core_Heatmap_Male.png", dpi=300)

p.femalecore <- plot_core(ps_female, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE)+ ggtitle("Core - Female subjects")
print(p.femalecore)
ggsave("./16S/output/Core_Heatmap_Female.png", dpi=300)


# names of core taxa
core_IDs = which(seq_to_tax$SeqID %in% taxa(ps_core_m) )
core_taxa_names_m = as.character(unique(seq_to_tax[core_IDs,"Tax"]))
core_IDs = which(seq_to_tax$SeqID %in% taxa(ps_core_f) )
core_taxa_names_f = as.character(unique(seq_to_tax[core_IDs,"Tax"]))

# write output to text file (names of core taxa)
sink("./16S/output/Core_Taxa_Names.txt")
print("Detection .1%, prevalence 10%",quote = FALSE)
print("For all samples:", quote = FALSE)
core_taxa_names
print("For males:", quote = FALSE)
core_taxa_names_m
print("For Females:", quote = FALSE)
core_taxa_names_f
sink(NULL)


# Core barplots
classGlommedcore = tax_glom(ps_core_samples, "Class")
classGlommedcore = subset_samples(classGlommedcore, Gender != "NA")
genders = levels(classGlommedcore@sam_data$Gender)


p=plot_bar(classGlommedcore, fill="Class")
p + facet_grid(~Class)

ggsave("./16S/output/barplot_class_core.png",dpi=300)

# core tree plot

plot_tree(ps2,color = "Phylum")
################################

# barplots
phylumGlommed = tax_glom(ps_ra, "Phylum")
phylumGlommed = subset_samples(phylumGlommed, Gender != "NA")
genders = levels(phylumGlommed@sam_data$Gender)

phylumGlommed_gender = merge_samples(phylumGlommed,group = "Gender")
phylumGlommed_gender@sam_data$Gender <- genders

p = plot_bar(phylumGlommed_gender, fill = "Phylum")
p + facet_wrap(~Phylum)
ggsave("./16S/output/barplot_phylum_by_gender_faceted.png",dpi=300)

phylumGlommed = tax_glom(ps_ra, "Class")
phylumGlommed = subset_samples(phylumGlommed, Gender != "NA")
genders = levels(phylumGlommed@sam_data$Gender)

phylumGlommed_gender = merge_samples(phylumGlommed,group = "Gender")
phylumGlommed_gender@sam_data$Gender <- genders

p = plot_bar(phylumGlommed_gender, fill = "Class")
p + facet_wrap(~Class)
ggsave("./16S/output/barplot_class_by_gender_faceted.png",dpi=300, width = 12, height = 8)



####################

dca = ordinate(ps_ra, method = "DCA")
plot_ordination(ps_ra,dca, type = "Samples")


phy_tree(ps2) <- fitGTR$tree
plot_tree(physeq = ps_ra,color = "Gender")
