#!/usr/bin/env Rscript
# Andrew W. Brooks
# Vanderbilt Genetics Institute
# 16_8_8

#########################################################################################
##### DESCRIPTION #####
# THIS IS A 'R' PIPELINE TO ANALYZE 16S rDNA MICROBIOME DATA GENERATED IN QIIME (OR MOTHUR)
# MICROBIOME COUNT DATA SHOULD BE IN BIOM FORMAT 

##### STRUCTURE AND FUNCTIONS #####
# ENTER THE DATA SOURCES INTO THE VARIABLES BELOW
# REQUIRED FUNCTIONS ARE NECESSARY TO LOAD IMPORTANT LIBRARIES (SHOULD NOT BE TOUCHED)
# OPTIONAL FUNCTIONS ARE DESIGNED FOR SPECIFIC ANALYSES
# THE MAIN SECTION IS WHERE REQUIRED FUNCTIONS ARE USED TO LOAD AND FORMAT DATA,
# AND OPTIONAL FUNCTIONS ARE USED TO PERFORM SPECIFIC PIPELINE ANALYSES

### CALLING THE SCRIPT ###
# ./phyloseq.R table.biom map.txt bacterial_phylogeny.tre
# ARGS[1] = BIOM TABLE
# ARGS[2] = MAP FILE
# ARGS[3] = BACTERIAL PHYLOGENY TREE

#########################################################################################
##### REQUIRED FUNCTIONS #####

###################################
##### INPUT ARGUMENTS #####
# IF YOU WANT COMMANDS PRINTED TO OUTPUT AS WELL = TRUE
options(echo=FALSE)

# INPUT #
args <- commandArgs(trailingOnly = TRUE)
print("INPUT ARGUMENTS:")

# BIOM TABLE #
biomPath <- args[1]
print("   BIOM TABLE:")
print(biomPath)

# MAPPING FILE #
mapPath <- args[2]
print("   MAPPING FILE:")
print(mapPath)

# TREE FILE #
treePath <- args[3]
print("   TREE FILE:")
print(treePath)

# OUTPUT DIRECTORY #
outDir <- args[4]
print("   OUTPUT DIRECTORY:")
print(treePath)

### MANUAL ENTRY (FOR TESTING) ###
#biomPath <- "/Users/brooks/Documents/brooks_utilities/microbiome/test_16S/nasonia/5_7_filter_wolbachia_from_otu_table.biom"
#mapPath <- "/Users/brooks/Documents/brooks_utilities/microbiome/test_16S/nasonia/map.txt"
#treePath <- "/Users/brooks/Documents/brooks_utilities/microbiome/test_16S/nasonia/4_5_make_phylogeny_gram_align.tre"
#outDir <- "/Users/brooks/Documents/brooks_utilities/microbiome/test_16S/nasonia/phyloseq_analysis/"
#getwd()

###################################
### MAKE AND NAVIGATE TO OUTPUT DIRECTORY ###
if(file.exists(outDir)){
   setwd(file.path(outDir))
} else{
   dir.create(file.path(outDir))
   setwd(file.path(outDir))
}

print("WORKING DIRECTORY:")
getwd()

###################################
### INSTALL LIBRARIES ###
install_libraries <- function(){
   install.packages("phyloseq")
   install.packages("gplots")
   install.packages("ggplot2")
   install.packages("scales")
   install.packages("grid")
   install.packages("GMD")
}

### IMPORT LIBRARIES ###
import_libraries <- function(){
   library("phyloseq")
   packageVersion("phyloseq")
   library("ggplot2")
   packageVersion("ggplot2")
   library("scales")
   packageVersion("scales")
   library("grid")
   packageVersion("grid")
   theme_set(theme_bw())
   library("gplots")
   require("GMD")
}
import_libraries()

###################################
### IMPORT DATA ###
import_data <- function(biomPath, mapPath, treePath){
   biomot <- import_biom(biomPath, treePath, parseFunction = parse_taxonomy_greengenes)
   bmsd <- import_qiime_sample_data(mapPath)
   class(bmsd)
   dim(bmsd)
   biomot
   # Merge Files (!)
   phylo_biom = merge_phyloseq(biomot, bmsd)
   return(phylo_biom)
}

phyloBiom <- import_data(biomPath, mapPath, treePath)

###################################
### BACKUP DATASET ###
backupBiom <- phyloBiom

# IF THINGS GO AWRY IN THE ANALYSIS AND THIS IS BEING RUN MANUALLY THEN YOU CAN RESET THE
# DATA TO AVOID REOADING WITH THIS CALL
# phyloBiom <- backupBiom

#########################################################################################
##### OPTIONAL FUNCTIONS #####

###################################
### DATA CONTROL FUNCTIONS ###

# Merge Samples
metaCat <- "Treatment"
phyloCategory <- merge_samples(phyloBiom, metaCat)

# Merge OTUs by Taxa
phylo_biom = tax_glom(phylo_biom, "Genus")

# Filter to Specific Taxonomy - subset $taxa
phylo_biom = subset_taxa(phylo_biom, Phylum=="Chlamydiae")

# Remove Samples with Less than min_reads
min_reads <- 20
phylo_biom = prune_samples(sampleSums(phylo_biom)>=min_reads, phylo_biom)

# Relative Abundance
phylo_biom = transform_sample_counts(phylo_biom, function(x) x/sum(x))

# Filter top N OTUs
phylo_biom <- prune_taxa(names(sort(taxa_sums(phylo_biom), TRUE)[1:1000]), phylo_biom)

# Relative Abundance Filter (.00001 = .001% - Below Removed)
phylo_biom = filter_taxa(phylo_biom, function(x) mean(x) > 1e-5, TRUE)


###################################
##### TAXONOMY PIPELINE #####

### PHYLUM ###
taxLevel <- "Phylum"

# COLLAPSE BY TAXONOMY #
phyloTax = tax_glom(phyloBiom, taxLevel)

# CONVERT TO RELATIVE ABUNDANCE #
phyloTaxRel = transform_sample_counts(phyloTax, function(x) x/sum(x))

# 
pdf(paste(taxLevel, '_tree.pdf', sep = ""), width=30, height=30)
plot_tree(phyloBiom, ladderize=TRUE, label.tips="taxa_names", size="abundance", color="Phylum", shape="Treatment") + coord_polar(theta="y")
dev.off()

# PLOT BAR RELATIVE #
pdf(paste(taxLevel, '_barchart.pdf', sep = ""), width=40, height=30)
plot_bar(phyloTaxRel, fill = taxLevel)
dev.off()

plot_net(phyloTaxRel, dist.fun="bray", color="Treatment")

GP.ord <- ordinate(phyloTaxRel, "NMDS", "bray")
plot_ordination(phyloTaxRel, GP.ord, type="taxa", color="Phylum", title="taxa")

phyloTaxRelPrune <- prune_taxa(names(sort(taxa_sums(phyloTaxRel),TRUE)[1:5]), phyloTaxRel)
plot_heatmap(phyloTaxRelPrune)
plot_heatmap(phyloRelative, "NMDS", "bray", "Taxonomy", "Family")


plot_heatmap(phyloTaxRel, "Phylum")

### CLASS ###
taxLevel <- "Class"

# COLLAPSE BY TAXONOMY #
phyloTax = tax_glom(phyloBiom, taxLevel)

# CONVERT TO RELATIVE ABUNDANCE #
phyloTaxRel = transform_sample_counts(phyloTax, function(x) x/sum(x))

# PLOT BAR RELATIVE #
pdf(paste(taxLevel, '.pdf', sep = ""), width=40, height=30)
plot_bar(phyloTaxRel, fill = taxLevel)
dev.off()


### ORDER ###
taxLevel <- "Order"

# COLLAPSE BY TAXONOMY #
phyloTax = tax_glom(phyloBiom, taxLevel)

# CONVERT TO RELATIVE ABUNDANCE #
phyloTaxRel = transform_sample_counts(phyloTax, function(x) x/sum(x))

# PLOT BAR RELATIVE #
pdf(paste(taxLevel, '.pdf', sep = ""), width=40, height=30)
plot_bar(phyloTaxRel, fill = taxLevel)
dev.off()


### FAMILY ###
taxLevel <- "Family"

# COLLAPSE BY TAXONOMY #
phyloTax = tax_glom(phyloBiom, taxLevel)

# CONVERT TO RELATIVE ABUNDANCE #
phyloTaxRel = transform_sample_counts(phyloTax, function(x) x/sum(x))

# PLOT BAR RELATIVE #
pdf(paste(taxLevel, '.pdf', sep = ""), width=40, height=30)
plot_bar(phyloTaxRel, fill = taxLevel)
dev.off()


### GENUS ###
taxLevel <- "Genus"

# COLLAPSE BY TAXONOMY #
phyloTax = tax_glom(phyloBiom, taxLevel)

# CONVERT TO RELATIVE ABUNDANCE #
phyloTaxRel = transform_sample_counts(phyloTax, function(x) x/sum(x))

# PLOT BAR RELATIVE #
pdf(paste(taxLevel, '.pdf', sep = ""), width=40, height=30)
plot_bar(phyloTaxRel, fill = taxLevel)
dev.off()

#########################################################################################