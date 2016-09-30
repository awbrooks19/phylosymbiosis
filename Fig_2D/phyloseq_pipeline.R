setwd("/Users/brooksaw/Desktop/phyloseq_analysis/nasonia")

# Import Files Individually
biom_file = "5_7_filter_wolbachia_from_otu_table.biom"
map_file = "map.txt"


biom_file = "6_5_summarize_otu_by_cat.biom"
map_file = "map_cat.txt"

tre_file = "4_5_make_phylogeny_gram_align.tre"
biomot = import_biom(biom_file, tre_file, parseFunction = parse_taxonomy_greengenes)
bmsd = import_qiime_sample_data(map_file)

#Get Info
class(bmsd)
dim(bmsd)
biomot

# Merge Files
phylo_biom = merge_phyloseq(biomot, bmsd)
phylo_biom

# Basic Plot of Alpha Diversity
plot_richness(phylo_biom, x = "Treatment", color = "Treatment") + geom_boxplot()

# Interacting with Sample Data
ntaxa(phylo_biom)
nsamples(phylo_biom)
sample_names(phylo_biom)[1:10]
taxa_names(phylo_biom)[1:10]
sample_variables(phylo_biom)[1:10]
length(sample_variables(phylo_biom))
get_variable(phylo_biom, sample_variables(phylo_biom)[5])[1:10]

# Interacting with Taxonomic Data
rank_names(phylo_biom)
get_taxa_unique(phylo_biom, "Phylum")

# Basic Plots
  # Alpha Diversity by Sample
plot_richness(phylo_biom)
  # Alpha Diversity by Group
(p = plot_richness(phylo_biom, x = "Treatment"))
  #  Alpha Diversity with Quartile Range
p + geom_boxplot(data = p$data, aes(x = Treatment, y = value, color = NULL), alpha = 0.1)

# Phylogenetic Tree Plots
get_taxa_unique(phylo_biom, "Phylum")
phylo_biom.phyla_tree = subset_taxa(phylo_biom, Phylum == "Proteobacteria")
plot_tree(phylo_biom.phyla_tree, color = "Treatment", label.tips = "Class", 
          size = "abundance", plot.margin = 1.5, ladderize = TRUE)
plot_tree(GP.chl, "treeonly", nodeplotblank)

# Abundance Plots
TopNOTUs = names(sort(taxa_sums(phylo_biom), TRUE)[1:20])
ent10 = prune_taxa(TopNOTUs, phylo_biom)
plot_bar(ent10, "Treatment", fill = "Sample",facet_grid = ~Order)

# prune OTUs that are not present in at least one sample
GP = prune_taxa(taxa_sums(phylo_biom) > 0, phylo_biom)
# Define a human-associated versus non-human categorical variable:
sample_data(GP)$human = factor(get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))

# Rarify Even Depth
rarefied_phylo = rarefy_even_depth(phylo_biom)
UniFrac(rarefied_phylo)

GP.subset = subset_taxa(phylo_biom, Phylum == "Proteobacteria")
# remove the samples that have less than 20 total reads from Chlamydiae
GP.subset = prune_samples(names(which(sample_sums(GP.subset) >= 20)), GP.subset)
# (p = plot_tree(GP.chl, color='SampleType', shape='Family',
# label.tips='Genus', size='abundance'))
GP.subset.r = rarefy_even_depth(GP.subset)
# Unweighted UNIFRAC
plot_ordination(GP.subset, ordinate(GP.subset , "MDS"), color = "Treatment") + geom_point(size = 5)
plot_ordination(GP.subset.r, ordinate(GP.subset.r, "MDS"), color = "Treatment") + geom_point(size = 5)

# UNIFRAC WHOLE SET
plot_ordination(phylo_biom, ordinate(phylo_biom , "MDS"), color = "Treatment") + geom_point(size = 5)
  # Rarefied
plot_ordination(rarefied_phylo, ordinate(rarefied_phylo, "MDS"), color = "Treatment") + geom_point(size = 5)


