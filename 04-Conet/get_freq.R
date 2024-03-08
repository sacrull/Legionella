rm(list = ls())
#library(devtools)
#library(spieceasi)
#library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
#library(RColorBrewer) # nice color options
#library(ggplot2) # publication quality figures, based on ggplot2
#library(dplyr) # data handling
#library(network)
#library(intergraph)
#devtools::install_github("briatte/ggnet")
#library(ggnet)
#library(igraph)
library(tidyverse)
#library(qiime2R)
#library(readr)
library(RColorBrewer)
#library(viridis)
#library(Matrix)
library(Hmisc)
#library(corrplot)
#library(reshape2)
#library(data.table)
library(circlize) #https://mran.microsoft.com/snapshot/2016-01-05/web/packages/circlize/vignettes/visualize_relations_by_chord_diagram.pdf
#library(reshape2)
library(chorddiag)
setwd("~/legionella/R/conet")
#reading sequence tables
seqtab_16S <- read.table("../feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_16S)
seqtab_18S <- read.table("../18S-ASV-renamed.txt", header=T, row.names=1)

seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)

#reading in taxonomy
tax_16S <- read.table("../taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)
tax_18S <- read.table("../18S-tax-renamed.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
#reeading in metadata
map_16S <- read.table("../metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)
map_18S <- read.table("../metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
#make 16S and 18S phyloseq object from qiime2
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
physeq_16S <- subset_taxa(physeq_16S, !V3=="Bacteria_unknown" & !V2=="Eukaryota")
physeq_16S1 <- prune_samples(sample_sums(physeq_16S) > 10000, physeq_16S) #remove less than 2000 reads
rare_16S <- rarefy_even_depth(physeq_16S1, rngseed=1, sample.size=0.99*min(sample_sums(physeq_16S1)), replace=F)
physeq_16S_filter2 = subset_taxa(rare_16S, V7=="Legionella")

physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
physeq_18S <- subset_taxa(physeq_18S, !V2=="Bacteria" & !V3=="Eukaryota_unknown" & !V3=="uncultured")
physeq_18S2 <- prune_samples(sample_sums(physeq_18S) > 5000, physeq_18S) #remove less than 2000 reads
rare_18S <- rarefy_even_depth(physeq_18S2, rngseed=1, sample.size=0.99*min(sample_sums(physeq_18S2)), replace=F)
#known hosts and predator of legionella
physeq_18S_filter1 = subset_taxa(rare_18S, V3=="Amoebozoa" | V3=="Heterolobosea" | V3=="Ciliophora" |V3=="Cercozoa")
#physeq_18S_filter1 <- microbiome::transform(physeq_18S_filter1,'compositional')
#merge phyloseq objects
physeq_merge1 <- merge_phyloseq(physeq_16S_filter2,physeq_18S_filter1)
physeq_merge2 <- subset_samples(physeq_merge1, month =="march" | month =="april")
asv_freq <- t(otu_table(physeq_merge2))
write.table(asv_freq, "mar_ap_freq.txt", sep = "\t")
asv_taxa <- tax_table(physeq_merge2)
write.table(asv_taxa, "mar_ap_taxa.txt", sep = "\t")
system("sed -i '1d' mar_ap_taxa.txt")

physeq_merge2 <- subset_samples(physeq_merge1, month =="may" | month =="june")
asv_freq <- t(otu_table(physeq_merge2))
write.table(asv_freq, "may_june_freq.txt", sep = "\t")
asv_taxa <- tax_table(physeq_merge2)
write.table(asv_taxa, "may_june_taxa.txt", sep = "\t")
system("sed -i '1d' may_june_taxa.txt")

physeq_merge2 <- subset_samples(physeq_merge1, month =="july" | month =="august")
asv_freq <- t(otu_table(physeq_merge2))
write.table(asv_freq, "july_aug_freq.txt", sep = "\t")
asv_taxa <- tax_table(physeq_merge2)
write.table(asv_taxa, "july_aug_taxa.txt", sep = "\t")
system("sed -i '1d' july_aug_taxa.txt")