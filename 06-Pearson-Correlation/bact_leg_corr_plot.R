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
library(microbiome)
library(tabletools)
library(metagMisc)


setwd("/home/suzanne/legionella/R/pearson")
#reading sequence tables
seqtab_16S <- read.table("../feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)

#reading in taxonomy
tax_16S <- read.table("../taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)

#reeading in metadata
map_16S <- read.table("../metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)

physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
physeq_16S <- subset_taxa(physeq_16S, !V7=="Bacteria_unknown")
physeq_16S_glom <- tax_glom(physeq_16S, taxrank=rank_names(physeq_16S)[6]) 
physeq_16S_glom1 <- subset_samples(physeq_16S_glom, month =="march" | month =="april")

#getting just genus level dataframe
names <- as.data.frame(tax_table(physeq_16S_glom1)[,"V7"])
names$names <- rownames(names)

ASV_freq <- t(as(otu_table(physeq_16S_glom1), "matrix"))
ASV_freq <- as.data.frame(ASV_freq)
ASV_freq$names <- rownames(ASV_freq)
genus_freq <- inner_join(names, ASV_freq, by=c('names'='names'))
rownames(genus_freq) <- genus_freq[,1]
genus_freq2 <- genus_freq %>% select( 3:ncol(.) )
genus_freq3 <- as.matrix(t(genus_freq2))

#running correlation
ASV_pear_corr1 <- rcorr(genus_freq3,type=c("pearson"))


cor_enviro_adjust <- rcorr_padjust(ASV_pear_corr1, method = "BH")
diag(cor_enviro_adjust$P) <- 0


pdf("leg_bact_corr.pdf", height = 100 , width = 100)
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off()