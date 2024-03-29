rm(list = ls())
library(devtools)
#library(spieceasi)
#library(microbiome) # data analysis and visualisation
library(tidyverse)
library(dplyr)
library(phyloseq)
library(microbiome)
library(scales)
library(ggplot2)
library(metagMisc)

setwd("~/legionella/R/rel_abudance")
#reading sequence tables
seqtab_16S <- read.table("../feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_16S)

#reading in taxonomy
tax_16S <- read.table("../taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)

#reeading in metadata
map_16S <- read.table("../metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)

physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
physeq_16S1 <- prune_samples(sample_sums(physeq_16S) > 10000, physeq_16S) #remove less than 2000 reads
rare_16S <- rarefy_even_depth(physeq_16S1, rngseed=1, sample.size=0.99*min(sample_sums(physeq_16S1)), replace=F)
sample_sums(rare_16S)
#remove unassigned bacter
physeq_16S_filter1 = subset_taxa(rare_16S, !V3=="Bacteria_unknown")

ps1 <- merge_samples(physeq_16S_filter1, "month") #combine by month

#get relative abundance 
rel.abund <- transform_sample_counts(ps1, function(x) x/sum(x))
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[2]) #collapse by genus level
data <- psmelt(glom)
data$V3 <- as.character(data$V3) #convert character
data$V3[data$Abundance < 0.01] <- "< 1% abund" # rename low freq phyla
#looking at overall median
medians <- plyr::ddply(data, ~V3, function(x) c(median=median(x$Abundance)))
medians
#get just useful columns
rel_16S <- data[,c(2,3,39)]

subset(rel_16S, V3 == "Bacteroidetes")

#sanity check
sum(data[data$Sample == 'april',]$Abundance)


pdf("16S_rel_abudance_month.pdf", width=30, height=18)
ggplot(rel_16S)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V3),stat="identity", position="stack")+
  theme_minimal() +
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+
  labs(fill='Phylum') +
  theme(axis.text = element_text(size = 23), axis.title = element_text(size = 30),legend.text = element_text(size = 25),legend.title = element_text(size = 30))
dev.off()




##not by month
rel.abund.all <- transform_sample_counts(rare_16S, function(x) x/sum(x))
glom_all <- tax_glom(rel.abund.all, taxrank=rank_names(rel.abund.all)[2], NArm=TRUE) #collapse by phlyum level
data_all <- psmelt(glom_all)
data_all$V3 <- as.character(data_all$V3) #convert character
data_all$V3[data_all$Abundance < 0.01] <- "< 1% abund" # rename low freq phyla
#looking at overall median
medians <- plyr::ddply(data_all, ~V3, function(x) c(median=median(x$Abundance)))
medians
#get just useful columns
rel_16S_all <- data_all[,c(2,3,39)]


#sanity check
sum(data_all[data_all$Sample == 'april',]$Abundance)
sample_order <- c("mar15_W_38", "mar24_E_46", "mar24_W_47", "mar24_H2O_48", "mar30_E_55", "mar30_W_56", "mar30_H2O_57", "apr5_E_64", "apr5_W_65", "apr5_H2O_66", "apr13_E_73", "apr13_W_74", "apr13_H2O_75", "apr20_E_82", "apr20_W_83", "apr20_H2O_84", "apr27_W_92", "apr27_H2O_93", "may5_E_100", "may5_W_101", "may10_E_109", "may10_W_110", "may10_H2O_111", "may18_E_118", "may18_W_119", "may18_H2O_120", "may25_E_127", "may25_W_128", "may25_H2O_129", "jun1_E_136", "jun1_W_137", "jun1_H2O_138", "jun7_E_148", "jun7_W_149", "jun7_H2O_150", "jun21_E_157", "jun21_W_158", "jun21_H2O_159", "jul7_E_166", "jul7_W_167", "jul7_H2O_168", "jul20_E_175", "jul20_H2O_177", "aug3_E_184", "aug3_W_185", "aug3_H2O_186")

pdf("16S_rel_abudance_all.pdf", width=30, height=18)
ggplot(rel_16S_all)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V3),stat="identity", position="stack")+
  theme_minimal() +
  scale_x_discrete(limits = sample_order,guide = guide_axis(angle = 90))+
  labs(fill='Phylum') +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20),legend.text = element_text(size = 18),legend.title = element_text(size = 20))
dev.off()



