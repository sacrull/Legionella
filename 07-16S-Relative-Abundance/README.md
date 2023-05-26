# Relative Abudance
## 16S Relative Abundance
#### Load Libraries
```R
library(devtools)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(microbiome)
library(scales)
library(ggplot2)
```
#### Create Phyloseq object
##### ASV Table
```R
seqtab_16S <- read.table("feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
```
##### Taxonomy
```R
tax_16S <- read.table("taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)
```
##### Metadata
map_16S <- read.table("metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)
##### Merging phyloseq objects into 1
```R
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
```
#### Filter
Remove unassigned bacteria
```R
physeq_16S_filter1 = subset_taxa(physeq_16S, !V3=="Bacteria_unknown")
```
### By Month
#### Combine by month
```R
ps16_month <- merge_samples(physeq_16S_filter1, "month")
```
#### Get Relative abundance and Collapse by Phyla
Relative Abundance
```R
rel.abund_16_month <- transform_sample_counts(ps16_month, function(x) x/sum(x))
```
Collapse by phyla
```R
glom_16_month <- tax_glom(rel.abund_16_month, taxrank=rank_names(rel.abund_16_month)[2])
```
#### Rename Low Abundance Phyla
Less than 1% phyla renamed
```R
data$V3[data$Abundance < 0.01] <- "< 1% abund"
```
#### Convert to dataframe
```R
data_16_month <- psmelt(glom_16_month)
data_16_month$V3 <- as.character(data_16_month$V3) #convert character
```
####
#### Look at overall Median
```R
medians_16_month <- plyr::ddply(data_16_month, ~V3, function(x) c(median=median(x$Abundance)))
medians_16S_month
```
#### Get just useful columns
```R
rel_16S_month <- data_16_month[,c(2,3,38)]
```
#### Graph by month
```R
pdf("16S_rel_abudance_month.pdf", width=30, height=18)
ggplot(rel_16S)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V3),stat="identity", position="stack")+
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()
```
### Not by Month
#### Remove empty samples
```R
to_remove <- c("jul20_W_176") 
physeq_16S_filter2 <- prune_samples(!(sample_names(physeq_16S_filter1) %in% to_remove), physeq_16S_filter1)
```
#### Get relative abundance at phyla level
Get relative abundance
```R
rel.abund.16S_all <- transform_sample_counts(physeq_16S_filter2, function(x) x/sum(x))
```
Collapase to phyla level
```R
glom_16S_all <- tax_glom(rel.abund.16S_all, taxrank=rank_names(rel.abund.16S_all)[2], NArm=TRUE)
```
#### Convert to data frame
```R
data_16S_all <- psmelt(glom_16S_ all)
data_16S_all$V3 <- as.character(data_16S_all$V3) #convert character
```
#### Rename Low Frequency Abundance
```R
data_all$V3[data_all$Abundance < 0.01] <- "< 1% abund" # rename low freq phyla
```
#### Overall median
```R
medians_16S_all <- plyr::ddply(data_16S_all, ~V3, function(x) c(median=median(x$Abundance)))
medians_16S_month
```
#### Get just useful columns
```R
rel_16S_all <- data_16S_all[,c(2,3,38)]
```
#### Graph Samples
Sample Order
```R
sample_order <- c("mar15_W_38", "mar24_E_46", "mar24_W_47", "mar24_H2O_48", "mar30_E_55", "mar30_W_56", "mar30_H2O_57", "apr5_E_64", "apr5_W_65", "apr5_H2O_66", "apr13_E_73", "apr13_W_74", "apr13_H2O_75", "apr20_E_82", "apr20_W_83", "apr20_H2O_84", "apr27_W_92", "apr27_H2O_93", "may5_E_100", "may5_W_101", "may10_E_109", "may10_W_110", "may10_H2O_111", "may18_E_118", "may18_W_119", "may18_H2O_120", "may25_E_127", "may25_W_128", "may25_H2O_129", "jun1_E_136", "jun1_W_137", "jun1_H2O_138", "jun7_E_148", "jun7_W_149", "jun7_H2O_150", "jun21_E_157", "jun21_W_158", "jun21_H2O_159", "jul7_E_166", "jul7_W_167", "jul7_H2O_168", "jul20_E_175", "jul20_H2O_177", "aug3_E_184", "aug3_W_185", "aug3_H2O_186")
```
Graphing
```R
pdf("16S_rel_abudance_all.pdf", width=30, height=18)
ggplot(rel_16S_all)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V3),stat="identity", position="stack")+
  scale_x_discrete(limits = sample_order)+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()
```





