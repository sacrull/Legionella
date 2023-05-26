# Relative Abudance
## 18S Relative Abundance
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
seqtab_18S <- read.table("18S-ASV-renamed.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
```
##### Taxonomy
```R
tax_18S <- read.table("18S-tax-renamed.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
```
##### Metadata
map_18S <- read.table("metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
##### Merging phyloseq objects into 1
```R
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
```
#### Filter
Remove unassigned Eukryotes and bacteria
```R
physeq_18S_filter1 = subset_taxa(physeq_18S, !V2=="Bacteria", !V3=="Eukaryota_unknown")
```
### By Month
#### Combine by month
```R
ps18_month <- merge_samples(physeq_18S_filter1, "month")
```
#### Get Relative abundance and Collapse by Phyla
Relative Abundance
```R
rel.abund_18_month <- transform_sample_counts(ps18_month, function(x) x/sum(x))
```
Collapse by phyla
```R
glom_18_month <- tax_glom(rel.abund_18_month, taxrank=rank_names(rel.abund_18_month)[2])
```
#### Rename Low Abundance Phyla
Less than 1% phyla renamed
```R
data$V3[data$Abundance < 0.01] <- "< 1% abund"
```
#### Convert to dataframe
```R
data_18_month <- psmelt(glom_18_month)
data_18_month$V3 <- as.character(data_18_month$V3) #convert character
```
####
#### Look at overall Median
```R
medians_18_month <- plyr::ddply(data_18_month, ~V3, function(x) c(median=median(x$Abundance)))
medians_18S_month
```
#### Get just useful columns
```R
rel_18S_month <- data_18_month[,c(2,3,38)]
```
#### Graph by month
```R
pdf("18S_rel_abudance_month.pdf", width=30, height=18)
ggplot(rel_18S)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V3),stat="identity", position="stack")+
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()
```
### Not by Month
#### Get relative abundance at phyla level
Get relative abundance
```R
rel.abund.18S_all <- transform_sample_counts(physeq_18S_filter1, function(x) x/sum(x))
```
Collapase to phyla level
```R
glom_18S_all <- tax_glom(rel.abund.18S_all, taxrank=rank_names(rel.abund.18S_all)[2], NArm=TRUE)
```
#### Convert to data frame
```R
data_18S_all <- psmelt(glom_18S_ all)
data_18S_all$V3 <- as.character(data_18S_all$V3) #convert character
```
#### Rename Low Frequency Abundance
```R
data_all$V3[data_all$Abundance < 0.01] <- "< 1% abund" # rename low freq phyla
```
#### Overall median
```R
medians_18S_all <- plyr::ddply(data_18S_all, ~V3, function(x) c(median=median(x$Abundance)))
medians_18S_month
```
#### Get just useful columns
```R
rel_18S_all <- data_18S_all[,c(2,3,38)]
```
#### Graph Samples
Sample Order
```R
sample_order <- c("mar15_W_38", "mar24_E_46", "mar24_W_47", "mar24_H2O_48", "mar30_E_55", "mar30_W_56", "mar30_H2O_57", "apr5_E_64", "apr5_W_65", "apr5_H2O_66", "apr13_E_73", "apr13_W_74", "apr13_H2O_75", "apr20_E_82", "apr20_W_83", "apr20_H2O_84", "apr27_W_92", "apr27_H2O_93", "may5_E_100", "may5_W_101", "may10_E_109", "may10_W_110", "may10_H2O_111", "may18_E_118", "may18_W_119", "may18_H2O_120", "may25_E_127", "may25_W_128", "may25_H2O_129", "jun1_E_136", "jun1_W_137", "jun1_H2O_138", "jun7_E_148", "jun7_W_149", "jun7_H2O_150", "jun21_E_157", "jun21_W_158", "jun21_H2O_159", "jul7_E_166", "jul7_W_167", "jul7_H2O_168", "jul20_E_175", "jul20_W_176", "jul20_H2O_177", "aug3_E_184", "aug3_W_185", "aug3_H2O_186")
```
Graphing
```R
pdf("18S_rel_abudance_all.pdf", width=30, height=18)
ggplot(rel_18S_all)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V3),stat="identity", position="stack")+
  scale_x_discrete(limits = sample_order)+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()
```





