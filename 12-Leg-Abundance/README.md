# Legionella Abundance

## 1. Load Libraries
```R
library(devtools)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ape)
library(microbiome)
library(ggpubr)
library(rstatix)
```
## 2. Getting 16S data into phyloseq
```R
#getting 16S frequency table
seqtab_16S <- read.table("feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_16S)

#reading in taxonomy
tax_16S <- read.table("taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)

#reeading in metadata
map_16S <- read.table("metadata_16S.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)

#make 16S and 18S phyloseq object from qiime2
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
```
## 3. Transform Abundance
```R
physeq_16S_bp_clr <- microbiome::transform(physeq_16S, "clr")
```
## 4. Collapse genera together
```R
glom_bp <- tax_glom(physeq_16S_bp_clr, taxrank=rank_names(physeq_16S_bp_clr)[6]) #collapse by genus level
```
## 5.Get just Legionella
```R
physeq_16S_clr_filter_bp = subset_taxa(glom_bp, V7=="Legionella")
```
## 6. Get abundance data into dataframe
```R
data_16S_bp <- psmelt(physeq_16S_clr_filter_bp)
data_pb_2 <- data_16S_bp[, c(2,6,3)] #get just important columns
```
## 7. Month order
```R
month_order=c("march","april","may","june", "july","august")
```
## 8. Make boxplots
```R
pdf("leg_asv_boxplot_clr.pdf")
ggplot(data_pb_2, aes(x=factor(month, levels=month_order),y=Abundance))+
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "none") +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size=1)+
  labs(x ="Month", y = "CLR Abundance")+
  theme_classic()
dev.off()
```