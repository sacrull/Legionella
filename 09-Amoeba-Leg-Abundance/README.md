# Relative Abundance PLot of Amoeba and 
## 1. Load in packages
```R
library(phyloseq)
library(tidyverse)
library(dplyr)
library(microbiome)
library(scales)
library(ggplot2)
```
## 2. Get 18S data into phyloseq
### Frequency table
```R
seqtab_18S <- read.table("18S-ASV-renamed.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
```
### Taxonomy
```R
tax_18S <- read.table("18S-tax-renamed.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
```
### Metadata
```R
map_18S <- read.table("metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
```
### Merging 18S object and Filtering
Getting only known predators and host phyla of legionella
```R
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
physeq_18S_filter1 = subset_taxa(physeq_18S, V3=="Amoebozoa" | V3=="Heterolobosea" | V3=="Ciliophora" |V3=="Cercozoa")
```
### Remove empty samples
```R
to_remove <- c("mar30_E_55", "mar30_H2O_57") #removed those samples since they didnt have any OTUs of interest
physeq_18S_filter2 <- prune_samples(!(sample_names(physeq_18S_filter1) %in% to_remove), physeq_18S_filter1)
```
### Merge Samples by Month
```R
ps1 <- merge_samples(physeq_18S_filter2, "month") #combine by month
```
## 3. 18S Relative Abundance
### Relative Abundance
```R
rel.abund <- transform_sample_counts(ps1, function(x) x/sum(x))
```
Collapse by Genus Level
```R
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[6]) #collapse by genus level
```
### Convert to Dataframe
```R
data <- psmelt(glom)
data$V7 <- as.character(data$V7) #convert character
```
## 4. Renaming Low Frequency ASVs
```R
data$V7[data$Abundance < 0.01] <- "< 1% abund" # rename low freq phyla
```
## 6. Overall Medians
```R
medians <- plyr::ddply(data, ~V7, function(x) c(median=median(x$Abundance)))
medians
```
## 7. Getting Just Useful 18S Columns
```R
rel_18S <- data[,c(2,3,42)]
```
## 8. Getting 16S Phyloseq Object
### Frequency Table
```R
seqtab_16S <- read.table("feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
```
### Taxonomy
```R
tax_16S <- read.table("taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
```
### Metadata
```R
map_16S <- read.table("metadata_16S.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)
```
### Merging 16S Objet
```R
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo
```
### Merge Samples by Month
```R
ps16 <- merge_samples(physeq_16S, "month") #combine by month
```
## 9. Get Relative Abundance
### Relative Abundance
```R
rel.abund.16S <- transform_sample_counts(ps16, function(x) x/sum(x))
```
### Collapase by Genus
```R
glom_16S <- tax_glom(rel.abund.16S, taxrank=rank_names(rel.abund.16S)[6]) #collapse by genus level
```
### Make into Dataframe
```R
data_16S <- psmelt(glom_16S)
data_16S$V7 <- as.character(data_16S$V7) #convert character
```
## 10. Get Legionella Relative Abundance
### Get just Legionella Relative Abundane
```R
leg <- data_16S[data_16S$V7 == 'Legionella',]
```
### Get just needed columns
```R
leg_rel <- leg[, c(2,3,42)]
```
### Rename Columns
```R
colnames(leg_rel)[2] <- "Abundance_leg"
colnames(leg_rel)[3] <- "Genus"
```
## 11. Combine 18S and Legionella Dataframes
```R
combined_leg_18S <- left_join(rel_18S, leg_rel, by = join_by(Sample == Sample))
```
## 12. Graph
### Colors for gtraph
```R
leg_pred_host_diff=c("Legionella"="plum3","Acanthamoeba" = "#D53E4F", "Echinamoebida"="#FDAE61", "Korotnevella" ="#FEE08B", "Naegleria"="#E6F598", "Tetrahymena"="#ABDDA4", "Vannella" = "#66C2A5","Vermamoeba"="#3288BD", 
      "Arcellinida_unknown" = "gray52", "BIO10-D10" = "gray52", "Dactylopodida" = "gray52", "Stygamoebida" = "gray52", "Euamoebida" = "gray52", "BOLA868" = "gray52", "Centramoebida" = "gray52", "Mycamoeba" = "gray52", "Vannella" = "gray52", "Euamoebida_unknown" = "gray52", "Vannellida" = "gray52", "Tubulinea_unknown" = "gray52", "Tubulinea" = "gray52", "uncultured" = "gray52", "Vannellida_unknown" = "gray52", "Amoebozoa_unknown" = "gray52", "Cryptodifflugia" = "gray52", "Amoebozoa" = "gray52", "Vermistella" = "gray52", "Protosteliopsis" = "gray52", "Arcellinida" = "gray52", "Arcella" = "gray52", "Gymnophrys" = "gray67", "Heteromita" = "gray67", "Cercozoa_unknown" = "gray67", "uncultured" = "gray67", "Paracercomonas" = "gray67", "Cercozoa" = "gray67", "Glissomonadida_unknown" = "gray67", "Vampyrellidae" = "gray67", "Tracheleuglypha" = "gray67", "Cercomonadidae" = "gray67", "Eocercomonas" = "gray67", "Kraken" = "gray67", "Euglypha" = "gray67", "Glissomonadida" = "gray67", "Trinema" = "gray67", "Thecofilosea_unknown" = "gray67", "Chilodonella" = "gray52", "Amphileptus" = "gray52", "Cyrtolophosis" = "gray52", "Leptopharynx" = "gray52", "Hymenostomatia" = "gray52", "Cyclidium" = "gray52", "Colpodea_unknown" = "gray52", "Hypotrichia_unknown" = "gray52", "Vorticella" = "gray52", "Protocyclidium" = "gray52", "Peritrichia" = "gray52", "Spirotrichea_unknown" = "gray52", "Oligohymenophorea" = "gray52", "Conthreep_unknown" = "gray52", "Nassophorea" = "gray52", "Colpodida" = "gray52", "Telotrochidium" = "gray52", "Nassophorea_unknown" = "gray52", "Oligohymenophorea_unknown" = "gray52", "Aspidisca" = "gray52", "Haptoria_unknown" = "gray52", "Ephelota" = "gray52", "Cyrtolophosidida" = "gray52", "Allovahlkampfia" = "gray52", "Tetramitia_unknown" = "gray52", "Vahlkampfia" = "gray52", "Neovahlkampfia" = "gray52", "Cercomonadidae_unknown" = "gray67", "Euglyphida_unknown" = "gray67", "Hypotrichia_unknown" = "gray52", "Oligohymenphorea_unknown"="gray52", "Oligohymenophorea_unkown"="gray52", "Vampyrellidae_unknown"="gray67")
```
### Graphinh
```R
pdf("pred_host_leg_rel_abund.pdf", width=30, height=18)
ggplot(combined_leg_18S)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V7),stat="identity", position="stack")+
  scale_fill_manual(values = leg_pred_host_diff)+ #"Vahlkampfia" = "#66C2A5" is a maybe
  geom_line(aes(x=Sample, y=100*Abundance_leg, group=1),stat="identity",color="red",size=5)+
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()
```


















