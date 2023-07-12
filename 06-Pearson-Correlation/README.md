# Pearson Correlation
## 1. Make a New Directory and Copy Files
Make a new directory for pearson correlation
```bash
cd ~/legionella/R/
mkdir pearson
```
Copy your metadata, frequency table, and taxonomy files for the 16S and 18S data into the new directory
```bash
cp *.txt ./pearson
```
## 2. Install the packages
Create a conda enviroment and open up R
```bash
conda create -n pearson
conda activate pearson
conda install -c conda-forge r-base
conda install r-sf
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install r-devtools
R
```
Install the packages in R
```R
library(devtools)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("tidyverse")
install.packages("RColorBrewer")
install.packages("Hmisc")
install.packages("circlize")
```
## 3. Load Packages
```R
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
library(Hmisc)
library(circlize)
```
## 4. Get data into phyloseq
### Frequency Tables
```R
seqtab_16S <- read.table("../feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_16S)
seqtab_18S <- read.table("../18S-ASV-renamed.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
```
### Taxonomy
```R
tax_16S <- read.table("../taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)
tax_18S <- read.table("../18S-tax-renamed.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
```
### Metadata
```R
#reeading in metadata
map_16S <- read.table("../metadata_16S.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)
map_18S <- read.table("../metadata_18S.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
```
### Combine Phyloseq Objects and Merge
Combine 16S data into 1 object and filter out everything, but was is assigned to Legionella at the genus level
```R
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
physeq_16S_filter1 = subset_taxa(physeq_16S, V7=="Legionella")
```
Combine 18S data into 1 phyloseq object and filter out everything that does not belong to phylums identified as haveing a relationship with Legionella (Amoebozoa, Heterolobosea, Ciliophora, Cercozoa)
```R
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
physeq_18S_filter1 = subset_taxa(physeq_18S, V3=="Amoebozoa" | V3=="Heterolobosea" | V3=="Ciliophora" |V3=="Cercozoa")
```
Merge phyloseq objects together. In this example I filter it by March and April and I would then change it out to May and June then July and August for those months. 
```R
physeq_merge1 <- merge_phyloseq(physeq_16S_filter1,physeq_18S_filter1)
physeq_merge2 <- subset_samples(physeq_merge1, month =="march" | month =="april")
```
## 5. Get information out of Phyloseq
### Get merged ASV frequency phyloseq object as a matrix
```R
ASV_freq <- as(otu_table(physeq_merge2), "matrix")
```
### Get taxonomy of ASVs
Getting the taxa from the genus level out of the phyloseq object
```R
taxa = as(tax_table(physeq_merge2), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V7)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
```
Order the chord diagram in alphabetical order for the amoebas. This for later in the chord diagram.
```R
list_all<-as.list(orderdf)
ordered_list <- unique(sort(list_all$V7))
ordered_list2 <- c("Legionella", ordered_list)
```
## 7. Pearson Correlation
Running the correlation from the matrix
```R
ASV_pear_corr1 <- rcorr(ASV_freq,type=c("pearson"))
```
Creating a function for the next step
```R
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
```
Converting matrix from wide format to long format with correlation values and P values
```R
flatten_ASV_pear1 <- flattenCorrMatrix(ASV_pear_corr1$r, ASV_pear_corr1$P)
```
## 8. Renaming ASVs to Genus Level
Getting the genus names for the first column
```R
flatten_ASV_pear2 <- left_join(flatten_ASV_pear1, orderdf, by=c('row'='ASV'))
colnames(flatten_ASV_pear2)[5] <- "origin"
```
Getting genus names for the second column
```R
flatten_ASV_pear3 <- left_join(flatten_ASV_pear2, orderdf, by=c('column'='ASV'))
colnames(flatten_ASV_pear3)[6] <- "destination"
```
Getting just the genus level names, correlation values, and p-values
```R
flatten_ASV_pear4 <- flatten_ASV_pear3[, c(5,6,3,4)]
colnames(flatten_ASV_pear4)[3] <- "value"
```
Order the p-values from smallest to largest
```R
flatten_ASV_pear4 = flatten_ASV_pear4[order(flatten_ASV_pear4$p),]
```
## 9. Benjamini Hochberg Correction
Adjust p-values using benjamini hochberg correction
```R
flatten_ASV_pear4$BH =
      p.adjust(flatten_ASV_pear4$p,
               method = "BH")
```
Get just the new p value
```R
flatten_ASV_pear5 <- flatten_ASV_pear4[, c(1,2,3,5)]
```
Get p-values only less than 0.05
```R
filter1_BH_pear <- subset(flatten_ASV_pear5, BH < .05) 
```
Get just the genus level names (origin and destination) and the correlation coefficent
```R
filter2_BH_pear <- filter1_BH_pear[, c(1:3)]
```
## 10. Plot Chorddiagram
Get colors. Host phyla are dark gray. Predator phyla are light gray. Genera with previously established host connection are colored along with Legionella.
```R
leg_pred_host_diff=c("Legionella"="plum3","Acanthamoeba" = "#D53E4F", "Echinamoeba"="#FDAE61", "Korotnevella" ="#FEE08B", "Naegleria"="#E6F598", "Tetrahymena"="#ABDDA4", "Vannella" = "#66C2A5","Vermamoeba"="#3288BD", 
      "Arcellinida_unknown" = "gray52", "BIO10-D10" = "gray52", "Dactylopodida" = "gray52", "Stygamoebida" = "gray52", "Euamoebida" = "gray52", "BOLA868" = "gray52", "Centramoebida" = "gray52", "Mycamoeba" = "gray52", "Vannella" = "gray52", "Euamoebida_unknown" = "gray52", "Vannellida" = "gray52", "Tubulinea_unknown" = "gray52", "Tubulinea" = "gray52", "uncultured" = "gray52", "Vannellida_unknown" = "gray52", "Amoebozoa_unknown" = "gray52", "Cryptodifflugia" = "gray52", "Amoebozoa" = "gray52", "Vermistella" = "gray52", "Protosteliopsis" = "gray52", "Arcellinida" = "gray52", "Arcella" = "gray52", "Gymnophrys" = "gray67", "Heteromita" = "gray67", "Cercozoa_unknown" = "gray67", "uncultured" = "gray67", "Paracercomonas" = "gray67", "Cercozoa" = "gray67", "Glissomonadida_unknown" = "gray67", "Vampyrellidae" = "gray67", "Tracheleuglypha" = "gray67", "Cercomonadidae" = "gray67", "Eocercomonas" = "gray67", "Kraken" = "gray67", "Euglypha" = "gray67", "Glissomonadida" = "gray67", "Trinema" = "gray67", "Thecofilosea_unknown" = "gray67", "Chilodonella" = "gray52", "Amphileptus" = "gray52", "Cyrtolophosis" = "gray52", "Leptopharynx" = "gray52", "Hymenostomatia" = "gray52", "Cyclidium" = "gray52", "Colpodea_unknown" = "gray52", "Hypotrichia_unknown" = "gray52", "Vorticella" = "gray52", "Protocyclidium" = "gray52", "Peritrichia" = "gray52", "Spirotrichea_unknown" = "gray52", "Oligohymenophorea" = "gray52", "Conthreep_unknown" = "gray52", "Nassophorea" = "gray52", "Colpodida" = "gray52", "Telotrochidium" = "gray52", "Nassophorea_unknown" = "gray52", "Oligohymenophorea_unknown" = "gray52", "Aspidisca" = "gray52", "Haptoria_unknown" = "gray52", "Ephelota" = "gray52", "Cyrtolophosidida" = "gray52", "Allovahlkampfia" = "gray52", "Tetramitia_unknown" = "gray52", "Vahlkampfia" = "gray52", "Neovahlkampfia" = "gray52", "Cercomonadidae_unknown" = "gray67", "Euglyphida_unknown" = "gray67", "Hypotrichia_unknown" = "gray52", "Oligohymenphorea_unknown"="gray52", "Oligohymenophorea_unkown"="gray52", "Vampyrellidae_unknown"="gray67", "Platyamoeba"="gray52")
#
```
Plot diagram
Chords that are connected to legionella are colored based on the correlation strength. It is ordered based on alphabetical order. The colors of the edges were assigned in the previous step. Can change pdf name depending on months. 
```R
pdf("BH_mar_ap_pear_ordered_clean.pdf")
chordDiagram(filter2_BH_pear, order = ordered_list2, transparency = 0.5, grid.col = leg_pred_host_diff, col = ifelse(filter2_BH_pear$origin == "Legionella", ifelse(filter2_BH_pear$value > 0, ifelse(filter2_BH_pear$value > 0.7, "#0487fb", ifelse(filter2_BH_pear$value > 0.5, "#04F0FB", "#04FB38")), "red") , "gray"),
	annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
	xlim = get.cell.meta.data("xlim")
	xplot = get.cell.meta.data("xplot")
	ylim = get.cell.meta.data("ylim")
	sector.name = get.cell.meta.data("sector.index")
	if(abs(xplot[2] - xplot[1]) < 20) {
	circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
	niceFacing = TRUE, adj = c(0, 0.5), cex=.29)
	} else {
	circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
	niceFacing = TRUE, adj = c(0.5, 0), cex=.29)
	}
}, bg.border = NA)
dev.off()
```