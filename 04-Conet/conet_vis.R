library(devtools)
#library(spieceasi)
#library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
#library(ggplot2) # publication quality figures, based on ggplot2
#library(dplyr) # data handling
#library(network)
#library(intergraph)
#library(ggnet)
#library(igraph)
library(tidyverse)
#library(qiime2R)
library(readr)
library(viridis)
library(Matrix)
library(circlize) #https://mran.microsoft.com/snapshot/2016-01-05/web/packages/circlize/vignettes/visualize_relations_by_chord_diagram.pdf
library(reshape2)
library(chorddiag)
#library(phylosmith)


setwd("/Users/suzannecrull/Downloads")
#colors 
leg_pred_host_diff=c("Legionella"="plum3","Acanthamoeba" = "#D53E4F", "Echinamoeba"="#FDAE61", "Korotnevella" ="#FEE08B", "Naegleria"="#E6F598", "Tetrahymena"="#ABDDA4", "Vannella" = "#66C2A5","Vermamoeba"="#3288BD", 
      "Arcellinida_unknown" = "gray52", "BIO10-D10" = "gray52", "Dactylopodida" = "gray52", "Stygamoebida" = "gray52", "Euamoebida" = "gray52", "BOLA868" = "gray52", "Centramoebida" = "gray52", "Mycamoeba" = "gray52", "Vannella" = "gray52", "Euamoebida_unknown" = "gray52", "Vannellida" = "gray52", "Tubulinea_unknown" = "gray52", "Tubulinea" = "gray52", "uncultured" = "gray52", "Vannellida_unknown" = "gray52", "Amoebozoa_unknown" = "gray52", "Cryptodifflugia" = "gray52", "Amoebozoa" = "gray52", "Vermistella" = "gray52", "Protosteliopsis" = "gray52", "Arcellinida" = "gray52", "Arcella" = "gray52", "Gymnophrys" = "gray67", "Heteromita" = "gray67", "Cercozoa_unknown" = "gray67", "uncultured" = "gray67", "Paracercomonas" = "gray67", "Cercozoa" = "gray67", "Glissomonadida_unknown" = "gray67", "Vampyrellidae" = "gray67", "Tracheleuglypha" = "gray67", "Cercomonadidae" = "gray67", "Eocercomonas" = "gray67", "Kraken" = "gray67", "Euglypha" = "gray67", "Glissomonadida" = "gray67", "Trinema" = "gray67", "Thecofilosea_unknown" = "gray67", "Chilodonella" = "gray52", "Amphileptus" = "gray52", "Cyrtolophosis" = "gray52", "Leptopharynx" = "gray52", "Hymenostomatia" = "gray52", "Cyclidium" = "gray52", "Colpodea_unknown" = "gray52", "Hypotrichia_unknown" = "gray52", "Vorticella" = "gray52", "Protocyclidium" = "gray52", "Peritrichia" = "gray52", "Spirotrichea_unknown" = "gray52", "Oligohymenophorea" = "gray52", "Conthreep_unknown" = "gray52", "Nassophorea" = "gray52", "Colpodida" = "gray52", "Telotrochidium" = "gray52", "Nassophorea_unknown" = "gray52", "Oligohymenophorea_unknown" = "gray52", "Aspidisca" = "gray52", "Haptoria_unknown" = "gray52", "Ephelota" = "gray52", "Cyrtolophosidida" = "gray52", "Allovahlkampfia" = "gray52", "Tetramitia_unknown" = "gray52", "Vahlkampfia" = "gray52", "Neovahlkampfia" = "gray52", "Cercomonadidae_unknown" = "gray67", "Euglyphida_unknown" = "gray67", "Hypotrichia_unknown" = "gray52", "Oligohymenphorea_unknown"="gray52", "Oligohymenophorea_unkown"="gray52", "Vampyrellidae_unknown"="gray67")
#reading sequence tables)
#reading sequence tables
seqtab_16S <- read.table("feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_16S)
seqtab_18S <- read.table("feature-table-18S.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
tax_16S <- read.table("taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)
tax_18S <- read.table("taxonomy-18S-clean2.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
#reeading in metadata
map_16S <- read.table("metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)
map_18S <- read.table("metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
#make 16S and 18S phyloseq object from qiime2
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
#physeq_16S_filter1 = subset_taxa(physeq_16S, V5=="Legionellales")
physeq_16S_filter2 = subset_taxa(physeq_16S, V7=="Legionella")
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
#known hosts and predator of legionella
physeq_18S_filter1 = subset_taxa(physeq_18S, V3=="Amoebozoa" | V3=="Heterolobosea" | V3=="Ciliophora" |V3=="Cercozoa")
#merge phyloseq objects
physeq_merge1 <- merge_phyloseq(physeq_16S_filter2,physeq_18S_filter1)
physeq_merge2 <- subset_samples(physeq_merge1, month =="july" | month =="august")

#get asv freqyuency matrix from 
ASV_freq <- as(otu_table(physeq_merge2), "matrix")


#get genus level names
taxa = as(tax_table(physeq_merge2), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V7)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
#read in file
df <- read.table("july_aug_edges.csv", header=TRUE, sep=",")
#renaming to genus level
df1 <- left_join(df, orderdf, by=c('to'='ASV'))
colnames(df1)[4] <- "origin"
df2 <- left_join(df1, orderdf, by=c('from'='ASV'))
colnames(df2)[5] <- "destination"
df3 <- df2[, c(4,5,3)]
colnames(df3)[3] <- "value"

#order list
list_all<-as.list(orderdf)
ordered_list <- unique(sort(list_all$V7))
ordered_list2 <- c("Legionella", ordered_list)
#making chord diagram

pdf("july_aug_ordered_conet_clean.pdf")
chordDiagram(df3, order = ordered_list2, transparency = 0.5, grid.col = leg_pred_host_diff, col = ifelse(df3$origin == "Legionella", ifelse(df3$value > 0, ifelse(df3$value > 0.7, "#0487fb", ifelse(df3$value > 0.5, "#04F0FB", "#04FB38")), "red") , 
	ifelse(df3$destination == "Legionella",  ifelse(df3$value > 0, ifelse(df3$value > 0.7, "#0487fb", ifelse(df3$value > 0.5, "#04F0FB", "#04FB38")), "red"), "gray")),
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


