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
library(microbiome)


setwd("/home/suzanne/legionella/R/pearson")
leg_pred_host_diff=c("Legionella"="plum3","Acanthamoeba" = "#D53E4F", "Echinamoeba"="#FDAE61", "Korotnevella" ="#FEE08B", "Naegleria"="#E6F598", "Tetrahymena"="#ABDDA4", "Vannella" = "#66C2A5","Vermamoeba"="#3288BD", 
      "Arcellinida_unknown" = "gray52", "BIO10-D10" = "gray52", "Dactylopodida" = "gray52", "Stygamoebida" = "gray52", "Euamoebida" = "gray52", "BOLA868" = "gray52", "Centramoebida" = "gray52", "Mycamoeba" = "gray52", "Vannella" = "gray52", "Euamoebida_unknown" = "gray52", "Vannellida" = "gray52", "Tubulinea_unknown" = "gray52", "Tubulinea" = "gray52", "uncultured" = "gray52", "Vannellida_unknown" = "gray52", "Amoebozoa_unknown" = "gray52", "Cryptodifflugia" = "gray52", "Amoebozoa" = "gray52", "Vermistella" = "gray52", "Protosteliopsis" = "gray52", "Arcellinida" = "gray52", "Arcella" = "gray52", "Gymnophrys" = "gray67", "Heteromita" = "gray67", "Cercozoa_unknown" = "gray67", "uncultured" = "gray67", "Paracercomonas" = "gray67", "Cercozoa" = "gray67", "Glissomonadida_unknown" = "gray67", "Vampyrellidae" = "gray67", "Tracheleuglypha" = "gray67", "Cercomonadidae" = "gray67", "Eocercomonas" = "gray67", "Kraken" = "gray67", "Euglypha" = "gray67", "Glissomonadida" = "gray67", "Trinema" = "gray67", "Thecofilosea_unknown" = "gray67", "Chilodonella" = "gray52", "Amphileptus" = "gray52", "Cyrtolophosis" = "gray52", "Leptopharynx" = "gray52", "Hymenostomatia" = "gray52", "Cyclidium" = "gray52", "Colpodea_unknown" = "gray52", "Hypotrichia_unknown" = "gray52", "Vorticella" = "gray52", "Protocyclidium" = "gray52", "Peritrichia" = "gray52", "Spirotrichea_unknown" = "gray52", "Oligohymenophorea" = "gray52", "Conthreep_unknown" = "gray52", "Nassophorea" = "gray52", "Colpodida" = "gray52", "Telotrochidium" = "gray52", "Nassophorea_unknown" = "gray52", "Oligohymenophorea_unknown" = "gray52", "Aspidisca" = "gray52", "Haptoria_unknown" = "gray52", "Ephelota" = "gray52", "Cyrtolophosidida" = "gray52", "Allovahlkampfia" = "gray52", "Tetramitia_unknown" = "gray52", "Vahlkampfia" = "gray52", "Neovahlkampfia" = "gray52", "Cercomonadidae_unknown" = "gray67", "Euglyphida_unknown" = "gray67", "Hypotrichia_unknown" = "gray52", "Oligohymenphorea_unknown"="gray52", "Oligohymenophorea_unkown"="gray52", "Vampyrellidae_unknown"="gray67", "Platyamoeba"="gray52")
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
#physeq_16S_filter1 = subset_taxa(physeq_16S, V5=="Legionellales")
physeq_16S_filter2 <- microbiome::transform(physeq_16S,'compositional')
physeq_16S_filter2 = subset_taxa(physeq_16S_filter2, V7=="Legionella")
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
#known hosts and predator of legionella
physeq_18S_filter1 = subset_taxa(physeq_18S, V3=="Amoebozoa" | V3=="Heterolobosea" | V3=="Ciliophora" |V3=="Cercozoa")
physeq_18S_filter1 <- microbiome::transform(physeq_18S_filter1,'compositional')
#merge phyloseq objects
physeq_merge1 <- merge_phyloseq(physeq_16S_filter2,physeq_18S_filter1)
physeq_merge2 <- subset_samples(physeq_merge1, month =="may" | month =="june")

#getting ASV frequency from physeq object
ASV_freq <- as(otu_table(physeq_merge2), "matrix")
#trans_ASV_freq <- t(ASV_freq)
#trans_ASV_freq_df1 <- as.data.frame(trans_ASV_freq)
#trans_ASV_freq_df2 <- rownames_to_column(trans_ASV_freq_df1, var = "ASV")

#write.table(OTU1, 'pred_host_freq.txt',sep = "\t",col.names=NA)

#getting taxa genus
taxa = as(tax_table(physeq_merge2), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V7)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")

#order taxa for chord diagram
list_all<-as.list(orderdf)
ordered_list <- unique(sort(list_all$V7))
ordered_list2 <- c("Legionella", ordered_list)

#renaming ASVs to genus level
	#renamed_ASV_freq1 <- left_join(trans_ASV_freq_df2, orderdf, by=c('ASV'='ASV'))
	#colnames(renamed_ASV_freq1)[49] <- "Genus"
	#renamed_ASV_freq2 <- renamed_ASV_freq1[, c(49,2:48)]
	#renamed_ASV_freq3 <- renamed_ASV_freq2 %>%
	#  group_by(Genus) %>%
	#  summarise(across(c(1:47), sum))
	#renamed_ASV_freq4 <- as.data.frame(renamed_ASV_freq3)
	#rownames(renamed_ASV_freq4) <- renamed_ASV_freq4[,1]
	#renamed_ASV_freq4[,1] <- NULL
	#renamed_ASV_freq5 <- t(as.matrix(renamed_ASV_freq4))

#running correlation
ASV_pear_corr1 <- rcorr(ASV_freq,type=c("pearson"))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
flatten_ASV_pear1 <- flattenCorrMatrix(ASV_pear_corr1$r, ASV_pear_corr1$P)

#renaming to genus level
flatten_ASV_pear2 <- left_join(flatten_ASV_pear1, orderdf, by=c('row'='ASV'))
colnames(flatten_ASV_pear2)[5] <- "origin"
flatten_ASV_pear3 <- left_join(flatten_ASV_pear2, orderdf, by=c('column'='ASV'))
colnames(flatten_ASV_pear3)[6] <- "destination"
flatten_ASV_pear4 <- flatten_ASV_pear3[, c(5,6,3,4)]
colnames(flatten_ASV_pear4)[3] <- "value"

flatten_ASV_pear4 = flatten_ASV_pear4[order(flatten_ASV_pear4$p),]
#headtail(flatten_ASV_pear4)


flatten_ASV_pear4$BH =
      p.adjust(flatten_ASV_pear4$p,
               method = "BH")

##sanity check
	#flatten_ASV_pear4[flatten_ASV_pear4$destination == "Vannella", ]

#false discovery rate 

flatten_ASV_pear5 <- flatten_ASV_pear4[, c(1,2,3,5)]
filter1_BH_pear <- subset(flatten_ASV_pear5, BH < .05) 
filter2_BH_pear <- filter1_BH_pear[, c(1:3)]


pdf("BH_may_june_pear_ordered_clean.pdf", width = 14, height= 14)
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