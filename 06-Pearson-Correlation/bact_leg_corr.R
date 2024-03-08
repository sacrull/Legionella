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
physeq_16S_glom <- tax_glom(physeq_16S, taxrank=rank_names(physeq_16S)[7]) 
physeq_16S_glom1 <- subset_samples(physeq_16S_glom, month =="july" | month =="august")

ASV_freq <- as(otu_table(physeq_16S_glom1), "matrix")
#trans_ASV_freq <- t(ASV_freq)
#trans_ASV_freq_df1 <- as.data.frame(trans_ASV_freq)
#trans_ASV_freq_df2 <- rownames_to_column(trans_ASV_freq_df1, var = "ASV")

#write.table(OTU1, 'pred_host_freq.txt',sep = "\t",col.names=NA)

#getting taxa genus
taxa = as(tax_table(physeq_16S_glom1), "matrix")
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
filter3_BH_pear <- subset(filter2_BH_pear, origin == "Legionella" | destination == "Legionella")
filter4_BH_pear <- subset(filter3_BH_pear, value > .9)

pdf("BH_bact_july_august.pdf", width = 14, height= 14)
chordDiagram(filter4_BH_pear, order = ordered_list2, transparency = 0.5, col = ifelse(filter4_BH_pear$origin == "Legionella", ifelse(filter4_BH_pear$value > 0, ifelse(filter4_BH_pear$value > 0.7, "#0487fb", ifelse(filter4_BH_pear$value > 0.5, "#04F0FB", "#04FB38")), "red") , "gray"),
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