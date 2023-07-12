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


setwd("/home/suzanne/legionella/silva_plots/spieceasi/pred_host_am")
#reading sequence tables
seqtab_18S <- read.table("feature-table-18S.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)

#reading in taxonomy
tax_18S <- read.table("taxonomy_18S_clean.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)

#reeading in metadata
map_18S <- read.table("metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)

#physeq_16S_filter1 = subset_taxa(physeq_16S, V5=="Legionellales")
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
#known hosts and predator of legionella
physeq_18S_filter1 = subset_taxa(physeq_18S, V3=="Amoebozoa" | V3=="Heterolobosea" | V3=="Ciliophora" |V3=="Cercozoa")

#remove empty samples
to_remove <- c("mar30_E_55", "mar30_H2O_57") #removed those samples since they didnt have any OTUs of interest
physeq_18S_filter2 <- prune_samples(!(sample_names(physeq_18S_filter1) %in% to_remove), physeq_18S_filter1)

ps1 <- merge_samples(physeq_18S_filter2, "month") #combine by month

#get relative abundance 
rel.abund <- transform_sample_counts(ps1, function(x) x/sum(x))
glom <- tax_glom(rel.abund, taxrank=rank_names(rel.abund)[6]) #collapse by genus level
data <- psmelt(glom)
data$V7 <- as.character(data$V7) #convert character
data$V7[data$Abundance < 0.01] <- "< 1% abund" # rename low freq phyla
#looking at overall median
medians <- plyr::ddply(data, ~V7, function(x) c(median=median(x$Abundance)))
medians
#get just useful columns
rel_18S <- data[,c(2,3,42)]


#sanity check
sum(data[data$Sample == 'april',]$Abundance)


pdf("pred_host_rel_abund_month_test4.pdf", width=30, height=18)
ggplot(data, aes(x=Sample, y=Abundance, fill=V7)) + geom_bar(aes(), stat="identity", position="stack") + 
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august")) + theme_minimal()
dev.off()


#getting legionella abudance data
seqtab_16S <- read.table("feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_16S)

#reading in taxonomy
tax_16S <- read.table("taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)

#reeading in metadata
map_16S <- read.table("metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)

#make 16S and 18S phyloseq object from qiime2
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)

ps16 <- merge_samples(physeq_16S, "month") #combine by month

#get relative abundance 
rel.abund.16S <- transform_sample_counts(ps16, function(x) x/sum(x))
glom_16S <- tax_glom(rel.abund.16S, taxrank=rank_names(rel.abund.16S)[6]) #collapse by genus level
data_16S <- psmelt(glom_16S)
data_16S$V7 <- as.character(data_16S$V7) #convert character

#get just legionella abundance
leg <- data_16S[data_16S$V7 == 'Legionella',]
leg_rel <- leg[, c(2,3,42)]
colnames(leg_rel)[2] <- "Abundance_leg"
colnames(leg_rel)[3] <- "Genus"


#combine 18S and legionella abudance data
combined_leg_18S <- left_join(rel_18S, leg_rel, by = join_by(Sample == Sample))

pdf("pred_host_leg.pdf", width=30, height=18)
ggplot(combined_leg_18S)  + 
  geom_bar(aes(x=Sample, y=Abundance, fill=V7),stat="identity", position="stack")+
  geom_line(aes(x=Sample, y=100*Abundance_leg, group=1),stat="identity",color="red",size=5)+
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()

import_am=c("Vermamoeba" = "gray", "Echinamoeba" = "gray", "Acanthamoeba" = "plum3", "Dactylopodida" = "gray", "Arcellinida_unknown" = "gray", "BIO10-D10" = "gray", "Stygamoebida" = "gray", "Euamoebida" = "gray", "BOLA868" = "gray", "Korotnevella" = "gray", "Centramoebida" = "gray", "Mycamoeba" = "gray", "Vannella" = "gray", "Euamoebida_unknown" = "gray", "Vannellida" = "gray", "Tubulinea_unknown" = "gray", "Tubulinea" = "gray", "Vermamoeba" = "gray", "uncultured" = "gray", "Vannellida_unknown" = "gray", "Amoebozoa_unknown" = "gray", "Cryptodifflugia" = "gray", "Amoebozoa" = "gray", "Vermistella" = "gray", "Protosteliopsis" = "gray", "Arcellinida" = "gray", "Arcella" = "gray", "Gymnophrys" = "gray", "Heteromita" = "gray", "Cercozoa_unknown" = "gray", "uncultured" = "gray", "Paracercomonas" = "gray", "Cercozoa" = "gray", "Glissomonadida_unknown" = "gray", "Vampyrellidae" = "gray", "Tracheleuglypha" = "gray", "Cercomonadidae" = "gray", "Eocercomonas" = "gray", "Kraken" = "gray", "Euglypha" = "gray", "Glissomonadida" = "gray", "Trinema" = "gray", "Thecofilosea_unknown" = "gray", "Chilodonella" = "gray", "Tetrahymena" = "gray", "Amphileptus" = "gray", "Cyrtolophosis" = "gray", "Leptopharynx" = "gray", "Hymenostomatia" = "gray", "Cyclidium" = "gray", "Colpodea_unknown" = "gray", "Hypotrichia_unknown" = "gray", "Vorticella" = "gray", "Protocyclidium" = "gray", "Peritrichia" = "gray", "Spirotrichea_unknown" = "gray", "Oligohymenophorea" = "gray", "Conthreep_unknown" = "gray", "Nassophorea" = "gray", "Colpodida" = "gray", "Telotrochidium" = "gray", "Nassophorea_unknown" = "gray", "Oligohymenophorea_unknown" = "gray", "Aspidisca" = "gray", "Haptoria_unknown" = "gray", "Ephelota" = "gray", "Cyrtolophosidida" = "gray", "Naegleria" = "gray", "Allovahlkampfia" = "gray", "Tetramitia_unknown" = "gray", "Vahlkampfia" = "gray", "Neovahlkampfia" = "gray")
leg_pred_host_diff=c("Legionella"="plum3","Acanthamoeba" = "#D53E4F", "Echinamoeba"="#FDAE61", "Korotnevella" ="#FEE08B", "Naegleria"="#E6F598", "Tetrahymena"="#ABDDA4", "Vannella" = "#66C2A5","Vermamoeba"="#3288BD", 
      "Arcellinida_unknown" = "gray52", "BIO10-D10" = "gray52", "Dactylopodida" = "gray52", "Stygamoebida" = "gray52", "Euamoebida" = "gray52", "BOLA868" = "gray52", "Centramoebida" = "gray52", "Mycamoeba" = "gray52", "Vannella" = "gray52", "Euamoebida_unknown" = "gray52", "Vannellida" = "gray52", "Tubulinea_unknown" = "gray52", "Tubulinea" = "gray52", "uncultured" = "gray52", "Vannellida_unknown" = "gray52", "Amoebozoa_unknown" = "gray52", "Cryptodifflugia" = "gray52", "Amoebozoa" = "gray52", "Vermistella" = "gray52", "Protosteliopsis" = "gray52", "Arcellinida" = "gray52", "Arcella" = "gray52", "Gymnophrys" = "gray67", "Heteromita" = "gray67", "Cercozoa_unknown" = "gray67", "uncultured" = "gray67", "Paracercomonas" = "gray67", "Cercozoa" = "gray67", "Glissomonadida_unknown" = "gray67", "Vampyrellidae" = "gray67", "Tracheleuglypha" = "gray67", "Cercomonadidae" = "gray67", "Eocercomonas" = "gray67", "Kraken" = "gray67", "Euglypha" = "gray67", "Glissomonadida" = "gray67", "Trinema" = "gray67", "Thecofilosea_unknown" = "gray67", "Chilodonella" = "gray52", "Amphileptus" = "gray52", "Cyrtolophosis" = "gray52", "Leptopharynx" = "gray52", "Hymenostomatia" = "gray52", "Cyclidium" = "gray52", "Colpodea_unknown" = "gray52", "Hypotrichia_unknown" = "gray52", "Vorticella" = "gray52", "Protocyclidium" = "gray52", "Peritrichia" = "gray52", "Spirotrichea_unknown" = "gray52", "Oligohymenophorea" = "gray52", "Conthreep_unknown" = "gray52", "Nassophorea" = "gray52", "Colpodida" = "gray52", "Telotrochidium" = "gray52", "Nassophorea_unknown" = "gray52", "Oligohymenophorea_unknown" = "gray52", "Aspidisca" = "gray52", "Haptoria_unknown" = "gray52", "Ephelota" = "gray52", "Cyrtolophosidida" = "gray52", "Allovahlkampfia" = "gray52", "Tetramitia_unknown" = "gray52", "Vahlkampfia" = "gray52", "Neovahlkampfia" = "gray52", "Cercomonadidae_unknown" = "gray67", "Euglyphida_unknown" = "gray67", "Hypotrichia_unknown" = "gray52", "Oligohymenphorea_unknown"="gray52", "Oligohymenophorea_unkown"="gray52", "Vampyrellidae_unknown"="gray67")
pdf("pred_host_leg_clean.pdf", width=30, height=18)
ggplot(combined_leg_18S)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V7),stat="identity", position="stack")+
  scale_fill_manual(values = c("Acanthamoeba" = "#D53E4F", "Dactylopodida" = "#F46D43", "Echinamoeba"="#FDAE61", "Korotnevella" ="#FEE08B", 
      "Naegleria"="#E6F598", "Tetrahymena"="#ABDDA4", "Vannella" = "#66C2A5","Vermamoeba"="#3288BD"))+ #"Vahlkampfia" = "#66C2A5" is a maybe
  geom_line(aes(x=Sample, y=100*Abundance_leg, group=1),stat="identity",color="red",size=5)+
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()

pdf("pred_host_leg_test2.pdf", width=30, height=18)
ggplot(combined_leg_18S)  + 
  geom_bar(aes(x=Sample, y=Abundance,fill=V7),stat="identity", position="stack")+
  scale_fill_manual(values = ifelse(combined_leg_18S$V7=="Acanthamoeba", "red", "gray"))+
  geom_line(aes(x=Sample, y=100*Abundance_leg, group=1),stat="identity",color="red",size=5)+
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+
  scale_y_continuous(sec.axis = sec_axis(~./1000)) + theme_minimal()
dev.off()

#confidence interval for mean legionella relative abudance 
rel.abund.16S <- transform_sample_counts(ps16, function(x) x/sum(x))
data_16S_nonglom <- psmelt(rel.abund.16S)
data_16S_nonglom$V7 <- as.character(data_16S_nonglom$V7) #convert character
leg_all <- data_16S_nonglom[data_16S_nonglom$V7 == 'Legionella',]
leg_all_rel <- leg_all[, c(2,3,42)]

dt <- leg_all_rel%>%
  dplyr::group_by(Sample)%>%
  dplyr::summarise(
    mean = mean(Abundance),
    lci = t.test(Abundance, conf.level = 0.95)$conf.int[1],
    uci = t.test(Abundance, conf.level = 0.95)$conf.int[2])

pdf("test1.pdf")
ggplot(data = dt) + geom_line(aes(x=Sample, y=mean, group=1),stat="identity",color="red") + 
  geom_point(aes(x=Sample, y=mean), color= "red") + 
  geom_errorbar(aes(x=Sample, ymin=lci, ymax= uci), width = 0.4, color ="red", size = 1) +
  geom_text(aes(x=Sample, y=lci, label = round(lci,1)), size= 2, vjust = 1) + 
  geom_text(aes(x=Sample, y=uci, label = round(uci,1)), size= 2, vjust = -1) +
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+ 
  theme_minimal()
dev.off()
