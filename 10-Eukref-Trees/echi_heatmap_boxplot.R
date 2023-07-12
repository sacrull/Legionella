#from eukref
	#sed -i '1 i\OTU apr13_E_73' temp5 
	#sed 's/\t.*/\t 0/g' temp5 | sed 's/ /\t/' | sed 's/ //g'> ref_ASV.txt
	#cp ref_ASV.txt ~/legionella/R/phyloseq_tree/vannella

	#awk '{print $1, $7}' 18S-tax-renamed.txt > genus
	#cat genus temp5 > taxa_all_vannella
	#sort taxa_all_vannella | uniq > taxa_vannella
	#sed 's/\t/\t /g' temp5 > ref_test.txt
library(devtools)
library(phyloseq)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ape)
library(microbiome)
#conda install -c conda-forge r-ggpubr
library(ggpubr)
library(rstatix)

#reading sequence tables
seqtab_18S <- read.table("18S-ASV-renamed.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
#reading in taxonomy
tax_18S <- read.table("18S-tax-renamed.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
#reeading in metadata
map_18S <- read.table("metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
#physeq_16S_filter1 = subset_taxa(physeq_16S, V5=="Legionellales")
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
physeq_18S_clr_filter = subset_taxa(physeq_18S, V7=="Echinamoebida")

#get monthly CLR
ps1 <- merge_samples(physeq_18S, "month") #combine by month

#get relative abundance 
physeq_18S_clr <- microbiome::transform(ps1, "clr")
physeq_18S_clr_filter = subset_taxa(physeq_18S_clr, V7=="Echinamoebida")
#getting ASVs in datafrane
data_18S <- psmelt(physeq_18S_clr_filter)
data_18S_2 <- data_18S[, c(1,2,3)]


#getting refrence files
tre <- read.tree("echi.tre")
#reference otu table
seqtab_18S_1 <- read.table("18S-ASV-renamed.txt", header=T)
seqtab_ref <- read.table("ref_ASV.txt", header=T)
#combing the otu tables
otu <- dplyr::full_join(seqtab_18S_1, seqtab_ref)
row.names(otu) <- otu$OTU
otu <- otu[,-1]
otu <- data.matrix(otu)
otu[is.na(otu)] <- as.numeric(0)
otu_df <- as.data.frame(otu)
otu_df <- cbind(rownames(otu_df), data.frame(otu_df, row.names=NULL))
colnames(otu_df)[1] <- "OTU"
filter_ref <- otu_df %>% 
  filter(!grepl('ASV', OTU))
row.names(filter_ref) <- filter_ref$OTU
filter_ref <- filter_ref[,-1]
filter_ref <- data.matrix(filter_ref)
otu_ref <-otu_table(filter_ref, taxa_are_rows=T)
#ref taxa
tax_other <- read.table("annotations.echi.txt", header=F, row.names=1, sep="\t")
tax_other_phylo <- tax_table(as.matrix(tax_other))
#combine ref data
physeq_ref <- merge_phyloseq(otu_ref,map_map_18S,tax_other_phylo, tre)
physeq_ref_month <- merge_samples(physeq_ref, "month") #combine by month
#getting ref genomes in order
data_ref <- psmelt(physeq_ref_month)
data_ref_2 <- data_ref[, c(1,2,3)]
#replace 0 with NA
data_ref_2[data_ref_2 == 0] <- NA

#join ref and ASVs together
otu <- dplyr::full_join(data_18S_2, data_ref_2)


#getting tip order and month order
month_order=c("march","april","may","june", "july","august")
tip_order <- rev(tre$tip.label)
tip_order2 <- gsub("'","",tip_order)

#heatmap
pdf("echinamoebida_asv_heatmap_clr.pdf")
ggplot(otu, aes(x=factor(Sample, levels=month_order),y=factor(OTU, levels=tip_order2))) +
  geom_tile((aes(fill = Abundance)))+
  scale_fill_gradient(low = "#FFC300", high = "#900C3F") +
  coord_fixed() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 2.5, vjust = 0.6, hjust = 0.0,), axis.text.y = element_text(size = 3),
        panel.background = element_rect(fill = 'grey67'), axis.title.x = element_text(size= 7), 
        axis.title.y = element_text(size= 6), plot.title = element_text(size = 10), 
        legend.title = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text = element_text(size = 4)) 
dev.off()


#boxplot
physeq_18S_bp_clr <- microbiome::transform(physeq_18S, "clr")
glom_bp <- tax_glom(physeq_18S_bp_clr, taxrank=rank_names(physeq_18S_bp_clr)[6]) #collapse by genus level
physeq_18S_clr_filter_bp = subset_taxa(glom_bp, V7=="Echinamoebida")
data_18S_bp <- psmelt(physeq_18S_clr_filter_bp)
data_pb_2 <- data_18S_bp[, c(2,6,3)]


pdf("echi_asv_boxplot_clr.pdf")
ggplot(data_pb_2, aes(x=factor(month, levels=month_order),y=Abundance))+
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "none") +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size=1)+
  labs(x ="Month", y = "CLR Abundance")+
  theme_classic()
dev.off()

