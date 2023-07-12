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


#getting 18S ASVs
seqtab_18S <- read.table("18S-ASV-renamed.txt", header=T)
#getting refrence files
seqtab_ref <- read.table("ref_ASV.txt", header=T)
#combing the otu tables
otu <- dplyr::full_join(seqtab_18S, seqtab_ref)
row.names(otu) <- otu$OTU
otu <- data.matrix(otu)
otu[is.na(otu)] <- as.numeric(0)
otutable <- otu_table(otu, taxa_are_rows=T)
#get all the taxa
tax_other <- read.table("annotations.koro.txt", header=F, sep="\t")
my_taxa <- read.table("18S-tax-renamed.txt", header=F, sep="\t")
taxa <- dplyr::full_join(my_taxa, tax_other, by=c('V7'='V2', 'V1'='V1'))
taxa[is.na(taxa)] <- as.character("bleep")
row.names(taxa) <- taxa$V1
taxa2 <- taxa[, c(2,7)]
taxa_phylo <- tax_table(as.matrix(taxa2))

#taxa_names(tax_18S_phylo)
#reeading in metadata
map_18S <- read.table("metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
#read in tree
tre <- read.tree("korotnevella.tre")
#combine object
physeq <- merge_phyloseq(otutable,map_map_18S, tre,taxa_phylo)
#make tree
ps18 <- merge_samples(physeq, "month") #combine by month
rel.abund.18S <- transform_sample_counts(ps18, function(x) x/sum(x))
data_18S <- psmelt(rel.abund.18S)
data_18S_2 <- data_18S[, c(1,2,3)]
data_18S_3 <- as.matrix(data_18S_2)

month_order=c("march","april","may","june", "july","august")
tip_order <- rev(tre$tip.label)
tip_order2 <- gsub("'","",tip_order)

pdf("korotnevella_asv_heatmap.pdf")
ggplot(data_18S_2, aes(x=factor(Sample, levels=month_order),y=factor(OTU, levels=tip_order2))) +
  geom_tile((aes(fill = log(Abundance))))+
  scale_fill_gradient(low = "#FFC300", high = "#900C3F") +
  coord_fixed() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 2.5, vjust = 0.6, hjust = 0.0,), axis.text.y = element_text(size = 3),
        panel.background = element_rect(fill = 'grey67'), axis.title.x = element_text(size= 7), 
        axis.title.y = element_text(size= 6), plot.title = element_text(size = 10), 
        legend.title = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text = element_text(size = 4)) 
dev.off()



