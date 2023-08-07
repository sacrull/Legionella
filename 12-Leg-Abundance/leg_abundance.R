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

#getting legionella abudance data
seqtab_16S <- read.table("../feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_16S)

#reading in taxonomy
tax_16S <- read.table("../taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)

#reeading in metadata
map_16S <- read.table("../../metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)

#make 16S and 18S phyloseq object from qiime2
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)


#boxplot
physeq_16S_bp_clr <- microbiome::transform(physeq_16S, "clr")
glom_bp <- tax_glom(physeq_16S_bp_clr, taxrank=rank_names(physeq_16S_bp_clr)[6]) #collapse by genus level
physeq_16S_clr_filter_bp = subset_taxa(glom_bp, V7=="Legionella")
data_16S_bp <- psmelt(physeq_16S_clr_filter_bp)
data_pb_2 <- data_16S_bp[, c(2,29,3,6)]
network_order=c("march_april","may_june","july_august")
month_colors=c("march"="darkorchid1","april"="#147BD1","may"="#2DC84D","june"="#F7EA48","july"="orange1","august"="#E03C31")

pdf("leg_asv_network_boxplot_clr.pdf")
ggplot(data_pb_2, aes(x=factor(network, levels=network_order),y=Abundance))+
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
  geom_boxplot() +
  geom_jitter(aes(color=month),shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = month_colors)+
  labs(x ="Month", y = "CLR Abundance")+
  theme_classic()
dev.off()

pdf("leg_asv_boxplot_clr.pdf")
ggplot(data_pb_2, aes(x=factor(month, levels=month_order),y=Abundance))+
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
  geom_boxplot() +
  geom_jitter(aes(color=month),shape=16, position=position_jitter(0.2), size=2.5)+
  scale_color_manual(values = month_colors)+
  labs(x ="Month", y = "CLR Abundance")+
  theme_classic()
dev.off()






physeq_16S_bp_rel <- microbiome::transform(physeq_16S, "compositional")
glom_bp <- tax_glom(physeq_16S_bp_rel, taxrank=rank_names(physeq_16S_bp_rel)[6]) #collapse by genus level
physeq_16S_clr_filter_bp = subset_taxa(glom_bp, V7=="Legionella")
data_16S_bp <- psmelt(physeq_16S_clr_filter_bp)
data_pb_2 <- data_16S_bp[, c(2,6,3)]
data_pb_3 <- subset(data_pb_2, Abundance > 0)

dt <- data_pb_3%>%
  dplyr::group_by(month)%>%
  dplyr::summarise(
    mean = mean(Abundance),
    lci = t.test(Abundance, conf.level = 0.95)$conf.int[1],
    uci = t.test(Abundance, conf.level = 0.95)$conf.int[2]) 


pdf("test1.pdf")
ggplot(data = dt) + geom_line(aes(x=month, y=mean, group=1),stat="identity",color="red") + 
  geom_point(aes(x=month, y=mean), color= "red") + 
  geom_errorbar(aes(x=month, ymin=lci, ymax= uci), width = 0.4, color ="red", size = 1) +
  geom_text(aes(x=month, y=lci, label = round(lci,1)), size= 2, vjust = 1) + 
  geom_text(aes(x=month, y=uci, label = round(uci,1)), size= 2, vjust = -1) +
  scale_x_discrete(limits = c("march", "april", "may", "june", "july","august"))+ 
  theme_minimal()
dev.off()


pdf("leg_asv_boxplot_rel.pdf")
ggplot(data_pb_3, aes(x=factor(month, levels=month_order),y=Abundance))+
  geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "none") +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), size=1)+
  labs(x ="Month", y = "Relative Abundance")+
  theme_classic()
dev.off()