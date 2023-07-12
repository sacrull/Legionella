library(phyloseq)
library(tidyverse)
library(microbiome)
library(Hmisc)
library(tabletools)
library(corrplot)
library(RColorBrewer)

#reading in envriomental metadata
enviro_data <- read.table("metadata_unsure_R.txt", 
                   header=TRUE, row.names=1, sep = "\t")

seqtab_18S <- read.table("../18S-ASV-renamed.txt", header=T, row.names=1)
seqtab_18S_trans <- t(seqtab_18S)
otu_18S <-otu_table(seqtab_18S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)
#reading in taxonomy
tax_18S <- read.table("../18S-tax-renamed.txt", header=F, row.names=1, sep="\t")
tax_18S_phylo <- tax_table(as.matrix(tax_18S))
#taxa_names(tax_18S_phylo)
#reeading in metadata
map_18S <- read.table("../metadata_unsure_18S_R.txt", sep="\t", header=T, row.names=1)
map_map_18S <- sample_data(map_18S)
#make 16S and 18S phyloseq object from qiime2
#physeq_16S_filter1 = subset_taxa(physeq_16S, V5=="Legionellales")
physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
physeq_18S_rel <- microbiome::transform(physeq_18S, "compositional")
#known hosts and predator of legionella
ameobas <- c("Acanthamoeba", "Echinamoeba", "Korotnevella", "Naegleria", "Tetrahymena","Vannella","Vermamoeba")
for(i in ameobas) {
physeq_18S_filter1 = subset_taxa(physeq_18S_rel, V7 == i)
#merge phyloseq objects
#combinging all legionella 
ASV_freq <- as(otu_table(physeq_18S_filter1), "matrix")
trans_ASV_freq <- t(ASV_freq)
trans_ASV_freq_df1 <- as.data.frame(trans_ASV_freq)
trans_ASV_freq_df2 <- rownames_to_column(trans_ASV_freq_df1, var = "ASV")
#getting taxa genus
taxa = as(tax_table(physeq_18S_filter1), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V7)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
#renmaing all to genus level
renamed_ASV_freq1 <- left_join(trans_ASV_freq_df2, orderdf, by=c('ASV'='ASV'))
colnames(renamed_ASV_freq1)[49] <- "Genus"
renamed_ASV_freq2 <- renamed_ASV_freq1[, c(49,2:48)]
renamed_ASV_freq3 <- renamed_ASV_freq2 %>%
  group_by(Genus) %>%
  summarise(across(c(1:47), sum))
renamed_ASV_freq4 <- as.data.frame(renamed_ASV_freq3)
rownames(renamed_ASV_freq4) <- renamed_ASV_freq4[,1]
renamed_ASV_freq4[,1] <- NULL
renamed_ASV_freq5 <- t(as.matrix(renamed_ASV_freq4))
#running correlation
ASV_freq_df <- as.data.frame(ASV_freq, row.names=NULL)
ASV_freq_df$sample <- row.names(ASV_freq_df)
enviro_data2<- enviro_data[,c(2,3,7:13,15:19)]

merdged_data_enviro1 <- merge(enviro_data2, renamed_ASV_freq5,
                          by = 'row.names', all = TRUE)
colnames(merdged_data_enviro1)[1] <- "sample"
merdged_data_enviro3 <- merge(merdged_data_enviro1, ASV_freq_df,
                          by = c('sample'='sample'), all = TRUE)

merdged_data_enviro3$water <- +(merdged_data_enviro3$type == "water" & !is.na(merdged_data_enviro3$type)) 
merdged_data_enviro3$months <- ifelse(merdged_data_enviro3$month == "march" , 3, 
    ifelse(merdged_data_enviro3$month == "april", 4 , 
    ifelse(merdged_data_enviro3$month == "may", 5, 
    ifelse(merdged_data_enviro3$month == "june", 6, 
    ifelse(merdged_data_enviro3$month == "july", 7, 8)))))
merdged_data_enviro4<- merdged_data_enviro3 %>%
  select(water, months, 4:ncol(.)) 

enviro_data_mat <- as.matrix(merdged_data_enviro4)
#running correlation for spearman
cor_enviro_spear <- rcorr(enviro_data_mat,type = c("spearman"))
#adjust p value
cor_enviro_adjust <- rcorr_padjust(cor_enviro_spear, method = "BH")
diag(cor_enviro_adjust$P) <- 0

pdf(paste0("enviro_spear_", i, ".pdf"))
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off()
}