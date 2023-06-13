# Enviromental Correlation
## 1. Load Packages
```R
library(phyloseq)
library(tidyverse)
library(microbiome)
library(Hmisc)
library(tabletools)
library(corrplot)
library(RColorBrewer)
```
## 2. Make Phyloseq Object
```R
enviro_data <- read.table("metadata.txt", 
                   header=TRUE, row.names=1, sep = "\t")
#reading in 16S data into phyloseq object and filtering to just legionella
#frequency table
seqtab_16S <- read.table("feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxonomy
tax_16S <- read.table("taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#metadata
map_16S <- read.table("metadata_16S.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)
#Combining phyloseq objects
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
```
## 3. Transforming the data into relative abudance
```R
physeq_16S_rel <- microbiome::transform(physeq_16S, "compositional")
```
## 4. Getting just Legionella data
```R
physeq_16S_filter2 = subset_taxa(physeq_16S_rel, V7=="Legionella")
```
## 5. Getting genus names
```R
taxa = as(tax_table(physeq_16S_filter2), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V7)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
```
## 6. Getting Freuqency Matrix out of Phyloseq
```R
ASV_freq <- as(otu_table(physeq_16S_filter2), "matrix")
write.table(ASV_freq,"leg_ASV_freq.txt", sep = "\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
trans_ASV_freq <- t(ASV_freq)
trans_ASV_freq_df1 <- as.data.frame(trans_ASV_freq)
trans_ASV_freq_df2 <- rownames_to_column(trans_ASV_freq_df1, var = "ASV")
```
## 7. Combining all Legionella 
```R
#rename to genus level
renamed_ASV_freq1 <- left_join(trans_ASV_freq_df2, orderdf, by=c('ASV'='ASV'))
colnames(renamed_ASV_freq1)[49] <- "Genus"
renamed_ASV_freq2 <- renamed_ASV_freq1[, c(49,2:48)]
#combine all legionella
renamed_ASV_freq3 <- renamed_ASV_freq2 %>%
  group_by(Genus) %>%
  summarise(across(c(1:47), sum))
#reformat into dataframe
renamed_ASV_freq4 <- as.data.frame(renamed_ASV_freq3)
rownames(renamed_ASV_freq4) <- renamed_ASV_freq4[,1]
renamed_ASV_freq4[,1] <- NULL
renamed_ASV_freq5 <- t(as.matrix(renamed_ASV_freq4))
```
## 8. Formatting for running Correlation
### 8a. Getting all Legionella ASVs and renaming them
```R
ASV_freq_df <- as.data.frame(ASV_freq, row.names=NULL)
names(ASV_freq_df) <- c("ASV1",  "ASV2", "ASV3", "ASV4", "ASV5", "ASV6", "ASV7", "ASV8", "ASV9", "ASV10",  "ASV11",  "ASV12",  "ASV13",  "ASV14",  "ASV15",  "ASV16",  "ASV17",  "ASV18",  "ASV19",  "ASV20",  "ASV21",  "ASV22",  "ASV23",  "ASV24",  "ASV25",  "ASV26",  "ASV27",  "ASV28",  "ASV29",  "ASV30",  "ASV31",  "ASV32",  "ASV33",  "ASV34",  "ASV35",  "ASV36",  "ASV37",  "ASV38",  "ASV39",  "ASV40",  "ASV41",  "ASV42",  "ASV43",  "ASV44",  "ASV45",  "ASV46",  "ASV47",  "ASV48",  "ASV49",  "ASV50",  "ASV51",  "ASV52",  "ASV53",  "ASV54",  "ASV55",  "ASV56",  "ASV57",  "ASV58",  "ASV59",  "ASV60",  "ASV61",  "ASV62",  "ASV63",  "ASV64",  "ASV65",  "ASV66",  "ASV67",  "ASV68",  "ASV69",  "ASV70")
ASV_freq_df$sample <- row.names(ASV_freq_df)
```
### 8b. Getting Enviromental Data
```R
enviro_data2<- enviro_data[,c(2,3,7:13,15:19)]

```
### 8c. Merging Data
```R
#combining enviromental data and all Legionella data
merdged_data_enviro1 <- merge(enviro_data2, renamed_ASV_freq5,
                          by = 'row.names', all = TRUE)
colnames(merdged_data_enviro1)[1] <- "sample"
#combining enviro and all Leg data with Leg ASV data
merdged_data_enviro3 <- merge(merdged_data_enviro1, ASV_freq_df,
                          by = c('sample'='sample'), all = TRUE)
#seperating water and non water samples
merdged_data_enviro3$water <- +(merdged_data_enviro3$type == "water" & !is.na(merdged_data_enviro3$type)) 
#renaming Months
merdged_data_enviro3$months <- ifelse(merdged_data_enviro3$month == "march" , 3, 
    ifelse(merdged_data_enviro3$month == "april", 4 , 
    ifelse(merdged_data_enviro3$month == "may", 5, 
    ifelse(merdged_data_enviro3$month == "june", 6, 
    ifelse(merdged_data_enviro3$month == "july", 7, 8)))))
merdged_data_enviro4<- merdged_data_enviro3[,c(87,88,4:86)]
#converting to a matrix for correlation
enviro_data_mat <- as.matrix(merdged_data_enviro4)
```
## 9. Running Correlation
```R
#running correlation for spearman
cor_enviro_spear <- rcorr(enviro_data_mat,type = c("spearman"))
```
## 10. Adjusting P-value
```R
cor_enviro_adjust <- rcorr_padjust(cor_enviro_spear, method = "BH")
diag(cor_enviro_adjust$P) <- 0
```
## 11. Visualize
```R
pdf("enviro_spear_BH.pdf")
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off()
``` 









