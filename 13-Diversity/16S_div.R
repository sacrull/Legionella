#remotes::install_github("phytomosaic/ecole")
##remotes::install_github("bcdudek/bcdstats")
#install.packages("ggpubr")
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(tidyverse)
#library(qiime2R)
#library(readr)
library(RColorBrewer)
library(Hmisc)
library(circlize) #https://mran.microsoft.com/snapshot/2016-01-05/web/packages/circlize/vignettes/visualize_relations_by_chord_diagram.pdf
#library(reshape2)
library(chorddiag)
library(vegan)
library(sinkr)
library(ecole)
library(ggpubr)

#reading sequence tables
seqtab_16S <- read.table("../feature-table-16S.txt", header=T, row.names=1)
seqtab_16S_trans <- t(seqtab_16S)
otu_16S <-otu_table(seqtab_16S_trans, taxa_are_rows=F)
#taxa_names(otu_18S)

#reading in taxonomy
tax_16S <- read.table("../taxonomy_16S.txt", header=F, row.names=1, sep="\t")
tax_16S_phylo <- tax_table(as.matrix(tax_16S))
#taxa_names(tax_16S_phylo)

#reeading in metadata
map_16S <- read.table("../metadata_unsure_16S_R.txt", sep="\t", header=T, row.names=1)
map_map_16S <- sample_data(map_16S)

#make 16S and 18S phyloseq object from qiime2
physeq_16S <- merge_phyloseq(otu_16S, map_map_16S, tax_16S_phylo)
physeq_16S <- subset_taxa(physeq_16S, !V3=="Bacteria_unknown", !V2=="Eukaryota")
clr_16S <- microbiome::transform(physeq_16S,'clr')

#beta diversity
pcoa_bc = ordinate(clr_16S, "PCoA", "euclidean") 
pdf("bd.16S.type.pdf")
plot_ordination(clr_16S, pcoa_bc, color = "type") + 
  geom_point(size = 3)
dev.off()
permanova_pairwise(otu_table(clr_16S), grp=sample_data(clr_16S)$type, method="euclidean") # check signficane
#test correlation between enviromental factors and betadiversity
types <- c("water", "bright", "dark")

env <- sample_data(clr_16S)
env.mat <- env[, c(8,9,10,11,12,13,14)]
all.sum.16S <- bioenv((otu_table(clr_16S)), env.mat, index = "euclidean")
bio_all.16S <- bioEnv(
  as.matrix(otu_table(clr_16S)),
  env.mat,
  fix.dist.method = "euclidean",
  var.dist.method = "euclidean",
  scale.fix = FALSE,
  scale.var = TRUE,
  output.best = 10,
  var.max = ncol(env.mat)
)

bioenv_16S <- vector("list", 10)
bioENV_16S <- vector("list", 10)

for(i in types) {
subset.16S <- subset_samples(clr_16S, type==i)
as.matrix(otu_table(subset.16S))
env <- sample_data(subset.16S)
env.mat <- env[, c(8,9,10,11,12,13,14)]
bioenv_16S[[i]] <- bioenv((otu_table(subset.16S)), env.mat, index = "euclidean")
bioENV_16S[[i]] <- bioEnv(
  as.matrix(otu_table(subset.16S)),
  env.mat,
  fix.dist.method = "euclidean",
  var.dist.method = "euclidean",
  scale.fix = FALSE,
  scale.var = TRUE,
  output.best = 10,
  var.max = ncol(env.mat)
)
}
bioenv_16S
bioENV_16S$water$best.model.vars
bioENV_16S$water$best.model.rho
bioENV_16S$bright$best.model.vars
bioENV_16S$bright$best.model.rho
bioENV_16S$dark$best.model.vars
bioENV_16S$dark$best.model.rho

#alpha diversity
to_remove <- c("jul20_W_176") #removed those samples since they didnt have any OTUs of interest
physeq_16S_filter1 <- prune_samples(!(sample_names(physeq_16S) %in% to_remove), physeq_16S)
rare_16S <- rarefy_even_depth(otu_table(physeq_16S_filter1), rngseed = TRUE, replace = FALSE)
data_otu_filt_rar = data.frame(otu_table(rare_16S)) # create a separated file
ps.rar.16S <- phyloseq(rare_16S, tax_16S_phylo, map_map_16S) # create a phyloseq object

permanova_pairwise(otu_table(ps.rar.16S), grp=sample_data(ps.rar.16S)$type, method="bray") #check beta diversity between samples

#plot differences in shannon diveristy across samples
pdf("./adiv.type.16S.pdf")
plot_richness(ps.rar.16S, measures=c("Observed", "Shannon"), x="type") + 
    theme_minimal() + 
    geom_jitter(alpha=0.25) +
    geom_pwc(label = "{p.format}{p.signif}", hide.ns =TRUE, p.adjust.method = "fdr") +
    geom_point() +
    stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 5,
    shape = 23,
    fill = "red"
  )
dev.off()
#get correlation between shannon diveristy for alls samples and enviroment
alpha <- estimate_richness(ps.rar.16S, split = TRUE, measures = NULL)
env <- sample_data(ps.rar.16S)
env.mat <- env[, c(8,9,10,11,12,13,14)]
env.mat$names <- rownames(env.mat)
alpha$names <- rownames(alpha)

names_df <- as.data.frame(alpha$names)
names_df$no <- rownames(names_df)
shannon <- as.data.frame(alpha$Shannon)
shannon$no <- rownames(shannon)
shannon_index <- inner_join(names_df, shannon, by=c('no'='no'))
shannon_index <- shannon_index[,c(1,3)]
colnames(shannon_index)[1] ="names"
colnames(shannon_index)[2] ="shannon"

meta_alpha <- inner_join(shannon_index, env.mat, by=c('names'='names'))
rownames(meta_alpha) <- meta_alpha[,1]
meta_alpha2 <- meta_alpha[,c(2:9)]
meta_alpha3 <-as.matrix(meta_alpha2)

#corr
cor_enviro_spear <- rcorr(meta_alpha3,type = c("spearman"))
cor_enviro_adjust <- rcorr_padjust(cor_enviro_spear, method = "BH")
diag(cor_enviro_adjust$P) <- 0

pdf("shannon16S_all_spear_BH.pdf")
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off()


#alpha corr for each type
for(i in types) {
subset.16S <- subset_samples(ps.rar.16S, type==i)
alpha <- estimate_richness(subset.16S, split = TRUE, measures = NULL)
env <- sample_data(subset.16S)
env.mat <- env[, c(8,9,10,11,12,13,14)]
env.mat$names <- rownames(env.mat)
alpha$names <- rownames(alpha)

names_df <- as.data.frame(alpha$names)
names_df$no <- rownames(names_df)
shannon <- as.data.frame(alpha$Shannon)
shannon$no <- rownames(shannon)
shannon_index <- inner_join(names_df, shannon, by=c('no'='no'))
shannon_index <- shannon_index[,c(1,3)]
colnames(shannon_index)[1] ="names"
colnames(shannon_index)[2] ="shannon"

meta_alpha <- inner_join(shannon_index, env.mat, by=c('names'='names'))
rownames(meta_alpha) <- meta_alpha[,1]
meta_alpha2 <- meta_alpha[,c(2:9)]
meta_alpha3 <-as.matrix(meta_alpha2)

#corr
cor_enviro_spear <- rcorr(meta_alpha3,type = c("spearman"))
cor_enviro_adjust <- rcorr_padjust(cor_enviro_spear, method = "BH")
diag(cor_enviro_adjust$P) <- 0

pdf(paste0("ab_enviro_spear16S_", i, ".pdf"))
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off() 
}
