#remotes::install_github("phytomosaic/ecole")
##remotes::install_github("bcdudek/bcdstats")
#install.packages("ggpubr")
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(tidyverse)
#library(qiime2R)
#library(readr)
library(RColorBrewer)
library(Hmisc)
library(circlize) #https://mran.microsoft.com/snapshot/2018-01-05/web/packages/circlize/vignettes/visualize_relations_by_chord_diagram.pdf
#library(reshape2)
library(chorddiag)
library(vegan)
library(sinkr)
library(ecole)
library(ranacapa)
library(bcdstats)
library(ggpubr)
library(tabletools)
setwd("~/legionella/R/diversity")
#reading sequence tables
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

physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)

#make 18S and 18S phyloseq object from qiime2
physeq_18S = subset_taxa(physeq_18S, !V2=="Bacteria" & !V3=="Eukaryota_unknown")
clr_18S <- microbiome::transform(physeq_18S,'clr')

#beta diversity
pcoa_bc = ordinate(clr_18S, "PCoA", "euclidean") 
pdf("bd.18S.type.pdf")
plot_ordination(clr_18S, pcoa_bc, color = "type") + 
  geom_point(size = 3)
dev.off()
permanova_pairwise(otu_table(clr_18S), grp=sample_data(clr_18S)$type, method="euclidean") # check signficane
#test correlation between enviromental factors and betadiversity
types <- c("water", "bright", "dark")

env <- sample_data(clr_18S)
env.mat <- env[, c(8,9,10,11,12,13,14)]
all.sum.18S <- bioenv((otu_table(clr_18S)), env.mat, index = "euclidean")
bio_all.18S <- bioEnv(
  as.matrix(otu_table(clr_18S)),
  env.mat,
  fix.dist.method = "euclidean",
  var.dist.method = "euclidean",
  scale.fix = FALSE,
  scale.var = TRUE,
  output.best = 10,
  var.max = ncol(env.mat)
)

bioenv_18S <- vector("list", 10)
bioENV_18S <- vector("list", 10)

for(i in types) {
subset.18S <- subset_samples(clr_18S, type==i)
as.matrix(otu_table(subset.18S))
env <- sample_data(subset.18S)
env.mat <- env[, c(8,9,10,11,12,13,14)]
bioenv_18S[[i]] <- bioenv((otu_table(subset.18S)), env.mat, index = "euclidean")
bioENV_18S[[i]] <- bioEnv(
  as.matrix(otu_table(subset.18S)),
  env.mat,
  fix.dist.method = "euclidean",
  var.dist.method = "euclidean",
  scale.fix = FALSE,
  scale.var = TRUE,
  output.best = 10,
  var.max = ncol(env.mat)
)
}
bioenv_18S
bioENV_18S$water$best.model.vars
bioENV_18S$water$best.model.rho
bioENV_18S$bright$best.model.vars
bioENV_18S$bright$best.model.rho
bioENV_18S$dark$best.model.vars
bioENV_18S$dark$best.model.rho

#alpha diversity
physeq_18S2 <- prune_samples(sample_sums(physeq_18S) > 5000, physeq_18S) #remove less than 2000 reads
rare_18S <- rarefy_even_depth(physeq_18S2, rngseed=1, sample.size=0.99*min(sample_sums(physeq_18S2)), replace=F)
sample_sums(rare_18S)

permanova_pairwise(otu_table(rare_18S), grp=sample_data(rare_18S)$type, method="bray") #check beta diversity between samples

#plot differences in shannon diveristy across samples
pdf("./adiv.type.18S.pdf")
plot_richness(rare_18S, measures=c("Observed", "Shannon"), x="type") + 
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
alpha <- estimate_richness(rare_18S, split = TRUE, measures = NULL)
env <- sample_data(rare_18S)
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

pdf("shannon18S_all_spear_BH.pdf")
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off()


#alpha corr for each type
for(i in types) {
subset.18S <- subset_samples(rare_18S, type==i)
alpha <- estimate_richness(subset.18S, split = TRUE, measures = NULL)
env <- sample_data(subset.18S)
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

pdf(paste0("ab_enviro_spear18S_", i, ".pdf"))
corrplot::corrplot(cor_enviro_adjust$r, type="upper", order="original", 
         p.mat = cor_enviro_adjust$P, sig.level = 0.05, insig = "blank",
         tl.cex = 0.4, na.label= " ")
dev.off() 
}

#rarefaction curve
options(warn=-1) # suppress warnings
p <- ggrare(physeq_18S, step = 1000, color = "month", se = TRUE)
#p <- p + facet_wrap(~Tooth_Classification)
pdf("./rarefaction_plots.18S.pdf")
p + theme_minimal() + scale_x_continuous(labels = scales::comma)
dev.off()
options(warn=0) # back on

options(warn=-1) # suppress warnings
p <- ggrare(rare_18S, step = 1000, color = "month", se = TRUE)
#p <- p + facet_wrap(~Tooth_Classification)
pdf("./raredcurve_18S.pdf")
p + theme_minimal() + scale_x_continuous(labels = scales::comma)
dev.off()
options(warn=0) # back on