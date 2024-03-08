#install.packages('randomForest')
#install.packages("caTools")
library(phyloseq)
library(randomForest)
library(caret)
library(e1071)
library(reshape2)
library(tidyverse)
library(dplyr)
library(caTools)
library(ranger)
library(fastshap)
library(kernelshap)
library(iml)
library(vip)
library(shapviz)
#set seed and load data
set.seed(12349)
setwd("/home/suzanne/legionella/R/pearson")


# format data
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
physeq_16S <- subset_taxa(physeq_16S, !V3=="Bacteria_unknown" & !V2=="Eukaryota")
physeq_16S1 <- prune_samples(sample_sums(physeq_16S) > 10000, physeq_16S) #remove less than 2000 reads
rare_16S <- rarefy_even_depth(physeq_16S1, rngseed=1, sample.size=0.99*min(sample_sums(physeq_16S1)), replace=F)
#rare_16S <- tax_glom(rare_16S, taxrank=rank_names(rare_16S)[6])
physeq_16S_filter2 <- subset_taxa(rare_16S, V7=="Legionella")

physeq_18S <- merge_phyloseq(otu_18S, map_map_18S, tax_18S_phylo)
physeq_18S <- subset_taxa(physeq_18S, !V2=="Bacteria" & !V3=="Eukaryota_unknown" & !V3=="uncultured")
physeq_18S2 <- prune_samples(sample_sums(physeq_18S) > 5000, physeq_18S) #remove less than 2000 reads
rare_18S <- rarefy_even_depth(physeq_18S2, rngseed=1, sample.size=0.99*min(sample_sums(physeq_18S2)), replace=F)
#known hosts and predator of legionella
physeq_18S_filter1 = subset_taxa(rare_18S, V3=="Amoebozoa" | V3=="Heterolobosea" | V3=="Ciliophora" |V3=="Cercozoa")

#merge phyloseq objects
physeq_merge1 <- merge_phyloseq(physeq_16S_filter2,physeq_18S_filter1)
glom <- tax_glom(physeq_merge1, taxrank=rank_names(physeq_merge1)[6]) #collapse by genus level
#rename to genus levels
#getting taxa genus
taxa = as(tax_table(physeq_merge1), "matrix")
taxadf = as.data.frame(taxa)
orderdf = select(taxadf, V7)
orderdf <- orderdf %>% 
  rownames_to_column(var = "ASV")
# get dataframe abudance
asv_tab_var <- as.data.frame(t(otu_table(glom)))
#rename to genus
asv_tab_var2 <- rownames_to_column(asv_tab_var, var = "ASV")
asv_tab_var3 <- left_join(asv_tab_var2, orderdf, by=c('ASV'='ASV'))
colnames(asv_tab_var3)[48] <- "genus"
rownames(asv_tab_var3) <- asv_tab_var3[,48]
asv_tab_var3 <- t(asv_tab_var3[,c(2:47)])
# remove nas
asv_tab_var <- asv_tab_var3[, colSums(is.na(asv_tab_var3)) == 0]
#get metadata
#claffiy into low and high leg
glom2 <- tax_glom(rare_16S, taxrank=rank_names(rare_16S)[6])
comp_16S <- microbiome::transform(glom2,'compositional')
leg <- as.data.frame(otu_table(subset_taxa(comp_16S, V7=="Legionella")))
leg <- rbind(names(leg), leg)
names(leg) <- names(leg)
leg$leg_val <- ifelse(leg$fa19522f4a04124c8aeb8b8247645879 > 0.003 , "high_leg", "low_leg")
leg <- rownames_to_column(leg, var = "sample")

meta <- microbiome::meta(glom)

meta$water <- +(meta$type == "water" & !is.na(meta$type))
meta$bright <- +(meta$type == "bright" & !is.na(meta$type))
meta$dark <- +(meta$type == "dark" & !is.na(meta$type))    
meta$months <- ifelse(meta$month == "march" , 3, 
    ifelse(meta$month == "april", 4 , 
    ifelse(meta$month == "may", 5, 
    ifelse(meta$month == "june", 6, 
    ifelse(meta$month == "july", 7, 8)))))
meta[, c(7:13)]
#meta2 <- meta %>% select(water, dark, bright, months, temp, dissolved_oxygen, conduct., pH, free_Cl, Br, turbidity) 
meta2 <- meta %>% select(months, temp, dissolved_oxygen, conduct., pH, free_Cl, Br, turbidity) 

# by high legionella location
rownames(leg) <- leg[,1]
#leg$leg_val <- factor(leg$leg_val)
asv_tab_var2 <- rownames_to_column(as.data.frame(asv_tab_var), var = "samples")
asv_tab_var3 <- left_join(asv_tab_var2, leg, by=c('samples'='sample'))
meta3 <- rownames_to_column(as.data.frame(meta2), var = "samples")
asv_tab_var3 <- left_join(asv_tab_var3, meta3, by=c('samples'='samples'))
rownames(asv_tab_var3) <- asv_tab_var3[,1]
asv_tab_var <- asv_tab_var3[, c(3:64,67:74,65)]
colnames(asv_tab_var)[71] <- "var"
asv_tab_var$var <- as.numeric(asv_tab_var$var)
asv_tab_var$var <- asv_tab_var$var *100


#make random forrest model
fit <- ranger(var ~ ., data = asv_tab_var,
	scale.permutation.importance = TRUE, importance = 'permutation')
#fit$confusion.matrix
#variable importance
pdf("./rf.leg_levels.importance.pdf")
vip(fit, title = "Variable Importance")
dev.off()
ranger::importance(fit)

#goodness of fit
p.ra <- predict(fit, data=asv_tab_var)
str(p.ra)
par(mfrow=c(1,2))
pdf("goodfit.pdf")
plot(asv_tab_var$var ~ p.ra$predictions, asp=1, pch=20, xlab="fitted", ylab="actual", xlim=c(0,3.3),
          ylim=c(0,3.7),main="log10(Zn), Meuse topsoils, Random Forest")
grid(); abline(0,1)
dev.off()

#fixed
x <-  asv_tab_var[,1:(ncol(asv_tab_var)-1)]
predictor <- Predictor$new(model =fit, data = x, y = asv_tab_var$var)
shapley <- iml::Shapley$new(predictor, x.interest = x[1, ])
pdf("sample1.pdf")
shapley$plot()
dev.off()


#shapley values for all samples
s <- kernelshap(fit, asv_tab_var[,1:(ncol(asv_tab_var)-1)], bg_X = asv_tab_var)
# Turn them into a shapviz object
sv <- shapviz(s)
# Step 3: Gain insights...
pdf("shaple_values_all.pdf", height = 20, width =20)
sv_importance(sv, kind = "bee")
dev.off()

pdf("shap_overall_import.pdf")
sv_importance(sv, show_numbers = TRUE)
dev.off()

import_feature <- c("Tracheleuglypha", "Vampyrellidae", "Mycamoeba", "Colpodea_unknown", "Cyclidium",
	"free_Cl", "Stygamoebida", "Arcellinida_unknown", "Cercozoa_unknown", "Paracercomonas")
pdf("feature_values.pdf", height =30 , width =30)
sv_dependence(sv, v = import_feature)
dev.off()
pdf("sample1_decomp.pdf")
sv_waterfall(sv, row_id = 1)
dev.off()

save.image("random_forest.RData")






























#shapley
library(tidyverse)
library(funModeling)
library(xgboost)
library(caret)
library(shapr)
#creat functions
# return matrix of shap score and mean ranked score list
shap.score.rank <- function(xgb_model = xgb_mod, shap_approx = TRUE, 
                            X_train = mydata$train_mm){
  require(xgboost)
  require(data.table)
  shap_contrib <- predict(xgb_model, X_train,
                          predcontrib = TRUE, approxcontrib = shap_approx)
  shap_contrib <- as.data.table(shap_contrib)
  shap_contrib[,BIAS:=NULL]
  cat('make SHAP score by decreasing order\n\n')
  mean_shap_score <- colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = T)]
  return(list(shap_score = shap_contrib,
              mean_shap_score = (mean_shap_score)))
}

# a function to standardize feature values into same range
std1 <- function(x){
  return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}


# prep shap data
shap.prep <- function(shap  = shap_result, X_train = mydata$train_mm, top_n){
  require(ggforce)
  # descending order
  if (missing(top_n)) top_n <- dim(X_train)[2] # by default, use all features
  if (!top_n%in%c(1:dim(X_train)[2])) stop('supply correct top_n')
  require(data.table)
  shap_score_sub <- as.data.table(shap$shap_score)
  shap_score_sub <- shap_score_sub[, names(shap$mean_shap_score)[1:top_n], with = F]
  shap_score_long <- melt.data.table(shap_score_sub, measure.vars = colnames(shap_score_sub))
  
  # feature values: the values in the original dataset
  fv_sub <- as.data.table(X_train)[, names(shap$mean_shap_score)[1:top_n], with = F]
  # standardize feature values
  fv_sub_long <- melt.data.table(fv_sub, measure.vars = colnames(fv_sub))
  fv_sub_long[, stdfvalue := std1(value), by = "variable"]
  # SHAP value: value
  # raw feature value: rfvalue; 
  # standarized: stdfvalue
  names(fv_sub_long) <- c("variable", "rfvalue", "stdfvalue" )
  shap_long2 <- cbind(shap_score_long, fv_sub_long[,c('rfvalue','stdfvalue')])
  shap_long2[, mean_value := mean(abs(value)), by = variable]
  setkey(shap_long2, variable)
  return(shap_long2) 
}

plot.shap.summary <- function(data_long){
  x_bound <- max(abs(data_long$value))
  require('ggforce') # for `geom_sina`
  plot1 <- ggplot(data = data_long)+
    coord_flip() + 
    # sina plot: 
    geom_sina(aes(x = variable, y = value, color = stdfvalue)) +
    # print the mean absolute value: 
    geom_text(data = unique(data_long[, c("variable", "mean_value"), with = F]),
              aes(x = variable, y=-Inf, label = sprintf("%.3f", mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2, 
              fontface = "bold") + # bold
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) + 
    scale_color_gradient(low="#FFCC33", high="#6600CC", 
                         breaks=c(0,1), labels=c("Low","High")) +
    theme_bw() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), # remove axis line
          legend.position="bottom") + 
    geom_hline(yintercept = 0) + # the vertical line
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    # reverse the order of features
    scale_x_discrete(limits = rev(levels(data_long$variable)) 
    ) + 
    labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") 
  return(plot1)
}


var_importance <- function(shap_result, top_n=10)
{
  var_importance=tibble(var=names(shap_result$mean_shap_score), importance=shap_result$mean_shap_score)
  
  var_importance=var_importance[1:top_n,]
  
  ggplot(var_importance, aes(x=reorder(var,importance), y=importance)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_light() + 
    theme(axis.title.y=element_blank()) 
}

#start analysis
asv_tab_var <- asv_tab_var3[, c(3:65,68:78,67)]
colnames(asv_tab_var)[75] <- "var"
asv_tab_var_2 <- asv_tab_var[, c(67:74)]
leg_shap = dummyVars(" ~ .", data = asv_tab_var_2, fullRank=T)
leg_x = predict(leg_shap, newdata = asv_tab_var_2)
target_var=ifelse(as.character(asv_tab_var$var)=="low_leg", 1,0)


## Create the xgboost model
model_leg = xgboost(data = leg_x, 
                   nround = 10, 
                   objective="reg:linear",
                   label= target_var)

pred_test = predict(model_leg, leg_x)
pred_test[(pred_test>0.5)] = 1
pred_test[(pred_test<0.5)] = 0

confusionMatrix(factor(pred),factor(test_set8$Employee_Turnover))


#calculate shap values
shap_result_leg = shap.score.rank(xgb_model = model_leg, 
                              X_train =leg_x,
                              shap_approx = F
                              ) 
#plot top 10 variables
pdf("var_predications.pdf")
var_importance(shap_result_leg, top_n=10)
dev.off()

## Prepare data for top N variables
shap_long_leg = shap.prep(shap = shap_result_leg,
                           X_train = leg_x , 
                           top_n = 10
                           )
#shap summary
pdf("shap_summary.pdf")
plot.shap.summary(data_long = shap_long_leg)
dev.off()

pdf("shap_plots.pdf")
xgb.plot.shap(data = leg_x, # input data
              model = model_leg, # xgboost model
              features = names(shap_result_leg$mean_shap_score[1:10]), # only top 10 var
              n_col = 3, # layout option
              plot_loess = T # add red line to plot
              )
dev.off()
