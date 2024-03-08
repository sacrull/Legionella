alpha <- estimate_richness(physeq_18S_filter1, split = TRUE, measures = NULL)
env <- sample_data(physeq_18S_filter1)
env.mat <- env[, c(8,9,10,11,12,13,14)]
env.mat$names <- rownames(env.mat)
alpha$names <- rownames(alpha)

meta_alpha <- inner_join(env.mat, alpha, by=c('names'='names'))

cor.test(meta_alpha$Shannon, meta_alpha$temp)

alpha <- estimate_richness(physeq_16S, split = TRUE, measures = NULL)
env <- sample_data(physeq_16S)
env.mat <- env[, c(8,9,10,11,12,13,14)]
env.mat$names <- rownames(env.mat)
alpha$names <- rownames(alpha)

meta_alpha <- inner_join(env.mat, alpha, by=c('names'='names'))

cor.test(meta_alpha$Shannon, meta_alpha$turbidity)




#16S rarefaction
to_remove <- c("jul20_W_176") #removed those samples since they didnt have any OTUs of interest
physeq_16S_filter1 <- prune_samples(!(sample_names(physeq_16S) %in% to_remove), physeq_16S)
options(warn=-1) # suppress warnings
p <- ggrare(physeq_16S_filter1, step = 1000, color = "type", se = TRUE)
#p <- p + facet_wrap(~aliquot_type)
pdf("./rarefaction_plots_16S.pdf")
p + theme_minimal()
dev.off()
options(warn=0) # back on

#18S rarefaction
to_remove <- c("jul20_W_176") #removed those samples since they didnt have any OTUs of interest
physeq_18S_filter1 <- prune_samples(!(sample_names(physeq_16S) %in% to_remove), physeq_16S)
options(warn=-1) # suppress warnings
p <- ggrare(physeq_18S_filter1, step = 1000, color = "type", se = TRUE)
#p <- p + facet_wrap(~aliquot_type)
pdf("./rarefaction_plots_18S.pdf")
p + theme_minimal()
dev.off()
options(warn=0)