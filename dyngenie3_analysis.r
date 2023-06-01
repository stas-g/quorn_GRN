library(reshape)
library(magrittr)
library(ggplot2)
library(data.table)
library(scales)

setwd("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/")

FILE <- "/home/ngrinber/quorn_grn/dynGENIE3_results"
##reading the data in-------->>

# base <- read.csv("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/aracne_all_conditions/vst.txt", header = TRUE, sep = "\t", row.names = 1)[, 1 : 5]

# regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]

# raw_dat <- lapply(1 : 4, FUN = function(i) {
#     read.csv(sprintf("aracne_all_conditions/vst_conditon_%s.txt", i), sep = '\t', row.names = 1)
# })
# degs <- lapply(1 : 4, FUN = function(i) {
#     read.csv(sprintf("l2FC_tables/CN_data_condition%s_DEG.txt", i), header = FALSE)[,1]
# })
# tfs <- lapply(1 : 4, FUN = function(i) {
#     read.csv(sprintf("aracne_all_conditions/transcription_factors_conditon_%s.txt", i), header = FALSE)[,1]
# })
##reading the data in-------->>

tri_regs <- c("g6430", "g6432")
tri_genes <- c("g6429", "g6431", "g6426", "g6434", "g6435", "g6436", "g6437")

TRI_GENES <- c("g6429" = "TRI4", "g6431" = "TRI5", "g6430" = "TRI6", "g6426" = "TRI8", "g6432" = "TRI10", "g6434" = "TRI11", "g6435" = "TRI12", "g6436" = "TRI13", "g6437" = "TRI14")


dyn_res0 <- lapply(1 : 4, FUN = function(i) {
    print(i)
    out <- readRDS(file.path(FILE, sprintf("cond%s_RF_sqrt.rds", i)))
    out$cond <- i
    out$mode <- "data"
    out
}) 

dyn_perm0 <- lapply(1 : 4, FUN = function(i) {
    print(i)
    out <- readRDS(file.path(FILE, sprintf("perm_cond%s_RF_sqrt.rds", i)))
    out$cond <- i
    out$mode <- "perm"
    out
}) 

dyn_res <- lapply(dyn_res0, FUN = function(x) x[x$weight > 0,]) %>% do.call(rbind, .)
dyn_perm <- lapply(dyn_perm0, FUN = function(x) x[x$weight > 0,]) %>% do.call(rbind, .)

dyn_all <- data.table(do.call(rbind, list(dyn_res, dyn_perm)))
dyn_all[, c("cond", "mode") := list(factor(dyn_all$cond), factor(dyn_all$mode))]

crit_val <- quantile(dyn_perm$weight, 0.95)
#what proportion of results pass the critical value
nrow(dyn_res[dyn_res$weight > crit_val,])/nrow(dyn_res)

# ggplot(dyn_all, aes(x = mode, y = weight)) + geom_boxplot() + facet_wrap(.~cond)
ggplot(dyn_all, aes(x = weight, fill = mode)) + geom_histogram(aes(y =..density..), alpha = 0.5) + facet_wrap(.~cond)

#how many of tri__regs in TFs and in the list of differentially expressed genes?
sapply(tfs, FUN = function(z) sum(c(tri_regs, tri_genes) %in% z))
sapply(degs, FUN = function(z) sum(c(tri_regs, tri_genes) %in% z))

#how many of the tri-gene connections are in the output? (per condition)
ids <- expand.grid(tri_regs, tri_genes)
ids <- paste(ids[,1], ids[,2], sep = "_")

#how many tri-gene connection pass the empirical cut-off?
dyn_res[dyn_res$weight > crit_val & dyn_res$id %in% ids,]
dyn_perm[dyn_perm$weight > crit_val & dyn_perm$id %in% ids,]

summary(dyn_res$weight)
summary(dyn_perm$weight)






##