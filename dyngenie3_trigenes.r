library(reshape)
library(magrittr)
library(ggplot2)
library(data.table)
library(scales)

setwd("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/")

##reading the data in-------->>

base <- read.csv("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/aracne_all_conditions/vst.txt", header = TRUE, sep = "\t", row.names = 1)[, 1 : 5]

regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]

raw_dat <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("aracne_all_conditions/vst_conditon_%s.txt", i), sep = '\t', row.names = 1)
})
degs <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("l2FC_tables/CN_data_condition%s_DEG.txt", i), header = FALSE)[,1]
})
tfs <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("aracne_all_conditions/transcription_factors_conditon_%s.txt", i), header = FALSE)[,1]
})
##reading the data in-------->>

tri_regs <- c("g6430", "g6432")
tri_genes <- c("g6429", "g6431", "g6426", "g6434", "g6435", "g6436", "g6437")

#--------------------------------------------------------------------
#separating conditions; creating expr dfs per rep (list of 4 conditions containing list of 5 time-series, 7 times points each)
TS_data <- lapply(1 : 4, FUN = function(j) {
  cond <- read.table(sprintf("/home/ngrinber/quorn_grn/L2FC_filtered_data/L2FC_filter_condition_%s.txt", j))
  out <- lapply(1 : 5, FUN = function(i) {
    ind <- grep(paste0("_", i), colnames(cond))
    x <-  cond[, ind]
    x[, 1] <- x[, 1] + 0.000001 #offset contant genes but nudging expr for the first time-point
    colnames(x) <- paste0("tp", 1 : ncol(x))
    as.matrix(x)
  })
  out
})
names(TS_data) <- paste("cond", 1 : 4)

#creating time-points
tvec <- 0 : 6
names(tvec) <- colnames(TS_data[[1]][[1]])
time_points <- list(tvec, tvec, tvec, tvec, tvec)

#indices of regulator genes
reg_ind <- lapply(1 : 4, FUN = function(i) which(rownames(TS_data[[i]][[1]]) %in% tfs[[i]]))





###