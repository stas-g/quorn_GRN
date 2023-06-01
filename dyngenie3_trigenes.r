library(reshape)
library(magrittr)
library(ggplot2)
library(data.table)
library(scales)

source("/home/ngrinber/projects/niab/gene_regulatory_network/install/dynGENIE3.R")
setwd("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/")

FILE <- "/home/ngrinber/quorn_grn/dynGENIE3_results"
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

TRI_GENES <- c("g6429" = "TRI4", "g6431" = "TRI5", "g6430" = "TRI6", "g6426" = "TRI8", "g6432" = "TRI10", "g6434" = "TRI11", "g6435" = "TRI12", "g6436" = "TRI13", "g6437" = "TRI14")

#--------------------------------------------------------------------
#separating conditions; creating expr dfs per rep (list of 4 conditions containing list of 5 (rep) time-series, 7 times points each)
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

#run condition cond
setwd("/home/ngrinber/projects/niab/gene_regulatory_network/install/")
# for(cond in 1 : 4) {
#   mod <- dynGENIE3(TS_data[[cond]], time_points, regulators = tfs[[cond]], tree.method = "RF", K = "sqrt", verbose = TRUE)
#   res <- get.link.list(mod$weight.matrix)
#   res$id <- paste(res$regulatory.gene, res$target.gene, sep = "_")
#   res <- res[order(res$id),]

#   saveRDS(res, file = file.path(FILE, sprintf("cond%s_RF_sqrt.rds", cond)))

#   message(sprintf("dynGENIE3 on condition %s: DONE!\n", cond))
# }

##permutation tests-------------->>>>>
TS_rand <- lapply(TS_data, FUN = function(x){
    lapply(x, FUN = function(z) {
      rn <- rownames(z)
      z <- lapply(1 : ncol(z), FUN = function(i) sample(z[,i])) %>% do.call(cbind, .)
      z <- lapply(1 : nrow(z), FUN = function(i) sample(z[i,])) %>% do.call(rbind, .)
      dimnames(z) <- list(rn, paste0("tp", 1 : ncol(z)))
      z
    })
  })

for(cond in 1 : 4) {
  mod <- dynGENIE3(TS_rand[[cond]], time_points, regulators = tfs[[cond]], tree.method = "RF", K = "sqrt", verbose = TRUE)
  res <- get.link.list(mod$weight.matrix)
  res$id <- paste(res$regulatory.gene, res$target.gene, sep = "_")
  res <- res[order(res$id),]

  saveRDS(res, file = file.path(FILE, sprintf("perm_cond%s_RF_sqrt.rds", cond)))

  message(sprintf("dynGENIE3 on permutated condition %s: DONE!\n", cond))
}


message('end of script')
q('no')



TRI_DAT <- lapply(TS_dat, FUN = function(z) z)