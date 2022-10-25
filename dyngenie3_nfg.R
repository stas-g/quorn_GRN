setwd("/home/ngrinber/projects/niab/gene_regulatory_network/all_data")
source("/home/ngrinber/projects/niab/gene_regulatory_network/install/dynGENIE3.R")
library(reshape2)
library(doRNG)
library(doParallel)
library(magrittr)

#reading in data and tfs
regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]
all_data <- read.table("/home/ngrinber/quorn_grn/vst.txt", header = TRUE)

#--------------------------------------------------------------------
#separating conditions; creating expr dfs per rep (list of conditions containing list of 5 time-series, 7 times points each)
TS_data <- lapply(1 : 4, FUN = function(j) {
  cond <- all_data[(6 + (j - 1) * 35) : (5 + 35 * j)]
  out <- lapply(1 : 5, FUN = function(i) {
    ind <- grep(paste0("_", i), colnames(cond))
    x <-  cond[, ind]
    x[, 1] <- x[, 1] + 0.000001
    colnames(x) <- paste0("tp", 1 : ncol(x))
    as.matrix(x)
  })
  out
})
names(TS_data) <- paste("cond", 1 : 4)

tvec <- 0 : 6
names(tvec) <- colnames(TS_data[[1]][[1]])
time_points <- list(tvec, tvec, tvec, tvec, tvec)

reg_ind <- which(rownames(TS_data[[1]][[1]]) %in% regs)
#run condition i
i <- 1 #i in (1,2,3,4)
res <- dynGENIE3(TS_data[[i]], time_points, regulators = reg_ind)
res2 <- dynGENIE3(TS_data[[i]], time_points, regulators = reg_ind, K = "all")
