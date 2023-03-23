setwd("/home/ngrinber/projects/niab/gene_regulatory_network/all_data")
source("/home/ngrinber/projects/niab/gene_regulatory_network/install/dynGENIE3.R")
library(reshape2)
library(doRNG)
library(doParallel)
library(magrittr)
library(randomFunctions) #devtools::install_github('chr1swallace/random-functions') or remotes::install_github('chr1swallace/random-functions')

#reading parameters of the experiment in
ARGS <- getArgs()
i <- as.numeric(ARGS$i)

#path to results folder
FILE <- "/home/ngrinber/quorn_grn/dynGENIE3_results"

#reading in data and tfs
regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]
#4 conditions, 7 time points, 5 reps each (4 x 7 x 5 = 140)
all_data <- read.table("/home/ngrinber/quorn_grn/vst.txt", header = TRUE)

#parameter grid:  DynGenie3 hyperparameters and condition number
params <- expand.grid(tree.method = c('RF', 'ET'), c('all', 'sqrt'), cond = 1 : 4)
tree.method <- params[i,]$tree.method
K <- params[i,]$K
cond <- params[i,]$cond

message(sprintf("running dynGENIE3 on condition %s with tree.method = %s and K = %s", cond, tree.method, K))
#--------------------------------------------------------------------
#separating conditions; creating expr dfs per rep (list of 4 conditions containing list of 5 time-series, 7 times points each)
TS_data <- lapply(1 : 4, FUN = function(j) {
  cond <- all_data[(6 + (j - 1) * 35) : (5 + 35 * j)]
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
reg_ind <- which(rownames(TS_data[[1]][[1]]) %in% regs)

#run condition cond
mod <- dynGENIE3(TS_data[[cond]], time_points, regulators = reg_ind, tree.method = tree.method, K = K, verbose = TRUE)
res <- get.link.list(mod$weight.matrix)
res$ID <- paste(res$regulatory.gene, res$target.gene, sep = ".")
res <- res[order(res$ID),]

message("Done! Saving...")
saveRDS(res, file = file.path(FILE, sprintf("cond%s_%s_%s.rds", cond, tree.method, K)))

message(sprintf("dynGENIE3 on condition %s with tree.method = %s and K = %s: DONE!", cond, tree.method, K))


q('no')
###
