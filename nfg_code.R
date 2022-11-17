setwd("/home/ngrinber/projects/niab/gene_regulatory_network/all_data")
source("/home/ngrinber/projects/niab/gene_regulatory_network/install/dynGENIE3.R")
library(reshape2)
library(doRNG)
library(doParallel)
library(magrittr)
library(ggplot2)
library(scales)

regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]
all_data <- read.table("/home/ngrinber/quorn_grn/vst.txt", header = TRUE)

#--------------------------------------------------------------------
#dynGENIE
#separating conditions; creating expr dfs
TS_data <- lapply(1 : 4, FUN = function(j) {
  cond <- all_data[(6 + (j - 1) * 35) : (5 + 35 * j)]
  out <- lapply(1 : 5, FUN = function(i) {
    ind <- grep(paste0("_", i), colnames(cond))
    x <-  cond[, ind]
    x[, 1] <- x[, 1] + 0.000001
    colnames(x) <- paste0("tp", 1 : ncol(x))
    # dimnames(x) <- list(paste0("gene_", 1 : nrow(x)), paste0("sample_", 1 : ncol(x)))
    as.matrix(x)
  })
  out
})
names(TS_data) <- paste("cond", 1 : 4)

setwd("/home/ngrinber/projects/niab/gene_regulatory_network/install/")
tvec <- 0 : 6
names(tvec) <- colnames(TS_data[[1]][[1]])
time_points <- list(tvec, tvec, tvec, tvec, tvec)
# k <- 300
# ind <- sample(1 : nrow(TS_data[[1]][[1]]), k, replace = FALSE)
#a subset of first 500 genes
ind <- 1 : 500
ts_test <- lapply(TS_data[[1]], FUN = function(x) x[ind, ])
reg_ind <- which(rownames(ts_test[[1]]) %in% regs)


params <- expand.grid(c('RF', 'ET'), c('all', 'sqrt'))
colnames(params) <- c("tree.method", "K")

res <- lapply(1 : nrow(params), FUN = function(i) {
  mod <- dynGENIE3(ts_test[1], time_points[1], regulators = reg_ind, tree.method = params[i,]$tree.method, K = params[i,]$K)
  res <- get.link.list(mod$weight.matrix)
  res$ID <- paste(res$regulatory.gene, res$target.gene, sep = ".")
  res <- res[order(res$ID),]
  res$tree.method <- params[i,]$tree.method
  res$K <- params[i,]$K
  res
})
names(res) <- paste(params[, 1], params[, 2], sep = ".")

res_all <- res[[1]][, c("regulatory.gene", "target.gene", "ID")]
for(i in 1 : length(res)) res_all[[names(res)[i]]] <- res[[i]]$weight
res_all2 <- do.call(rbind, res)
res_all2$method <- rep(names(res), each = nrow(res[[1]]))
res_all2$exp <- "test"

#distribution of weights for various parameter combinations
par(mfrow = c(2, 2))
for(i in 1 : 4) hist(res[[i]]$weight[res[[i]]$weight > 0], main = names(res)[i], xlab = "weight")
#pairwise comparison
par(mfrow = c(2, 2))
plot(res[['RF.all']]$weight, res[['RF.sqrt']]$weight, xlab = "RF.all", ylab = "RF.sqrt"); abline(0, 1, col = 'firebrick1', lty = 2)
plot(res[['ET.all']]$weight, res[['ET.sqrt']]$weight, xlab = "ET.all", ylab = "ET.sqrt"); abline(0, 1, col = 'firebrick1', lty = 2)
plot(res[['RF.all']]$weight, res[['ET.all']]$weight, xlab = "RF.all", ylab = "ET.all"); abline(0, 1, col = 'firebrick1', lty = 2)
plot(res[['RF.sqrt']]$weight, res[['ET.sqrt']]$weight, xlab = "RF.sqrt", ylab = "ET.sqrt"); abline(0, 1, col = 'firebrick1', lty = 2)

boxplot(weight ~ method, data = res_all2[res_all2$weight > 0, ])

#permutation tests
ts_rand <- lapply(ts_test, FUN = function(x){
    rn <- rownames(x)
    x <- lapply(1 : ncol(x), FUN = function(i) sample(x[,i])) %>% do.call(cbind, .)
    rownames(x) <- rn
    x
  })

res_rand <- lapply(1 : nrow(params), FUN = function(i) {
  mod <- dynGENIE3(ts_rand[1], time_points[1], regulators = reg_ind, tree.method = params[i,]$tree.method, K = params[i,]$K)
  res <- get.link.list(mod$weight.matrix)
  res$ID <- paste(res$regulatory.gene, res$target.gene, sep = ".")
  res <- res[order(res$ID),]
  res$tree.method <- params[i,]$tree.method
  res$K <- params[i,]$K
  res
})
names(res_rand) <- paste(params[, 1], params[, 2], sep = ".")
res_rand2 <- do.call(rbind, res_rand)
res_rand2$method <- rep(names(res_rand), each = nrow(res_rand[[1]]))
res_rand2$exp <- "random"

#comparing real run with permutation run
table(res_all2$weight == 0, res_all2$method)
table(res_rand2$weight == 0, res_rand2$method)

res_combined <- rbind(res_all2[res_all2$weight > 0,], res_rand2[res_rand2$weight > 0, ])
ggplot(res_combined, aes(x = weight, fill = as.factor(exp))) + geom_histogram(alpha = 0.5, position = 'identity') + facet_wrap(. ~ K + tree.method)
ggplot(res_combined, aes(x = exp, y = weight, fill = as.factor(exp))) + geom_boxplot(alpha = 0.5) + facet_wrap(. ~ K + tree.method)

par(mfrow = c(2, 2))
for(i in 1 : length(res)) {
  plot(res[[i]]$weight, res_rand[[i]]$weight, xlab = 'test', ylab = 'random')
  abline(0, 1, col = 'firebrick1', lty = 2)
}

saveRDS(res_combined, file = "/home/ngrinber/quorn_grn/res_combined_test.rds")
#--------------------------------------------------------------------
#ARACNE
i <- j <- 1
cond <- all_data[(6 + (j - 1) * 35) : (5 + 35 * j)]
# x <- cond[1 : 300, (1 + (i - 1) * 5) : (i * 5)]
x <- cond[1 : 500, ]
rr <- regs[which(regs %in% rownames(x))]

write.table(x, file = "/home/ngrinber/quorn_grn/test_expr.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(rr, file = "/home/ngrinber/quorn_grn/test_tfs.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(do.call(cbind, ts_rand), file = "/home/ngrinber/quorn_grn/rand_expr.txt", sep = "\t", quote = FALSE, col.names = NA)

#Run aracne Condition1
#Gerate bootstrap values
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/test_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap_rand \
--pvalue 1E-8 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/test_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap_rand \
--pvalue 1E-8 \
-t /home/ngrinber/quorn_grn/test_tfs.txt \
--seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap_rand \
--consolidate

#------------------------------------------
#Gerate bootstrap values
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/rand_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap_rand \
--pvalue 1E-8 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/rand_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap_rand \
--pvalue 1E-8 \
-t /home/ngrinber/quorn_grn/test_tfs.txt \
--seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap_rand \
--consolidate

aracne <- read.table("/home/ngrinber/quorn_grn/aracne_test/bootstrap/network.txt", header = TRUE)
arand <- read.table("/home/ngrinber/quorn_grn/aracne_test/bootstrap_rand/network.txt", header = TRUE)
aracne$ID <- paste(aracne$Regulator, aracne$Target, sep = ".")
arand$ID <- paste(arand$Regulator, arand$Target, sep = ".")
aracne$exp <- "aracne"
arand$exp <- "random"

#compare random vs aracne
par(mfrow = c(1, 2))
r <- range(c(aracne$MI, arand$MI))
hist(aracne$MI, col = alpha('firebrick1', 0.5), prob = TRUE)
hist(arand$MI, col = alpha('dodgerblue', 0.5), add = TRUE, prob = TRUE)
boxplot(MI ~ exp, data = rbind(aracne, arand))

#comparing to dng3
dng3 <- res[["RF.all"]]
dng3$ID <- paste(dng3$regulatory.gene, dng3$target.gene, sep = ".")
ind <- intersect(aracne$ID, dng3$ID)

par(mfrow = c(1, 2))
hist(dng3$weight)
hist(aracne$MI)

dng3$in_aracne <- dng3$ID %in% ind
aracne_ <- aracne[aracne$ID %in% ind, ]
aracne_ <- aracne_[order(aracne_$ID),]
dng3_ <- dng3[dng3$ID %in% ind, ]
dng3_ <- dng3_[order(dng3_$ID),]

par(mfrow = c(1, 2))
plot(dng3_$weight, aracne_$MI)
boxplot(dng3$weight ~ dng3$in_aracne)

hist(dng3[!(dng3$id %in% ind),]$weight)
hist(dng3_$weight)







##RUBBISH--------------------------------------------------------------------------------------------------->>>

aracne <- read.csv("/home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/C1.1_results/network.txt", sep = "\t")
# aracne <- read.csv("all_data_aracne/condition_1_aracne.txt", sep = "\t")
dng3 <- read.csv("all_data_dyngenie/condition_1_dyngenie.txt", sep = "\t")

aracne$id <- paste(aracne$Regulator, aracne$Target, sep = ".")
dng3$id <- paste(dng3$regulatory.gene, dng3$target.gene, sep = ".")
ind <- intersect(aracne$id, dng3$id)
aracne_ <- aracne[aracne$id %in% ind, ]
aracne_ <- aracne_[order(aracne_$id),]
dng3_ <- dng3[dng3$id %in% ind, ]
dng3_ <- dng3_[order(dng3_$id),]

plot(dng3_$weight, aracne_$MI)
par(mfrow = c(1, 2))
hist(-log10(dng3_$weight))
hist(-log10(aracne_$weight))




/home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar
/home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP

C1.1_result












java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/test/matrix.txt \
-o /home/ngrinber/quorn_grn/aracne_test \
--pvalue 0 \
--seed 1 \
--calculateThreshold

for i in {1..100}
do
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/test/matrix.txt \
-o /home/ngrinber/quorn_grn/aracne_test \
--pvalue 0 \
-t /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/test/tfs.txt
--seed $i
done

java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o /home/ngrinber/quorn_grn/aracne_test \
--consolidate


java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C1_t1_7_uncollapsed.txt \
-o /home/ngrinber/quorn_grn/aracne_test \
--pvalue 0 \
--seed 1 \
--calculateThreshold












source("/home/ngrinber/projects/niab/gene_regulatory_network/install/dynGENIE3.R")
setwd("/home/ngrinber/projects/niab/gene_regulatory_network/install")

TS1 <- read.expr.matrix("time_series_1.txt", form="rows.are.samples")
TS2 <- read.expr.matrix("time_series_2.txt", form="rows.are.samples")
TS3 <- read.expr.matrix("time_series_3.txt", form="rows.are.samples")
time.points <- list(TS1[1,], TS2[1,], TS3[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),])
TS.data_r <- lapply(TS.data, FUN = function(x) {
  rn <- rownames(x)
  cn <- colnames(x)
  x <- lapply(1 : ncol(x), FUN = function(i) sample(x[,i])) %>% do.call(cbind, .)
  rownames(x) <- rn
  colnames(x) <- cn
  x
})
res <- dynGENIE3(TS.data,time.points)
res <- get.link.list(res$weight.matrix)

ar_res <- read.table("/home/ngrinber/quorn_grn/aracne_test/bootstrap3/network.txt", header = TRUE)

res2 <- dynGENIE3(TS.data,time.points)
res2 <- get.link.list(res2$weight.matrix)
res_rand <- dynGENIE3(TS.data_r,time.points)
res_rand <- get.link.list(res_rand$weight.matrix)
par(mfrow = c(1, 2)); hist(res$weight); hist(res_rand$weight)

res <- res[order(res$regulatory.gene, res$target.gene),]
res2 <- res2[order(res2$regulatory.gene, res2$target.gene),]
res_rand <- res_rand[order(res_rand$regulatory.gene, res_rand$target.gene),]


par(mfrow = c(1, 2))
hist(res$weight); hist(res_rand$weight)

dat <- do.call(cbind, TS.data)
write.table(dat, file = "/home/ngrinber/quorn_grn/dyn_expr.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(rownames(dat), file = "/home/ngrinber/quorn_grn/dyn_tfs.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


#Gerate bootstrap values
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/dyn_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap3 \
--pvalue 1E-8 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/dyn_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap3 \
--pvalue 1E-8 \
-t /home/ngrinber/quorn_grn/dyn_tfs.txt \
--seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap3 \
--consolidate





##
