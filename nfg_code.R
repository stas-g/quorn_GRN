setwd("/home/ngrinber/projects/niab/gene_regulatory_network/all_data")
source("/home/ngrinber/projects/niab/gene_regulatory_network/install/dynGENIE3.R")
library(reshape2)
library(doRNG)
library(doParallel)
library(magrittr)

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
ind <- 1 : 500
ts_test <- lapply(TS_data[[1]], FUN = function(x) x[ind, ])
reg_ind <- which(rownames(ts_test[[1]]) %in% regs)

res_test <- dynGENIE3(ts_test[1], time_points[1], regulators = reg_ind, tree.method = "RF", K = "all")
res_test0 <- get.link.list(res_test$weight.matrix)
dng3 <- res_test0[res_test0$weight > 0, ]
dng3 <- dng3[order(dng3$regulatory.gene, dng3$target.gene),]

res_test <- dynGENIE3(ts_test[1][1], time_points[1][1], regulators = reg_ind, tree.method = "ET")
res_test0 <- get.link.list(res_test$weight.matrix)
dng3_2 <- res_test0[res_test0$weight > 0, ]
dng3_2 <- dng3_2[order(dng3_2$regulatory.gene, dng3_2$target.gene),]

saveRDS(dng3, file = "/home/ngrinber/quorn_grn/dyng3_test.rds")

res <- dynGENIE3(TS_data[[1]], time_points, verbose = TRUE)

##permutation test
ts_rand <- lapply(TS_data[[1]], FUN = function(x){
    x <- x[1:500,]
    rn <- rownames(x)
    x <- lapply(1 : ncol(x), FUN = function(i) sample(x[,i])) %>% do.call(cbind, .)
    rownames(x) <- rn
    x
  })
res_rand <- dynGENIE3(ts_rand[1], time_points[1], regulators = reg_ind, K = 'all')
res_rand <- get.link.list(res_rand$weight.matrix)
res_rand <- res_rand[res_rand$weight > 0, ]
res_rand <- res_rand[order(res_rand$regulatory.gene, res_rand$target.gene),]

res_rand2 <- dynGENIE3(ts_rand[1], time_points[1], regulators = reg_ind)
res_rand2 <- get.link.list(res_rand2$weight.matrix)
res_rand2 <- res_rand2[res_rand2$weight > 0, ]
res_rand2 <- res_rand2[order(res_rand2$regulatory.gene, res_rand2$target.gene),]

plot(res_rand$weight, res_rand2$weight)

#--------------------------------------------------------------------
#ARACNE
i <- j <- 1
cond <- all_data[(6 + (j - 1) * 35) : (5 + 35 * j)]
# x <- cond[1 : 300, (1 + (i - 1) * 5) : (i * 5)]
x <- cond[1 : 200, ]
rr <- regs[which(regs %in% rownames(TS_data[[1]][[1]][1:2000,]))]
# rn <- rownames(x)
# x <- lapply(1 : ncol(x), FUN = function(i) sample(x[,i])) %>% do.call(cbind, .)
# rownames(x) <- rn

write.table(x, file = "/home/ngrinber/quorn_grn/test_expr.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(rr, file = "/home/ngrinber/quorn_grn/test_tfs.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(do.call(cbind, ts_rand), file = "/home/ngrinber/quorn_grn/rand_expr.txt", sep = "\t", quote = FALSE, col.names = NA)

#Run aracne Condition1
#Gerate bootstrap values
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/test_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap \
--pvalue 1E-8 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/test_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap \
--pvalue 1E-8 \
-t /home/ngrinber/quorn_grn/test_tfs.txt \
--seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap \
--consolidate

#------------------------------------------
#Gerate bootstrap values
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/rand_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap2 \
--pvalue 1E-8 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/rand_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap2 \
--pvalue 1E-8 \
-t /home/ngrinber/quorn_grn/test_tfs.txt \
--seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o /home/ngrinber/quorn_grn/aracne_test/bootstrap2 \
--consolidate

aracne <- read.table("/home/ngrinber/quorn_grn/aracne_test/bootstrap/network.txt", header = TRUE)
arand <- read.table("/home/ngrinber/quorn_grn/aracne_test/bootstrap2/network.txt", header = TRUE)

aracne$id <- paste(aracne$Regulator, aracne$Target, sep = ".")
dng3$id <- paste(dng3$regulatory.gene, dng3$target.gene, sep = ".")
ind <- intersect(aracne$id, dng3$id)
dng3$in_aracne <- dng3$id %in% ind
aracne_ <- aracne[aracne$id %in% ind, ]
aracne_ <- aracne_[order(aracne_$id),]
dng3_ <- dng3[dng3$id %in% ind, ]
dng3_ <- dng3_[order(dng3_$id),]

plot(dng3_$weight, aracne_$MI)
plot(dng3_$weight, -log10(aracne_$pvalue))

par(mfrow = c(1, 2))
hist(dng3[!(dng3$id %in% ind),]$weight)
hist(dng3_$weight)
boxplot(dng3$weight ~ dng3$in_aracne)







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
