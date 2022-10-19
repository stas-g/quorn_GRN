setwd("/home/ngrinber/projects/niab/gene_regulatory_network/all_data")
source("/home/ngrinber/projects/niab/gene_regulatory_network/install/dynGENIE3.R")
library(reshape2)
library(doRNG)
library(doParallel)

regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]
all_data <- read.table("/home/ngrinber/quorn_grn/vst.txt", header = TRUE)

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

#--------------------------------------------------------------------
#dynGENIE
#separating conditions; creating expr dfs
TS_data <- lapply(1 : 4, FUN = function(j) {
  cond <- all_data[(6 + (j - 1) * 35) : (5 + 35 * j)]
  out <- lapply(1 : 7, FUN = function(i) {
    x <- cond[(1 + (i - 1) * 5) : (i * 5)]
    x[, 1] <- x[, 1] + 0.000001
    colnames(x) <- paste0("tp", 1 : ncol(x))
    # dimnames(x) <- list(paste0("gene_", 1 : nrow(x)), paste0("sample_", 1 : ncol(x)))
    as.matrix(x)
  })
  out
})
names(TS_data) <- paste("cond", 1 : 4)

setwd("/home/ngrinber/projects/niab/gene_regulatory_network/install/")
tvec <- 1 : 5
names(tvec) <- colnames(TS_data[[1]][[1]])
time_points <- list(tvec, tvec, tvec, tvec, tvec, tvec, tvec)
ts_test <- lapply(TS_data[[1]], FUN = function(x) x[1:100, ])
reg_ind <- which(rownames(ts_test[[1]]) %in% regs)
res_test <- dynGENIE3(ts_test, time_points, regulators = reg_ind)
res_test0 <- get.link.list(res_test$weight.matrix)


res <- dynGENIE3(TS_data[[1]], time_points)

#--------------------------------------------------------------------
#ARACNE
i <- j <- 1
cond <- all_data[(6 + (j - 1) * 35) : (5 + 35 * j)]
x <- cond[1 : 100, (1 + (i - 1) * 5) : (i * 5)]
write.table(x, file = "/home/ngrinber/quorn_grn/test_expr.txt", sep = "\t", quote = FALSE, col.names = NA)

#Run aracne Condition1
#Gerate bootstrap values
java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/quorn_grn/test_expr.txt \
-o /home/ngrinber/quorn_grn/aracne_test \
--pvalue 0 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C1_t1_7_uncollapsed.txt \
-o C1.1_results \
--pvalue 0 \
-t WT_minion_TF_expression_only_head.txt \
--seed $i
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o C1.1_results \
--consolidate





/home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar
/home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP

C1.1_result












java -Xmx5G -jar /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e /home/ngrinber/projects/niab/gene_regulatory_network/ARACNe-AP/test/matrix.txt \
-o /home/ngrinber/quorn_grn/aracne_test \
--pvalue 0 \
--seed 1 \
--calculateThreshold

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
res <- dynGENIE3(TS.data,time.points)
##
