library(reshape)
library(magrittr)
library(ggplot2)
library(data.table)
library(scales)

setwd("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/")

tri_regs <- c("g6430", "g6432")
tri_genes <- c("g6429", "g6431", "g6426", "g6434", "g6435", "g6436", "g6437")

TRI_GENES <- c("g6429" = "TRI4", "g6431" = "TRI5", "g6430" = "TRI6", "g6426" = "TRI8", "g6432" = "TRI10", "g6434" = "TRI11", "g6435" = "TRI12", "g6436" = "TRI13", "g6437" = "TRI14")

##reading the data in-------->>
regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]
#differentially expressed genes, per condition 
degs <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("l2FC_tables/CN_data_condition%s_DEG.txt", i), header = FALSE)[,1]
})
#adding tri-genes to the list of differentially expressed genes (they might or might not be there already)
degs <- lapply(degs, FUN = function(z) unique(c(tri_regs, tri_genes, z)))
#regulators, per condition (differentially expressed regulators + the two tri-gene regulators)
tfs <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("/home/ngrinber/quorn_grn/L2FC_filtered_data/transcription_factors_conditon_%s.txt", i), header = FALSE)[,1]
})
##reading the data in-------->>

#NOTE: only tri_gene analysis is the same for old and new
##new output
aracne_output <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps/cond_%s/network.txt", i), sep = "\t")
    tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp$cond <- i
    tmp$mode <- "new"
    tmp
}) %>% do.call(rbind,.)

#permuted output
aracne_perm <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps/cond_perm_%s/network.txt", i), sep = "\t")
    tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp$cond <- i
    tmp$mode <- "new_perm"
    tmp
}) %>% do.call(rbind,.)

aracne_tri <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps/tri_genes_cond_%s/network.txt", i), sep = "\t")
    if(nrow(tmp) > 0) {
        tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
        tmp <- tmp[order(tmp$id),]
        tmp$cond <- i
        tmp$mode <- "new_tri_genes"
        return(tmp)
    }
}) %>% do.call(rbind,.)

aracne_vst <- read.csv("/home/ngrinber/quorn_grn/aracne_bootstraps/vst_nobase/network.txt", sep = "\t")
aracne_vst$id <- paste(aracne_vst$Regulator, aracne_vst$Target, sep = "_")
aracne_vst <- aracne_vst[order(aracne_vst$id),]
aracne_vst$cond <- "all"
aracne_vst$mode <- "vst_new"

all_output <- rbind(aracne_output, aracne_perm, aracne_tri)

##old output
output_old <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps_old/cond_%s/network.txt", i), sep = "\t")
    tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp$cond <- i
    tmp$mode <- "old"
    tmp
}) %>% do.call(rbind,.)

#permuted output
perm_old <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps_old/cond_perm_%s/network.txt", i), sep = "\t")
    tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp$cond <- i
    tmp$mode <- "old_perm"
    tmp
}) %>% do.call(rbind,.)

tri_old <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps_old/tri_genes_cond_%s/network.txt", i), sep = "\t")
    if(nrow(tmp) > 0) {
        tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
        tmp <- tmp[order(tmp$id),]
        tmp$cond <- i
        tmp$mode <- "old_tri_genes"
        return(tmp)
    }
}) %>% do.call(rbind,.)

vst_old <- read.csv("/home/ngrinber/quorn_grn/aracne_bootstraps_old/vst_nobase/network.txt", sep = "\t")
vst_old$id <- paste(vst_old$Regulator, vst_old$Target, sep = "_")
vst_old <- vst_old[order(vst_old$id),]
vst_old$cond <- "all"
vst_old$mode <- "vst_old"

output_old <- rbind(output_old, perm_old, tri_old)

##reading the data in-------->>

#how many of the tri-gene connections are in the output? (per condition)
ids <- expand.grid(tri_regs, c(tri_regs, tri_genes))
ids <- paste(ids[,1], ids[,2], sep = "_")

table(all_output$cond, all_output$id %in% ids, all_output$mode)
table(output_old$cond, output_old$id %in% ids, output_old$mode)

all_output[all_output$id %in% ids,]
output_old[output_old$id %in% ids,]

dd <- all_output[all_output$id %in% ids & all_output$mode == "new_tri_genes" & all_output$cond == 2,]
plot(graph_from_data_frame(dd, directed = TRUE, vertices = NULL))

#concordance between differentially expressed genes between different conditions: Jaccard index
mat_id <- matrix(NA, nrow = 4, ncol = 4)
for(i in 1 : 3) {
    for(j in (i + 1) : 4) {
        mat_id[i, j] <- sum(aracne_output[aracne_output$cond == i,]$id %in% aracne_output[aracne_output$cond == j,]$id)/min(nrow(aracne_output[aracne_output$cond == i,]), nrow(aracne_output[aracne_output$cond == j,]))
    }
}

#comparing the new run with the old
comb_output <- rbind(all_output, output_old)
comb_output$mode <- factor(comb_output$mode, c("new", "old", "new_tri_genes", "old_tri_genes", "new_perm", "old_perm"))
ggplot(all_output, aes(y = MI, x = mode, fill = mode)) + geom_boxplot() + facet_wrap(cond ~ .)
ggplot(comb_output, aes(y = MI, x = mode, fill = mode)) + geom_boxplot() + facet_wrap(cond ~ .)

ggplot(rbind(all_output, aracne_vst), aes(y = MI, x = mode, fill = mode)) + geom_boxplot()
ggplot(rbind(output_old, vst_old), aes(y = MI, x = mode, fill = mode)) + geom_boxplot()    

#pairwise comaprison of MI between differentn conditions (for those pairs that exist for both)
par(mfrow = c(2, 3))
for(i in 1 : 3) {
    for(j in (i + 1) : 4) {
        ind <- intersect(aracne_output[[i]]$id, aracne_output[[j]]$id)
        mi1 <- aracne_output[[i]][aracne_output[[i]]$id %in% ind, ]$MI
        mi2 <- aracne_output[[j]][aracne_output[[j]]$id %in% ind, ]$MI
        plot(mi1, mi2, xlab = paste0("condition ", i), ylab = paste0("condition ", j))
        abline(0, 1, lty = 2, col = "firebrick2")
    }
}


##





























##################
dat <- read.csv("aracne_all_conditions/vst_conditon_1.txt", sep = '\t', row.names = 1)
z <- read.csv("l2FC_tables/CN_data_condition1_DEG.txt", header = FALSE)
dat_f <- read.csv("aracne_all_conditions/L2FC_filter_condition_1.txt", sep = "\t", row.names = 1)

dd <- dat[rownames(dat) %in% z$V1,]
dd <- dd[order(rownames(dd)),]
dat_f <- dat_f[order(rownames(dat_f)),]
# L2FC_filter_condition_1.txt -- filtered data according to l2FC_tables/CN_data_condition1_DEG.txt
# aracne_all_conditions/vst_conditon_1.txt -- raw data
# transcription_factors_conditon_1.txt --- WT_minion_TF_expression_only_head.txt subsetted by l2FC_tables/CN_data_condition1_DEG.txt, i.e. only use differentially expressed transcription factors

w <- read.csv("aracne_all_conditions/transcription_factors_conditon_1.txt", sep = "\t", header = FALSE)

