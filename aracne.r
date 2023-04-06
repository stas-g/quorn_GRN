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

# plotting gene expression for the tri-genes in the 4 conditions (usinng raw data)
tri_dat <- lapply(1 : 4, FUN = function(i) {
    tmp <- raw_dat[[i]][rownames(raw_dat[[i]]) %in% c(tri_regs, tri_genes), ]
    tmp <- lapply(1 : 7, FUN = function(k) {
        rowMeans(tmp[, (1 + 5 * (k - 1)) : (5 * k)])
    }) %>% do.call(cbind, .) %>% data.frame()
    colnames(tmp) <- 1 : ncol(tmp)
    tmp$gene <- rownames(tmp)
    tmp$cond <- i
    tmp
}) %>% do.call(rbind, .)
tri_dat0 <- base[rownames(base) %in% c(tri_regs, tri_genes),]
tri_dat0 <- rowMeans(tri_dat0)

tri_dat <- melt(tri_dat, id = c("gene", "cond"))
setnames(tri_dat, c("variable", "value"), c("tp", "expr"))
tri_dat$cond <- factor(tri_dat$cond)
ggplot(data = tri_dat, aes(x = as.numeric(tp), y = expr, colour = gene)) + geom_line() + facet_wrap(cond ~ .)

#run on averaged data
#run just on tri group

#John's output
aracne_output0 <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("aracne_all_conditions/link_list_condition_%s/network.txt", i), sep = "\t")
    tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp$cond <- i
    tmp$human <- "john"
    tmp
}) %>% do.call(rbind,.)

#my output
aracne_output <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps/cond_%s/network.txt", i), sep = "\t")
    tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp$cond <- i
    tmp$human <- "stas"
    tmp
}) %>% do.call(rbind,.)

#permuted output
aracne_perm <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps/cond_perm_%s/network.txt", i), sep = "\t")
    tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp$cond <- i
    tmp$human <- "stas_perm"
    tmp
}) %>% do.call(rbind,.)

aracne_tri <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps/tri_genes_cond_%s/network.txt", i), sep = "\t")
    if(nrow(tmp) > 0) {
        tmp$id <- paste(tmp$Regulator, tmp$Target, sep = "_")
        tmp <- tmp[order(tmp$id),]
        tmp$cond <- i
        tmp$human <- "stas_tri_genes"
        return(tmp)
    }
}) %>% do.call(rbind,.)

all_output <- rbind(aracne_output0, aracne_output, aracne_perm, aracne_tri)
##reading the data in-------->>

#how many of the tri-gene connections are in the output? (per condition)
ids <- expand.grid(tri_regs, tri_genes)
ids <- paste(ids[,1], ids[,2], sep = "_")

table(all_output$cond, all_output$id %in% ids, all_output$human)

#how many of tri__regs in TFs and in the list of differentially expressed genes?
sapply(tfs, FUN = function(z) sum(c(tri_regs, tri_genes) %in% z))
sapply(degs, FUN = function(z) sum(c(tri_regs, tri_genes) %in% z))


#concordance between differentially expressed genes between different conditions: Jaccard index
mat_deg <- matrix(NA, nrow = 4, ncol = 4)
mat_tf <- matrix(NA, nrow = 4, ncol = 4)
mat_id <- matrix(NA, nrow = 4, ncol = 4)
for(i in 1 : 3) {
    for(j in (i + 1) : 4) {
        mat_deg[i, j] <- sum(degs[[i]] %in% degs[[j]])/min(length(degs[[i]]), length(degs[[j]]))
        mat_tf[i, j] <- sum(tfs[[i]] %in% tfs[[j]])/min(length(tfs[[i]]), length(tfs[[j]]))
        mat_id[i, j] <- sum(aracne_output[[i]]$id %in% aracne_output[[j]]$id)/min(nrow(aracne_output[[i]]), nrow(aracne_output[[j]]))
    }
}

#comparing my run with John's
par(mfrow = c(2, 2))
for(i in 1 : 4) hist(-log10(aracne_output[[i]]$pvalue), main = paste0("condition ", i))

ggplot(all_output, aes(y = MI, x = human, fill = human)) + geom_boxplot() + facet_wrap(cond ~ .)
ggplot(all_output, aes(y = ..density.., x = MI, fill = human)) + geom_histogram(alpha = 0.5) + facet_wrap(cond ~ .)

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







#creating permuted datasets-------->>
setwd("/home/ngrinber/quorn_grn/")

for(i in 1 : 4) {
    z <- read.csv(sprintf("L2FC_filtered_data/L2FC_filter_condition_%s.txt", i), sep = "\t")
    for(j in 2 : ncol(z)) {
        z[,j] <- sample(z[,j])
    }
    colnames(z)[1] <- ""
    write.table(z, file = sprintf("L2FC_filtered_data/L2FC_filter_condition_perm_%s.txt", i), sep = "\t", row.names = FALSE, quote = FALSE)
}
#creating permuted datasets-------->>

#combining all conditions
deg <- unique(unlist(degs))

vst_base <- read.table("/home/ngrinber/quorn_grn/vst.txt", header = TRUE)
vst_base <- vst_base[rownames(vst_base) %in% deg, ]
vst <- vst_base[, -(1 : 5)]
vst <- data.frame(rownames(vst), vst)
vst_base <- data.frame(rownames(vst_base), vst_base)
colnames(vst)[1] <- ''
colnames(vst_base)[1] <- ''

write.table(vst_base, file = "L2FC_filtered_data/vst_L2FC_filtered.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(vst, file = "L2FC_filtered_data/vst_L2FC_filtered_nobase.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(regs[regs %in% deg], file = "L2FC_filtered_data/tfs_combined.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

#creating data just for tri-genes
tri_dat <- lapply(1 : 4, FUN = function(i) {
    tmp <- raw_dat[[i]][rownames(raw_dat[[i]]) %in% c(tri_regs, tri_genes), ]
    tmp <- data.frame(gene = rownames(tmp), tmp)
    write.table(tmp, file = sprintf("L2FC_filtered_data/L2FC_filter_condition_tri_%s.txt", i), sep = "\t", row.names = FALSE, quote = FALSE)
})
write.table(tri_regs, file = "L2FC_filtered_data/tri_regs.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)




































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

