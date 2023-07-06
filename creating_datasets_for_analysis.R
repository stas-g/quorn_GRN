setwd("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/")
MY_FOLDER = "/home/ngrinber/quorn_grn/L2FC_filtered_data"

##reading the data in-------->>
#base expression
base <- read.csv("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/aracne_all_conditions/vst.txt", header = TRUE, sep = "\t", row.names = 1)[, 1 : 5]
raw_dat <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("aracne_all_conditions/vst_conditon_%s.txt", i), sep = '\t', row.names = 1)
})
#list of regulator genes (un-filtered)
regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]
#list of genes with signal per condition
degs <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("l2FC_tables/CN_data_condition%s_DEG.txt", i), header = FALSE)[,1]
})
##reading the data in-------->>

#list of tri-genes; those expected to be regulators (2) and the rest (7)
tri_regs <- c("g6430", "g6432")
tri_genes <- c("g6429", "g6431", "g6426", "g6434", "g6435", "g6436", "g6437")

#adding tri-genes to differentially expressed genes
degs <- lapply(degs, FUN = function(z) unique(c(tri_regs, tri_genes, z)))
#creating subsets of expression data with just the differentially expressed genes + tri-genes 
for(i in 1 : 4) {
    zz <- read.csv(sprintf("aracne_all_conditions/vst_conditon_%s.txt", i), sep = '\t')
    zz <- zz[order(zz[,1]),]
    zz <- zz[zz[,1] %in% degs[[i]],]
    colnames(zz)[1] <- ''

    write.table(zz, file = file.path(MY_FOLDER, sprintf("L2FC_filter_condition_%s.txt", i)), sep = "\t", row.names = FALSE, quote = FALSE)
}

#adding tri-genes to the list of regulators
regs <- unique(c(tri_regs, regs))
#subsetting regulators according to which are expressed in each condition
for(i in 1 : 4) {
    z <- regs[regs %in% degs[[i]]]
    write.table(z, file = file.path(MY_FOLDER, sprintf("transcription_factors_conditon_%s.txt", i)), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


#######CREATING DATASETS FOR DYNGENIE3
#creating permuted datasets-------->>
for(i in 1 : 4) {
    z <- read.csv(file.path(MY_FOLDER, sprintf("L2FC_filter_condition_%s.txt", i)), sep = "\t")
    for(j in 2 : ncol(z)) {
        z[,j] <- sample(z[,j])
    }
    colnames(z)[1] <- ""
    write.table(z, file = file.path(MY_FOLDER, sprintf("L2FC_filter_condition_perm_%s.txt", i)), sep = "\t", row.names = FALSE, quote = FALSE)
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

write.table(vst_base, file = file.path(MY_FOLDER, "vst_L2FC_filtered.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(vst, file = file.path(MY_FOLDER, "vst_L2FC_filtered_nobase.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(regs[regs %in% deg], file = file.path(MY_FOLDER, "tfs_combined.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)

#creating data just for tri-genes
tri_dat <- lapply(1 : 4, FUN = function(i) {
    tmp <- raw_dat[[i]][rownames(raw_dat[[i]]) %in% c(tri_regs, tri_genes), ]
    tmp <- data.frame(gene = rownames(tmp), tmp)
    write.table(tmp, file = file.path(MY_FOLDER, sprintf("L2FC_filter_condition_tri_%s.txt", i)), sep = "\t", row.names = FALSE, quote = FALSE)
})
write.table(tri_regs, file = file.path(MY_FOLDER, "tri_regs.txt"), row.names = FALSE, quote = FALSE, col.names = FALSE)






