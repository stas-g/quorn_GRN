library(reshape)
library(magrittr)
library(ggplot2)
library(data.table)
library(scales)
library(igraph)

setwd("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/")

tri_regs <- c("g6430", "g6432")
tri_genes <- c("g6429", "g6431", "g6426", "g6434", "g6435", "g6436", "g6437")

TRI_GENES <- c("g6429" = "TRI4", "g6431" = "TRI5", "g6430" = "TRI6", "g6426" = "TRI8", "g6432" = "TRI10", "g6434" = "TRI11", "g6435" = "TRI12", "g6436" = "TRI13", "g6437" = "TRI14")

##critical value for dyngenir3--------->>
dyn_perm <- lapply(1 : 4, FUN = function(i) {
    print(i)
    out <- readRDS(sprintf("/home/ngrinber/quorn_grn/dynGENIE3_results/perm_cond%s_RF_sqrt.rds", i))
    out <- out[out$weight > 0,]
    out
}) %>% do.call(rbind, .)

crit_val <- quantile(dyn_perm$weight, 0.95)
##critical value for dyngenir3--------->>

##reading the data & results in-------->>
raw_dat <- lapply(1 : 4, FUN = function(i) {
    read.csv(sprintf("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/aracne_all_conditions/vst_conditon_%s.txt", i), sep = '\t', row.names = 1)
})

aracne_output <- lapply(1 : 4, FUN = function(i) {
    tmp <- read.csv(sprintf("/home/ngrinber/quorn_grn/aracne_bootstraps/cond_%s/network.txt", i), sep = "\t")
    tmp$cond <- i
    tmp$method <- "aracne"
    tmp$id <- paste(tmp$Regulator, tmp$Target, tmp$cond, sep = "_")
    tmp <- tmp[order(tmp$id),]
    tmp
}) %>% do.call(rbind,.) %>% data.table

dyn_output <- lapply(1 : 4, FUN = function(i) {
    print(i)
    out <- readRDS(sprintf("/home/ngrinber/quorn_grn/dynGENIE3_results/cond%s_RF_sqrt.rds", i))
    out <- out[out$weight > 0,]
    out$regulatory.gene <- as.character(out$regulatory.gene)
    out$target.gene <- as.character(out$target.gene)
    out$cond <- i
    out$method <- "dyngenie3"
    out$id <- paste(out$regulatory.gene, out$target.gene, out$cond, sep = "_")
    out <- out[order(out$id),]
    out
}) %>% do.call(rbind, .) %>% data.table

#subset of results only featuring tri_genes
aracne_trigene <- aracne_output[Regulator %in% c(tri_regs, tri_genes) & Target %in% c(tri_regs, tri_genes)]
dyn3_trigene <- dyn_output[regulatory.gene %in% c(tri_regs, tri_genes) & target.gene %in% c(tri_regs, tri_genes)]
dyn3_trigene <- dyn3_trigene[weight > crit_val]

aracne_trigene[, c("Regulator", "Target") := list(TRI_GENES[Regulator], TRI_GENES[Target])]
aracne_trigene[, "id" := paste(Regulator, Target, cond, sep = "_")]
dyn3_trigene[, c("regulatory.gene", "target.gene") := list(TRI_GENES[regulatory.gene], TRI_GENES[target.gene])]
dyn3_trigene[, "id" := paste(regulatory.gene, target.gene, cond, sep = "_")]
##reading the data & results in-------->>

col_ <- ifelse(dyn3_trigene$cond == 1, "#E69F00", ifelse(dyn3_trigene$cond == 2, "#56B4E9", "#CC79A7"))
plot(graph_from_data_frame(aracne_trigene), vertex.size = 23, edge.width = aracne_trigene$MI, edge.color = as.factor(aracne_trigene$cond))
plot(graph_from_data_frame(dyn3_trigene), vertex.size = 23, edge.width = aracne_trigene$weight, edge.color = as.factor(dyn3_trigene$cond))
plot(graph_from_data_frame(dyn3_trigene), vertex.size = 23, edge.width = aracne_trigene$weight, edge.color = col_)

table(dyn3_trigene$id %in% aracne_trigene$id)
#concordance between dyngenie3 and aracne
dyn_output_signif <- dyn_output[weight > crit_val,]
table(aracne_output$id %in% dyn_output_signif$id)

aracne_output[, "in_dyng3" := id %in% dyn_output_signif$id] 
ggplot(aracne_output, aes(x = factor(in_dyng3), y = MI, colour = in_dyng3)) + geom_boxplot() + xlab("is gene pair significant in DynGenie3?") + ylab("Mutual Information score (Aracne)")

compare_df <- data.frame(aracne_MI = aracne_output[id %in% dyn_output_signif$id]$MI, dyng3_weight = dyn_output_signif[id %in% aracne_output$id]$weight, cond = aracne_output[id %in% dyn_output_signif$id]$cond)
ggplot(compare_df, aes(x = aracne_MI, y = dyng3_weight)) + geom_point()
ggplot(compare_df, aes(x = aracne_MI, y = dyng3_weight)) + geom_point() + facet_wrap(cond ~ .)


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
tri_dat$mode <- factor(ifelse(tri_dat$gene %in% tri_genes, "target", "reg"), levels = c("target", "reg"))
tri_dat$gene <- TRI_GENES[tri_dat$gene]

tri_dat <- melt(tri_dat, id = c("gene", "cond", "mode"))
setnames(tri_dat, c("variable", "value"), c("tp", "expr"))
tri_dat$cond <- factor(tri_dat$cond)

ggplot(data = tri_dat, aes(x = as.numeric(tp), y = expr, colour = gene, linetype = mode)) + scale_colour_brewer(palette = "Set1") + xlab("time point") + ylab("expression") + geom_line() + facet_wrap(cond ~ .)


