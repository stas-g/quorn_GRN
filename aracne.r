list.files("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/aracne_all_conditions")

list.files("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/l2FC_tables")

setwd("/home/ngrinber/projects/niab/gene_regulatory_network/carbon_nitrogen_data/")
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

regs <- read.table("/home/ngrinber/projects/niab/gene_regulatory_network/install/WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")[, 1]
