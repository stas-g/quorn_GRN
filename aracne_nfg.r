setwd("/home/ngrinber/quorn_grn/")

#creating permuted datasets
for(i in 1 : 4) {
    z <- read.csv(sprintf("L2FC_filtered_data/L2FC_filter_condition_%s.txt", i), sep = "\t")
    for(j in 2 : ncol(z)) {
        z[,j] <- sample(z[,j])
    }
    colnames(z)[1] <- ""
    write.table(z, file = sprintf("L2FC_filtered_data/L2FC_filter_condition_perm_%s.txt", i), sep = "\t", row.names = FALSE, quote = FALSE)
}