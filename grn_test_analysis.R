setwd("/home/ngrinber/projects/niab/gene_regulatory_network/compare_swing_dyngenie_aracne")

aracne <- read.table("aracne_compare_result.txt", header = TRUE)
dyng3 <- read.table("dyn_compare_result.txt")
swing <- read.table("SWING_compare_result.txt", header = TRUE)

dyng3$id <- paste(dyng3$regulatory.gene, dyng3$target.gene, sep = ".")
swing$id <- paste(swing$regulatory.gene, swing$target.gene, sep = ".")
aracne$id <- paste(aracne$Regulator, aracne$Target, sep = ".")

par(mfrow = c(1, 3))
hist(dyng3$weight)
hist(aracne$MI)
hist(swing$weight)

dyng3 <- dyng3[order(dyng3$id),]
swing <- swing[order(swing$id),]
aracne <- aracne[order(aracne$id),]

par(mfrow = c(2, 2))
plot(dyng3$weight, swing$weight, xlab = "DYNG3", ylab = "SWING")
plot(aracne$MI, swing$weight[swing$id %in% aracne$id], xlab = "ARACNE", ylab = "SWING")
plot(aracne$MI, dyng3$weight[dyng3$id %in% aracne$id], XLAB = "ARACNE", ylab = "DYN3")
plot(swing$weight[swing$id %in% aracne$id], dyng3$weight[dyng3$id %in% aracne$id], xlab = "SWING", ylab = "DYNG3")





































##
