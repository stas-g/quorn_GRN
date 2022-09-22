#Part1 
#Run DESeq2 with C/N count files reps uncollapsed 
library(DESeq2)
library(tximport)
setwd("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/deseq2")

#Create tx2gene table 
gene.t <- paste0("g", 1:13809, ".t", 1)
gene.t <- as.data.frame(gene.t)
gene <- paste0("g", 1:13809)
gene <- as.data.frame(gene)
tx2gene <- cbind(gene.t, gene)
colnames(tx2gene) <- c("Name", "Name")

#Import data
txi.reps <- tximport(paste(list.dirs("./", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)
mysamples <- list.dirs("./",full.names=F,recursive=F)
txi.genes <- summarizeToGene(txi.reps,tx2gene)
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

#Create experiment design table 
designtable <- mysamples
designtable <- as.data.frame(designtable)
designtable$timepoint <- sapply(strsplit(designtable$designtable, split = "FvC[0,1,2,3,4]"),'[',2)
designtable$timepoint <- sapply(strsplit(designtable$timepoint, split = "_"),'[',1)
designtable$rep <- sapply(strsplit(designtable$designtable, split = "_"),'[',2)
designtable$ID <- designtable[,1]

condition <- function(timepoint){
  if(timepoint == "T0")
    return ("control")
   else 
    return ("experimental")
}

designtable$condition <- lapply(X = designtable$timepoint, FUN = condition)
designtable$ID[designtable$condition == "control"] <- "control"

#Set up experiment 
colData <- designtable
colData$indrep <- paste0(colData$ID)

#Define the GLM parameters 
design <- ~ indrep
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

#Set rowsums threshold 
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

#Library normalisation
dds <- estimateSizeFactors(dds)

#Run Deseq
dds <- DESeq(dds)
resultsNames(dds)

#Calculate variance stabilizing transformation from the fitted dispersion-mean relation, for machine learning 
vst <- varianceStabilizingTransformation(dds)
vst <- vst(dds,blind=FALSE)
write.table(assay(vst), file="vst.txt", sep = '\t', quote = FALSE)


#Part2 Run dynGENIE reverse network reconstruction 
library(GENIE3)
library (reshape2)
library (doRNG)
library (doParallel)
source ("dynGENIE3.R")

all_data <- read.table("vst.txt", header = T)

#Condition1 
condition1 <- all_data[6:40]
C1TP1 <- condition1[1:5]
C1TP2 <- condition1[6:10]
C1TP3 <- condition1[11:15]
C1TP4 <- condition1[16:20]
C1TP5 <- condition1[21:25]
C1TP6 <- condition1[26:30]
C1TP7 <- condition1[31:35]

C1TP1$FvC1T1_1 <- C1TP1$FvC1T1_1 + 0.000001
C1TP2$FvC1T2_1 <- C1TP2$FvC1T2_1 + 0.000001
C1TP3$FvC1T3_1 <- C1TP3$FvC1T3_1 + 0.000001
C1TP4$FvC1T4_1 <- C1TP4$FvC1T4_1 + 0.000001
C1TP5$FvC1T5_1 <- C1TP5$FvC1T5_1 + 0.000001
C1TP6$FvC1T6_1 <- C1TP6$FvC1T6_1 + 0.000001
C1TP7$FvC1T7_1 <- C1TP7$FvC1T7_1 + 0.000001

colnames(C1TP1) <- c("1", "2", "3", "4", "5")
colnames(C1TP2) <- c("1", "2", "3", "4", "5")
colnames(C1TP3) <- c("1", "2", "3", "4", "5")
colnames(C1TP4) <- c("1", "2", "3", "4", "5")
colnames(C1TP5) <- c("1", "2", "3", "4", "5")
colnames(C1TP6) <- c("1", "2", "3", "4", "5")
colnames(C1TP7) <- c("1", "2", "3", "4", "5")

C1TP1 <- t(C1TP1)
C1TP2 <- t(C1TP2)
C1TP3 <- t(C1TP3)
C1TP4 <- t(C1TP4)
C1TP5 <- t(C1TP5)
C1TP6 <- t(C1TP6)
C1TP7 <- t(C1TP7)

C1TP1 <- data.frame("time_points"=rownames(C1TP1), C1TP1)
C1TP2 <- data.frame("time_points"=rownames(C1TP2), C1TP2)
C1TP3 <- data.frame("time_points"=rownames(C1TP3), C1TP3)
C1TP4 <- data.frame("time_points"=rownames(C1TP4), C1TP4)
C1TP5 <- data.frame("time_points"=rownames(C1TP5), C1TP5)
C1TP6 <- data.frame("time_points"=rownames(C1TP6), C1TP6)
C1TP7 <- data.frame("time_points"=rownames(C1TP7), C1TP7)

dir.create("../VST_C1")
write.table(C1TP1, file = "../C1TP1.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C1TP2, file = "../C1TP2.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C1TP3, file = "../C1TP3.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C1TP4, file = "../C1TP4.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C1TP5, file = "../C1TP5.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C1TP6, file = "../C1TP6.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C1TP7, file = "../C1TP7.txt", sep = '\t', quote = FALSE, col.names = NA) 

TS1 <- read.expr.matrix("C1TP1.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C1TP2.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C1TP3.txt",form="rows.are.samples")                            
TS4 <- read.expr.matrix("C1TP4.txt",form="rows.are.samples")
TS5 <- read.expr.matrix("C1TP5.txt",form="rows.are.samples")
TS6 <- read.expr.matrix("C1TP6.txt",form="rows.are.samples")
TS7 <- read.expr.matrix("C1TP7.txt",form="rows.are.samples")                            


#setwd("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C1")
#C1 <- read.table("C1TP1.txt", header = T)
#C2 <- read.table("C1TP2.txt", header = T)
#C3 <- read.table("C1TP3.txt", header = T)
#C4 <- read.table("C1TP4.txt", header = T)
#C5 <- read.table("C1TP5.txt", header = T)
#C6 <- read.table("C1TP6.txt", header = T)
#C7 <- read.table("C1TP7.txt", header = T)


#HG_list <- list(C1, C2, C3, C4, C5, C6, C7)
#merged_data <- Reduce(function(x, y) merge(x, y, all=TRUE), HG_list)
#write.table(merged_data, file = "merged_tp1_con1.txt", sep = '\t', quote = FALSE, row.names = FALSE)


time.points <- list (TS1[1,], TS2[1,], TS3[1,], TS4[1,],TS5[1,], TS6[1,], TS7[1,])
TS.data <- list (TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),], TS6[2:nrow(TS6),], TS7[2:nrow(TS7),]) 
regs <- read.table("WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")
reg <- regs[,1]
reg <- as.vector(reg)
res <- dynGENIE3 (TS.data,time.points, regulators=reg)
link.list <- get.link.list(res$weight.matrix, threshold=0.00)
write.table(link.list, "Time_series_con1_regs.txt", sep="\t")

#Condition2 
condition2 <- all_data[41:75]
C2TP1 <- condition2[1:5]
C2TP2 <- condition2[6:10]
C2TP3 <- condition2[11:15]
C2TP4 <- condition2[16:20]
C2TP5 <- condition2[21:25]
C2TP6 <- condition2[26:30]
C2TP7 <- condition2[31:35]

C2TP1$FvC2T1_1 <- C2TP1$FvC2T1_1 + 0.000001
C2TP2$FvC2T2_1 <- C2TP2$FvC2T2_1 + 0.000001
C2TP3$FvC2T3_1 <- C2TP3$FvC2T3_1 + 0.000001
C2TP4$FvC2T4_1 <- C2TP4$FvC2T4_1 + 0.000001
C2TP5$FvC2T5_1 <- C2TP5$FvC2T5_1 + 0.000001
C2TP6$FvC2T6_1 <- C2TP6$FvC2T6_1 + 0.000001
C2TP7$FvC2T7_1 <- C2TP7$FvC2T7_1 + 0.000001

colnames(C2TP1) <- c("1", "2", "3", "4", "5")
colnames(C2TP2) <- c("1", "2", "3", "4", "5")
colnames(C2TP3) <- c("1", "2", "3", "4", "5")
colnames(C2TP4) <- c("1", "2", "3", "4", "5")
colnames(C2TP5) <- c("1", "2", "3", "4", "5")
colnames(C2TP6) <- c("1", "2", "3", "4", "5")
colnames(C2TP7) <- c("1", "2", "3", "4", "5")

C2TP1 <- t(C2TP1)
C2TP2 <- t(C2TP2)
C2TP3 <- t(C2TP3)
C2TP4 <- t(C2TP4)
C2TP5 <- t(C2TP5)
C2TP6 <- t(C2TP6)
C2TP7 <- t(C2TP7)

C2TP1 <- data.frame("time_points"=rownames(C2TP1), C2TP1)
C2TP2 <- data.frame("time_points"=rownames(C2TP2), C2TP2)
C2TP3 <- data.frame("time_points"=rownames(C2TP3), C2TP3)
C2TP4 <- data.frame("time_points"=rownames(C2TP4), C2TP4)
C2TP5 <- data.frame("time_points"=rownames(C2TP5), C2TP5)
C2TP6 <- data.frame("time_points"=rownames(C2TP6), C2TP6)
C2TP7 <- data.frame("time_points"=rownames(C2TP7), C2TP7)

dir.create("../VST_C2")
write.table(C2TP1, file = "../C2TP1.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C2TP2, file = "../C2TP2.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C2TP3, file = "../C2TP3.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C2TP4, file = "../C2TP4.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C2TP5, file = "../C2TP5.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C2TP6, file = "../C2TP6.txt", sep = '\t', quote = FALSE, col.names = NA)
write.table(C2TP7, file = "../C2TP7.txt", sep = '\t', quote = FALSE, col.names = NA)


#setwd("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C2")
#C1 <- read.table("C2TP1.txt", header = T)
#C2 <- read.table("C2TP2.txt", header = T)
#C3 <- read.table("C2TP3.txt", header = T)
#C4 <- read.table("C2TP4.txt", header = T)
#C5 <- read.table("C2TP5.txt", header = T)
#C6 <- read.table("C2TP6.txt", header = T)
#C7 <- read.table("C2TP7.txt", header = T)


#HG_list <- list(C1, C2, C3, C4, C5, C6, C7)
#merged_data <- Reduce(function(x, y) merge(x, y, all=TRUE), HG_list)
#write.table(merged_data, file = "merged_tp1_con1.txt", sep = '\t', quote = FALSE, row.names = FALSE)


TS1 <- read.expr.matrix("C2TP1.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C2TP2.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C2TP3.txt",form="rows.are.samples")                            
TS4 <- read.expr.matrix("C2TP4.txt",form="rows.are.samples")
TS5 <- read.expr.matrix("C2TP5.txt",form="rows.are.samples")
TS6 <- read.expr.matrix("C2TP6.txt",form="rows.are.samples")
TS7 <- read.expr.matrix("C2TP7.txt",form="rows.are.samples")

time.points <- list (TS1[1,], TS2[1,], TS3[1,], TS4[1,],TS5[1,], TS6[1,], TS7[1,])
TS.data <- list (TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),], TS6[2:nrow(TS6),], TS7[2:nrow(TS7),]) 
regs <- read.table("WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")
reg <- regs[,1]
reg <- as.vector(reg)
res <- dynGENIE3 (TS.data,time.points, regulators=reg)
link.list <- get.link.list(res$weight.matrix, threshold=0.00)
write.table(link.list, "Time_seriescon2_th0.txt", sep="\t")

#Condition3
condition3 <- all_data[76:110]
C3TP1 <- condition3[1:5]
C3TP2 <- condition3[6:10]
C3TP3 <- condition3[11:15]
C3TP4 <- condition3[16:20]
C3TP5 <- condition3[21:25]
C3TP6 <- condition3[26:30]
C3TP7 <- condition3[31:35]

C3TP1$FvC3T1_1 <- C3TP1$FvC3T1_1 + 0.000001
C3TP2$FvC3T2_1 <- C3TP2$FvC3T2_1 + 0.000001
C3TP3$FvC3T3_1 <- C3TP3$FvC3T3_1 + 0.000001
C3TP4$FvC3T4_1 <- C3TP4$FvC3T4_1 + 0.000001
C3TP5$FvC3T5_1 <- C3TP5$FvC3T5_1 + 0.000001
C3TP6$FvC3T6_1 <- C3TP6$FvC3T6_1 + 0.000001
C3TP7$FvC3T7_1 <- C3TP7$FvC3T7_1 + 0.000001

colnames(C3TP1) <- c("1", "2", "3", "4", "5")
colnames(C3TP2) <- c("1", "2", "3", "4", "5")
colnames(C3TP3) <- c("1", "2", "3", "4", "5")
colnames(C3TP4) <- c("1", "2", "3", "4", "5")
colnames(C3TP5) <- c("1", "2", "3", "4", "5")
colnames(C3TP6) <- c("1", "2", "3", "4", "5")
colnames(C3TP7) <- c("1", "2", "3", "4", "5")

C3TP1 <- t(C3TP1)
C3TP2 <- t(C3TP2)
C3TP3 <- t(C3TP3)
C3TP4 <- t(C3TP4)
C3TP5 <- t(C3TP5)
C3TP6 <- t(C3TP6)
C3TP7 <- t(C3TP7)

C3TP1 <- data.frame("time_points"=rownames(C3TP1), C3TP1)
C3TP2 <- data.frame("time_points"=rownames(C3TP2), C3TP2)
C3TP3 <- data.frame("time_points"=rownames(C3TP3), C3TP3)
C3TP4 <- data.frame("time_points"=rownames(C3TP4), C3TP4)
C3TP5 <- data.frame("time_points"=rownames(C3TP5), C3TP5)
C3TP6 <- data.frame("time_points"=rownames(C3TP6), C3TP6)
C3TP7 <- data.frame("time_points"=rownames(C3TP7), C3TP7)

dir.create("../VST_C3")
write.table(C3TP1, file = "../VST_C3/C3TP1.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C3TP2, file = "../VST_C3/C3TP2.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C3TP3, file = "../VST_C3/C3TP3.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C3TP4, file = "../VST_C3/C3TP4.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C3TP5, file = "../VST_C3/C3TP5.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C3TP6, file = "../VST_C3/C3TP6.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C3TP7, file = "../VST_C3/C3TP7.txt", sep = '\t', quote = FALSE, row.names = FALSE)

TS1 <- read.expr.matrix("C3TP1.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C3TP2.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C3TP3.txt",form="rows.are.samples")                            
TS4 <- read.expr.matrix("C3TP4.txt",form="rows.are.samples")
TS5 <- read.expr.matrix("C3TP5.txt",form="rows.are.samples")
TS6 <- read.expr.matrix("C3TP6.txt",form="rows.are.samples")
TS7 <- read.expr.matrix("C3TP7.txt",form="rows.are.samples")

time.points <- list (TS1[1,], TS2[1,], TS3[1,], TS4[1,],TS5[1,], TS6[1,], TS7[1,])
TS.data <- list (TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),], TS6[2:nrow(TS6),], TS7[2:nrow(TS7),]) 
regs <- read.table("WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")
reg <- regs[,1]
reg <- as.vector(reg)
res <- dynGENIE3 (TS.data,time.points, regulators=reg)
link.list <- get.link.list(res$weight.matrix, threshold=0.00)
write.table(link.list, "Time_seriescon3_th0.txt", sep="\t")


#Condition4 
condition4 <- all_data[111:145]
C4TP1 <- condition4[1:5]
C4TP2 <- condition4[6:10]
C4TP3 <- condition4[11:15]
C4TP4 <- condition4[16:20]
C4TP5 <- condition4[21:25]
C4TP6 <- condition4[26:30]
C4TP7 <- condition4[31:35]

C4TP1$FvC4T1_1 <- C4TP1$FvC4T1_1 + 0.000001
C4TP2$FvC4T2_1 <- C4TP2$FvC4T2_1 + 0.000001
C4TP3$FvC4T3_1 <- C4TP3$FvC4T3_1 + 0.000001
C4TP4$FvC4T4_1 <- C4TP4$FvC4T4_1 + 0.000001
C4TP5$FvC4T5_1 <- C4TP5$FvC4T5_1 + 0.000001
C4TP6$FvC4T6_1 <- C4TP6$FvC4T6_1 + 0.000001
C4TP7$FvC4T7_1 <- C4TP7$FvC4T7_1 + 0.000001

colnames(C4TP1) <- c("1", "2", "3", "4", "5")
colnames(C4TP2) <- c("1", "2", "3", "4", "5")
colnames(C4TP3) <- c("1", "2", "3", "4", "5")
colnames(C4TP4) <- c("1", "2", "3", "4", "5")
colnames(C4TP5) <- c("1", "2", "3", "4", "5")
colnames(C4TP6) <- c("1", "2", "3", "4", "5")
colnames(C4TP7) <- c("1", "2", "3", "4", "5")

C4TP1 <- t(C4TP1)
C4TP2 <- t(C4TP2)
C4TP3 <- t(C4TP3)
C4TP4 <- t(C4TP4)
C4TP5 <- t(C4TP5)
C4TP6 <- t(C4TP6)
C4TP7 <- t(C4TP7)

C4TP1 <- data.frame("time_points"=rownames(C4TP1), C4TP1)
C4TP2 <- data.frame("time_points"=rownames(C4TP2), C4TP2)
C4TP3 <- data.frame("time_points"=rownames(C4TP3), C4TP3)
C4TP4 <- data.frame("time_points"=rownames(C4TP4), C4TP4)
C4TP5 <- data.frame("time_points"=rownames(C4TP5), C4TP5)
C4TP6 <- data.frame("time_points"=rownames(C4TP6), C4TP6)
C4TP7 <- data.frame("time_points"=rownames(C4TP7), C4TP7)

dir.create("../VST_C4")
write.table(C4TP1, file = "../VST_C4/C4TP1.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C4TP2, file = "../VST_C4/C4TP2.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C4TP3, file = "../VST_C4/C4TP3.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C4TP4, file = "../VST_C4/C4TP4.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C4TP5, file = "../VST_C4/C4TP5.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C4TP6, file = "../VST_C4/C4TP6.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(C4TP7, file = "../VST_C4/C4TP7.txt", sep = '\t', quote = FALSE, row.names = FALSE)

TS1 <- read.expr.matrix("C4TP1.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("C4TP2.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("C4TP3.txt",form="rows.are.samples")                            
TS4 <- read.expr.matrix("C4TP4.txt",form="rows.are.samples")
TS5 <- read.expr.matrix("C4TP5.txt",form="rows.are.samples")
TS6 <- read.expr.matrix("C4TP6.txt",form="rows.are.samples")
TS7 <- read.expr.matrix("C4TP7.txt",form="rows.are.samples")

time.points <- list (TS1[1,], TS2[1,], TS3[1,], TS4[1,],TS5[1,], TS6[1,], TS7[1,])
TS.data <- list (TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),], TS4[2:nrow(TS4),], TS5[2:nrow(TS5),], TS6[2:nrow(TS6),], TS7[2:nrow(TS7),]) 
regs <- read.table("WT_minion_TF_expression_only_head.txt", header = TRUE, sep = "\t")
reg <- regs[,1]
reg <- as.vector(reg)
res <- dynGENIE3 (TS.data,time.points, regulators=reg)
link.list <- get.link.list(res$weight.matrix, threshold=0.00)
write.table(link.list, "Time_seriescon4_th0.txt", sep="\t")


#Part3 run aracne reverse network reconstruction 
#Condition1
aracne_data_C1_t1_7_uncollapsed <- condition1
write.table(aracne_data_C1_t1_7_uncollapsed, file = "../aracne_data_C1_t1_7_uncollapsed.txt", quote = FALSE, sep = '\t', col.names = NA)
#Run aracne Condition1
#Gerate bootstrap values
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C1_t1_7_uncollapsed.txt \
-o C1.1_results \
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
#Filter results to take top scores in relation to dyngenie results 
#Filter aracne condition1 
filter <- read.table("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C1/Time_series1_7_con1_arac.txt", header = T)
ordered <- filter[order(filter$pvalue),]
fo <- ordered[ordered[,4]>0,]
araccompare <- tail(fo, n =1042)
araccompare <- araccompare[,c('Regulator', 'Target', 'pvalue')]
colnames(araccompare) <- c('regulatory.gene', 'target.gene', 'weight')

#Condition2
aracne_data_C2_t1_7_uncollapsed <- condition2
write.table(aracne_data_C2_t1_7_uncollapsed, file = "../aracne_data_C2_t1_7_uncollapsed.txt", quote = FALSE, sep = '\t', col.names = NA)
#Run aracne Condition2
#Gerate bootstrap values
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C2_t1_7_uncollapsed.txt \
-o C2.1_results \
--pvalue 0 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C2_t1_7_uncollapsed.txt \
-o C2.1_results \
--pvalue 0 \
-t WT_minion_TF_expression_only_head.txt \
--seed $i 
done
# Consolidate, i.e. combine the bootstraps into a final network file

#Filter aracne condition2
filter <- read.table("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C2/Time_series1_7_con2_arac.txt", header = T)
ordered <- filter[order(filter$pvalue),]
fo <- ordered[ordered[,4]>0,]
araccompare <- tail(fo, n =2200)
araccompare <- araccompare[,c('Regulator', 'Target', 'pvalue')]
colnames(araccompare) <- c('regulatory.gene', 'target.gene', 'weight')
write.table(araccompare, file = "C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C2/Time_series1_7_con2_arac_edited.txt", quote = FALSE, row.names = FALSE, col.names = FALSE )

#Condition3
aracne_data_C3_t1_7_uncollapsed <- condition3
write.table(aracne_data_C3_t1_7_uncollapsed, file = "../VST_C3/aracne_data_C3_t1_7_uncollapsed.txt", quote = FALSE, sep = '\t', col.names = NA)
#Run aracne Condition2
#Gerate bootstrap values
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C3_t1_7_uncollapsed.txt \
-o C3.1_results \
--pvalue 0 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C3_t1_7_uncollapsed.txt \
-o C3.1_results \
--pvalue 0 \
-t WT_minion_TF_expression_only_head.txt \
--seed $i 
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar /home/jconnell/projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o C3.1_results \
--consolidate
#Filter aracne condition3
filter <- read.table("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C3/Time_series_1_7_con3_arac.txt", header = T)
ordered <- filter[order(filter$pvalue),]
fo <- ordered[ordered[,4]>0,]
araccompare <- tail(fo, n =10000)
araccompare <- araccompare[,c('Regulator', 'Target', 'pvalue')]
colnames(araccompare) <- c('regulatory.gene', 'target.gene', 'weight')
write.table(araccompare, file = "C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C3/Time_series1_7_con3_arac_edited.txt", quote = FALSE, row.names = FALSE, col.names = FALSE )

#Condition4
aracne_data_C4_t1_7_uncollapsed <- condition4
write.table(aracne_data_C4_t1_7_uncollapsed, file = "../VST_C4/aracne_data_C4_t1_7_uncollapsed.txt", quote = FALSE, sep = '\t', col.names = NA)
#Run aracne Condition2
#Gerate bootstrap values
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C4_t1_7_uncollapsed.txt \
-o C4.1_results \
--pvalue 0 \
--seed 1 \
--calculateThreshold
#Run aracne on bootstraps of the input matrix
for i in {1..100}
do
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-e projects/niab/gene_regulatory_network/ARACNe-AP/aracne_data_C4_t1_7_uncollapsed.txt \
-o C4.1_results \
--pvalue 0 \
-t WT_minion_TF_expression_only_head.txt \
--seed $i 
done
# Consolidate, i.e. combine the bootstraps into a final network file
java -Xmx5G -jar projects/niab/gene_regulatory_network/ARACNe-AP/dist/aracne.jar \
-o C4.1_results \
--consolidate
#Filter aracne condition4
filter <- read.table("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C4/Time_series_1_7_con4_arac.txt", header = T)
ordered <- filter[order(filter$pvalue),]
fo <- ordered[ordered[,4]>0,]
araccompare <- tail(fo, n =5000)
araccompare <- araccompare[,c('Regulator', 'Target', 'pvalue')]
colnames(araccompare) <- c('regulatory.gene', 'target.gene', 'weight')
write.table(araccompare, file = "C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C4/Time_series1_7_con4_arac_edited.txt", quote = FALSE, row.names = FALSE, col.names = FALSE )

#Part4 create data for SWING
all_data <- read.table("vst.txt", header = T)

#Condition1 
#condition1 <- all_data[6:40]
#C1TP1 <- condition1[1:5]
#C1TP2 <- condition1[6:10]
#C1TP3 <- condition1[11:15]
#C1TP4 <- condition1[16:20]
#C1TP5 <- condition1[21:25]
#C1TP6 <- condition1[26:30]
#C1TP7 <- condition1[31:35]

#C1TP1$FvC1T1_1 <- C1TP1$FvC1T1_1 + 0.000001
#C1TP2$FvC1T2_1 <- C1TP2$FvC1T2_1 + 0.000001
#C1TP3$FvC1T3_1 <- C1TP3$FvC1T3_1 + 0.000001
#C1TP4$FvC1T4_1 <- C1TP4$FvC1T4_1 + 0.000001
#C1TP5$FvC1T5_1 <- C1TP5$FvC1T5_1 + 0.000001
#C1TP6$FvC1T6_1 <- C1TP6$FvC1T6_1 + 0.000001
#C1TP7$FvC1T7_1 <- C1TP7$FvC1T7_1 + 0.000001

#C1TP1$rows <- rowMeans(C1TP1)
#C1TP2$rows <- rowMeans(C1TP2)
#C1TP3$rows <- rowMeans(C1TP3)
#C1TP4$rows <- rowMeans(C1TP4)
#C1TP5$rows <- rowMeans(C1TP5)
#C1TP6$rows <- rowMeans(C1TP6)
#C1TP7$rows <- rowMeans(C1TP7)

#C1TA <- data.frame(C1TP1$rows, C1TP2$rows, C1TP3$rows, C1TP4$rows, C1TP5$rows, C1TP6$rows, C1TP7$rows)
#row.names(C1TA) <- row.names(C1TP1)
#colnames(C1TA) <- c("1", "2", "3", "4", "5", "6", "7")
#C1TA <- t(C1TA)
#C1TA <- data.frame("Time"=rownames(C1TA), C1TA)
#dir.create("../SWING")
#write.table(C1TA, file = "../SWING/C1TA.txt", sep = '\t', row.names = FALSE, quote = FALSE)
library ("berryFunctions")

condition1 <- all_data[6:40]
C1TP1 <- condition1[1:5]
C1TP2 <- condition1[6:10]
C1TP3 <- condition1[11:15]
C1TP4 <- condition1[16:20]
C1TP5 <- condition1[21:25]
C1TP6 <- condition1[26:30]
C1TP7 <- condition1[31:35]

C1TP1$FvC1T1_1 <- C1TP1$FvC1T1_1 + 0.000001
C1TP2$FvC1T2_1 <- C1TP2$FvC1T2_1 + 0.000001
C1TP3$FvC1T3_1 <- C1TP3$FvC1T3_1 + 0.000001
C1TP4$FvC1T4_1 <- C1TP4$FvC1T4_1 + 0.000001
C1TP5$FvC1T5_1 <- C1TP5$FvC1T5_1 + 0.000001
C1TP6$FvC1T6_1 <- C1TP6$FvC1T6_1 + 0.000001
C1TP7$FvC1T7_1 <- C1TP7$FvC1T7_1 + 0.000001

colnames(C1TP1) <- c("1", "2", "3", "4", "5")
colnames(C1TP2) <- c("1", "2", "3", "4", "5")
colnames(C1TP3) <- c("1", "2", "3", "4", "5")
colnames(C1TP4) <- c("1", "2", "3", "4", "5")
colnames(C1TP5) <- c("1", "2", "3", "4", "5")
colnames(C1TP6) <- c("1", "2", "3", "4", "5")
colnames(C1TP7) <- c("1", "2", "3", "4", "5")

C1TP1 <- t(C1TP1)
C1TP2 <- t(C1TP2)
C1TP3 <- t(C1TP3)
C1TP4 <- t(C1TP4)
C1TP5 <- t(C1TP5)
C1TP6 <- t(C1TP6)
C1TP7 <- t(C1TP7)

bindedC1 <- rbind(C1TP1, C1TP2, C1TP3, C1TP4, C1TP5, C1TP6, C1TP7) 
bindedC1 <- as.data.frame(bindedC1)
Time <- rep(c(1,2,3,4,5), times = 7)
Time <- as.data.frame(Time)
bindedC1 <- cbind(bindedC1, Time)
bindedC1 <- bindedC1[,c(12805, 1:12804)]
bindedC1 <- insertRows(bindedC1, c(6,12,18,24,30,36) , new = "")
write.table(bindedC1, file = "../SWING/C1_swing.txt", sep = '\t', quote = FALSE, row.names = FALSE)

#Try with rowmeans 

C1TP1$r <- rowMeans(C1TP1)
C1TP2$r <- rowMeans(C1TP2)
C1TP3$r <- rowMeans(C1TP3)
C1TP4$r <- rowMeans(C1TP4)
C1TP5$r <- rowMeans(C1TP5)
C1TP6$r <- rowMeans(C1TP6)
C1TP7$r <- rowMeans(C1TP7)

bindedC1 <- cbind(C1TP1$r, C1TP2$r, C1TP3$r, C1TP4$r, C1TP5$r, C1TP6$r, C1TP7$r) 
row.names(bindedC1) <- rownames(C1TP1)
colnames(bindedC1) <- c("1", "2", "3", "4", "5", "6", "7")
bindedC1 <- as.data.frame(bindedC1)
bindedC1$`1`<- bindedC1$`1` + 0.000001
bindedC1 <- t(bindedC1)
bindedC1 <- data.frame("Time" =  rownames(bindedC1), bindedC1)
write.table(bindedC1, file = "../SWING/bindedc1.txt", quote = FALSE, row.names = FALSE, sep = '\t')


#SWING result C1

res1 <- read.table("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/SWING/Results/SWING_link_list.txt", header = TRUE)
ordered <- res1[order(res1$weight),]
fo <- ordered[ordered[,3]>0,]
swingcompare <- head(fo, n =10000)

#Condition2
condition2 <- all_data[41:75]
C2TP1 <- condition2[1:5]
C2TP2 <- condition2[6:10]
C2TP3 <- condition2[11:15]
C2TP4 <- condition2[16:20]
C2TP5 <- condition2[21:25]
C2TP6 <- condition2[26:30]
C2TP7 <- condition2[31:35]

C2TP1$FvC2T1_1 <- C2TP1$FvC2T1_1 + 0.000001
C2TP2$FvC2T2_1 <- C2TP2$FvC2T2_1 + 0.000001
C2TP3$FvC2T3_1 <- C2TP3$FvC2T3_1 + 0.000001
C2TP4$FvC2T4_1 <- C2TP4$FvC2T4_1 + 0.000001
C2TP5$FvC2T5_1 <- C2TP5$FvC2T5_1 + 0.000001
C2TP6$FvC2T6_1 <- C2TP6$FvC2T6_1 + 0.000001
C2TP7$FvC2T7_1 <- C2TP7$FvC2T7_1 + 0.000001

colnames(C2TP1) <- c("1", "2", "3", "4", "5")
colnames(C2TP2) <- c("1", "2", "3", "4", "5")
colnames(C2TP3) <- c("1", "2", "3", "4", "5")
colnames(C2TP4) <- c("1", "2", "3", "4", "5")
colnames(C2TP5) <- c("1", "2", "3", "4", "5")
colnames(C2TP6) <- c("1", "2", "3", "4", "5")
colnames(C2TP7) <- c("1", "2", "3", "4", "5")

C2TP1 <- t(C2TP1)
C2TP2 <- t(C2TP2)
C2TP3 <- t(C2TP3)
C2TP4 <- t(C2TP4)
C2TP5 <- t(C2TP5)
C2TP6 <- t(C2TP6)
C2TP7 <- t(C2TP7)

bindedC2 <- rbind(C2TP1, C2TP2, C2TP3, C2TP4, C2TP5, C2TP6, C2TP7) 
bindedC2 <- as.data.frame(bindedC2)
Time <- rep(c(1,2,3,4,5), times = 7)
Time <- as.data.frame(Time)
bindedC2 <- cbind(bindedC2, Time)
bindedC2 <- bindedC2[,c(12805, 1:12804)]
bindedC2 <- insertRows(bindedC2, c(6,12,18,24,30,36) , new = "")
write.table(bindedC2, file = "../SWING/C2_swing.txt", sep = '\t', quote = FALSE, row.names = FALSE)

#Condition3
condition3 <- all_data[76:110]
C3TP1 <- condition3[1:5]
C3TP2 <- condition3[6:10]
C3TP3 <- condition3[11:15]
C3TP4 <- condition3[16:20]
C3TP5 <- condition3[21:25]
C3TP6 <- condition3[26:30]
C3TP7 <- condition3[31:35]

C3TP1$FvC3T1_1 <- C3TP1$FvC3T1_1 + 0.000001
C3TP2$FvC3T2_1 <- C3TP2$FvC3T2_1 + 0.000001
C3TP3$FvC3T3_1 <- C3TP3$FvC3T3_1 + 0.000001
C3TP4$FvC3T4_1 <- C3TP4$FvC3T4_1 + 0.000001
C3TP5$FvC3T5_1 <- C3TP5$FvC3T5_1 + 0.000001
C3TP6$FvC3T6_1 <- C3TP6$FvC3T6_1 + 0.000001
C3TP7$FvC3T7_1 <- C3TP7$FvC3T7_1 + 0.000001

colnames(C3TP1) <- c("1", "2", "3", "4", "5")
colnames(C3TP2) <- c("1", "2", "3", "4", "5")
colnames(C3TP3) <- c("1", "2", "3", "4", "5")
colnames(C3TP4) <- c("1", "2", "3", "4", "5")
colnames(C3TP5) <- c("1", "2", "3", "4", "5")
colnames(C3TP6) <- c("1", "2", "3", "4", "5")
colnames(C3TP7) <- c("1", "2", "3", "4", "5")

C3TP1 <- t(C3TP1)
C3TP2 <- t(C3TP2)
C3TP3 <- t(C3TP3)
C3TP4 <- t(C3TP4)
C3TP5 <- t(C3TP5)
C3TP6 <- t(C3TP6)
C3TP7 <- t(C3TP7)

bindedC3 <- rbind(C3TP1, C3TP2, C3TP3, C3TP4, C3TP5, C3TP6, C3TP7) 
bindedC3 <- as.data.frame(bindedC3)
Time <- rep(c(1,2,3,4,5), times = 7)
Time <- as.data.frame(Time)
bindedC3 <- cbind(bindedC3, Time)
bindedC3 <- bindedC3[,c(12805, 1:12804)]
bindedC3 <- insertRows(bindedC3, c(6,12,18,24,30,36) , new = "")
write.table(bindedC3, file = "../SWING/C3_swing.txt", sep = '\t', quote = FALSE, row.names = FALSE)                

#Condition4
condition4 <- all_data[111:145]
C4TP1 <- condition4[1:5]
C4TP2 <- condition4[6:10]
C4TP3 <- condition4[11:15]
C4TP4 <- condition4[16:20]
C4TP5 <- condition4[21:25]
C4TP6 <- condition4[26:30]
C4TP7 <- condition4[31:35]

C4TP1$FvC4T1_1 <- C4TP1$FvC4T1_1 + 0.000001
C4TP2$FvC4T2_1 <- C4TP2$FvC4T2_1 + 0.000001
C4TP3$FvC4T3_1 <- C4TP3$FvC4T3_1 + 0.000001
C4TP4$FvC4T4_1 <- C4TP4$FvC4T4_1 + 0.000001
C4TP5$FvC4T5_1 <- C4TP5$FvC4T5_1 + 0.000001
C4TP6$FvC4T6_1 <- C4TP6$FvC4T6_1 + 0.000001
C4TP7$FvC4T7_1 <- C4TP7$FvC4T7_1 + 0.000001

colnames(C4TP1) <- c("1", "2", "3", "4", "5")
colnames(C4TP2) <- c("1", "2", "3", "4", "5")
colnames(C4TP3) <- c("1", "2", "3", "4", "5")
colnames(C4TP4) <- c("1", "2", "3", "4", "5")
colnames(C4TP5) <- c("1", "2", "3", "4", "5")
colnames(C4TP6) <- c("1", "2", "3", "4", "5")
colnames(C4TP7) <- c("1", "2", "3", "4", "5")

C4TP1 <- t(C4TP1)
C4TP2 <- t(C4TP2)
C4TP3 <- t(C4TP3)
C4TP4 <- t(C4TP4)
C4TP5 <- t(C4TP5)
C4TP6 <- t(C4TP6)
C4TP7 <- t(C4TP7)

bindedC4 <- rbind(C4TP1, C4TP2, C4TP3, C4TP4, C4TP5, C4TP6, C4TP7) 
bindedC4 <- as.data.frame(bindedC4)
Time <- rep(c(1,2,3,4,5), times = 7)
Time <- as.data.frame(Time)
bindedC4 <- cbind(bindedC4, Time)
bindedC4 <- bindedC4[,c(12805, 1:12804)]
bindedC4 <- insertRows(bindedC4, c(6,12,18,24,30,36) , new = "")
write.table(bindedC4, file = "../SWING/C4_swing.txt", sep = '\t', quote = FALSE, row.names = FALSE)


#Part5
#Plot linklist in cytoscape
library(igraph)
library(RCy3)

input_1 <- "C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C1/filtered_con1_dyn_unstrict_2500.txt"
data <- read.table(input_1)
edge_listsi <- data
Gsi <- graph.data.frame(edge_listsi,directed = F)
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)
cw = createNetworkFromGraph(graph=g.cyto)

input_2 <- "C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C2/filtered_con2_dyn_unstrict_5000.txt"
data <- read.table(input_2)
edge_listsi <- data
Gsi <- graph.data.frame(edge_listsi,directed = F)
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)
cw = createNetworkFromGraph(graph=g.cyto)

input_3 <- "C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C3/Time_series3_7_con1_dyn.txt"
data <- read.table(input_3)
edge_listsi <- data
Gsi <- graph.data.frame(edge_listsi,directed = F)
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)
cw = createNetworkFromGraph(graph=g.cyto)

input_4 <- "C:/Users/john.connell/Documents/Bioinformatics_2022/projects/GRN/VST_C4/Time_series4_7_con1_dyn.txt"
data <- read.table(input_4)
edge_listsi <- data
Gsi <- graph.data.frame(edge_listsi,directed = F)
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)
cw = createNetworkFromGraph(graph=g.cyto)


edge_listsi <- araccompare
Gsi <- graph.data.frame(edge_listsi,directed = F)
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)
cw = createNetworkFromGraph(graph=g.cyto)

#Plot swing data

edge_listsi <- swingcompare
Gsi <- graph.data.frame(edge_listsi,directed = F)
Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
g_arasi <- graph.adjacency(Asi,mode = "undirected",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)
cw = createNetworkFromGraph(graph=g.cyto)


 