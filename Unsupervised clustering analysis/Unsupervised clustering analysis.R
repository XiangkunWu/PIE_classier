###Unsupervised clustering analysis by CNMF algorithm of “CancerSubtype” R package
library(stringr)
library(CancerSubtypes)
library(Biobase)
# set working path
setwd("") 

# read Omics_data1 #colunm: sample; row:genes
Omics_data1 <- read.table("Omics_data1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "",row.names = 1)
Omics_data1[1:3, 1:3]

# read Omics_data2 #colunm: sample; row:genes
Omics_data2 <- read.table("Omics_data2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "",row.names = 1)
Omics_data2[1:3, 1:3]

# clinical information
surv <- read.table("survival.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,quote = "",row.names = 1)
surv[1:3, ]

#check sample names
if(!(all(colnames(Omics_data1)==colnames(Omics_data2)) & all(colnames(Omics_data1)==rownames(surv)))) {
  cat("Samples mismatched! Processing now...\n")
  commonSam <- intersect(intersect(colnames(Omics_data1),colnames(Omics_data2)),rownames(surv))
  Omics_data1 <- Omics_data1[,commonSam]
  Omics_data2 <- Omics_data1[,commonSam]
  surv <- surv[commonSam,]
}

if((all(colnames(Omics_data1)==colnames(Omics_data2)) & all(colnames(Omics_data1)==rownames(surv)))) {
  cat("Samples successfully matched!\n") }

time <- as.numeric(surv$futime)
status <- as.numeric(surv$fustat) # 1: dead 0: alive

Omics_data1 <- as.matrix(Omics_data1)
Omics_data2 <- as.matrix(Omics_data2)

###data distribution
data.checkDistribution(Omics_data1)
data.checkDistribution(Omics_data2)

#merge
data1 <- Omics_data1
data2 <- Omics_data2
TCGA_target <- list(GeneExp=data1,miRNAExp=data2) 

#run CNMF clustering
source("ExcuteCNMF.R")
result_NMF <- ExecuteCNMF(TCGA_target, clusterNum=2,nrun=100) ##There is randomness in the results, and more cycles (nrun) are set to get more stable results
NMF_group <- result_NMF$group
NMF_distanceMatrix <- result_NMF$distanceMatrix
p_value <- survAnalysis(mainTitle="CNMF",time,status,NMF_group, NMF_distanceMatrix,similarity=TRUE)

#output KM curve
pdf("CNMF.pdf")
p_value <- survAnalysis(mainTitle="CNMF",time,status,NMF_group, NMF_distanceMatrix,similarity=TRUE)
invisible(dev.off())

#output cluster file for NTP analysis
subtype <- data.frame(samID = rownames(surv),
                      Subtype = paste0("C",NMF_group),
                      row.names = rownames(surv),
                      stringsAsFactors = F)
write.table(subtype,"CNMF_subtype.txt",row.names = T,col.names = NA,quote = F,sep = "\t")



