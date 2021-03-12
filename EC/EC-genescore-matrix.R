rm(list = ls())
getwd()
setwd("/public/home/Chenzh275/Data/scATAC/ArchR-Endothelial_Cell")


library(ArchR)
ArchR.EC <- loadArchRProject(path = "~/Data/scATAC/ArchR-Endothelial_Cell/Save-Proj2-EC3")
GSM<-getMatrixFromProject(ArchRProj = ArchR.EC, useMatrix = "GeneScoreMatrix")
head(GSM@assays@data$GeneScoreMatrix)[,1:10]
genename<-rowData(GSM)$name ####共24333个基因
table(ArchR.EC@cellColData$Subtype)
EC1 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Subtype == "EC1" ),]
table(ArchR.EC@cellColData$Subtype)
EC2 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Subtype == "EC2" ),]
EC3 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Subtype == "EC3" ),]
EC4 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Subtype == "EC4" ),]
EC5 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Subtype == "EC5" ),]

EC1count<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(EC1))]
class(EC1count)
EC2count<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(EC2))]
EC3count<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(EC3))]
EC4count<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(EC4))]
EC5count<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(EC5))]

EC1count<-as.matrix(EC1count)
EC2count<-as.matrix(EC2count)
EC3count<-as.matrix(EC3count)
EC4count<-as.matrix(EC4count)
EC5count<-as.matrix(EC5count)

#EC1<-data.frame(mean = rowMeans(EC1count))
#EC2<-data.frame(mean = rowMeans(EC2count))
#EC3<-data.frame(mean = rowMeans(EC3count))
#EC4<-data.frame(mean = rowMeans(EC4count))
#EC5<-data.frame(mean = rowMeans(EC5count))

SUM<-cbind(EC1count,EC2count)
SUM<-cbind(SUM,EC3count)
SUM<-cbind(SUM,EC4count)
SUM<-cbind(SUM,EC5count) 

#dim(EC1)
#SUM<-cbind(EC1,EC2)
#SUM<-cbind(SUM,EC3)
#SUM<-cbind(SUM,EC4)
#SUM<-cbind(SUM,EC5)

SUM<-cbind(genename,SUM)
#SUM<-as.data.frame(SUM)
#colnames(SUM) <- c("gene","EC1","EC2","EC3","EC4","EC5")
#dim(SUM)
write.table(SUM,"All-cell-geneScore-ECs.txt")
