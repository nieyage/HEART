setwd("/public/home/nieyg/scATAC-ArchR/")

library(ArchR)
ArchR.EC <- loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")

#getAvailableMatrices(ArchR.EC)

table(ArchR.EC@cellColData$Clusters_3)
table(ArchR.EC$Clusters_3)

C1 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C1" ),]
C2 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C2" ),]
C3 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C3" ),]
C4 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C4" ),]
C5 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C5" ),]
C6 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C6" ),]
C7 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C7" ),]
C8 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C8" ),]
C9 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C9" ),]
C10 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C10" ),]
C11 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C11" ),]
C12 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C12" ),]
C13 <- ArchR.EC[which(ArchR.EC@cellColData$Clusters_3 == "C13" ),]

table(C1$Sample)
table(C2$Sample)
table(C3$Sample)
table(C4$Sample)
table(C5$Sample)
table(C6$Sample)
table(C7$Sample)
table(C8$Sample)
table(C9$Sample)
table(C10$Sample)
table(C11$Sample)
table(C12$Sample)
table(C13$Sample)

####Extract geneScore Matrix####
GSM<-getMatrixFromProject(ArchRProj = ArchR.EC, useMatrix = "GeneScoreMatrix")
#head(GSM@assays@data$GeneScoreMatrix)[,1:10]
genename<-rowData(GSM)$name ####共24333个基因

gene<-GSM@assays@data$GeneScoreMatrix
gene<-as.data.frame(gene)
class(gene)
dim(gene)
gene<-cbind(genename,gene)
rownames(gene)<-gene[,1]
head(gene[,1:5])
gene<-gene[,-1]

C1 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C1" ),]
C2 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C2" ),]
C3 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C3" ),]
C4 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C4" ),]
C5 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C5" ),]
C6 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C6" ),]
C7 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C7" ),]
C8 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C8" ),]
C9 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C9" ),]
C10 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C10" ),]
C11 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C11" ),]
C12 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C12" ),]
C13 <- ArchR.EC@cellColData[which(ArchR.EC@cellColData$Clusters_3 == "C13" ),]

count<-GSM@assays@data$GeneScoreMatrix   
C1<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C1",])
C2<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C2",])
C3<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C3",])
C4<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C4",])
C5<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C5",])
C6<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C6",])
C7<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C7",])
C8<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C8",])
C9<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C9",])
C10<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C10",])
C11<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C11",])
C12<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C12",])
C13<-rownames(colData(GSM)[colData(GSM)$Clusters_3=="C13",])

C1count<-count[,which(colnames(count)%in%rownames(C1))]
C2count<-count[,which(colnames(count)%in%rownames(C2))]
C3count<-count[,which(colnames(count)%in%rownames(C3))]
C4count<-count[,which(colnames(count)%in%rownames(C4))]
C5count<-count[,which(colnames(count)%in%rownames(C5))]
C6count<-count[,which(colnames(count)%in%rownames(C6))]
C7count<-count[,which(colnames(count)%in%rownames(C7))]
C8count<-count[,which(colnames(count)%in%rownames(C8))]
C9count<-count[,which(colnames(count)%in%rownames(C9))]
C10count<-count[,which(colnames(count)%in%rownames(C10))]
C11count<-count[,which(colnames(count)%in%rownames(C11))]
C12count<-count[,which(colnames(count)%in%rownames(C12))]
C13count<-count[,which(colnames(count)%in%rownames(C13))]

C1count<-as.matrix(C1count)
C2count<-as.matrix(C2count)
C3count<-as.matrix(C3count)
C4count<-as.matrix(C4count)
C5count<-as.matrix(C5count)
C6count<-as.matrix(C6count)
C7count<-as.matrix(C7count)
C8count<-as.matrix(C8count)
C9count<-as.matrix(C9count)
C10count<-as.matrix(C10count)
C11count<-as.matrix(C11count)
C12count<-as.matrix(C12count)
C13count<-as.matrix(C13count)

C1<-rowMeans(C1count)
C2<-rowMeans(C2count)
C3<-rowMeans(C3count)
C4<-rowMeans(C4count)
C5<-rowMeans(C5count)
C6<-rowMeans(C6count)
C7<-rowMeans(C7count)
C8<-rowMeans(C8count)
C9<-rowMeans(C9count)
C10<-rowMeans(C10count)
C11<-rowMeans(C11count)
C12<-rowMeans(C12count)
C13<-rowMeans(C13count)

subcount<-cbind(C1,C2)
subcount<-cbind(subcount,C3)
subcount<-cbind(subcount,C4)
subcount<-cbind(subcount,C5)
subcount<-cbind(subcount,C6)
subcount<-cbind(subcount,C7)
subcount<-cbind(subcount,C8)
subcount<-cbind(subcount,C9)
subcount<-cbind(subcount,C10)
subcount<-cbind(subcount,C11)
subcount<-cbind(subcount,C12)
subcount<-cbind(subcount,C13)

subcount<-cbind(genename,subcount)
subcount<-as.data.frame(subcount)
#write.table(subcount,"Endhothelial-cell-GenescoreMatrix-cluster_3.txt")

head(subcount[,1:2])
rownames(subcount)= subcount$genename
subcount =subcount[,-1]
tsubcount = t(subcount)
tdist=dist(tsubcount,method="euclidean")
hc=hclust(tdist,method="complete")
plclust(hc)
plot(hc, hang=-1, xlab="");

# 注释：在聚类中求两点的距离有：
# 1，绝对距离：manhattan
# 2，欧氏距离：euclidean 默认
# 3，闵科夫斯基距离：minkowski
# 4，切比雪夫距离：chebyshev
# 5，马氏距离：mahalanobis
# 6，蓝氏距离：canberra

# 1，类平均法：average
# 2，重心法：centroid
# 3，中间距离法:median
# 4，最长距离法：complete 默认
# 5，最短距离法：single
# 6，离差平方和法：ward
# 7，密度估计法：density

####time point
table(ArchR.EC$Clusters_3)
table(ArchR.EC$Sample)
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("sample","Timepoint");
meta[,1] = ArchR.EC$Sample
meta<-as.data.frame(meta)
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "AR3" |cl[i] == "NAR3" |cl[i] == "P3")
    meta[i,2] = "Day3"
  else if (cl[i] == "AR7" |cl[i] == "NAR7" |cl[i] == "P7")
    meta[i,2]  = "Day7"
  else
    meta[i,2]  = "Day14"
}
head(meta[3525:3664,])
table(ArchR.EC$Timepoint)
meta = data.frame(meta)
ArchR.EC$Timepoint  = meta$Timepoint

ArchR.EC <- addUMAP(ArchRProj = ArchR.EC, reducedDims = "IterativeLSI", nNeighbors = 10, minDist = 0.3, metric = "cosine", force = TRUE)
p5 <- plotEmbedding(ArchRProj = ArchR.EC, colorBy = "cellColData", name = "Timepoint", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p5, name = "UMAP-EC-timepoint-Clusters_3-iter2-01.pdf", ArchRProj = ArchR.EC, addDOC = FALSE, width = 5, height = 5)

####add subtype
table(ArchR.EC$Clusters_3)
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("cluster","subtype");
meta[,1] = ArchR.EC$Clusters_3
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "C4" |cl[i] == "C5" |cl[i] == "C6" |cl[i] == "C7" |cl[i] == "C10"  )
    meta[i,2] = "EC1"
  else if (cl[i] == "C2"|cl[i] == "C3" )
    meta[i,2]  = "EC2"
  else if (cl[i] == "C8" )
    meta[i,2]  = "EC3"
  else if (cl[i] == "C1" | cl[i] == "C12" |cl[i] == "C13")
    meta[i,2]  = "EC4"
  else
    meta[i,2]  = "EC5"
}
head(meta)
meta = data.frame(meta)
ArchR.EC$Subtype  = meta$subtype
table(ArchR.EC$Subtype)
p6 <- plotEmbedding(ArchRProj = ArchR.EC, colorBy = "cellColData", name = "Subtype", embedding = "UMAP")
plotPDF(p6, name = "UMAP-EC-Subtype-Clusters_3-iter2.pdf", ArchRProj = ArchR.EC, addDOC = FALSE, width = 5, height = 5)

EC1<-ArchR.EC[which(ArchR.EC$Subtype == "EC1" ),]
EC2<-ArchR.EC[which(ArchR.EC$Subtype == "EC2" ),]
EC3<-ArchR.EC[which(ArchR.EC$Subtype == "EC3" ),]
EC4<-ArchR.EC[which(ArchR.EC$Subtype == "EC4" ),]
EC5<-ArchR.EC[which(ArchR.EC$Subtype == "EC5" ),]

table(EC1$Sample)
table(EC2$Sample)
table(EC3$Sample)
table(EC4$Sample)
table(EC5$Sample)
