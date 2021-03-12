########gene score matrix####
GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)

genename<-rowData(GM)$name ####共24333个基因
gene<-GM@assays@data$GeneScoreMatrix
gene<-as.data.frame(gene)
gene<-cbind(genename,gene)
rownames(gene)<-gene[,1]
gene<-gene[,-1]
########distance tree########

coldata<-nproj@cellColData
C1<-coldata[which(coldata$Clusters2=="C1"),]
C1gene<-gene[,which(colnames(gene)%in% rownames(C1))]
C1<-rowMeans(C1gene)


C2<-coldata[which(coldata$Clusters2=="C2"),]
C2gene<-gene[,which(colnames(gene)%in% rownames(C2))]
C2<-rowMeans(C2gene)

C3<-coldata[which(coldata$Clusters2=="C3"),]
C3gene<-gene[,which(colnames(gene)%in% rownames(C3))]
C3<-rowMeans(C3gene)


C4<-coldata[which(coldata$Clusters2=="C4"),]
C4gene<-gene[,which(colnames(gene)%in% rownames(C4))]
C4<-rowMeans(C4gene)

C5<-coldata[which(coldata$Clusters2=="C5"),]
C5gene<-gene[,which(colnames(gene)%in% rownames(C5))]
C5<-rowMeans(C5gene)


C6<-coldata[which(coldata$Clusters2=="C6"),]
C6gene<-gene[,which(colnames(gene)%in% rownames(C6))]
C6<-rowMeans(C6gene)

C7<-coldata[which(coldata$Clusters2=="C7"),]
C7gene<-gene[,which(colnames(gene)%in% rownames(C7))]
C7<-rowMeans(C7gene)

C8<-coldata[which(coldata$Clusters2=="C8"),]
C8gene<-gene[,which(colnames(gene)%in% rownames(C8))]
C8<-rowMeans(C8gene)

C9<-coldata[which(coldata$Clusters2=="C9"),]
C9gene<-gene[,which(colnames(gene)%in% rownames(C9))]
C9<-rowMeans(C9gene)

C10<-coldata[which(coldata$Clusters2=="C10"),]
C10gene<-gene[,which(colnames(gene)%in% rownames(C10))]
C10<-rowMeans(C10gene)

C11<-coldata[which(coldata$Clusters2=="C11"),]
C11gene<-gene[,which(colnames(gene)%in% rownames(C11))]
C11<-rowMeans(C11gene)


C12<-coldata[which(coldata$Clusters2=="C12"),]
C12gene<-gene[,which(colnames(gene)%in% rownames(C12))]
C12<-rowMeans(C12gene)

all<-cbind(C1,C2)
all<-cbind(all,C3)
all<-cbind(all,C4)
all<-cbind(all,C5)
all<-cbind(all,C6)
all<-cbind(all,C7)
all<-cbind(all,C8)
all<-cbind(all,C9)
all<-cbind(all,C10)
all<-cbind(all,C11)
all<-cbind(all,C12)


tall<-t(all)
tdist=dist(tall,method="euclidean")
hc=hclust(tdist,method="complete")
plclust(hc)
plot(hc, hang=-1, xlab="");