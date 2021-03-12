library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")

#####调取各个亚群的基因分数矩阵######
GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)

genename<-rowData(GM)$name ####共24333个基因
count<-GM@assays@data$GeneScoreMatrix   #####
FB1<-rownames(colData(GM)[colData(GM)$subtype=="FB1",])
FB2<-rownames(colData(GM)[colData(GM)$subtype=="FB2",])
FB3<-rownames(colData(GM)[colData(GM)$subtype=="FB3",])
FB4<-rownames(colData(GM)[colData(GM)$subtype=="FB4",])

FB1count<-count[,which(colnames(count)%in%FB1)]
FB2count<-count[,which(colnames(count)%in%FB2)]
FB3count<-count[,which(colnames(count)%in%FB3)]
FB4count<-count[,which(colnames(count)%in%FB4)]
FB1count<-as.matrix(FB1count)
FB2count<-as.matrix(FB2count)  
FB3count<-as.matrix(FB3count)
FB4count<-as.matrix(FB4count)
FB1count<-rowMeans(FB1count)
FB2count<-rowMeans(FB2count)
FB3count<-rowMeans(FB3count)
FB4count<-rowMeans(FB4count)
subcount<-cbind(FB1count,FB2count)
subcount<-cbind(subcount,FB3count)
subcount<-cbind(subcount,FB4count)
subcount<-cbind(genename,subcount)
subcount<-as.data.frame(subcount)
write.table(subcount,"fibroblast-allgenescorematrix-subcluster.txt")
######从mart中寻找ECM相关基因绘制violin plot#####
library(ggplot2)
library(reshape2)

testgo<-mart[which(mart$GO.term.name=="negative regulation of extracellular matrix constituent secretion"),]
testgo<-testgo$Gene.name
testgo<-subcount[which(subcount[,1]%in%testgo),]
testgo<-as.data.frame(testgo)
###准备完成矩阵，开始画小提琴图####
pdf("negative regulation of extracellular matrix constituent secretion.pdf")
#####多个基因的小提琴图#####
a<-melt(testgo, id.vars="genename",variable.name = "subtype", value.name = "genescore")
a$genescore<-as.numeric(a$genescore)	
P2<- ggplot(a, aes(x=subtype, y=genescore,fill=subtype)) + 
  geom_violin(trim=FALSE,color="white") +   geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  theme_bw()+ labs(title="GO term:negative regulation of extracellular matrix constituent secretion")+  ylab("GeneScore")+xlab("subtype")
P2

dev.off();

######单个基因的小提琴图######
p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    colorBy = "GeneScoreMatrix", 
	name = aa,######genename
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
aa<-c("Rgcc","Smad3","Sox9","Emilin1","Rgcc","Cpb2","Agt","Tmem38b","Eng","Col5a1","Col5a2","Col5a3")

library(gplots)
subcount$FB1count<-as.numeric(subcount$FB1count)
subcount$FB2count<-as.numeric(subcount$FB2count)
subcount$FB3count<-as.numeric(subcount$FB3count)
subcount$FB4count<-as.numeric(subcount$FB4count)
subcount=t(scale(t(subcount),scale = T,center = T))
mart<-read.csv("mart_export.csv")
testgo<-mart[which(mart$GO.term.name=="extracellular matrix"),]
testgo<-testgo$Gene.name
testgo<-subcount[which(rownames(subcount)%in%testgo),]
testgo<-as.data.frame(testgo)
pdf("testgo.pdf")
list<-pheatmap(testgo, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               color = colorRampPalette(colors = c("blue","white","red"))(100),
               cluster_cols=F,cutree_rows =7,
               show_colnames = T,fontsize_row = 7.5)
			   
newOrder=testgo[list$tree_row$order,]
tail(newOrder)
#####用cutree()给每一行加上cluster信息，提取出符合条件的cluster的行######
row_cluster=cutree(list$tree_row,k=7)

pheatmap(testgo, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               show_rownames=T,border_color = NA,
               color = colorRampPalette(colors = c("blue","white","red"))(100),
               cluster_cols=F,cutree_rows =7,
               show_colnames = T)


Cebpe
Cebpd
Cebpz
Cebpb
Cebpg
Cebpa
