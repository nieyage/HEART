#########Zonation for SMC#############
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
ArchR.SMC<-loadArchRProject(path = "Save-SMC-merge12")


library(monocle)
library(ggplot2)
library(patchwork)
library(magrittr)
#Extract data, phenotype data, and feature data from the SeuratObject
GM<-getMatrixFromProject(
  ArchRProj =ArchR.SMC,
  useMatrix = "GeneScoreMatrix"
)

data<-assays(GM)$GeneScoreMatrix
pd <- new('AnnotatedDataFrame', data = as.data.frame(ArchR.SMC@cellColData[,22:23]))
##########data 与 pd 行名列名顺序要一致############
colnames(data)[1:5]
row.names(pd)[1:5]
order<-row.names(pd)
#data<-data[,order]
data<-as(as.matrix(data), 'sparseMatrix')
gene<-rowData(GM)$name
rownames(data)<-gene
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
rownames(data)[1:5]
row.names(fd)[1:5]

library(reshape2)
GSM<-getMatrixFromProject(ArchRProj = ArchR.SMC, useMatrix = "GeneScoreMatrix")
head(GSM@assays@data$GeneScoreMatrix)[,1:10]
head(assays(GM)$GeneScoreMatrix)[,1:10]

genename<-rowData(GSM)$name 
count<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(ArchR.SMC))]
count<-as.matrix(count)
data <- as(as.matrix(count), 'sparseMatrix')
count2 <- count
rownames(count2) = genename
count2 <- as.data.frame(count2)
rownames(data) <- rownames(count2)

metadata<- data.frame(ArchR.SMC@cellColData[,c(1,21:23)])
head(metadata[,1:4])
metadata <- metadata[colnames(data),]
pd <- new('AnnotatedDataFrame', data = metadata)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily=VGAM:::tobit(Lower=0.1))



# 归一化 
 monocle_cds<- estimateSizeFactors(monocle_cds)
 monocle_cds<- estimateDispersions(monocle_cds)
###Filtering low-quality cells
 monocle_cds <- monocle::detectGenes(monocle_cds, min_expr = 3 )
head(featureData(monocle_cds)@data)
expressed_genes <- row.names(subset(featureData(monocle_cds)@data,
                                     num_cells_expressed >= 10))
features<-c(
  "Hspa1a","Hspa8","Dnajb14","Hspa1l","Hspa1b","Atf3","Fos",##SMC1
  "Myh11","Acta2","Cnn1","Cnn2","Cnn3","Des","Vim","Smtn","Vcl","Cald1","Tagln",#####SMC2
  "Spp1","Ereg","Eln","Mgp","Thbs1","Thbs2"####SMC3
  )

disp_table <- monocle::dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

f<-union(expressed_genes,features)
monocle_cds <- monocle::setOrderingFilter(monocle_cds, f)
#pdf("monocle_test.pdf")
#plot_ordering_genes(monocle_cds)

#Trajectory step 2: reduce data dimensionality  
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                            method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)

pdf("SMC-trajectory.pdf")
plot_cell_trajectory(monocle_cds,  color_by = "subtype",cell_size=1)+ 
 scale_color_manual(breaks = c("SMC1","SMC2","SMC3"), 
  values=c("#206A5D","#FFCC29","#BE0000")) + theme(legend.position = "right")

plot_cell_trajectory(monocle_cds,  color_by = "subtype",cell_size=0.5) +
    facet_wrap(~subtype, nrow =4)+ 
  scale_color_manual(breaks = c("SMC1","SMC2","SMC3"), 
  values=c("#206A5D","#FFCC29","#BE0000")) + theme(legend.position = "right")
plot_cell_trajectory(monocle_cds,  color_by = "State",cell_size=1)
plot_cell_trajectory(monocle_cds,  color_by = "Pseudotime",cell_size=1)
plot_cell_trajectory(monocle_cds,  color_by = "seurat_clusters",cell_size=1)
dev.off()

