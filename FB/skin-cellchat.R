library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
fetal<-readRDS("fetalbrain.rds")
library(devtools)
install_github("arc85/celltalker")
library(celltalker)

#鉴定差异表达的配体和受体

#在我们的数据集中识别配体和受体

ligs <- as.character(unique(ramilowski_pairs$ligand))
recs <- as.character(unique(ramilowski_pairs$receptor))

ligs.present <- rownames(fetal)[rownames(fetal) %in% ligs]
recs.present <- rownames(fetal)[rownames(fetal) %in% recs]

genes.to.use <- union(ligs.present,recs.present) # union用于合并子集

# 使用FindAllMarkers区分组之间差异表达的配体和受体

Idents(fetal) <- "type"
markers <- FindAllMarkers(fetal, assay="RNA",features=genes.to.use,only.pos=TRUE)
ligs.recs.use <- unique(markers$gene)
length(ligs.recs.use)


interactions.forward1 <- ramilowski_pairs[as.character(ramilowski_pairs$ligand) %in% ligs.recs.use,]
interactions.forward2 <- ramilowski_pairs[as.character(ramilowski_pairs$receptor) %in% ligs.recs.use,]
interact.for <- rbind(interactions.forward1, interactions.forward2)
dim(interact.for)

#产生3238个配体和受体相互作用

#Create data for celltalker

expr.mat <- GetAssayData(fetal,slot="counts")
rownames(expr.mat)<-gsub('\\.[0-9]*', '', rownames(expr.mat))

defined.type <- fetal@meta.data$type
defined.subtype <- fetal@meta.data$type
defined.pos <- fetal@meta.data$Pos

reshaped.matrices <- reshape_matrices(count.matrix=expr.mat,clusters=defined.type,groups=defined.type,replicates=defined.pos,ligands.and.receptors=interact.for)

# freq.group.in.cluster: 只对包含细胞数大于总细胞数5%的簇进行互作分析
put.int <- putative_interactions(ligand.receptor.tibble=reshaped.matrices,clusters=defined.type,groups=defined.subtype,freq.group.in.cluster=0.05,ligands.and.receptors=interact.for)

put.int <- putative_interactions(ligand.receptor.tibble=reshaped.matrices,clusters=defined.type,groups=defined.subtype,
	ligands.and.receptors=interact.for)


#Identify unique ligand/receptor interactions present in each sample
unique.ints <- unique_interactions(put.int,group1="pbmc",group2="tonsil",interact.for)

#Get data to plot circos for PBMC
pbmc.to.plot <- pull(unique.ints[1,2])[[1]]
for.circos.pbmc <- pull(put.int[1,2])[[1]][pbmc.to.plot]

circos_plot(interactions=for.circos.pbmc,clusters=defined.clusters)




library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)

options(stringsAsFactors = FALSE)
fetal<-readRDS("fetalbrain.rds")
data.input  <- fetal@assays$RNA@data
identity = data.frame(group =fetal$type   , row.names = names(fetal$type)) # create a dataframe consisting of the cell labels
unique(identity$group)
cellchat <- createCellChat(fetal,group.by ="type")
cellchat
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use # set the used database in the object

>  unique(CellChatDB$interaction$annotation)
[1] "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
 cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  

cellchat <- computeCommunProb(cellchat)  #注意这个函数如果你可以用就用，这个是作者的。
mycomputeCommunProb <-edit(computeCommunProb)  # computeCommunProb内部似乎有一些bug，同一套数据在window10上没事，到了Linux上有报错。发现是computeExpr_antagonist这个函数有问题，(matrix(1, nrow = 1, ncol = length((group))))，中应为(matrix(1, nrow = 1, ncol = length(unique(group))))？ 不然矩阵返回的不对。de了它。
environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- mycomputeCommunProb(cellchat)  # 这儿是我de过的。


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
pdf("test-cellchat.pdf")
vertex.receiver = seq(1,4)
pathways.show <- "EGF"
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)   # 原函数
netAnalysis_contribution(cellchat, signaling = pathways.show)

groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



###########skin looom file #####
library(Seurat)
library(loomR)
library(dplyr)
skin<-read.table("GSE162183_Raw_gene_counts_matrix.tab.gz",header=TRUE,row.names=1)
skin  <-  CreateSeuratObject(counts = skin, min.cells=3,project="skin", assay = "RNA")

skin_loom<-create(filename = "GSE162183_Raw_gene_counts_matrix_LoomFile.loom.gz",data = skin@assays$RNA@counts,overwrite = T)

skin_loom$get.attribute.df(MARGIN = 1,attributes = "Gene")[1:6,]
skin_loom$get.attribute.df(MARGIN = 2,attributes = "CellID")[1:6,]

skin_seurat<-as.Seurat(skin_loom)

skin <- NormalizeData(skin, normalization.method = "LogNormalize", scale.factor = 10000)
skin <- FindVariableFeatures(skin, selection.method = "vst", nfeatures = 2000)
skin <- ScaleData(skin)
skin <- RunPCA(skin, features = VariableFeatures(object = skin))

# 对于较大的数据集 JackStraw 比较耗时
skin <- JackStraw(skin, num.replicate = 100)
skin <- ScoreJackStraw(skin, dims = 1:30)

ElbowPlot(skin)


skin <- FindNeighbors(skin, dims = 1:20)
skin <- FindClusters(skin, resolution = 0.5)
skin <- RunUMAP(skin, dims = 1:20)

# 可用 `label = TRUE` 参数标识每个 cluster
pdf("skin.pdf")
DimPlot(skin, reduction = "umap",label = TRUE)

dev.off()
gene<-c("CDH5","PECAM1","CD34","KRT10","KRT14","KRT1","PTPRC","ITGAM","CD3E","LUM","PDGFRB","COL1A1","SOX10","MLANA","S100B")




DefaultAssay(skin) <- "RNA"

pdf("skin_cluster_annotation-all.pdf")
DotPlot(skin, features = gene,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()

new.cluster.ids <- c("EpD","Endo","Mes","Mes","EpD","Mes","Mes","IM","EpD","EpD",
	"IM","EpD","Endo","EpD","unknown","Mes","EpD","SchM","Endo","SchM","IM")

names(new.cluster.ids) <- levels(skin)
skin <- RenameIdents(skin, new.cluster.ids)

pdf("skin-umap.pdf")
DimPlot(skin, reduction = "umap",label = TRUE)
dev.off()
skin@active.ident<-skin$orig.ident
skin<-subset(skin, idents = c("Ctrl1","Ctrl2","Ctrl3"))
skin@active.ident<-skin$seurat_clusters
skin$celltype<-skin@active.ident
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)

skin<-readRDS("skin.rds")
options(stringsAsFactors = FALSE)
data.input  <- skin@assays$RNA@data
identity = data.frame(group =skin$celltype   , row.names = names(skin$celltype)) # create a dataframe consisting of the cell labels
unique(identity$group)
cellchat <- createCellChat(skin,group.by ="celltype")
cellchat
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  

cellchat <- computeCommunProb(cellchat)  #注意这个函数如果你可以用就用，这个是作者的。
mycomputeCommunProb <-edit(computeCommunProb)  # computeCommunProb内部似乎有一些bug，同一套数据在window10上没事，到了Linux上有报错。
#发现是computeExpr_antagonist这个函数有问题，(matrix(1, nrow = 1, ncol = length((group))))，
#中应为(matrix(1, nrow = 1, ncol = length(unique(group))))？ 不然矩阵返回的不对。de了它。
environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- mycomputeCommunProb(cellchat)  # 这儿是我de过的。


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
pdf("skin-cellchat.pdf")
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
