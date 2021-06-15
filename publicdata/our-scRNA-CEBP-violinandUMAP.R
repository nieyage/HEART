library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
FB1<-markerList$FB1$name
FB2<-markerList$FB2$name
FB3<-markerList$FB3$name
FB4<-markerList$FB4$name

#######给全部cluster添加注释名，并根据细胞类型分类 #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")

####RNA#########
library(Seurat)
scRNA <- readRDS("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/integrative_analysis/cmc_sct.rds")
selRNA <- scRNA[,which(scRNA@meta.data$orig.ident == "A3" |
                         scRNA@meta.data$orig.ident == "A7" | scRNA@meta.data$orig.ident == "A14" |
                         scRNA@meta.data$orig.ident == "P7" | scRNA@meta.data$orig.ident == "P14" )]

selRNA<-selRNA[,which(selRNA@meta.data$cell_types == "Fibroblasts" )]
DefaultAssay(selRNA) <- "RNA" # Create dummy new assay to demo switching default assays
selRNA <- NormalizeData(selRNA) ###Normalizing the data
#Identification of highly variable features (feature selection)
selRNA <- FindVariableFeatures(selRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(selRNA)
selRNA <- ScaleData(selRNA, features = all.genes)
#Perform linear dimensional reduction
selRNA <- RunPCA(selRNA, features = VariableFeatures(object = selRNA))
#Cluster the cells
selRNA.cluster.integrated <- FindNeighbors(selRNA, dims = 1:25)
selRNA.cluster.integrated <- FindClusters(selRNA.cluster.integrated, resolution = 0.5)
nproj<- addGeneIntegrationMatrix(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  seRNA = selRNA.cluster.integrated,
  addToArrow = TRUE,
   reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
  groupRNA = "cell_types",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  force = TRUE,
  useImputation = FALSE,
  nameScore = "predictedScore_Un"
)

CEBP<-c("Cebpa","Cebpb","Cebpg","Cebpe","Cebpz")

####features 中不存在（其他几个成员均存在）

library(reshape2)
GM<-getMatrixFromProject(
  ArchRProj =nproj,
  useMatrix = "GeneIntegrationMatrix"
)
data<-assays(GM)$GeneIntegrationMatrix
data<-as(as.matrix(data), 'sparseMatrix')
gene<-rowData(GM)$name
rownames(data)<-gene
data<-as.matrix(data)
vln.df=as.data.frame(data[CEBP,])
vln.df$gene=rownames(vln.df)

vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=data.frame(CB = rownames(nproj@cellColData),celltype=nproj$subtype)
library(dplyr)
vln.df=inner_join(vln.df,anno,by="CB")

pdf("ourscRNA-Vlnplot-FBCEBP-noimpute.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+geom_boxplot(width=0.1,outlier.size=0,
    notch=TRUE,position=position_dodge(width=0.9))+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 0,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()



nproj <-addImputeWeights(
       ArchRProj = nproj,
       reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
       dimsToUse = 1:15
     )
	 
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneIntegrationMatrix", 
    name = CEBP, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)
#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
    name = "fibroblast-CEBP-scRNA-umap.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)
