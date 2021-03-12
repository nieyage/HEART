library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast")

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
"Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
"Tie1","Fabp4","Esam","Kdr","Tek","Eng",###endothelial
"Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
"Dclk1","Fkbp5","Akt3",##fibroblast
"Rgs5","Ano1","Acta2",####smoothmuscle
"Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8"#####epicardial
)
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "test-fibroblast-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)


##########Visualizing Marker Genes on an Embedding###########

nproj <-addImputeWeights(
       ArchRProj = nproj,
       reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
       dimsToUse = 1:15
     )
	 
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
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
    name = "fibroblast-test-Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)
	


p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)
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
    name = "test-umap-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)
	
	
	
	
	
	
##########去掉C1，C2，C5重新cluster
idxSample <- BiocGenerics::which(nproj$Clusters2 %in% c("C1","C2","C5"))
cellsSample <- nproj$cellNames[-idxSample]
nproj=nproj[cellsSample, ]
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-fibroblast2", load = T)
nproj<-loadArchRProject(path = "Save-fibroblast2")


nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-dim15-fibroblast-0.8-40000", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(0.8), 
        sampleCells = 10000, 
        n.start = 10
    ), 
	force = TRUE,
    varFeatures = 40000, 
    dimsToUse = 1:15
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-dim15-fibroblast-0.8-40000",
    method = "Seurat",
    name = "Clusters3",
	force = TRUE,
    resolution = 0.8
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-dim15-fibroblast-0.8-40000", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.05, 
	force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "treatment", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters3", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1, p2,p3,p4,p5,p6, type = "h")
plotPDF(p1,p2,p3,p4,p5,p6, name = "fibroblast2-5-Plot-UMAP-2.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters3",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)
head(GM@assays@data$GeneScoreMatrix)[,1:10]



tdata<-MM@assays@data$z
sample_ann<-nproj@cellColData
sample_ann<-as.data.frame(sample_ann)
tdata<-as.data.frame(tdata)
sample_ann<-sample_ann[order(rownames(sample_ann)),]
tdata<-tdata[,order(colnames(tdata))]


saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-fibroblast2", load = T)
nproj<-addGeneScoreMatrix(input =nproj,force = TRUE)

pathToMacs2 <- findMacs2()
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "Clusters3")

nproj <- addReproduciblePeakSet(
    ArchRProj = nproj, 
    groupBy = "Clusters3", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(nproj)
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-fibroblast2", load = T)
nproj <- addPeakMatrix(nproj)
PM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "PeakMatrix"
)


########distance tree########

C1<-coldata[which(coldata$Clusters3=="C1"),]
rownames(C1)<-gsub("#", ".",rownames(C1))
rownames(C1)<-gsub("-", ".",rownames(C1))
C1gene<-gene[,which(colnames(gene)%in% rownames(C1))]
C1<-rowMeans(C1gene)


C2<-coldata[which(coldata$Clusters3=="C2"),]
rownames(C2)<-gsub("#", ".",rownames(C2))
rownames(C2)<-gsub("-", ".",rownames(C2))
C2gene<-gene[,which(colnames(gene)%in% rownames(C2))]
C2<-rowMeans(C2gene)

C10<-coldata[which(coldata$Clusters3=="C10"),]
rownames(C10)<-gsub("#", ".",rownames(C10))
rownames(C10)<-gsub("-", ".",rownames(C10))
C10gene<-gene[,which(colnames(gene)%in% rownames(C10))]
C10<-rowMeans(C10gene)


tall<-t(all)
tdist=dist(tall,method="euclidean")
hc=hclust(tdist,method="complete")

# 1，类平均法：average
# 2，重心法：centroid
# 3，中间距离法:median
# 4，最长距离法：complete 默认
# 5，最短距离法：single
# 6，离差平方和法：ward
# 7，密度估计法：density
y
plclust(hc)
plot(hc, hang=-1, xlab="");



#####add timepoint/treatment########
cluster<-nproj$Clusters3
table(nproj$Clusters3)
for (i in 1:20831){
if(cluster[i]=="C1"){cluster[i]="FB3"}
if(cluster[i]=="C2"){cluster[i]="FB3"}
if(cluster[i]=="C3"){cluster[i]="FB4"}
if(cluster[i]=="C4"){cluster[i]="FB4"}
if(cluster[i]=="C5"){cluster[i]="FB1"}
if(cluster[i]=="C6"){cluster[i]="FB2"}
if(cluster[i]=="C7"){cluster[i]="FB1"}
if(cluster[i]=="C8"){cluster[i]="FB1"}
if(cluster[i]=="C9"){cluster[i]="FB1"}
if(cluster[i]=="C10"){cluster[i]="FB1"}
}

cellNames <- nproj$cellNames
nproj$subtype<-cluster
nproj <- addCellColData(ArchRProj = nproj, data = paste0(cluster),
    cells=cellNames, name = "subtype",force = TRUE)




