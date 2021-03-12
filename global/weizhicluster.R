#######提取未知cluster #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")

idxSample <- BiocGenerics::which(proj5$Clusters %in% c("C3","C4","C5","C6","C9","C10","C11"))
cellsSample <- proj5$cellNames[idxSample]
nproj=proj5[cellsSample, ]
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-nProj-weizhicelltype", load = FALSE)
nproj<-loadArchRProject(path = "Save-nProj-weizhicelltype")
nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-dim8-nproj", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:8
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-dim8-nproj",
    method = "Seurat",
    name = "Clusters2",
    resolution = 0.8
)
nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-dim8-nproj", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.1, 
	force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1, p2,p3, type = "h")
plotPDF(p1,p2,p3, name = "Plot-UMAP-Clusters-nproj.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
marker<-read.csv("somemarkergene.csv",header=TRUE)
c<-as.character(marker$gene)
library(Hmisc)
c = tolower(c)
markergene<-capitalize(c)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markergene,
  transpose = TRUE
)


ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)












