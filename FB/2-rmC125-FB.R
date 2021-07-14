library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast")

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
