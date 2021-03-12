library("devtools")
library("BiocManager")
devtools::install_local("/md01/nieyg/ArchR-master.zip")


library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 
addArchRGenome("mm10")

inputFiles <- c("AR14_fragments.tsv.gz","AR3_fragments.tsv.gz","AR7_fragments.tsv.gz",
             "NAR14_fragments.tsv.gz","NAR3_fragments.tsv.gz","NAR7_fragments.tsv.gz",
			 "P14_fragments.tsv.gz","P3_fragments.tsv.gz","P7_fragments.tsv.gz")
names(inputFiles)<-c("AR14","AR3","AR7","NAR14","NAR3","NAR7","P14","P3","P7")

inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,   ####fragment files or bam files 
  sampleNames = names(inputFiles),
  filterTSS = 6, #Dont set this too high because you can always increase later
  filterFrags = 2800, 
  addTileMat = TRUE,
  force=TRUE,
  addGeneScoreMat = TRUE
)
####for bam file ###
ArrowFiles <- createArrowFiles(
   inputFiles = "../subbam/fibroblast.all.bam.reformat.bam",
   sampleNames ="fibroblast",
   filterTSS = 4, 
   filterFrags = 1000, 
   addTileMat = TRUE,
   addGeneScoreMat = TRUE,
   minFrags=1,
   maxFrags=1e+100
 )


doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,
	dimsToUse=1:8
)



proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "downstream",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
getAvailableMatrices(proj)
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj-allsample", load = FALSE)

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj2-allsample", load = FALSE)
proj2 <- filterDoublets(proj)

#####projHemeTmp <- filterDoublets(projHeme1, filterRatio = 1.5)
proj2 <- addIterativeLSI(
    ArchRProj = proj2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-dim8", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:8
)

####cluster###
proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI-dim8",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)
###效果看上去不错
table(proj2$Clusters)

cM <- confusionMatrix(paste0(pro2$Clusters), paste0(proj2$Sample))
cM
proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI-dim8", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI-dim8", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.1, 
	force = TRUE,
    metric = "cosine"
)
####minDist越小，分的越开####


p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters-all-3.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)


proj3 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI-dim8", 
    name = "UMAP", 
    nNeighbors = 50, 
    minDist = 0.5, 
	force = TRUE,
    metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "Plot-UMAP-Sample-Clusters-all-proj3.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-Proj3-allsample", load = FALSE)
dev.off();


proj <- loadArchRProject(path = "Save-Proj-allsample")
proj2 <- loadArchRProject(path = "Save-Proj2-allsample")
proj3 <- loadArchRProject(path = "Save-Proj3-allsample")
proj5<-loadArchRProject(path = "Save-Proj5")
saveArchRProject(ArchRProj = proj5, outputDirectory = "Save-Proj5", load = FALSE)

###downstream####
markersGS <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
"Itgam","Lgals3","Lyzl1","Mrc1",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd"74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1","Tie1","Fabp4","Esam","Kdr","Tek",###endothelial
"Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra",##fibroblast
"Abcc9","Rgs5","Ano1","Acta2",####smoothmuscle
"Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8"
#####epicardial
)
"
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)

proj2 <- addImputeWeights(proj2)
p <- plotEmbedding(
    ArchRProj = proj2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj2)
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
    name = "Plot-UMAP-Marker-Genes-Imputation.pdf", 
    ArchRProj = proj2, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = proj2, 
    addDOC = FALSE, width = 5, height = 5)
######jiehe RNA-chucuo####
seRNA <- readRDS("cmc_sct.rds")
seRNA
colnames(colData(seRNA))
table(colData(seRNA)$BioClassification)
proj2 <- addGeneIntegrationMatrix(
    ArchRProj = proj2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)
###############
proj4 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters")
pathToMacs2 <- findMacs2()
proj4 <- addReproduciblePeakSet(
    ArchRProj = proj4, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(proj4)
saveArchRProject(ArchRProj = proj4, outputDirectory = "Save-Proj4", load = FALSE)
proj5 <- addPeakMatrix(proj4)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)

#pma <- markerPlot(seMarker = markersPeaks, name = "Erythroid", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
#pma

proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif")

#######Not Run###########
motifsUp <- peakAnnoEnrichment(
    seMarker = markerList,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  
motifsUp
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
  df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo
plotPDF(ggUp, ggDo, name = "-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

###########
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
  
  if("Motif" %ni% names(proj5@peakAnnotation)){
    proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif")
}
proj5 <- addBgdPeaks(proj5)

proj5 <- addDeviationsMatrix(
  ArchRProj = proj5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj5, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj5, addDOC = FALSE)
p1 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1,  type = "h")

trajectory <- c("C1", "C8", "C21")
trajectory


proj5 <- addTrajectory(
    ArchRProj = proj5, 
    name = "MyeloidU", 
    groupBy = "Clusters",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)

p <- plotTrajectory(proj5, trajectory = "MyeloidU", colorBy = "cellColData", name = "MyeloidU")

plotPDF(p, name = "Plot-MyeloidU-Traj-UMAP.pdf", ArchRProj = proj5, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = proj5, outputDirectory = "Save-Proj5", load = FALSE)






####change plot color#####
###for condition 
for (i in 1:22160){
proj2$Sample[i]=="AR"}
for (i in 22161:39940){
proj2$Sample[i]=="NAR"}

for (i in 39941:60420){
proj2$Sample[i]=="P"}
p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
#####for timepoint####
###Day3
for (i in 1:9699){
proj2$Sample[i]=="Day14"}
for (i in 9699:14701){
proj2$Sample[i]=="Day3"}
for (i in 14702:22160){
proj2$Sample[i]=="Day7"}

for (i in 22161:29619){
proj2$Sample[i]=="Day14"}
for (i in 29620:36816){
proj2$Sample[i]=="Day3"}
for (i in 36816:41656){
proj2$Sample[i]=="Day7"}
for (i in 41657:49249){
proj2$Sample[i]=="Day14"}
for (i in 49250:54656){
proj2$Sample[i]=="Day3"}
for (i in 54657:22160){
proj2$Sample[i]=="Day7"}






p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "/public/home/nieyg/scATAC-ArchR/downstream/Plot-UMAP-Sample-time-condition.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)



motifs <- c("LINE5785", "ELF1","FOS","ASCL1","SP1","GATA4")
markerMotifs <- getFeatures(proj5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
p <- plotGroups(ArchRProj = proj5, 
  groupBy = "Clusters", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(proj5)
)
p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
plotPDF(p2, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = proj5, addDOC = FALSE)