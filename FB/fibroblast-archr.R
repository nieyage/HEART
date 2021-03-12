#######提取fibroblast-cluster #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")

idxSample <- BiocGenerics::which(proj5$Clusters %in% c("C3","C4","C5","C6","C2","C9","C10","C11"))
cellsSample <- proj5$cellNames[idxSample]
nproj=proj5[cellsSample, ]
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-fibroblast", load = T)
nproj<-loadArchRProject(path = "Save-fibroblast")
#####add timepoint/treatment########
sample<-nproj$Sample
table(nproj$Sample)
for (i in 1:21621){
if(sample[i]=="AR3"){sample[i]="Day3"}
if(sample[i]=="AR7"){sample[i]="Day7"}
if(sample[i]=="AR14"){sample[i]="Day14"}
if(sample[i]=="NAR3"){sample[i]="Day3"}
if(sample[i]=="NAR7"){sample[i]="Day7"}
if(sample[i]=="NAR14"){sample[i]="Day14"}
if(sample[i]=="P3"){sample[i]="Day3"}
if(sample[i]=="P7"){sample[i]="Day7"}
if(sample[i]=="P14"){sample[i]="Day14"}
}

cellNames <- nproj$cellNames
nproj$timepoint<-sample2
nproj <- addCellColData(ArchRProj = nproj, data = paste0(sample),
    cells=cellNames, name = "timepoint",force = TRUE)

sample2<-nproj$Sample
table(nproj$Sample)
for (i in 1:21621){
if(sample2[i]=="AR3"){sample2[i]="AR"}
if(sample2[i]=="AR7"){sample2[i]="AR"}
if(sample2[i]=="AR14"){sample2[i]="AR"}
if(sample2[i]=="NAR3"){sample2[i]="NAR"}
if(sample2[i]=="NAR7"){sample2[i]="NAR"}
if(sample2[i]=="NAR14"){sample2[i]="NAR"}
if(sample2[i]=="P3"){sample2[i]="P"}
if(sample2[i]=="P7"){sample2[i]="P"}
if(sample2[i]=="P14"){sample2[i]="P"}
}
nproj$treatment<-sample2
nproj <- addCellColData(ArchRProj = nproj, data = paste0(sample2),
    cells=cellNames, name = "treatment",force = TRUE)
	table(nproj$treatment)


nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-dim15-fibroblast-1.2-40000", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(1.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 50000, 
    dimsToUse = 1:15
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
    method = "Seurat",
    name = "Clusters2",
	force = TRUE,
    resolution = 0.8
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000", 
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

ggAlignPlots(p1, p2,p3,p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "fibroblast-Plot-UMAP-last.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("AR14")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "fibroblast-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)


#######
pathToMacs2 <- findMacs2()
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "Clusters2")

nproj <- addReproduciblePeakSet(
    ArchRProj = nproj, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(nproj)
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-fibroblast", load = T)
nproj <- addPeakMatrix(nproj)


markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "fibroblast-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

nproj <- addMotifAnnotations(ArchRProj = nproj, force = TRUE,motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 1"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "fibroblast-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
  
####去掉同一个基因家族的TF，仅保留一个，添加尽可能多的基因家族进来

  idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% c("Klf4_143","Klf1_214","Klf3_804","Klf2_819","Klf8_194","Klf12_236",
           "Sp3_802","Smad5_883","Sp2_150","Sp5_238","Sp6_191","Sp4_167","Klf7_171","Klf5_145","E2f4_257","Klf14_237","Klf11_798","Klf13_817","Bach2_119",
		   "Tcf3_31","Tcf12_59","Ascl2_23","Klf15_806","Klf16_874","Sp7_222","Sp8_207","Sp9_231","Zfp219_816","Gata5_385","Gata4_386","Gata1_387","Mesp2_57",
		   "Gata2_383","Sox9_725","Egr2_188","Egr3_183","Zfp161_209","Zic5_195","Tcfap2c_3","Gata3_384",""))
cellsSample <- enrichMotifs@NAMES[-idxSample]
enrichMotifs1=enrichMotifs[cellsSample, ]


heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 30, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "fibroblast-Motifs-Enriched-Marker-Heatmap-5",width = 12, height = 6,ArchRProj = nproj, addDOC = FALSE)

####轨迹分析#####


saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-fibroblast", load = FALSE)

trajectory <- c("NAR3", "NAR7", "NAR14")
trajectory
nproj <- addTrajectory(
    ArchRProj = nproj, 
    name = "NAR", 
    groupBy = "Sample",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)

p <- plotTrajectory(nproj, trajectory = "NAR", colorBy = "cellColData", name = "NAR")

plotPDF(p, name = "Plot-NAR-Traj-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)
####test######
####test######
trajectory <- c("P3", "P7", "P14")
trajectory
nproj <- addTrajectory(
    ArchRProj = nproj, 
    name = "P", 
    groupBy = "Sample",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)

p <- plotTrajectory(nproj, trajectory = "P", colorBy = "cellColData", name = "P")

plotPDF(p, name = "Plot-P-Traj-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)


p5 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
ggAlignPlots(p5, type = "h")
plotPDF(p5, name = "fibroblast-Sample-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-fibroblast", load = T)


###################
######monocle######
###################
nproj <- addBgdPeaks(nproj)
nproj <- addDeviationsMatrix(
  ArchRProj = nproj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(nproj, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = nproj, addDOC = FALSE)
getAvailableMatrices(nproj)
nproj@cellColData

markerMotifs <- getFeatures(projHeme5, useMatrix = "MotifMatrix")
markerMotifs


MM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "MotifMatrix"
)

tdata<-MM@assays@data$z
sample_ann<-nproj@cellColData
sample_ann<-as.data.frame(sample_ann)
tdata<-as.data.frame(tdata)
sample_ann<-sample_ann[order(rownames(sample_ann)),]
tdata<-tdata[,order(colnames(tdata))]

pd <- new('AnnotatedDataFrame', data = sample_ann)
fData <- data.frame(gene_short_name = row.names(tdata), row.names = row.names(tdata))
head(fData)
fd <- new('AnnotatedDataFrame', data = fData)
head(fd)

monocle_cds <- newCellDataSet(as.matrix(tdata),
                phenoData = pd,
                featureData = fd,
                lowerDetectionLimit = 0.1,
                expressionFamily = tobit(Lower = 0.1))

#check phenodata、featuredata
head(pData(monocle_cds))
head(fData(monocle_cds))
monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
head(pData(monocle_cds))
#monocle_cds$cluster <- as.character(monocle_cds$cluster)
cluster<-nproj$Clusters2
table(nproj$Clusters2)

for (i in 1:17479){
if(cluster[i]=="C1"){cluster[i]="sub1"}
if(cluster[i]=="C2"){cluster[i]="sub2"}
if(cluster[i]=="C3"){cluster[i]="sub5"}
if(cluster[i]=="C4"){cluster[i]="sub5"}
if(cluster[i]=="C5"){cluster[i]="sub3"}
if(cluster[i]=="C6"){cluster[i]="sub4"}
if(cluster[i]=="C7"){cluster[i]="sub4"}
if(cluster[i]=="C8"){cluster[i]="sub6"}
if(cluster[i]=="C9"){cluster[i]="sub6"}
if(cluster[i]=="C10"){cluster[i]="sub5"}
if(cluster[i]=="C11"){cluster[i]="sub5"}
if(cluster[i]=="C12"){cluster[i]="sub5"}
if(cluster[i]=="C13"){cluster[i]="sub5"}
if(cluster[i]=="C14"){cluster[i]="sub5"}


}
table(cluster)
nproj$sub<-cluster
cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(cluster),
    cells=cellNames, name = "sub",force = TRUE)
	
pdf("monocle2-4.pdf")
plot_cell_trajectory(monocle_cds,color_by = "Sample",cell_size = 0.01)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.01)
plot_cell_trajectory(monocle_cds,color_by = "timepoint",cell_size = 0.01)
plot_cell_trajectory(monocle_cds, color_by = "treatment",cell_size = 0.01)


plot_cell_trajectory(monocle_cds, color_by = "Sample",cell_size = 0.01) + facet_wrap(~Sample, nrow = 3)


pdf("monocle2-5.pdf")
plot_cell_trajectory(monocle_cds, color_by = "Clusters2",cell_size = 0.01) + facet_wrap(~Clusters2, nrow = 3)

dev.off();
pdf("monocle2-6.pdf")
plot_cell_trajectory(monocle_cds, color_by = "sub",cell_size = 0.01) + facet_wrap(~sub, nrow = 2)

dev.off();

#计算pseudotime值
pdf("monocle2-pseudotime-comp3.pdf")
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")






##############################################
######compute every cluster cell numbers######
##############################################
idxSample <- BiocGenerics::which(nproj$Clusters2 %in% c("C14"))
cellsSample <- nproj$cellNames[idxSample]
C=nproj[cellsSample, ]
table(C$Sample)













######添加celltype信息###
cluster<-nproj$Clusters2
table(nproj$Clusters2)

for (i in 1:17479){
if(cluster[i]=="C1"){cluster[i]="sub1"}
if(cluster[i]=="C2"){cluster[i]="sub3"}
if(cluster[i]=="C3"){cluster[i]="sub1"}
if(cluster[i]=="C4"){cluster[i]="sub3"}
if(cluster[i]=="C5"){cluster[i]="sub3"}
if(cluster[i]=="C6"){cluster[i]="sub3"}
if(cluster[i]=="C7"){cluster[i]="sub3"}
if(cluster[i]=="C8"){cluster[i]="sub2"}
if(cluster[i]=="C9"){cluster[i]="sub3"}

}
table(cluster)
nproj$sub<-cluster
cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(cluster),
    cells=cellNames, name = "sub",force = TRUE)
p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "sub", embedding = "UMAP")

ggAlignPlots(p1, type = "h")
plotPDF(p1, name = "fibroblast-Plot-UMAP-sub.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "sub",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
sub1:
	Hr, Ngf, Platr22, Ifi214, Calhm4, Glipr1, Mir546, Avil, Gm12295, Msgn1, 4930448C13Rik, Serpina3g, Serpina3h, Mir1983, D730050B12Rik, Edil3, 4930474H20Rik, Cdh9, Acot10, Rspo2
sub2:
	Platr22, Ifi214, Calhm4, Glipr1, Mir546, Avil, Gm12295, 4930448C13Rik, Serpina3g, Serpina3h, Mir1983, D730050B12Rik, Edil3, 4930474H20Rik, Cdh9, Acot10, Rspo2, Trat1, Vgll3, Vmn2r95
sub3:
	Msgn1, Tvp23a, Platr22, Ifi214, Calhm4, Glipr1, Mir546, Avil, Gm12295, 4930448C13Rik, Serpina3g, Serpina3h, Mir1983, D730050B12Rik, Edil3, Hr, 4930474H20Rik, Cdh9, Acot10, Rspo2
	
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "fibroblast-GeneScores-Marker-Heatmap-sub", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "sub",
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
plotPDF(heatmapPeaks, name = "fibroblast-Peak-Marker-Heatmap-sub", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

nproj <- addMotifAnnotations(ArchRProj = nproj, force = TRUE,motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap-sub", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
  

trajectory <- c("sub2", "sub3","sub1" )
trajectory
nproj <- addTrajectory(
    ArchRProj = nproj, 
    name = "MyeloidU", 
    groupBy = "sub",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)

p <- plotTrajectory(nproj, trajectory = "MyeloidU", colorBy = "cellColData", name = "MyeloidU")

plotPDF(p, name = "Plot-MyeloidU-Traj-UMAP-sub.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

