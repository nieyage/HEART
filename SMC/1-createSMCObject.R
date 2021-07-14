#######提取Smooth muscle cells-cluster #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")

idxSample <- BiocGenerics::which(proj5$Clusters %in% c("C7","C8"))
cellsSample <- proj5$cellNames[idxSample]
nproj=proj5[cellsSample, ]
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-SMC", load = T)
nproj<-loadArchRProject(path = "Save-SMC")
#####add timepoint/treatment########
sample<-nproj$Sample
table(nproj$Sample)
for (i in 1:4109){
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
nproj$timepoint<-sample

sample2<-nproj$Sample
table(nproj$Sample)
for (i in 1:4109){
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

nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-dim15-SMC-1.2-80000-4", 
    iterations = 4, 
    clusterParams = list( 
        resolution = c(1.2), 
        sampleCells = 3000, 
        n.start = 10
    ), 
    varFeatures = 80000, 
    dimsToUse = 1:15
)

nproj<- addGeneIntegrationMatrix(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI-dim15-SMC-1.2-80000-4",
    seRNA = RNA,
    addToArrow = FALSE,
    groupRNA = "seurat_clusters",##seurat_clusters
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-dim15-SMC-1.2-80000-4",
    method = "Seurat",
    name = "Clusters2",
	force = TRUE,
    resolution = 1.2
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-dim15-SMC-1.2-80000-4", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.01, 
	force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "treatment", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")

ggAlignPlots(p1, p2,p3,p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "SMC-Plot-UMAP-6.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("Rgs5","Ano1","Acta2")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "SMC-Cluster-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

#####根据以上结果分亚型#########
cluster<-nproj$Clusters2
table(nproj$Clusters2)

for (i in 1:4109){
if(cluster[i]=="C1"){cluster[i]="SMC2"}
if(cluster[i]=="C2"){cluster[i]="SMC2"}
if(cluster[i]=="C3"){cluster[i]="SMC2"}
if(cluster[i]=="C4"){cluster[i]="SMC2"}
if(cluster[i]=="C5"){cluster[i]="SMC1"}
if(cluster[i]=="C6"){cluster[i]="SMC1"}
if(cluster[i]=="C7"){cluster[i]="SMC3"}
if(cluster[i]=="C8"){cluster[i]="SMC3"}
if(cluster[i]=="C9"){cluster[i]="SMC3"}
if(cluster[i]=="C10"){cluster[i]="SMC1"}
if(cluster[i]=="C11"){cluster[i]="SMC1"}
if(cluster[i]=="C12"){cluster[i]="SMC1"}
}
table(cluster)
nproj$subtype<-cluster
cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(cluster),
    cells=cellNames, name = "subtype",force = TRUE)

#####根据subtype寻找特异性基因#########
markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes<-c("Rgs5","Ano1","Acta2","Talgln","Mustn1","Myh11","Mylk")####smoothmuscle
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "SMC-subtype-Marker-Heatmap", width = 8, height = 4, ArchRProj = nproj, addDOC = FALSE)

#######find peak#######
pathToMacs2 <- findMacs2()
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "subtype")
nproj <- addReproduciblePeakSet(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(nproj)
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-SMC", load = T)
nproj <- addPeakMatrix(nproj)
markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", returnGR = TRUE)
markerList
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "SMC-Peak-Marker-Heatmap", width = 6, height = 4, ArchRProj = nproj, addDOC = FALSE)
nproj <- addMotifAnnotations(ArchRProj = nproj, force = TRUE,motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 30, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "SMC-Motifs-Enriched-Marker-Heatmap", width = 12, height = 6, ArchRProj = nproj, addDOC = FALSE)
  
####去掉同一个基因家族的TF，仅保留一个，添加尽可能多的基因家族进来
  idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% c("Klf4_143","Klf1_214","Klf3_804","Klf2_819","Klf8_194","Klf12_236",
           "Sp3_802","Smad5_883","Sp2_150","Sp5_238","Sp6_191","Sp4_167","Klf7_171","Klf5_145","E2f4_257","Klf14_237","Klf11_798","Klf13_817","Bach2_119",
		   "Tcf3_31","Tcf12_59","Ascl2_23","Klf15_806","Klf16_874","Sp7_222","Sp8_207","Sp9_231","Zfp219_816","Gata5_385","Gata4_386","Gata1_387","Mesp2_57",
		   "Gata2_383","Sox9_725","Egr2_188","Egr3_183","Zfp161_209","Zic5_195","Tcfap2c_3","Gata3_384","Zfp281_193",
		   "Zfp202_169","LINE3883_876","Ctcf_146","Zbtb1_182","E2f6_264","Id4_35","Fosb_98","Junb_98","Fosl1_107","Batf3_789",
		   "Nfe2l2_101","Bcl11b_814","LINE5785_293","Pou4f2_840","Pou4f1_841","LINE10076_380","LINE9919_344","LINE10042_370",
		   "Pou3f2_555","Pou3f4_561","Foxc2_316","Mef2c_638","Mef2d_842","Mef2b_882","Nfix_722","Nfib_859","Nfia_862","LINE9979_373","Tbpl2_773","Hic2_211",
		   "Junb_127","Jun_126","Sox12_741","LINE9832_335","LINE9950_357","LINE10046_367","LINE10077_381","LINE9891_343","LINE10016_369","LINE15940_622","Pou6f1_416",
		   "LINE6279_607","Foxs1_329","Irf2_636","Zfp148_800"))
cellsSample <- enrichMotifs@NAMES[-idxSample]
enrichMotifs1=enrichMotifs[cellsSample, ]
heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 30, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "SMC-Motifs-rmrep",width = 12, height = 6,ArchRProj = nproj, addDOC = FALSE)
