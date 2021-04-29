
library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 
addArchRGenome("mm10")

inputFiles <- c("../fragment_tbi_arrow/AR14_fragments.tsv.gz","../fragment_tbi_arrow/AR3_fragments.tsv.gz","../fragment_tbi_arrow/AR7_fragments.tsv.gz",
             "../fragment_tbi_arrow/NAR14_fragments.tsv.gz","../fragment_tbi_arrow/NAR3_fragments.tsv.gz","../fragment_tbi_arrow/NAR7_fragments.tsv.gz",
			 "../fragment_tbi_arrow/P14_fragments.tsv.gz","../fragment_tbi_arrow/P3_fragments.tsv.gz","../fragment_tbi_arrow/P7_fragments.tsv.gz","../publish_data/GSE153479_fragments.tsv.gz")
names(inputFiles)<-c("AR14","AR3","AR7","NAR14","NAR3","NAR7","P14","P3","P7","publish")

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
#########QC plot##########
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "combine-TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1,p2, name = "combine-QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

getAvailableMatrices(proj)

Name            Sample        Sample Index    Cell Barcode Suffix
P1D3MI          P1+3 dpi        SI-NA-G1              -1
P1D3Sham        P1+3 dps        SI-NA-H1              -2
P8D3MI          P8+3 dpi        SI-NA-A2              -3
P8D3Sham        P8+3 dps        SI-NA-B2              -4

#####subtract sample
sample<-rownames(proj@cellColData[which(proj$Sample=="publish"),])

P8D3MI<-grep(".*-3",sample)
P8D3MI<-sample[P8D3MI]
other<-setdiff(sample,P8D3MI)
######remove other sample
nproj<-proj[-which(rownames(proj@cellColData) %in% other ),]
df <- getCellColData(nproj, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "combine-reother-TSS-vs-Frags.pdf", ArchRProj = nproj, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = nproj)
p2 <- plotTSSEnrichment(ArchRProj = nproj)
plotPDF(p1,p2, name = "combine-reother-QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

#proj5 <- addImputeWeights(proj5)

saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-combine-all", load = FALSE)
proj<-loadArchRProject(path = "Save-combine-all")

proj <- filterDoublets(ArchRProj = proj)


proj2 <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(0.1), 
        sampleCells = 50000, 
        n.start = 10), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

proj3<- addIterativeLSI(
    ArchRProj = proj2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI3", 
    iterations = 5, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.1, 0.2,0.4,0.6), 
        sampleCells = 50000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:15
)


proj3 <- addHarmony(
    ArchRProj = proj3,
    reducedDims = "IterativeLSI3",
    name = "Harmony",
    force = TRUE,
    groupBy = "Sample"
)

####cluster###
proj3 <- addClusters(
    input = proj3,
    reducedDims = "IterativeLSI3",
    method = "Seurat",
    name = "Clusters",
    force = TRUE,
    resolution = 0.4
)

proj3 <- addUMAP(
    ArchRProj = proj3, 
    reducedDims = "IterativeLSI3", 
    name = "UMAP", 
    nNeighbors = 30, 
    force = TRUE,
    minDist = 0.5, 
    metric = "cosine"
)


p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Combine-Plot-UMAP-Sample-Clusters-all-proj3.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

#####choose proj3进行后续分析###########
###downstream####
markersGS <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####CM
"Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
"Itgam","Lgals3","Lyzl1","Mrc1",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1","Tie1","Fabp4","Esam","Kdr","Tek",###endothelial
"Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra",##fibroblast
"Abcc9","Rgs5","Ano1","Acta2",####smoothmuscle
"Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8"#####epicardial
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Combine-proj3-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)


cluster<-proj3$Clusters
table(proj3$Clusters)

for (i in 1:77506){
if(cluster[i]=="C1"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C2"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C3"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C4"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C5"){cluster[i]="Endothelial"}
if(cluster[i]=="C6"){cluster[i]="Endothelial"}
if(cluster[i]=="C7"){cluster[i]="Endothelial"}
if(cluster[i]=="C8"){cluster[i]="Endothelial"}
if(cluster[i]=="C9"){cluster[i]="Endothelial"}
if(cluster[i]=="C10"){cluster[i]="Fibroblast"}
if(cluster[i]=="C11"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C12"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C13"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C14"){cluster[i]="Fibroblast"}
if(cluster[i]=="C15"){cluster[i]="Fibroblast"}
if(cluster[i]=="C16"){cluster[i]="Fibroblast"}
if(cluster[i]=="C17"){cluster[i]="Fibroblast"}
if(cluster[i]=="C18"){cluster[i]="Fibroblast"}
}
table(cluster)
proj3$celltype<-cluster
p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
ggAlignPlots(p3, p4,p5, type = "h")
plotPDF(p3,p4,p5, name = "combine-addcelltype-UMAP-all.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj3, overwrite = F,outputDirectory = "Save-combine-all", load = FALSE)

#######提取fibroblast-cluster #####

idxSample <- BiocGenerics::which(proj3$celltype %in% c("Fibroblast"))
cellsSample <- proj3$cellNames[idxSample]
nproj=proj3[cellsSample, ]
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-combine-fibroblast", load = T)
nproj<-loadArchRProject(path = "Save-combine-fibroblast")
table(nproj$Sample)

nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-combine-fibroblast", 
    iterations = 5, 
    clusterParams = list( 
        resolution = c(0.1, 0.2,0.4,0.6), 
        sampleCells = 20000, 
        n.start = 10
    ), 
    varFeatures = 50000, 
    dimsToUse = 1:15,
    force = TRUE
)
#nproj <- addIterativeLSI(
#    ArchRProj = nproj,
#    useMatrix = "TileMatrix", 
#    name = "IterativeLSI-combine-fibroblast", 
#    iterations = 5, 
#    clusterParams = list( 
#        resolution = c(0.1, 0.2,0.4,0.6), 
#        sampleCells = 20000, 
#        n.start = 10
#    ), 
#    varFeatures = 25000, 
#    dimsToUse = 1:15,
#    force = TRUE
#)
nproj <- addHarmony(
    ArchRProj = nproj,
    reducedDims = "IterativeLSI-combine-fibroblast",
    name = "Harmony",
    force = TRUE,
    groupBy = "Sample"
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-combine-fibroblast",
    method = "Seurat",
    name = "Clusters2",
    force = TRUE,
    resolution = 0.5
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-combine-fibroblast", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.05, 
    force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1, p3,p4, type = "h")
plotPDF(p1, p3,p4, name = "combine-fibroblast-Plot-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)


######rm C1 C2 #######

idxSample <- BiocGenerics::which(nproj$Clusters2 %in% c("C1","C2","C3"))
cellsSample <- nproj$cellNames[-idxSample]
proj5=nproj[cellsSample, ]
#proj5 <- addImputeWeights(proj5)
#saveArchRProject(ArchRProj = proj5, overwrite = F,outputDirectory = "Save-cellreport-fibroblast", load = FALSE)

#####re cluster########
nproj <- addIterativeLSI(
    ArchRProj = proj5,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-combine-fibroblast2", 
    iterations = 5, 
    clusterParams = list( 
        resolution = c(0.1, 0.2,0.4,0.6), 
        sampleCells = 20000, 
        n.start = 10
    ), 
    varFeatures = 50000, 
    dimsToUse = 1:15,
    force = TRUE
)
#nproj <- addIterativeLSI(
#    ArchRProj = nproj,
#    useMatrix = "TileMatrix", 
#    name = "IterativeLSI-combine-fibroblast", 
#    iterations = 5, 
#    clusterParams = list( 
#        resolution = c(0.1, 0.2,0.4,0.6), 
#        sampleCells = 20000, 
#        n.start = 10
#    ), 
#    varFeatures = 25000, 
#    dimsToUse = 1:15,
#    force = TRUE
#)
nproj <- addHarmony(
    ArchRProj = nproj,
    reducedDims = "IterativeLSI-combine-fibroblast2",
    name = "Harmony",
    force = TRUE,
    groupBy = "Sample"
)


nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-combine-fibroblast2",
    method = "Seurat",
    name = "Clusters4",
    force = TRUE,
    resolution = 0.8
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-combine-fibroblast2", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.05, 
  force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters4", embedding = "UMAP")

ggAlignPlots(p1, p2,p4, type = "h")
plotPDF(p1, p2,p4, name = "combine-fibroblast-Plot-UMAP-rmC123.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)


cluster<-nproj$Clusters3
table(nproj$Clusters3)

for (i in 1:16015){
if(cluster[i]=="C1"){cluster[i]="FB3"}
if(cluster[i]=="C2"){cluster[i]="FB3"}
if(cluster[i]=="C3"){cluster[i]="FB4"}
if(cluster[i]=="C4"){cluster[i]="FB4"}
if(cluster[i]=="C5"){cluster[i]="FB4"}
if(cluster[i]=="C6"){cluster[i]="FB2"}
if(cluster[i]=="C7"){cluster[i]="FB2"}
if(cluster[i]=="C8"){cluster[i]="FB1"}
if(cluster[i]=="C9"){cluster[i]="FB1"}
if(cluster[i]=="C10"){cluster[i]="FB1"}
if(cluster[i]=="C11"){cluster[i]="FB1"}
}
table(cluster)
nproj$subtype<-cluster
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters3", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "subtype", embedding = "UMAP")
ggAlignPlots(p3, p4,p5, type = "h")
plotPDF(p3,p4,p5, name = "combine-addsubtype-UMAP-all.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-combine-fibroblast", load = FALSE)


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
markerGenes  <- c( "Mdk","Mfap2","Col1a1","Ptn","Col3a1","Fabp4","Gpihbp1","Cav1","Egfl7","Cd36",#####early FB1
 "Mfap5","Pi16","Sema3c","Dcn","Gfpt2","Tcf21","Pdgfra","Ddr2",####Muture   FB4
 "Pim1","Brinp1","Cdh18","Grid2",####FB2
 "Col10a1","Rell2","Hspa8","Hspa1a","Hspa1b","Hspb6"###FB3
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "combine-fibroblast-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

    AR14  AR3  AR7  P14   P3   P7
FB1   70   45 2410   37 1493 3153
FB2    8 1344    9    2   41   22
FB3  580  120  585  391  248  502
FB4 2636    8  108 2143   10   50


#######
pathToMacs2 <- findMacs2()
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "subtype")

nproj <- addReproduciblePeakSet(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(nproj)
nproj <- addPeakMatrix(nproj)


markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
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
plotPDF(heatmapPeaks, name = "combine-fibroblast-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

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


