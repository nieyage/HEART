library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 
addArchRGenome("mm10")

inputFiles <- c("GSE153479_fragments.tsv.gz")
names(inputFiles)<-c("sample")

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
    LSIMethod = 1
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

plotPDF(p, name = "cellreport-TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1,p2, name = "cellreport-QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

getAvailableMatrices(proj)

Name            Sample        Sample Index    Cell Barcode Suffix
P1D3MI          P1+3 dpi        SI-NA-G1              -1
P1D3Sham        P1+3 dps        SI-NA-H1              -2
P8D3MI          P8+3 dpi        SI-NA-A2              -3
P8D3Sham        P8+3 dps        SI-NA-B2              -4

P1D3MI<-grep(".*-1",rownames(proj@cellColData))
P1D3MI<-rownames(proj@cellColData)[P1D3MI]
P1D3Sham<-grep(".*-2",rownames(proj@cellColData))
P1D3Sham<-rownames(proj@cellColData)[P1D3Sham]

P8D3MI<-grep(".*-3",rownames(proj@cellColData))
P8D3MI<-rownames(proj@cellColData)[P8D3MI]
P8D3Sham<-grep(".*-4",rownames(proj@cellColData))
P8D3Sham<-rownames(proj@cellColData)[P8D3Sham]

sample<-proj@cellColData$Sample
sample<-as.character(sample)

i=1
for (i in 1:37737){
  if(rownames(proj@cellColData)[i]%in%P1D3MI ){sample[i]="P1D3MI"}
  if(rownames(proj@cellColData)[i]%in%P1D3Sham ){sample[i]="P1D3Sham"}
  if(rownames(proj@cellColData)[i]%in%P8D3MI ){sample[i]="P8D3MI"}
  if(rownames(proj@cellColData)[i]%in%P8D3Sham ){sample[i]="P8D3Sham"}
  i=i+1
}
  P1D3MI P1D3Sham   P8D3MI P8D3Sham 
    8590     9504     8024    11619 
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-cellreport-data", load = FALSE)
proj<-loadArchRProject(path = "Save-cellreport-data")


proj <- filterDoublets(ArchRProj = proj)






proj2 <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10), 
    varFeatures = 25000, 
    dimsToUse = 1:8
)

####cluster###
proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

table(proj2$Clusters)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.1, 
  force = TRUE,
    metric = "cosine"
)
####minDist越小，分的越开####

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "CellReport-Plot-UMAP-Sample-Clusters-all.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

proj3 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 50, 
    minDist = 0.5, 
  force = TRUE,
    metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "CellReport-Plot-UMAP-Sample-Clusters-all.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-cellreport-all-filter", load = FALSE)
proj3<-loadArchRProject("Save-cellreport-all-filter")
#dev.off();
###downstream####
markersGS <- getMarkerFeatures(
    ArchRProj = proj3, 
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
plotPDF(heatmapGS, name = "CellReport-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)




idxSample <- BiocGenerics::which(proj3$Clusters %in% c("C1","C2","C3","C4","C5","C6",
  "C7","C9","C10","C11","C18","C12","C13","C14","C15","C16","C17"))
cellsSample <- proj3$cellNames[idxSample]
proj4=proj3[cellsSample, ]
proj4 <- addImputeWeights(proj4)

cluster<-proj4$Clusters
table(proj4$Clusters)

for (i in 1:33854){
if(cluster[i]=="C1"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C2"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C3"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C4"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C5"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C6"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C7"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C8"){cluster[i]="C8"}
if(cluster[i]=="C9"){cluster[i]="Epicardial"}
if(cluster[i]=="C10"){cluster[i]="Endothelial"}
if(cluster[i]=="C11"){cluster[i]="Endothelial"}
if(cluster[i]=="C12"){cluster[i]="Endothelial"}
if(cluster[i]=="C13"){cluster[i]="Endothelial"}
if(cluster[i]=="C14"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C15"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C16"){cluster[i]="Fibroblast"}
if(cluster[i]=="C17"){cluster[i]="Fibroblast"}
if(cluster[i]=="C18"){cluster[i]="Fibroblast"}
}
table(cluster)
proj4$celltype<-cluster
p3 <- plotEmbedding(ArchRProj = proj4, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj4, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj4, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
ggAlignPlots(p3, p4,p5, type = "h")
plotPDF(p3,p4,p5 name = "CellReport-addcelltype-UMAP-all.pdf", ArchRProj = proj4, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj4, overwrite = F,outputDirectory = "Save-cellreport-all-filter-addcelltype", load = FALSE)

##########extract fibroblast##########
idxSample <- BiocGenerics::which(proj4$celltype %in% c("Fibroblast"))
cellsSample <- proj4$cellNames[idxSample]
proj5=proj4[cellsSample, ]
proj5 <- addImputeWeights(proj5)
saveArchRProject(ArchRProj = proj5, overwrite = F,outputDirectory = "Save-cellreport-fibroblast", load = FALSE)

nproj <- addIterativeLSI(
    ArchRProj = proj5,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-FB", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(1.2), 
        sampleCells = 6000, 
        n.start = 10
    ), force = TRUE,
    varFeatures = 50000, 
    dimsToUse = 1:15
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-FB",
    method = "Seurat",
    #force = TRUE,
    name = "ClustersFB2",
  force = TRUE,
    resolution = 0.5)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-FB", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.05, 
  force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "ClustersFB2", embedding = "UMAP")

ggAlignPlots(p1, p2,p3, type = "h")
plotPDF(p1,p2,p3, name = "CellReport-fibroblast-UMAP-2.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

######rm C1 C2 #######

idxSample <- BiocGenerics::which(nproj$ClustersFB %in% c("C1","C2"))
cellsSample <- nproj$cellNames[-idxSample]
proj5=nproj[cellsSample, ]
#proj5 <- addImputeWeights(proj5)
saveArchRProject(ArchRProj = proj5, overwrite = F,outputDirectory = "Save-cellreport-fibroblast", load = FALSE)

######our top marker gene 在heatmap上的分布
ourFBproj<-loadArchRProject(path = "../Save-fibroblast2")

markersGS <- getMarkerFeatures(
    ArchRProj = ourFBproj, 
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
FBgene<-c(FB1,FB2,FB3,FB4)


markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "ClustersFB2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = FBgene,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "CellReport-FB-topmarker-GeneScores-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

cluster<-nproj$ClustersFB2
table(nproj$ClustersFB2)
for (i in 1:6913){
if(cluster[i]=="C1"){cluster[i]="F1"}
if(cluster[i]=="C2"){cluster[i]="F1"}
if(cluster[i]=="C3"){cluster[i]="F2"}
if(cluster[i]=="C4"){cluster[i]="F3"}
if(cluster[i]=="C5"){cluster[i]="F4"}
if(cluster[i]=="C6"){cluster[i]="F5"}
if(cluster[i]=="C7"){cluster[i]="F3"}
if(cluster[i]=="C8"){cluster[i]="F3"}
if(cluster[i]=="C9"){cluster[i]="F2"}
}
nproj$subtype<-cluster


markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = FBgene,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "CellReport-FB-topmarker-subtype", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

F1<-markerList$F1$name
F2<-markerList$F2$name
F3<-markerList$F3$name
F4<-markerList$F4$name
F5<-markerList$F5$name


#####GO analysis for F##########
library(clusterProfiler)
library(org.Mm.eg.db)

pdf("F-GO.pdf")
ego_geneID <- bitr(F5, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
###GO analysis####
ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"F5-BP.csv")





