

#######提取SMC-cluster #####

library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")

proj2<-loadArchRProject(path = "Save-rmNAR")

idxSample <- BiocGenerics::which(proj2$celltype %in% c("Smoothmuscle"))
cellsSample <- proj2$cellNames[idxSample]
nproj=proj2[cellsSample, ]
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-rmNAR-SMC", load = T)
nproj<-loadArchRProject(path = "Save-rmNAR-SMC")
table(nproj$Sample)

nproj <- addIterativeLSI(
    ArchRProj = nproj,
    name = "IterativeLSI-rmNAR-SMC",
    useMatrix = "TileMatrix",
    iterations = 2, 
    clusterParams = list(resolution = c(1.2),
    sampleCells = 2000, n.start = 10), 
    varFeatures = 80000, 
    force = TRUE,
    dimsToUse = 1:15
)


####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-rmNAR-SMC",
    name = "Clusters2", force = TRUE, resolution = 1.5, dimsToUse = 1:15
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-rmNAR-SMC", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.1, 
    force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, size = 1,colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, size = 1,colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")

ggAlignPlots(p1, p4, type = "h")
plotPDF(p1, p4, name = "rmNAR-SMC-Plot-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "rmNAR-SMC-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)



cluster<-nproj$Clusters2
table(nproj$Clusters2)

for (i in 1:3206){
if(cluster[i]=="C1"){cluster[i]="SMC2"}
if(cluster[i]=="C2"){cluster[i]="SMC2"}
if(cluster[i]=="C3"){cluster[i]="SMC2"}
if(cluster[i]=="C4"){cluster[i]="SMC2"}
if(cluster[i]=="C5"){cluster[i]="SMC1"}
if(cluster[i]=="C6"){cluster[i]="SMC1"}
if(cluster[i]=="C7"){cluster[i]="SMC1"}
if(cluster[i]=="C8"){cluster[i]="SMC3"}
if(cluster[i]=="C9"){cluster[i]="SMC3"}
if(cluster[i]=="C10"){cluster[i]="SMC3"}
if(cluster[i]=="C11"){cluster[i]="SMC3"}
if(cluster[i]=="C12"){cluster[i]="SMC3"}
if(cluster[i]=="C13"){cluster[i]="SMC1"}
}
table(cluster)
nproj$subtype<-cluster
p3 <- plotEmbedding(ArchRProj = nproj, size=1,colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, size=1,colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = nproj, size=1,colorBy = "cellColData", name = "subtype", embedding = "UMAP")
ggAlignPlots(p3, p4,p5, type = "h")
plotPDF(p3,p4,p5, name = "rmNAR-addsubtype-UMAP-SMC.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-rmNAR-SMC", load = FALSE)
nproj<-loadArchRProject("Save-rmNAR-SMC")
table(nproj$subtype,nproj$Sample)

      AR14 AR3 AR7 P14  P3  P7
  SMC1  180  64 132 175 102 260
  SMC2  200  69 190 245  87 285
  SMC3  180 145 228 227 168 269
       AR14 AR3 AR7 P14  P3  P7
  SMC1  170  54  82 157  19  96
  SMC2  200  69 190 245  87 285
  SMC3  190 155 278 245 251 433


markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
Contraction<-c("Myh11","Acta2","Cnn1","Cnn2","Cnn3","Des","Vim","Smtn","Vcl","Cald1","Tagln")
Secreted<-c("Spp1","Ereg","Eln","Mgp","Thbs1","Thbs2")

markerGenes<-c(
  "Hspa1a","Hspa8","Dnajb14","Hspa1l","Hspa1b","Atf3","Fos",##SMC1
  "Myh11","Acta2","Cnn1","Cnn2","Cnn3","Des","Vim","Smtn","Vcl","Cald1","Tagln",#####SMC2
  "Spp1","Ereg","Eln","Mgp","Thbs1","Thbs2"####SMC3
  )

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "SMC-Marker-Heatmap", width = 8, height = 6, ArchRProj =nproj, addDOC = FALSE)


markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
SMC1<-markerList$SMC1$name
SMC2<-markerList$SMC2$name
SMC3<-markerList$SMC3$name




library(clusterProfiler)
library(org.Mm.eg.db)

pdf("/md01/nieyg/scATAC-ArchR/Save-rmNAR-SMC/Plots/SMC-GO.pdf")
ego_geneID <- bitr(SMC3, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
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
write.csv(ego,"/md01/nieyg/scATAC-ArchR/Save-rmNAR-SMC/Plots/SMC3-BP.csv")


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
plotPDF(heatmapPeaks, name = "rmNAR-SMC-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

nproj <- addMotifAnnotations(ArchRProj = nproj, force = TRUE,motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "SMC-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
  
