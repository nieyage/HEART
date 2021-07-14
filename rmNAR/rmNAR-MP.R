

#######提取MP-cluster #####

library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")

proj2<-loadArchRProject(path = "Save-rmNAR")

idxSample <- BiocGenerics::which(proj2$celltype %in% c("ImmuneCells"))
cellsSample <- proj2$cellNames[idxSample]
nproj=proj2[cellsSample, ]
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-rmNAR-MP", load = T)
nproj<-loadArchRProject(path = "Save-rmNAR-MP")
table(nproj$Sample)

nproj <- addIterativeLSI(
    ArchRProj = nproj,
    name = "IterativeLSI-rmNAR-MP",
    useMatrix = "TileMatrix",
    iterations = 2, 
    clusterParams = list(resolution = c(1),
    sampleCells = 1500, n.start = 10), 
    varFeatures = 50000, 
    force = TRUE,
    dimsToUse = 1:15
)


####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-rmNAR-MP",
    name = "Clusters2", force = TRUE, resolution = 1, dimsToUse = 1:15
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-rmNAR-MP", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, size = 1,colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, size = 1,colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")

ggAlignPlots(p1, p4, type = "h")
plotPDF(p1, p4, name = "rmNAR-MP-Plot-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

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
plotPDF(heatmapGS, name = "rmNAR-MP-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

####macrophage ameker
markerGenes1  <- c("FCGR1", "ADGRE1", "CD14", "CSF1R","CD74", "CD163",
                  "CD68","ITGAM","LGALS3", "LYZ1","MRC1") ####Macrophages
markerGenes2 <- c("CD80","CD86","FCGR1","FCGR3","TLR2") ####M1 macrophages
markerGenes3 <- c("MRC1", "CD163","ARG1","CLEC7A","RETNLB", "CHIL3","CD68")####M2 macrophages

###other immune cell
markerGenes4 <- c("CD3E","CD3D", "CD8A", "CD8B1","NKG7","IGFBP4","LAT")### T cells
markerGenes5 <- c("CD79A","CD79B", "MZB1","LY6D")### B cells
markerGenes6 <- c("S100A9", "S100A8") #Granulocyte
markerGenes7 <- c("CD74","CD83","CD86","FLT3","CD209A")#####DC

markerGenes<-c(markerGenes1,markerGenes2,markerGenes3,markerGenes4,markerGenes5,markerGenes6,markerGenes7)
####Heatmap for mark genes######
heatmapGS <-plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = markerGenes, transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-ICs-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)



idxSample <- BiocGenerics::which(nproj$Clusters2 %in% c("C4","C5","C6","C7","C8"))
cellsSample <- nproj$cellNames[idxSample]
nproj=nproj[cellsSample, ]
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-rmNAR-MP", load = T)
nproj<-loadArchRProject(path = "Save-rmNAR-MP")









cluster<-nproj$Clusters2
table(nproj$Clusters2)

for (i in 1:3206){
if(cluster[i]=="C1"){cluster[i]="MP2"}
if(cluster[i]=="C2"){cluster[i]="MP2"}
if(cluster[i]=="C3"){cluster[i]="MP2"}
if(cluster[i]=="C4"){cluster[i]="MP2"}
if(cluster[i]=="C5"){cluster[i]="MP1"}
if(cluster[i]=="C6"){cluster[i]="MP1"}
if(cluster[i]=="C7"){cluster[i]="MP1"}
if(cluster[i]=="C8"){cluster[i]="MP3"}
if(cluster[i]=="C9"){cluster[i]="MP3"}
if(cluster[i]=="C10"){cluster[i]="MP3"}
if(cluster[i]=="C11"){cluster[i]="MP3"}
if(cluster[i]=="C12"){cluster[i]="MP3"}
if(cluster[i]=="C13"){cluster[i]="MP1"}
}
table(cluster)
nproj$subtype<-cluster
p3 <- plotEmbedding(ArchRProj = nproj, size=1,colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, size=1,colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = nproj, size=1,colorBy = "cellColData", name = "subtype", embedding = "UMAP")
ggAlignPlots(p3, p4,p5, type = "h")
plotPDF(p3,p4,p5, name = "rmNAR-addsubtype-UMAP-MP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-rmNAR-MP", load = FALSE)
nproj<-loadArchRProject("Save-rmNAR-MP")
table(nproj$subtype,nproj$Sample)

      AR14 AR3 AR7 P14  P3  P7
  MP1  180  64 132 175 102 260
  MP2  200  69 190 245  87 285
  MP3  180 145 228 227 168 269
       AR14 AR3 AR7 P14  P3  P7
  MP1  170  54  82 157  19  96
  MP2  200  69 190 245  87 285
  MP3  190 155 278 245 251 433


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
  "Hspa1a","Hspa8","Dnajb14","Hspa1l","Hspa1b","Atf3","Fos",##MP1
  "Myh11","Acta2","Cnn1","Cnn2","Cnn3","Des","Vim","Smtn","Vcl","Cald1","Tagln",#####MP2
  "Spp1","Ereg","Eln","Mgp","Thbs1","Thbs2"####MP3
  )

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "MP-Marker-Heatmap", width = 8, height = 6, ArchRProj =nproj, addDOC = FALSE)


markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
MP1<-markerList$MP1$name
MP2<-markerList$MP2$name
MP3<-markerList$MP3$name




library(clusterProfiler)
library(org.Mm.eg.db)

pdf("/md01/nieyg/scATAC-ArchR/Save-rmNAR-MP/Plots/MP-GO.pdf")
ego_geneID <- bitr(MP3, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
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
write.csv(ego,"/md01/nieyg/scATAC-ArchR/Save-rmNAR-MP/Plots/MP3-BP.csv")


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
plotPDF(heatmapPeaks, name = "rmNAR-MP-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

nproj <- addMotifAnnotations(ArchRProj = nproj, force = TRUE,motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "MP-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
  
