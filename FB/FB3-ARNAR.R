
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")


idxSample <- BiocGenerics::which(nproj$subtype %in% c("FB3"))
cellsSample <- nproj$cellNames[idxSample]
FB3proj=nproj[cellsSample, ]
saveArchRProject(ArchRProj = FB3proj, outputDirectory = "Save-fibroblast-FB3", load = T)
FB3proj<-loadArchRProject(path = "Save-fibroblast-FB3")
#####根据AR NAR P分组########

markersGS <- getMarkerFeatures(
    ArchRProj = FB3proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "treatment",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
p<-plotMarkerHeatmap(
       seMarker = markersGS,
       cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
       transpose = TRUE,
       labelRows = TRUE
     )

plotPDF(p, name = "FB3-treatment-FDR005-Heatmap", width = 12, height = 6, ArchRProj = FB3proj, addDOC = FALSE)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.5 & Log2FC >= 0.3")

library(clusterProfiler)
library(org.Mm.eg.db)
ego_geneID <- bitr(nar, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
pdf("FB3-NAR-GOandKEGG.pdf")
ego <- enrichKEGG(
  gene = ego_geneID$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"FB3-NAR-KEGG.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"FB3-NAR-BP.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"FB3-NAR-MF.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"FB3-NAR-CC.csv")


#######提取AR，P/NAR，P寻找差异######

idxSample <- BiocGenerics::which(FB3proj$treatment %in% c("AR","P"))
cellsSample <- FB3proj$cellNames[idxSample]
ARproj=FB3proj[cellsSample, ]


ARmarkersGS <- getMarkerFeatures(
    ArchRProj = ARproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "treatment",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
p<-plotMarkerHeatmap(
       seMarker = ARmarkersGS,
       cutOff = "FDR <= 0.1 & Log2FC >= 0",
       labelRows = TRUE
     )
plotPDF(p, name = "FB3-ARvsP-FDR01-Heatmap", width = 4, height = 12, ArchRProj = ARproj, addDOC = FALSE)

idxSample <- BiocGenerics::which(FB3proj$treatment %in% c("NAR","P"))
cellsSample <- FB3proj$cellNames[idxSample]
NARproj=FB3proj[cellsSample, ]



NARmarkersGS <- getMarkerFeatures(
    ArchRProj = NARproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "treatment",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
p<-plotMarkerHeatmap(
       seMarker = NARmarkersGS,
       cutOff = "FDR <= 0.1 & Log2FC >= 0",
       labelRows = TRUE
     )
plotPDF(p, name = "FB3-NARvsP-FDR01-Heatmap", width = 4, height = 12, ArchRProj = NARproj, addDOC = FALSE)
