
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
nproj <- addImputeWeights(nproj)

meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("Sample","Genotype");
meta[,1] = nproj$Sample
for ( i in 1:nrow(meta))
{
  if (meta[i,1]  == "AR3" |meta[i,1]  == "AR7" |meta[i,1]  == "AR14" )
    meta[i,2] = "AR"
  else if (meta[i,1]  == "NAR3"|meta[i,1]  == "NAR7" |meta[i,1]  == "NAR14")
    meta[i,2]  = "NAR"
  else
    meta[i,2]  = "P"
}
head(meta)
meta = data.frame(meta)
nproj$Genotype  = meta$Genotype
table(nproj$Genotype)

cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(subtype),
    cells=cellNames, name = "Subtype",force = TRUE)
nproj <- addCellColData(ArchRProj = nproj, data = paste0(meta$Genotype),
    cells=cellNames, name = "Genotype",force = TRUE)
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-Proj2-EC5-subtype", load = T)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
vein<-c("Mgp","Cfh","Apoe","Cpe","Cytl1","Bgn","Dcn","Plvap","Ctsh","Vwf","Rbp1","Fabp5","Npr3","H19","Vcam1","Tmem108","Tm4sf1","Id2",
	"Tmem176b","Hmcn1","Ptgs1","Ccdc80","F2r","Igfbp4","Il6st","Ctla2a","Emcn","Cd59a","Clu","Ahsg","Spint2","Smoc1","Vim","Comp","Tmem176a",
	"Mgst1","Plagl1","Lum","Igf2","Gpr182","Colec11","Cxcl16","Tmem2","Tmem100","Cdh11","Ablim1","Bmx","Ramp2","Pdlim3","Ednrb")

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = vein,
  clusterCols = FALSE,
  transpose = FALSE
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.3")

heatmapGS1 <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = marker,
  transpose = FALSE
)
ComplexHeatmap::draw(heatmapGS1, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS1, name = "Endothelial-Heatmap-forlabel",width=6,height=12,ArchRProj = nproj, addDOC = FALSE)

##########提取 marker gene matrix ############
heatmapGS<-plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  labelMarkers = NULL,
  returnMat = TRUE,
  transpose = FALSE
)
pdf("EC-heatmap-blueYellow.pdf",width=6,height=12)
pheatmap(heatmapGS,cluster_cols = F,cluster_rows = F,
         color = paletteContinuous(set = "blueYellow", n = 256, reverse = FALSE),
         #cellwidth = 15, cellheight = 12,
         show_rownames=F,show_colnames=T)
dev.off()


#####peak matrix######
paletteContinuous(set = "solarExtra", n = 256, reverse = FALSE)
#####motif matrix#####
paletteContinuous(set = "comet", n = 100)
#####gene matrix######
paletteContinuous(set = "blueYellow", n = 256, reverse = FALSE)
######################
####针对EC各个subtype进行富集分析#######
library(clusterProfiler)
library(org.Mm.eg.db)
#####for EC1/2/3/4/5############
ego_geneID <- bitr(ec5, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
pdf("EC5-GOandKEGG.pdf")
ego <- enrichKEGG(
  gene = ego_geneID$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"EC5-KEGG.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"EC5-BP.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"EC5-MF.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
barplot(ego, showCategory=20)
write.table(ego,"EC5-CC.csv")

dev.off()


###########find specific TFs#######


markersPeaks <- getMarkerFeatures(ArchRProj = ArchR.EC, useMatrix = "PeakMatrix", 
                                  groupBy = "Subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList

heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "EC-Peak-Marker-Heatmap-subtype", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

nproj <- addMotifAnnotations(ArchRProj = nproj, motifSet = "cisbp", name = "Motif", force = TRUE)
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,  ArchRProj = nproj, 
                                   peakAnnotation = "Motif",  cutOff = "FDR <= 0.5 & Log2FC >= 0.1")

enrichMotifs
df <- data.frame(TF = rownames(enrichMotifs), mlog10Padj = assay(enrichMotifs)[,3])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = FALSE,returnMatrix = TRUE,clusterCols = FALSE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "EC-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
                          EC1         EC2         EC3 EC4        EC5
Fosl1_107 (22)     100.000000  96.7975573 61.79774024   0   0.000000
Smarcc1_843 (146)  100.000000  23.4282762 16.86956006   0   2.605330
Bach1_108 (88)     100.000000  33.2041436 24.40305300   0   0.985043
Fosb_98 (76)       100.000000  22.7001481 23.93248757   0   2.056988
Jund_135 (58)      100.000000  42.6564663 33.07381465   0   0.000000
Nfkb1_701 (44)     100.000000   1.6896638  0.00000000   0   0.000000
Rel_697 (43)       100.000000   2.3425148  0.00000000   0   0.000000
Nfe2_132 (38)      100.000000  26.7617292 18.62431940   0   5.588061
Batf_790 (25)      100.000000  45.1838727 21.18939555   0   3.389119
Tcf21_79 (60)        0.000000 100.0000000  0.00000000   0   0.000000
Ascl2_23 (44)        0.000000 100.0000000  0.00000000   0   0.000000
Zfp238_227 (43)      0.000000 100.0000000  0.00000000   0   0.000000
Myod1_24 (38)        0.000000 100.0000000  0.00000000   0   0.000000
Elf1_285 (38)       41.528391 100.0000000  0.00000000   0   0.000000
Tcfap4_20 (32)       0.000000 100.0000000  0.00000000   0   0.000000
Etv2_270 (29)        2.321461 100.0000000  0.00000000   0   0.000000
Sox7_746 (29)       27.042936 100.0000000 19.38503496   0   0.000000
Nhlh2_84 (28)        0.000000 100.0000000  0.00000000   0   0.000000
Sp2_150 (245)        0.000000   0.0000000  0.00000000 100   0.000000
Wt1_872 (217)        0.000000   0.0000000  0.00000000 100   0.000000
Klf7_171 (181)       0.000000   0.0000000  0.00000000 100   0.000000
Zfp148_800 (175)     0.000000   0.0000000  0.00000000 100   0.000000
Tcfap2d_5 (162)      0.000000   0.0000000  0.00000000 100   0.000000
Smad1_863 (154)      0.000000   0.0000000  0.00000000 100   0.000000
Klf12_236 (144)      0.000000   0.0000000  0.00000000 100   0.000000
E2f4_257 (136)       0.000000   0.0000000  0.00000000 100   0.000000
Mbd2_641 (136)       0.000000   0.0000000  0.00000000 100   0.000000
Gata2_383 (141)      0.000000   0.0000000  2.65255504   0 100.000000
Foxd1_330 (92)       0.000000   0.2135617  0.00000000   0 100.000000
LINE9878_334 (92)    0.000000   0.2135617  0.00000000   0 100.000000
LINE10003_359 (92)   0.000000   0.2135617  0.00000000   0 100.000000
Nfatc2_700 (88)      0.000000   2.0824957  3.55071087   0 100.000000
Nfic_724 (78)        0.000000   0.0000000  0.00000000   0 100.000000
LINE9932_332 (54)    0.000000   5.3290367  0.00000000   0 100.000000
LINE10066_355 (54)   0.000000   5.3290367  0.00000000   0 100.000000
Foxc2_316 (53)       0.000000   3.9581782  0.00000000   0 100.000000
Foxa1_306 (51)       0.000000   6.9682013  0.06770915   0 100.000000

pdf("EC-motif-enrich.pdf")
pheatmap(o,cluster_cols = F,cluster_rows = F,
         color = paletteContinuous(set = "comet", n = 100),
         cellwidth = 12, cellheight = 12,
         show_rownames=T,show_colnames=T)
dev.off()

#########TF-Gene-Score##########

markerGenes  <- c("Fosl1","Smarcc1","Bach1","Fosb","Jund","Nfe2","Batf","Sox7")

p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "EC-TFgeneScore-UMAP-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)



if("Motif" %ni% names(nproj@peakAnnotation)){
    nproj <- addMotifAnnotations(ArchRProj = nproj, motifSet = "cisbp", name = "Motif")
}
nproj <- addBgdPeaks(nproj)
nproj <- addDeviationsMatrix(
  ArchRProj = nproj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(nproj, name = "MotifMatrix", plot = TRUE)
motifs <- c("Fosl1","Smarcc1","Bach1","Fosb","Jund","Nfe2","Batf","Sox7")
markerMotifs <- getFeatures(nproj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs


p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)
plotPDF(plotList = p, 
    name = "EC3-chromVAR-devscore-umap.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 10, height = 3)



motifPositions <- getPositions(nproj)
motifs <- c("Fosl1","Smarcc1","Bach1","Fosb","Jund","Nfe2","Batf","Sox7")

markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

markerMotifs
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "Subtype")
seFoot <- getFootprints(
  ArchRProj = nproj, 
  positions =motifPositions[markerMotifs], 
  groupBy = "Subtype"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = nproj, 
  normMethod = "Subtract",
  plotName = "EC3-motif-Footprints-Subtract-Bias",
  addDOC = FALSE,
  baseSize = 10,
  smoothWindow = 5
)