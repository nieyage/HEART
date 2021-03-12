
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")
nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI-dim15-fibroblast-1.2-40000")

motifPositions <- getPositions(nproj)
markerMotifs <- "Cebpd_97"
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "subtype")
seFoot <- getFootprints(
  ArchRProj = nproj, 
  positions =motifPositions[markerMotifs], 
  groupBy = "subtype"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = nproj, 
  normMethod = "Subtract",
  plotName = "Cebpd-Footprints-Subtract-Bias",
  addDOC = FALSE,
  baseSize = 10,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = nproj, 
  normMethod = "Divide",
  plotName = "Cebpd-Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = nproj, 
  normMethod = "None",
  plotName = "Cebpd-Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5
)

GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)

PM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "PeakMatrix"
)

getAvailableMatrices(nproj)
TM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "TileMatrix",
  binarize = TRUE
)

m<-c(
"FOXA1","FOXA2","FOXA3","FOXB1","FOXC1","FOXC2","FOXD1","FOXD2","FOXD3","FOXD4","FOXD5","FOXD6","FOXE1","FOXE2","FOXE3","FOXF1","FOXF2","FOXG1","FOXH1",
"FOXI1","FOXJ1","FOXJ2","FOXJ3","FOXK1","FOXK2","FOXL1","FOXL2","FOXM1","FOXN1","FOXN2","FOXN3","FOXN4","FOXN5","FOXN6","FOXO1","FOXO2","FOXO3","FOXO4","FOXP1","FOXP2","FOXP3","FOXP4","FOXQ1")
 
p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = m,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )
plotPDF(plotList = p2, 
    name = "Foxl1motifgene-violinplot.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = m, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "Foxl1motifgene-UMAP-Imputation.pdf", 
    ArchRProj = nproj,
    addDOC = FALSE, width = 5, height = 5)






projHemeTmp <- addReproduciblePeakSet(ArchRProj = nproj,groupBy = "subtype",peakMethod = "Tiles",method = "p")

######peak ID#####
length(unique(paste(seqnames(peakSet), peakSet$idx, sep=":")))







findMotifsGenome.pl FB2-2-2.bed mm10 ./ -find foxl2.motif > ./foxl2.txt;
annotatePeaks.pl motif1-1.bed mm10 > motif1_peak.annotation.xls





library(clusterProfiler)
gene<-c("Pim1","Auh","Rhob","Cdkl4","Mdm4","Hspa5","Uhrf1bp1l","Pag1","Gm5079","Celf2","Phactr1","Dram1","Slc44a3","Nfkbia","Fabp12","Dusp6",
	"4930563F08Rik","D030024E09Rik","Plac8","1700016H13Rik","4930549C15Rik","Atp6v1b2","Ociad1","Mir21a")


gene.df <- bitr(gene, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Mm.eg.db)
head(gene.df)
ggo <- groupGO(gene=gene.df$ENTREZID,OrgDb= org.Mm.eg.db, ont= "BP",level= 3,readable = TRUE)

head(ggo)

mkk2 <- gseMKEGG(geneList = l, organism = 'mmu')
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 1)


m<-"Mki67"
p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = m,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )
plotPDF(plotList = p2, 
    name = "Mki67-violinplot.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = m, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "Mki67-UMAP-Imputation.pdf", 
    ArchRProj = nproj,
    addDOC = FALSE, width = 5, height = 5)

m<-getGroupSE(
       ArchRProj = nproj,
       useMatrix = "GeneScoreMatrix",
       groupBy = "genotype",
       divideN = TRUE,
       scaleTo = NULL,
       threads = getArchRThreads(),
       verbose = TRUE,
       logFile = createLogFile("getGroupSE")
     )


Dab2ip
Edem2
Herpud1
Herpud2
Faf2
Derl1
Hspa5
Dnajc3
Thbs4
Parp16
Yod1
Stt3b
Ficd
Creb3l4
Atf6
Erp27
Hspd1
Ern2
Hsph1
Creb3l3
Hsp90aa1
Hspa8
Xbp1
Edem1
Dnajb5
Hspa2
Crebrf
Fbxo6
Hspa9
Ern1
Upf3a
Comp
Derl2
Tmbim6
Upf3b
Upf2
Mfn2
Ddit3
Dnajb9
Tmem129
4931417E11Rik
Jkamp
Hspa4l
Erp44
Ube2j2
Hspa1b
Hspa1a
Hspa1l
Hspa13
Chac1
Hspb1
Hspa14
Edem3
Creb3l2
Hspa4l
