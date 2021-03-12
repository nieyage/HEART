#######################SMC#############################
##################monocle##############################
#######提取Smooth muscle cells-cluster #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
ArchR.SMC<-loadArchRProject(path = "Save-SMC")


library(monocle)
library(ggplot2)
library(patchwork)
library(magrittr)
#Extract data, phenotype data, and feature data from the SeuratObject
GM<-getMatrixFromProject(
  ArchRProj =ArchR.SMC,
  useMatrix = "GeneScoreMatrix"
)
data<-assays(GM)$GeneScoreMatrix
data<-as(as.matrix(data), 'sparseMatrix')


gene<-rowData(GM)$name
rownames(data)<-gene
pd <- new('AnnotatedDataFrame', data = as.data.frame(ArchR.SMC@cellColData[,22:23]))


fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd)

################SMC subtype gene umap ##########
ArchR.SMC <- addImputeWeights(ArchR.SMC,reducedDims="IterativeLSI-dim15-SMC-1.2-50000-4")
Contraction<-c("Myh11","Acta2","Cnn1","Cnn2","Cnn3","Des","Vim","Smtn","Vcl","Cald1","Tagln")
Secreted<-c("Spp1","Ereg","Eln","Mgp","Thbs1","Thbs2")
p <- plotEmbedding(
    ArchRProj = ArchR.SMC, 
    colorBy = "GeneScoreMatrix", 
    name = Contraction, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(ArchR.SMC)
)
plotPDF(plotList = p, 
    name = "Contraction-Genes-umap.pdf", 
    ArchRProj = ArchR.SMC, 
    addDOC = FALSE, width = 5, height = 5)
#Rearrange for grid plotting
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
p3<-do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p3, 
    name = "Contraction-Genes-umap.pdf", 
    ArchRProj = ArchR.SMC, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotEmbedding(
    ArchRProj = ArchR.SMC, 
    colorBy = "GeneScoreMatrix", 
    name = Secreted, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(ArchR.SMC)
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
p3<-do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
pdf("Secreted-Genes-umap.pdf")
p3
dev.off()

################提取SMC１/2寻找差异##############


idxSample <- BiocGenerics::which(ArchR.SMC$subtype %in% c("SMC1","SMC2"))
cellsSample <- ArchR.SMC$cellNames[idxSample]
ArchR.SMC=ArchR.SMC[cellsSample, ]
saveArchRProject(ArchRProj = ArchR.SMC, outputDirectory = "Save-SMC12", load = T)
ArchR.SMC<-loadArchRProject(path = "Save-SMC12")


markersGS <- getMarkerFeatures(
    ArchRProj =ArchR.SMC, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.3")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.01,
                readable = TRUE)

dev.off()



#######重新把SMC12合在一起分析
#####根据以上结果分亚型#########
subtype<-ArchR.SMC$subtype

for (i in 1:4109){
if(subtype[i]=="SMC1"){subtype[i]="SMC1"}
if(subtype[i]=="SMC2"){subtype[i]="SMC1"}
if(subtype[i]=="SMC3"){subtype[i]="SMC3"}
if(subtype[i]=="SMC4"){subtype[i]="SMC2"}
}
table(subtype)
ArchR.SMC$subtype<-subtype
cellNames <- ArchR.SMC$cellNames
ArchR.SMC <- addCellColData(ArchRProj = ArchR.SMC, data = paste0(subtype),
    cells=cellNames, name = "subtype",force = TRUE)
saveArchRProject(ArchRProj =ArchR.SMC, outputDirectory = "Save-SMC-merge12", load = T)

ArchR.SMC<-loadArchRProject(path = "Save-SMC-merge12")

markersGS <- getMarkerFeatures(
    ArchRProj = ArchR.SMC, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
gene<-markerList$SMC1$name
library(clusterProfiler)
library(org.Mm.eg.db)
pdf("SMC3-GO-KEGG.pdf")
ego_geneID <- bitr(gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
###GO analysis####
ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.01,
                readable = TRUE)

ego
barplot(ego, showCategory=20,cex.axis=1)
#write.csv(ego,"SMC1-new-BP.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.01,
                readable = TRUE)

ego
#write.csv(ego,"SMC1-new-MF.csv")
barplot(ego, showCategory=20)
ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.01,
                readable = TRUE)

ego
#write.csv(ego,"SMC1-new-CC.csv")

barplot(ego, showCategory=20)

####KEGG analysis#######

ego <- enrichKEGG(
  gene = ego_geneID$ENTREZID,
  keyType = "kegg",
  organism  = 'mmu',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego

write.csv(ego,"SMC1-new-KEGG.csv")
barplot(ego, showCategory=20)
