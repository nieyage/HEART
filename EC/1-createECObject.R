setwd("/public/home/Chenzh275/Data/scATAC/ArchR-Endothelial_Cell")
library(ArchR)
ArchR <- loadArchRProject(path = "/public/home/nieyg/scATAC-ArchR/Save-Proj5")
ArchR.EC <- ArchR[which(ArchR@cellColData$Clusters == "C14" |ArchR@cellColData$Clusters == "C15" |ArchR@cellColData$Clusters == "C16" 
                         |ArchR@cellColData$Clusters == "C17" |ArchR@cellColData$Clusters == "C18" |ArchR@cellColData$Clusters == "C19" |
                           ArchR@cellColData$Clusters == "C20" |ArchR@cellColData$Clusters == "C21" ),]
ArchR.EC <- addIterativeLSI( ArchRProj = ArchR.EC,useMatrix = "TileMatrix",name = "IterativeLSI", iterations = 2, 
  clusterParams = list(resolution = c(1),sampleCells = 10000, n.start = 10), 
  varFeatures = 30000, dimsToUse = 1:15)

####add cluster####
ArchR.EC <- addClusters(input = ArchR.EC, reducedDims = "IterativeLSI",name = "Clusters_3", force = TRUE, resolution = 1, dimsToUse = 1:15)
ArchR.EC <- addUMAP(ArchRProj = ArchR.EC, reducedDims = "IterativeLSI", nNeighbors = 10, minDist = 0.3, metric = "cosine", force = TRUE)
#####define genotype####
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("Sample","Genotype");
meta[,1] = ArchR.EC$Sample
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
ArchR.EC$Genotype  = meta$Genotype
table(ArchR.EC$Genotype)

markersGS <- getMarkerFeatures(ArchRProj = ArchR.EC, useMatrix = "GeneScoreMatrix", groupBy = "Clusters_3",
                               bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersGS
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
markerGenes  <- c("CDH5","PECAM1","EDNRB","EGFL7","EMCN","EPAS1","FLT1","FABP4","ESAM","KDR","TEK","CD34")
####Heatmap for marker genes######
heatmapGS <-plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", labelMarkers = markerGenes, transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-ECs-Marker-Heatmap-log2FC05-umap02-cluster_3", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)

####Visualization of marker genes######
ArchR.EC <- addImputeWeights(ArchR.EC)
p <- plotEmbedding(ArchRProj = ArchR.EC, colorBy = "GeneScoreMatrix", name = markerGenes, 
                   embedding = "UMAP",imputeWeights = getImputeWeights(ArchR.EC))
p3 <- lapply(p, function(x){x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p3))
plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-EC-W-Imputation-umap2-cluster_3.pdf", ArchRProj = ArchR.EC, addDOC = FALSE, width = 5, height = 5)

####Making Pseudo-bulk Replicates######
ArchR.EC <- addGroupCoverages(ArchRProj = ArchR.EC, groupBy = "Clusters_3")
####Call peak#####
pathToMacs2 <- "/public/home/Chenzh275//miniconda3/bin/macs2"
ArchR.EC <- addReproduciblePeakSet(ArchRProj = ArchR.EC, groupBy = "Clusters_3", pathToMacs2 = pathToMacs2, cutOff = 0.01)
getPeakSet(ArchR.EC)

####add peak matrix####
ArchR.EC2 <- addPeakMatrix(ArchR.EC)
getAvailableMatrices(ArchR.EC2)
saveArchRProject(ArchRProj = ArchR.EC2, outputDirectory = "Save-Proj2-EC2", load = FALSE)
markersPeaks <- getMarkerFeatures(ArchRProj = ArchR.EC, useMatrix = "PeakMatrix", 
                                  groupBy = "Clusters_3", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

####marker peak#####
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "EC-Peak-Marker-Heatmap-umap02-cluster_3", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)
#####define subtype####
table(ArchR.EC$Clusters_3)
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("cluster","subtype");
meta[,1] = ArchR.EC$Clusters_3
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "C4" |cl[i] == "C5" |cl[i] == "C6" |cl[i] == "C7" |cl[i] == "C10"  )
    meta[i,2] = "EC1"
  else if (cl[i] == "C2"|cl[i] == "C3" )
    meta[i,2]  = "EC2"
  else if (cl[i] == "C8" )
    meta[i,2]  = "EC3"
  else if (cl[i] == "C1" | cl[i] == "C12" |cl[i] == "C13")
    meta[i,2]  = "EC4"
  else if (cl[i] == "C9" )
    meta[i,2]  = "EC5"
  else 
    meta[i,2]  = "EC6"
}
meta = data.frame(meta)
ArchR.EC$Subtype  = meta$subtype
p6 <- plotEmbedding(ArchRProj = ArchR.EC, colorBy = "cellColData", name = "Subtype", embedding = "UMAP")
plotPDF(p6, name = "UMAP-EC-Subtype-Clusters_3-iter2.pdf", ArchRProj = ArchR.EC, addDOC = FALSE, width = 5, height = 5)

####define EC5##
markerGenes <- c("CD34","FLT1","KDR","PROM1","PTPRC","ENG","VCAM1","KIT",
                 "KDR","MCAM","PECAM1","CDH5","NOS3" ,"VWF","PLVAP" )
####Visualization of EPC EC marker genes######
ArchR.EC <- addImputeWeights(ArchR.EC)
p <- plotEmbedding(ArchRProj = ArchR.EC, colorBy = "GeneScoreMatrix", name = markerGenes, 
                   embedding = "UMAP",imputeWeights = getImputeWeights(ArchR.EC))
p3 <- lapply(p, function(x){x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 3),p3))
plotPDF(plotList = p, name = "Plot-UMAP-EPC-EC-Marker-Genes-EC-umap02-cluster_3.pdf", 
        ArchRProj = ArchR.EC, addDOC = FALSE, width = 5, height = 5)
plotPDF(plotList = p, name = "Plot-UMAP-EPC-EC-Marker-Genes-EC-W-Imputation-umap2-cluster_3.pdf", ArchRProj = ArchR.EC, addDOC = FALSE, width = 5, height = 5)

#####define genotype####
table(ArchR.EC$Sample)
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("sample","genotype");
meta[,1] = ArchR.EC$Sample
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "AR3" |cl[i] == "AR7" |cl[i] == "AR14"  )
    meta[i,2] = "AR"
  else if (cl[i] == "NAR3"|cl[i] == "NAR7" |cl[i] == "NAR14")
    meta[i,2]  = "NAR"
  else
    meta[i,2]  = "P"
}
head(meta)
meta = data.frame(meta)
ArchR.EC$Genotype  = meta$genotype
table(ArchR.EC$Genotype)
p6 <- plotEmbedding(ArchRProj = ArchR.EC, colorBy = "cellColData", name = "Genotype", embedding = "UMAP")
#####define timepoint####
table(ArchR.EC$Sample)
meta<-matrix(data=NA,nrow=28055,ncol=2);
colnames(meta)=c("sample","timepoint");
meta[,1] = ArchR.EC$Sample
cl <- meta[,1]
for ( i in 1:length(cl))
{
  if (cl[i] == "AR3" |cl[i] == "NAR3" |cl[i] == "P3"  )
    meta[i,2] = "Day3"
  else if (cl[i] == "AR7"|cl[i] == "NAR7" |cl[i] == "P7")
    meta[i,2]  = "Day7"
  else
    meta[i,2]  = "Day14"
}
head(meta)
meta = data.frame(meta)
ArchR.EC$Timepoint  = meta$timepoint
table(ArchR.EC$Timepoint)
p6 <- plotEmbedding(ArchRProj = ArchR.EC, colorBy = "cellColData", name = "Timepoint", embedding = "UMAP")

#####re draw umap 
df = data.frame(ArchR.EC@embeddings$UMAP$df, ArchR.EC$Genotype)
df = data.frame(ArchR.EC@embeddings$UMAP$df, ArchR.EC$Timepoint)
df = data.frame(ArchR.EC@embeddings$UMAP$df, ArchR.EC$Subtype)
Timepoint = ArchR.EC$Timepoint
Subtype = ArchR.EC$Subtype
Sample = ArchR.EC$Genotype
###timepoint
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = Timepoint)) +
  geom_point(size=0.8,alpha=0.5) + 
  theme_bw() + 
  scale_color_manual(values = c("#FFA500","#0073C2FF", "#20854EFF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Timepoint")

##subtype umap
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = Subtype)) +
  geom_point(size=0.8,alpha=0.5) + 
  theme_bw() + 
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF", "#E18727FF","#20854EFF","#7876B1FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Subtype")

###sample
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2,color = Sample)) +
  geom_point(size=0.8,alpha=0.5) + 
  theme_bw() + 
  scale_color_manual(values = c("#E64B35B2","#00468BB2","#42B540B2" )) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Genotype")





