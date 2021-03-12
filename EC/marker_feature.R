rm(list = ls())
getwd()
setwd("/public/home/Chenzh275/Data/scATAC/ArchR-Endothelial_Cell")

library(ArchR)
#ArchR.EC <- loadArchRProject(path = "~/Data/scATAC/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
#####exclude EC5 investegate difference
#ArchR.EC <- ArchR.EC[which(ArchR.EC$Subtype == "EC1" |ArchR.EC$Subtype == "EC2"|
#                             ArchR.EC$Subtype == "EC3"|ArchR.EC$Subtype == "EC4"),]
#saveArchRProject(ArchRProj = ArchR.EC , outputDirectory = "Save-Proj2-EC4-noEC5", load = T)
#ArchR.EC2 <- loadArchRProject(path = "~/Data/scATAC/ArchR-Endothelial_Cell/Save-Proj2-EC4-noEC5")
#table(ArchR.EC2$Subtype)
#table(ArchR.EC2$Clusters_3)
####Making Pseudo-bulk Replicates######
#ArchR.EC2 <- addGroupCoverages(ArchRProj = ArchR.EC2, groupBy = "ArchR.EC2$Subtype")

####Call peak#####
#table(ArchR.EC2@cellColData$Clusters_3)
#pathToMacs2 <- "/public/home/Chenzh275//miniconda3/bin/macs2"
#ArchR.EC2 <- addReproduciblePeakSet(ArchRProj = ArchR.EC2, groupBy = "ArchR.EC2$Subtype", pathToMacs2 = pathToMacs2, cutOff = 0.01)
#getPeakSet(ArchR.EC2)

#projHemeTmp2 <- addReproduciblePeakSet(ArchRProj = ArchR.EC2, groupBy = "ArchR.EC2$Subtype", peakMethod = "Tiles", method = "p", cutOff = 0.01)
#getPeakSet(projHemeTmp2)

######Comparing the two peak calling methods####
#length(subsetByOverlaps(getPeakSet(ArchR.EC2), getPeakSet(projHemeTmp2))) / length(getPeakSet(ArchR.EC2))
#length(subsetByOverlaps(getPeakSet(projHemeTmp2), getPeakSet(ArchR.EC2))) / length(getPeakSet(projHemeTmp2))
#length(subsetByOverlaps(resize(getPeakSet(ArchR.EC2), 1000, "center"), getPeakSet(projHemeTmp2))) / length(getPeakSet(ArchR.EC2))
#length(subsetByOverlaps(getPeakSet(projHemeTmp2), resize(getPeakSet(ArchR.EC2), 1000, "center"))) / length(getPeakSet(projHemeTmp2))

####add peak matrix####
#ArchR.EC3 <- addPeakMatrix(ArchR.EC2)
#getAvailableMatrices(ArchR.EC3)

#saveArchRProject(ArchRProj = ArchR.EC3 , outputDirectory = "Save-Proj2-EC4-noEC5", load = T)
ArchR.EC <- loadArchRProject(path = "~/Data/scATAC/ArchR-Endothelial_Cell/Save-Proj2-EC4-noEC5")
table(ArchR.EC$Subtype)
markersGS <- getMarkerFeatures(ArchRProj = ArchR.EC, useMatrix = "GeneScoreMatrix", groupBy = "ArchR.EC$Subtype",
                               bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersGS
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
#write.table(as.data.frame(markerList),"markersList-EC-gene-subtype-cluster_3-FC1_noEC5.txt",sep="\t",quote = FALSE)

markerGenes  <- c("CDH5","PECAM1","EDNRB","EGFL7","EMCN","EPAS1","FLT1","FABP4","ESAM","KDR","TEK","CD34")

####Heatmap for mark genes######
heatmapGS <-plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", labelMarkers = markerGenes, transpose = TRUE)
draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-ECs-subtype-Marker-Heatmap-FDR001-log2FC05-cluster_3_noEC5", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)

getAvailableMatrices(ArchR.EC)
markersPeaks <- getMarkerFeatures(ArchRProj = ArchR.EC, useMatrix = "PeakMatrix", 
                                  groupBy = "ArchR.EC$Subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerList
write.table(as.data.frame(markerList),"markersList-EC-peak-subtype--FDR005-log2FC05_noEC5.txt",sep="\t",quote = FALSE)

####marker peak#####
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5", transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "EC-Peak-Marker-Heatmap-subtype-cluster_3-FDR005_noEC5", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)

#####marker motif
enrichMotifs <- peakAnnoEnrichment( seMarker = markersPeaks, ArchRProj = ArchR.EC, peakAnnotation = "Motif",
                                    cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)