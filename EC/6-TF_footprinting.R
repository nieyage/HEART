
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
#NonEC3 VSEC3



setwd("/public/home/Chenzh275/Data/scATAC/ArchR-Endothelial_Cell/new-analysis-20211221/EC3vsNonEC3")
####Making Pseudo-bulk Replicates######



EC4proj2 <- loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Endothelial_Cell/Save-Proj3-EC5-s
	ubtype");

EC3_NonEC3 <- EC4proj2$Subtype

for (i in 1:25532){
	if(EC3_NonEC3[i]%in%c("EC1","EC2","EC4")){EC3_NonEC3[i]="NonEC3"}
}
EC4proj2$EC3_NonEC3<-EC3_NonEC3
library(BSgenome.Mmusculus.UCSC.mm10)


####Call peak#####
pathToMacs2 <- "/public/home/Chenzh275/miniconda3/bin/macs2"
EC4proj2 <- addReproduciblePeakSet(ArchRProj = EC4proj2, groupBy = "EC4proj2$EC3_NonEC3",peaksPerCell = 1000, 
 pathToMacs2 = pathToMacs2, cutOff = 0.01)
getPeakSet(EC4proj2)

EC4proj2 <- addPeakMatrix(EC4proj2)
getAvailableMatrices(EC4proj2)

markersPeaksN3 <- getMarkerFeatures(ArchRProj = EC4proj2, useMatrix = "PeakMatrix", groupBy = "EC4proj2$EC3_NonEC3", 
 bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "binomial",binarize = TRUE)
markersPeaksN3
markerListN3 <- getMarkers(markersPeaksN3, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerListN3$EC3

write.table(markerListN3$EC3[,-2],"EC3-FDR005log2FC05.bed",quote = F,sep = "\t",rownames=F)

findMotifsGenome.pl EC3-FDR005log2FC05.bed mm10 EC3-FDR005log2FC05-peak_motif -len 8,10,12

EC3_peak<-makeGRangesFromDataFrame(markerListN3$EC3)



motifPositions <- getPositions(EC4proj2)
motifs <- c("AP1", "BATF", "BACH", "ATF","FOS","Fra")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), ignore.case = TRUE, value = TRUE)))

markerMotifs
EC4proj2 <- addGroupCoverages(ArchRProj = EC4proj2, maxCells = 1000, minCells = 40, 
 minReplicates = 3, groupBy = "EC4proj2$EC3_NonEC3",force = TRUE)

motif_peak<-intersect(EC3_peak,motifPositions[markerMotifs[2]][[1]],ignore.strand=T);


seFoot <- getFootprints(
  ArchRProj = EC4proj2, 
  positions = motifPositions[markerMotifs], 
  groupBy = "EC4proj2$EC3_NonEC3"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = EC4proj2, 
  normMethod = "subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
);

plotFootprints(
  seFoot = seFoot,
  ArchRProj = EC4proj2, 
  normMethod = "divide",
  plotName = "Footprints-divide-Bias-10",
  addDOC = FALSE,
  smoothWindow = 10
);

plotFootprints(
  seFoot = seFoot,
  ArchRProj = EC4proj2, 
  normMethod = "none",
  plotName = "Footprints-none-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

markersPeaksN3 <- getMarkerFeatures(ArchRProj = EC4proj2, useMatrix = "PeakMatrix", groupBy = "EC4proj2$EC3_nonEC3", 
 bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "binomial",binarize = TRUE)
markersPeaksN3
markerListN3 <- getMarkers(markersPeaksN3, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerListN3$EC3
write.table(markerListN3,"markerList-Peak-NonEC3-vs-EC3-FDR005log2FC05.txt",quote = F,sep = "\t")
heatmapEMN3 <- plotMarkerHeatmap(markersPeaksN3, transpose = F,returnMat = T, cutOff = "FDR <= 0.05 & Log2FC >= 0.5", plotLog2FC = TRUE)
head(heatmapEMN3)




####marker motif 
EC4proj2 <- addMotifAnnotations(ArchRProj = EC4proj2, motifSet = "homer", name = "Motif_N3")
enrichMotifsN3 <- peakAnnoEnrichment(seMarker = markersPeaksN3,  ArchRProj = EC4proj2,
                                   peakAnnotation = "Motif_N3",  cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
enrichMotifsN3
heatmapEMN3 <- plotEnrichHeatmap(enrichMotifsN3, transpose = F,returnMat = T)






