#######给全部cluster添加注释名，并根据细胞类型分类 #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")
######添加celltype信息###
cluster<-proj5$Clusters
table(proj5$Clusters)

for (i in 1:60420){
if(cluster[i]=="C1"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C2"){cluster[i]="Fibroblast"}
if(cluster[i]=="C3"){cluster[i]="Fibroblast"}
if(cluster[i]=="C4"){cluster[i]="Fibroblast"}
if(cluster[i]=="C5"){cluster[i]="Fibroblast"}
if(cluster[i]=="C6"){cluster[i]="Fibroblast"}
if(cluster[i]=="C7"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C8"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C9"){cluster[i]="Fibroblast"}
if(cluster[i]=="C10"){cluster[i]="Fibroblast"}
if(cluster[i]=="C11"){cluster[i]="Fibroblast"}
if(cluster[i]=="C12"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C13"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C14"){cluster[i]="Endothelial"}
if(cluster[i]=="C15"){cluster[i]="Endothelial"}
if(cluster[i]=="C16"){cluster[i]="Endothelial"}
if(cluster[i]=="C17"){cluster[i]="Endothelial"}
if(cluster[i]=="C18"){cluster[i]="Endothelial"}
if(cluster[i]=="C19"){cluster[i]="Endothelial"}
if(cluster[i]=="C20"){cluster[i]="Endothelial"}
if(cluster[i]=="C21"){cluster[i]="Endothelial"}
}
proj5$celltype<-cluster
cellNames <- proj5$cellNames
proj5 <- addCellColData(ArchRProj = proj5, data = paste0(cluster),
    cells=cellNames, name = "celltype",force = TRUE)

#########add timepoint/treatment#####
sample<-proj5$Sample
table(proj5$Sample)
for (i in 1:60420){
if(sample[i]=="AR3"){sample[i]="Day3"}
if(sample[i]=="AR7"){sample[i]="Day7"}
if(sample[i]=="AR14"){sample[i]="Day14"}
if(sample[i]=="NAR3"){sample[i]="Day3"}
if(sample[i]=="NAR7"){sample[i]="Day7"}
if(sample[i]=="NAR14"){sample[i]="Day14"}
if(sample[i]=="P3"){sample[i]="Day3"}
if(sample[i]=="P7"){sample[i]="Day7"}
if(sample[i]=="P14"){sample[i]="Day14"}
}
proj5$timepoint<-sample
proj5 <- addCellColData(ArchRProj = proj5, data = paste0(sample),
    cells=cellNames, name = "timepoint",force = TRUE)

sample2<-proj5$Sample
table(proj5$Sample)
for (i in 1:60420){
if(sample2[i]=="AR3"){sample2[i]="AR"}
if(sample2[i]=="AR7"){sample2[i]="AR"}
if(sample2[i]=="AR14"){sample2[i]="AR"}
if(sample2[i]=="NAR3"){sample2[i]="NAR"}
if(sample2[i]=="NAR7"){sample2[i]="NAR"}
if(sample2[i]=="NAR14"){sample2[i]="NAR"}
if(sample2[i]=="P3"){sample2[i]="P"}
if(sample2[i]=="P7"){sample2[i]="P"}
if(sample2[i]=="P14"){sample2[i]="P"}
}
proj5$treatment<-sample2
proj5 <- addCellColData(ArchRProj = proj5, data = paste0(sample2),
    cells=cellNames, name = "treatment",force = TRUE)
	table(proj5$treatment)

p1 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "celltype2", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "treatment", embedding = "UMAP")

ggAlignPlots(p1, p2, p3, type = "h")
plotPDF(p1,p2,p3, name = "New-UMAP-celltype-timepoint-treatment.pdf", ArchRProj = proj5, addDOC = FALSE, width = 5, height = 5)

#####add marker gene######
####have other genes###
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
"Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
"Tie1","Fabp4","Esam","Kdr","Tek","Eng",###endothelial
"Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
"Dclk1","Fkbp5","Akt3",##fibroblast
"Rgs5","Ano1","Acta2",####smoothmuscle
"Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8",#####epicardial
"Pdgfrb","Abcc9",#Pericyte
"Pde3b","Ghr","Sik2","Pparg","Pde8b","Mgst1","Mapk10","Mme","Lama4",
"Fasn","Igf1","Dhrs3","Adra1a","Adh1b","Rtn4","Ptger3",#####Adipocyte
"Mbnl1","Fn1","Pde4d",####VSMC
"Scn7a","Bcl2",####Neuronal
"Ccnd3","Ptprc"#####Lymphocyte
)

markersGS <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "New-GeneScores-Heatmap-haveother", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
####save marker list to GO ananlysis####
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
write.csv(markerList,"allsample-markerList.csv")

###find marker peak###
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "celltype-Peak-Heatmap", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif")

####marker TF#####

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 8, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "celltype-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
  
if("Motif" %ni% names(proj5@peakAnnotation)){
    proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif")
}
proj5 <- addBgdPeaks(proj5)
proj5 <- addDeviationsMatrix(
  ArchRProj = proj5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj5, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "celltype-Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj5, addDOC = FALSE)
