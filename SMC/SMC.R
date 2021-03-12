#######提取Smooth muscle cells-cluster #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")

idxSample <- BiocGenerics::which(proj5$Clusters %in% c("C7","C8"))
cellsSample <- proj5$cellNames[idxSample]
nproj=proj5[cellsSample, ]
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-SMC", load = T)
nproj<-loadArchRProject(path = "Save-SMC")
#####add timepoint/treatment########
sample<-nproj$Sample
table(nproj$Sample)
for (i in 1:4109){
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

cellNames <- nproj$cellNames
nproj$timepoint<-sample

sample2<-nproj$Sample
table(nproj$Sample)
for (i in 1:4109){
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
nproj$treatment<-sample2


nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-dim15-SMC-1.2-80000-4", 
    iterations = 4, 
    clusterParams = list( 
        resolution = c(1.2), 
        sampleCells = 3000, 
        n.start = 10
    ), 
    varFeatures = 80000, 
    dimsToUse = 1:15
)

nproj<- addGeneIntegrationMatrix(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI-dim15-SMC-1.2-80000-4",
    seRNA = RNA,
    addToArrow = FALSE,
    groupRNA = "seurat_clusters",##seurat_clusters
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-dim15-SMC-1.2-80000-4",
    method = "Seurat",
    name = "Clusters2",
	force = TRUE,
    resolution = 1.2
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-dim15-SMC-1.2-80000-4", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.01, 
	force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "treatment", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")

ggAlignPlots(p1, p2,p3,p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "SMC-Plot-UMAP-6.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)




markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
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
"Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "SMC-Cluster-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

########gene score matrix####
GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)

genename<-rowData(GM)$name ####共24333个基因
gene<-GM@assays@data$GeneScoreMatrix
gene<-as.data.frame(gene)
gene<-cbind(genename,gene)
rownames(gene)<-gene[,1]
gene<-gene[,-1]
########distance tree########

coldata<-nproj@cellColData
C1<-coldata[which(coldata$Clusters2=="C1"),]
C1gene<-gene[,which(colnames(gene)%in% rownames(C1))]
C1<-rowMeans(C1gene)


C2<-coldata[which(coldata$Clusters2=="C2"),]
C2gene<-gene[,which(colnames(gene)%in% rownames(C2))]
C2<-rowMeans(C2gene)

C3<-coldata[which(coldata$Clusters2=="C3"),]
C3gene<-gene[,which(colnames(gene)%in% rownames(C3))]
C3<-rowMeans(C3gene)


C4<-coldata[which(coldata$Clusters2=="C4"),]
C4gene<-gene[,which(colnames(gene)%in% rownames(C4))]
C4<-rowMeans(C4gene)

C5<-coldata[which(coldata$Clusters2=="C5"),]
C5gene<-gene[,which(colnames(gene)%in% rownames(C5))]
C5<-rowMeans(C5gene)


C6<-coldata[which(coldata$Clusters2=="C6"),]
C6gene<-gene[,which(colnames(gene)%in% rownames(C6))]
C6<-rowMeans(C6gene)

C7<-coldata[which(coldata$Clusters2=="C7"),]
C7gene<-gene[,which(colnames(gene)%in% rownames(C7))]
C7<-rowMeans(C7gene)

C8<-coldata[which(coldata$Clusters2=="C8"),]
C8gene<-gene[,which(colnames(gene)%in% rownames(C8))]
C8<-rowMeans(C8gene)

C9<-coldata[which(coldata$Clusters2=="C9"),]
C9gene<-gene[,which(colnames(gene)%in% rownames(C9))]
C9<-rowMeans(C9gene)

C10<-coldata[which(coldata$Clusters2=="C10"),]
C10gene<-gene[,which(colnames(gene)%in% rownames(C10))]
C10<-rowMeans(C10gene)

C11<-coldata[which(coldata$Clusters2=="C11"),]
C11gene<-gene[,which(colnames(gene)%in% rownames(C11))]
C11<-rowMeans(C11gene)


C12<-coldata[which(coldata$Clusters2=="C12"),]
C12gene<-gene[,which(colnames(gene)%in% rownames(C12))]
C12<-rowMeans(C12gene)

all<-cbind(C1,C2)
all<-cbind(all,C3)
all<-cbind(all,C4)
all<-cbind(all,C5)
all<-cbind(all,C6)
all<-cbind(all,C7)
all<-cbind(all,C8)
all<-cbind(all,C9)
all<-cbind(all,C10)
all<-cbind(all,C11)
all<-cbind(all,C12)


tall<-t(all)
tdist=dist(tall,method="euclidean")
hc=hclust(tdist,method="complete")
plclust(hc)
plot(hc, hang=-1, xlab="");

#####根据以上结果分亚型#########
cluster<-nproj$Clusters2
table(nproj$Clusters2)

for (i in 1:4109){
if(cluster[i]=="C1"){cluster[i]="SMC4"}
if(cluster[i]=="C2"){cluster[i]="SMC4"}
if(cluster[i]=="C3"){cluster[i]="SMC4"}
if(cluster[i]=="C4"){cluster[i]="SMC4"}
if(cluster[i]=="C5"){cluster[i]="SMC2"}
if(cluster[i]=="C6"){cluster[i]="SMC2"}
if(cluster[i]=="C7"){cluster[i]="SMC3"}
if(cluster[i]=="C8"){cluster[i]="SMC3"}
if(cluster[i]=="C9"){cluster[i]="SMC3"}
if(cluster[i]=="C10"){cluster[i]="SMC1"}
if(cluster[i]=="C11"){cluster[i]="SMC1"}
if(cluster[i]=="C12"){cluster[i]="SMC1"}

}
table(cluster)
nproj$subtype<-cluster
cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(cluster),
    cells=cellNames, name = "subtype",force = TRUE)
#####calculate the cell number of CMC#######
sample<-nproj$Sample
subtype<-nproj$subtype
info<-cbind(sample,subtype)
info<-as.data.frame(info)
AR3<-info[which(info$sample=="AR3"),]
table(AR3$subtype)



#####根据subtype寻找特异性基因#########

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
"Tie1","Fabp4","Esam","Kdr","Tek","Eng",###endothelial
"Mmp2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
"Dclk1","Fkbp5","Akt3",##fibroblast
"Rgs5","Ano1","Acta2","Talgln","Mustn1","Myh11","Mylk")####smoothmuscle
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "SMC-subtype-Marker-Heatmap-0.5", width = 8, height = 4, ArchRProj = nproj, addDOC = FALSE)


####change the treatment color######

sample<-nproj$treatment
sample[1]=""
sample[2]="  "
sample[3]=" "
sample[4]="   "
sample[5]="    "
table(nproj$sample)
nproj$treatment2<-sample
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "treatment2", embedding = "UMAP")
ggAlignPlots(p3, type = "h")
plotPDF(p3, name = "SMC-UMAP-treatment-changecolor.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

#######find peak#######
pathToMacs2 <- findMacs2()
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "subtype")

nproj <- addReproduciblePeakSet(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(nproj)
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-SMC", load = T)
nproj <- addPeakMatrix(nproj)


markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", returnGR = TRUE)
markerList
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "SMC-Peak-Marker-Heatmap", width = 6, height = 4, ArchRProj = nproj, addDOC = FALSE)

nproj <- addMotifAnnotations(ArchRProj = nproj, force = TRUE,motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 30, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "SMC-Motifs-Enriched-Marker-Heatmap", width = 12, height = 6, ArchRProj = nproj, addDOC = FALSE)
  
####去掉同一个基因家族的TF，仅保留一个，添加尽可能多的基因家族进来

  idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% c("Klf4_143","Klf1_214","Klf3_804","Klf2_819","Klf8_194","Klf12_236",
           "Sp3_802","Smad5_883","Sp2_150","Sp5_238","Sp6_191","Sp4_167","Klf7_171","Klf5_145","E2f4_257","Klf14_237","Klf11_798","Klf13_817","Bach2_119",
		   "Tcf3_31","Tcf12_59","Ascl2_23","Klf15_806","Klf16_874","Sp7_222","Sp8_207","Sp9_231","Zfp219_816","Gata5_385","Gata4_386","Gata1_387","Mesp2_57",
		   "Gata2_383","Sox9_725","Egr2_188","Egr3_183","Zfp161_209","Zic5_195","Tcfap2c_3","Gata3_384","Zfp281_193",
		   "Zfp202_169","LINE3883_876","Ctcf_146","Zbtb1_182","E2f6_264","Id4_35","Fosb_98","Junb_98","Fosl1_107","Batf3_789",
		   "Nfe2l2_101","Bcl11b_814","LINE5785_293","Pou4f2_840","Pou4f1_841","LINE10076_380","LINE9919_344","LINE10042_370",
		   "Pou3f2_555","Pou3f4_561","Foxc2_316","Mef2c_638","Mef2d_842","Mef2b_882","Nfix_722","Nfib_859","Nfia_862","LINE9979_373","Tbpl2_773","Hic2_211",
		   "Junb_127","Jun_126","Sox12_741","LINE9832_335","LINE9950_357","LINE10046_367","LINE10077_381","LINE9891_343","LINE10016_369","LINE15940_622","Pou6f1_416",
		   "LINE6279_607","Foxs1_329","Irf2_636","Zfp148_800"))
cellsSample <- enrichMotifs@NAMES[-idxSample]
enrichMotifs1=enrichMotifs[cellsSample, ]


heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 30, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "SMC-Motifs-rmrep",width = 12, height = 6,ArchRProj = nproj, addDOC = FALSE)


#####GO analysis for SMC##########
library(clusterProfiler)
library(org.Mm.eg.db)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

SMC1<-markerList$SMC1$name
SMC2<-markerList$SMC2$name
SMC3<-markerList$SMC3$name
SMC4<-markerList$SMC4$name

#更改gene id###
pdf("SMC4-GO-KEGG.pdf")
ego_geneID <- bitr(SMC3, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
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
barplot(ego, showCategory=20)
write.csv(ego,"SMC3-BP.csv")

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.01,
                readable = TRUE)

ego
write.csv(ego,"SMC3-MF.csv")
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
write.csv(ego,"SMC3-CC.csv")

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

write.csv(ego,"SMC3-KEGG.csv")


barplot(ego, showCategory=20)


###SMC3/SMC4  AR NAR P###

idxSample <- BiocGenerics::which(nproj$subtype %in% c("SMC4"))
cellsSample <- nproj$cellNames[idxSample]
SMC4=nproj[cellsSample, ]
saveArchRProject(ArchRProj = SMC4, outputDirectory = "Save-SMC4", load = T)
SMC4<-loadArchRProject(path = "Save-SMC4")
#########treatment#######
markersGS <- getMarkerFeatures(
    ArchRProj = SMC4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "treatment",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
"Tie1","Fabp4","Esam","Kdr","Tek","Eng",###endothelial
"Mmp2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
"Dclk1","Fkbp5","Akt3",##fibroblast
"Rgs5","Ano1","Acta2","Talgln","Mustn1","Myh11","Mylk")####smoothmuscle
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.1")
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "SMC4-treatment-Marker-Heatmap-0.5", width = 8, height = 4, ArchRProj = SMC3, addDOC = FALSE)


#########treatment#######
markersGS <- getMarkerFeatures(
    ArchRProj = SMC4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "timepoint",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Mmp2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
"Dclk1","Fkbp5","Akt3",##fibroblast
"Rgs5","Ano1","Acta2","Talgln","Mustn1","Myh11","Mylk",
"Dst")####smoothmuscle
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.2", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.2")
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "SMC4-timepoint-Marker-Heatmap", width = 8, height = 4, ArchRProj = SMC4, addDOC = FALSE)
##

Day3<-markerList$Day3$name
Day7<-markerList$Day7$name
Day14<-markerList$Day14$name

#更改gene id###
pdf("SMC4-Day17-GO-KEGG.pdf")
ego_geneID <- bitr(Day7, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
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
barplot(ego, showCategory=20)

ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.01,
                readable = TRUE)

ego
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
barplot(ego, showCategory=20)





