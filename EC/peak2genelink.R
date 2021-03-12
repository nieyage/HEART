#######peak2genelink
#######Creating Low-Overlapping Aggregates of Cells###############
#######Co-accessibility with ArchR################################
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
ECproj<-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
#####Peak2Peak corealation##########
#ECproj <- addCoAccessibility(
#    ArchRProj = ECproj,
#    reducedDims = "IterativeLSI"
#)
#cA <- getCoAccessibility(
#    ArchRProj = ECproj,
#    corCutOff = 0.5,
#    resolution = 1000,
#    returnLoops = TRUE
#)

#######Peak2GeneLinkage with ArchR
ECproj <- addPeak2GeneLinks(
    ArchRProj = ECproj,corCutOff = 0.1,dimsToUse = 1:20,useMatrix = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI"
)

#addPeak2GeneLinks(
#  ArchRProj = NULL,
#  reducedDims = "IterativeLSI",
#  scaleDims = NULL,
#  cellsToUse = NULL,
#  k = 100,
#  knnIteration = 500,
#  overlapCutoff = 0.8,
#  maxDist = 250000,
#  scaleTo = 10^4,
#  log2Norm = TRUE,
#  predictionCutoff = 0.4,
#  addEmpiricalPval = FALSE,
#  threads = max(floor(getArchRThreads()/2), 1),
#  verbose = TRUE
#)
p2g <- getPeak2GeneLinks(
    ArchRProj = ECproj,
    corCutOff = 0.1,
    resolution = 1,
    FDRCutOff = 1e-02,
    returnLoops = FALSE,
    varCutOffATAC =0,
    varCutOffRNA = 0
)

#getPeak2GeneLinks(
#  ArchRProj = NULL,
#  corCutOff = 0.45,
#  varCutOffATAC = 0.25,
#  varCutOffRNA = 0.25,
#  resolution = 1,
#  returnLoops = TRUE
#)


markerGenes  <- c("Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Tcf4","Col4a2")

p <- plotBrowserTrack(
    ArchRProj = ECproj, 
    groupBy = "Subtype", 
    geneSymbol ="Col4a2", 
    upstream = 10000,
    downstream = 150000,
    loops = getPeak2GeneLinks(ECproj)
)
grid::grid.newpage()
grid::grid.draw(p$Col4a2)

plotPDF(plotList = p, 
    name = "Anti-angio-Genes-with-Peak2GeneLinks-1.pdf", 
    ArchRProj = ECproj, 
    addDOC = FALSE, width = 5, height = 5)
####Col4a2#########
chr8  11303945  11304445  + Col4a2
chr8  11222454  11222954  + Col4a2
chr8  11299552  11300052  + Col4a2
chr8  11283213  11283713  + Col4a2
chr8  11309166  11309666  + Col4a2
chr8  11200337  11200837  + Col4a2
chr16 86094325  86094825  + Adamts1
chr17 73918809  73919309  + Xdh
chr6  122811621 122812121 + Foxj2
chr18 5510015   5510515   + Zeb1
chr12 36980995  36981495  + Meox2
chr5  17879814  17880314  + Cd36
chr5  17774100  17774600  + Cd36
chr18 69388080  69388580  + Tcf4

name<-c("Col4a2","Col4a2","Col4a2","Col4a2","Col4a2","Col4a2","Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Cd36","Tcf4")
end<-c(11304445,11222954,11300052,11283713,11309666,11200837,86094825,73919309,122812121,
  5510515,36981495,17880314,17774600,69388580)
start<-c(11303945,11222454,11299552,11283213,11309166,11200337,86094325,73918809,122811621,
  5510015,36980995,17879814,17774100,69388080)

peak<-GRanges(
  seqnames=c("chr8","chr8","chr8","chr8","chr8","chr8","chr16","chr17","chr6","chr18","chr12","chr5","chr5","chr18"),
  ranges=IRanges(start=start,end=end,names = name ))
TFid
 [1]  98 101 104 107 108 113 117 119 126 127 132 135 789 790 871

findOverlaps(position[98],peak)

Fosb/Fos/Fosl1/Bach1/Nfe2/Jund
Nfe2l2/Fos

Fosb/Fos/Fosl1/Bach1/Bach2/Jun/Junb/Jund
Nfe2l2/Fos/Nfe2/Jund/Batf3/Batf

Fosb/Nfe2l2/Fos/Fosl1/Bach1/Bach2/Jun/Junb/Nfe2/Jund






################EC4¡¢FB3¡¢MP3¡¢SMC1 distance tree ###############


p<-getGroupSE(
  ArchRProj = ECproj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "ECproj$Subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = ECproj,
  useMatrix = "GeneScoreMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
tsubcount = t(count)
tdist=dist(tsubcount,method="euclidean")
hc=hclust(tdist,method="complete")
plclust(hc)
pdf("EC-subtype-distance.pdf")
plot(hc, hang=-1, xlab="");

FBproj<-loadArchRProject(path = "Save-fibroblast2")
MPproj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")
SMCproj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")



p<-getGroupSE(
  ArchRProj = FBproj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = FBproj,
  useMatrix = "GeneScoreMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
tsubcount = t(count)
tdist=dist(tsubcount,method="euclidean")
hc=hclust(tdist,method="complete")
plclust(hc)
pdf("FB-subtype-distance.pdf")
plot(hc, hang=-1, xlab="");
dev.off()


p<-getGroupSE(
  ArchRProj = MPproj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "MPproj$Subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = MPproj,
  useMatrix = "GeneScoreMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
tsubcount = t(count)
tdist=dist(tsubcount,method="euclidean")
hc=hclust(tdist,method="complete")
plclust(hc)
pdf("MP-subtype-distance.pdf")
plot(hc, hang=-1, xlab="");
dev.off()


p<-getGroupSE(
  ArchRProj = SMCproj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = SMCproj,
  useMatrix = "GeneScoreMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
tsubcount = t(count)
tdist=dist(tsubcount,method="euclidean")
hc=hclust(tdist,method="complete")
plclust(hc)
pdf("SMC-subtype-distance.pdf")
plot(hc, hang=-1, xlab="");
dev.off()

########ÌáÈ¡promoter-proximalÇøÓòµÄpeak##############################


p<-getGroupSE(
  ArchRProj = ECproj,
  useMatrix = "PeakMatrix",
  groupBy = "ECproj$Subtype"
)
PM<-getMatrixFromProject(
  ArchRProj = ECproj,
  useMatrix = "PeakMatrix"
)


peakinfo<-rowData(p)
peak<-assays(p)$PeakMatrix

EC<-read.csv("EC4_promoter-proximal.CSV")

ECpeak<-GRanges(
	seqnames=EC$Chr,
	ranges=IRanges(start=EC$Start,end=EC$End,strand=EC$Strand),names = EC$Gene.Name)

allpeak<-GRanges(
	seqnames=peakinfo$seqnames,
	ranges=IRanges(start=peakinfo$start,end=peakinfo$end),names=rownames(peakinfo)
	)
hits<-findOverlaps(allpeak,ECpeak)@ from
ECpromoterid<-allpeak[hits,]$names
matrix<-peak[ECpromoterid,]

pdf("EC-promoter-proximal-peak-change-with-EC2.pdf")
plot(density(matrix[,2]), col="blue", lty=1, 
     xlab = "Peak openness", main = "EC2 and EC4 promoter-proximal peak openness distribution  ")
lines(density(matrix[,4]), col="red", lty=1)

FC<-matrix[,4]/matrix[,2]
plot(density(FC), col="green", lty=1, 
     xlab = "Peak openness Fold Change (EC4/EC2)", main = "Peak openness")

EC2<-as.numeric(matrix[,2])
EC4<-as.numeric(matrix[,4])
EC2<-ecdf(EC2)
EC4<-ecdf(EC4)

plot(EC2,col="blue")
lines(EC4, col = "red")
ks.test(EC2,FB1)

dev.off()
