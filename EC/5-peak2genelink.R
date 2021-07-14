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
p2g <- getPeak2GeneLinks(
    ArchRProj = ECproj,
    corCutOff = 0.1,
    resolution = 1,
    FDRCutOff = 1e-02,
    returnLoops = FALSE,
    varCutOffATAC =0,
    varCutOffRNA = 0
)
markerGenes  <- c("Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Tcf4","Col4a2")
markerGenes<-c("Plcg1","Acvrl1","Flt1","Itgb3","Itgb1bp1")
Hsp<-c("Hspa8","Hsp90aa1","Hspa1b","Hspa1a")
p <- plotBrowserTrack(
    ArchRProj = ECproj, 
    groupBy = "Subtype", 
    geneSymbol =Hsp, 
    upstream = 10000,
    downstream = 10000,
    loops = getPeak2GeneLinks(ECproj)
)
grid::grid.newpage()
grid::grid.draw(p$Col4a2)

plotPDF(plotList = p, 
    name = "Hsp-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = ECproj, 
    addDOC = FALSE, width = 5, height = 5)
name<-c("Col4a2","Col4a2","Col4a2","Col4a2","Col4a2","Col4a2","Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Cd36","Tcf4")
end<-c(11304445,11222954,11300052,11283713,11309666,11200837,86094825,73919309,122812121,
  5510515,36981495,17880314,17774600,69388580)
start<-c(11303945,11222454,11299552,11283213,11309166,11200337,86094325,73918809,122811621,
  5510015,36980995,17879814,17774100,69388080)

peak<-GRanges(
  seqnames=c("chr8","chr8","chr8","chr8","chr8","chr8","chr16","chr17","chr6","chr18","chr12","chr5","chr5","chr18"),
  ranges=IRanges(start=start,end=end,names = name ))
findOverlaps(position[98],peak)
