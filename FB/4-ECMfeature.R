#####find special ECM gene for subtype ####
####use archR to find overlap ####

library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.5 & Log2FC >= 0.3")
mart<-read.csv("mart_export.csv")
ECMgene<-mart[which(mart$GO.term.name=="extracellular matrix"),]
ECMgene<-ECMgene$Gene.name

FB1<-markerList$FB1$name
FB1ECM<-intersect(FB1,ECMgene)
FB1ECM

FB2<-markerList$FB2$name
FB2ECM<-intersect(FB2,ECMgene)
FB2ECM

FB3<-markerList$FB3$name
FB3ECM<-intersect(FB3,ECMgene)
FB3ECM

FB4<-markerList$FB4$name
FB4ECM<-intersect(FB4,ECMgene)
FB4ECM
g<-c(FB1ECM,FB2ECM,FB3ECM,FB4ECM)
count1<-subcount[which(rownames(subcount)%in%g),]
count1=t(scale(t(count1),scale = T,center = T))
pdf("ECM-overlap-ArchR.pdf",width=6,height=12)
list<-pheatmap(c, cluster_rows=F, 
               clustering_method="ward.D2",
               show_rownames=T,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)
dev.off()
