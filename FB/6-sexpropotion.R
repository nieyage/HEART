
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")
GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)

genename<-rowData(GM)$name ####共24333个基因
count<-GM@assays@data$GeneScoreMatrix
count<-as.matrix(count)
#count<-cbind(genename,count)
rownames(count)<-genename
Xist<-count[which(rownames(count)=="Xist"),]
AR14 =nproj@cellColData[which(nproj$Sample=="AR14"),]
NAR14=nproj@cellColData[which(nproj$Sample=="NAR14"),]
P14  =nproj@cellColData[which(nproj$Sample=="P14"),]
AR7  =nproj@cellColData[which(nproj$Sample=="AR7"),]
NAR7 =nproj@cellColData[which(nproj$Sample=="NAR7"),]
P7   =nproj@cellColData[which(nproj$Sample=="P7"),]
AR3  =nproj@cellColData[which(nproj$Sample=="AR3"),]
NAR3 =nproj@cellColData[which(nproj$Sample=="NAR3"),]
P3=nproj@cellColData[which(nproj$Sample=="P3"),]
xAR14<-Xist[which(names(Xist)%in%rownames(AR14))]
xNAR14<-Xist[which(names(Xist)%in%rownames(NAR14))]
xP14<-Xist[which(names(Xist)%in%rownames(P14))]
xAR7<-Xist[which(names(Xist)%in%rownames(AR7))]
xNAR7<-Xist[which(names(Xist)%in%rownames(NAR7))]
xP7<-Xist[which(names(Xist)%in%rownames(P7))]
xAR3<-Xist[which(names(Xist)%in%rownames(AR3))]
xNAR3<-Xist[which(names(Xist)%in%rownames(NAR3))]
xP3<-Xist[which(names(Xist)%in%rownames(P3))]
str(xAR3[which(xAR3==0.000)])
str(xAR3[which(xAR3>0.000)])
str(xAR3)


fAR14 <-Firre[which(names(Firre)%in%rownames(AR14))]
fNAR14<-Firre[which(names(Firre)%in%rownames(NAR14))]
fP14  <-Firre[which(names(Firre)%in%rownames(P14))]
fAR7  <-Firre[which(names(Firre)%in%rownames(AR7))]
fNAR7 <-Firre[which(names(Firre)%in%rownames(NAR7))]
fP7   <-Firre[which(names(Firre)%in%rownames(P7))]
fAR3  <-Firre[which(names(Firre)%in%rownames(AR3))]
fNAR3 <-Firre[which(names(Firre)%in%rownames(NAR3))]
fP3   <-Firre[which(names(Firre)%in%rownames(P3))]
str(fAR3[which(fAR3==0.000)])
str(fAR3[which(fAR3>0.000)])


tAR14 <-Tsix[which(names(Tsix)%in%rownames(AR14))]
tNAR14<-Tsix[which(names(Tsix)%in%rownames(NAR14))]
tP14  <-Tsix[which(names(Tsix)%in%rownames(P14))]
tAR7  <-Tsix[which(names(Tsix)%in%rownames(AR7))]
tNAR7 <-Tsix[which(names(Tsix)%in%rownames(NAR7))]
tP7   <-Tsix[which(names(Tsix)%in%rownames(P7))]
tAR3  <-Tsix[which(names(Tsix)%in%rownames(AR3))]
tNAR3 <-Tsix[which(names(Tsix)%in%rownames(NAR3))]
tP3   <-Tsix[which(names(Tsix)%in%rownames(P3))]
str(tP3[which(tP3==0.000)])
str(tP3[which(tP3>0.000)])
str(tAR3[which(tAR3==0.000)])
str(tAR3[which(tAR3>0.000)])
str(tNAR3[which(tNAR3==0.000)])
str(tNAR3[which(tNAR3>0.000)])

str(tP7[which(tP7==0.000)])
str(tP7[which(tP7>0.000)])
str(tAR7[which(tAR7==0.000)])
str(tAR7[which(tAR7>0.000)])
str(tNAR7[which(tNAR7==0.000)])
str(tNAR7[which(tNAR7>0.000)])

str(tP14[which(tP14==0.000)])
str(tP14[which(tP14>0.000)])
str(tAR14[which(tAR14==0.000)])
str(tAR14[which(tAR14>0.000)])
str(tNAR14[which(tNAR14==0.000)])
str(tNAR14[which(tNAR14>0.000)])

tAR14 <-u[which(names(u)%in%rownames(AR14))]
tNAR14<-u[which(names(u)%in%rownames(NAR14))]
tP14  <-u[which(names(u)%in%rownames(P14))]
tAR7  <-u[which(names(u)%in%rownames(AR7))]
tNAR7 <-u[which(names(u)%in%rownames(NAR7))]
tP7   <-u[which(names(u)%in%rownames(P7))]
tAR3  <-u[which(names(u)%in%rownames(AR3))]
tNAR3 <-u[which(names(u)%in%rownames(NAR3))]
tP3   <-u[which(names(u)%in%rownames(P3))]


####################CEBPD 调控通路 Nfkbia##########
#####Nfkbia的genescore以及 track#############
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = "Nfkbia", 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "Nfkbia-UMAP-Imputation.pdf", 
    ArchRProj = nproj,
    addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = "Nfkbia",
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )
plotPDF(plotList = p2, 
    name = "Nfkbia-violinplot.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

p3 <- plotBrowserTrack(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    geneSymbol = "Nfkbia", 
    upstream = 1000,
    downstream = 1000
)

plotPDF(plotList = p3, 
    name = "Nfkbia-Tracks-Marker-Genes-1000.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)
