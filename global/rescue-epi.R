####rescue-epi######
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast")
proj5<-loadArchRProject(path = "Save-allcelltype")
nproj2<-loadArchRProject(path = "Save-fibroblast2")
all<-rownames(proj5@cellColData)[which(proj5$celltype=="Fibroblast")]
FB<-rownames(nproj2@cellColData)
EP<-setdiff(all,FB)
nproj=proj5[EP, ]
saveArchRProject(ArchRProj = nproj, overwrite = FALSE,outputDirectory = "Save-EPi", load = T)
nproj<-loadArchRProject(path = "Save-EPi")

nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-Epi", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(1.2), 
        sampleCells = 300, 
        n.start = 10
    ), 
    varFeatures = 30000, 
    dimsToUse = 1:15
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-Epi",
    method = "Seurat",
    name = "Clusters2",
	force = TRUE,
    resolution = 1.2
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-Epi", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
	force = TRUE,
    metric = "cosine"
)
nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI-Epi")

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1,p4,p3, type = "h")
plotPDF(p1,p4, p3,name = "Epi-Plot-UMAP-2.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

Epimarker<-c("Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8","Msln")
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = Epimarker, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)
plotPDF(plotList = p, 
    name = "Epimarker-Genes-umap.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)


labelOld <- levels(nproj$Clusters2)
labelOld
#  [1] "Cluster11" "Cluster2"  "Cluster12" "Cluster1"  "Cluster8"  "Cluster4" 
#  [7] "Cluster9"  "Cluster5"  "Cluster7"  "Cluster14" "Cluster3"  "Cluster10"
# [13] "Cluster6"  "Cluster13"
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
remapClust <- c(
    "C1" = "sub1",
    "C2" = "sub1",
    "C3" = "sub2",
    "C4" = "sub3",
    "C5" = "sub3",
    "C6" = "sub3"
)
nproj$Clusters3 <- mapLabels(nproj$Clusters2, newLabels = remapClust, oldLabels = names(remapClust))





marker<-c("Kcnj8","Adgre1","Itgal",####Pericyte
	"Plp1")###Glial
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = marker, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)
plotPDF(plotList = p, 
    name = "othermarker-Genes-umap.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)


Glial<-rownames(nproj@cellColData)[which(nproj$Clusters3=="sub1")]
Pericyte<-rownames(nproj@cellColData)[which(nproj$Clusters3=="sub2")]
Epi<-rownames(nproj@cellColData)[which(nproj$Clusters3=="sub3")]

ArchR.IC <- loadArchRProject(path = "/public/home/Chenzh275/Data/scATAC/ArchR-Immune_Cell/Save-Proj1-IC3")
meta<-matrix(data=NA,nrow=4114,ncol=2);
colnames(meta)=c("Subtype","subtypename");
meta[,1] = ArchR.IC$Subtype
cl <- meta[,1]
for ( i in 1:nrow(meta))
{
  if (meta[i,1]  == "IC1" |meta[i,1]  == "IC2" |meta[i,1]  == "IC3" )
    meta[i,2] = "Macrophages"
  else if (meta[i,1]  == "IC4")
    meta[i,2]  = "granulocyte"
  else if (meta[i,1]  == "IC5")
    meta[i,2]  = "T_cell"

  else
    meta[i,2]  = "B_cell"
}
head(meta)
meta = data.frame(meta)
ArchR.IC$subtypename  = meta$subtypename
table(ArchR.IC$subtypename)

Gra_Barcode<- rownames(ArchR.IC[which(ArchR.IC$subtypename =="granulocyte"),])
T_Barcode<- rownames(ArchR.IC[which(ArchR.IC$subtypename =="T_cell"),])
B_Barcode<- rownames(ArchR.IC[which(ArchR.IC$subtypename =="B_cell"),])
MP_Barcode<- rownames(ArchR.IC[which(ArchR.IC$subtypename =="Macrophages"),])
IC_Barcode <- c(MP_Barcode,Gra_Barcode,T_Barcode,B_Barcode)


meta<-matrix(data=NA,nrow=60420,ncol=3);
colnames(meta)=c("Subtype","cellbarcode","subtypename");
meta[,1] = proj5$celltype
meta[,2] = rownames(proj5@cellColData)
#cl <- meta[,1]
for ( i in 1:nrow(meta))
{
  if (meta[i,2] %in% Gra_Barcode  )
    meta[i,1] = "granulocyte"
  if (meta[i,2] %in% T_Barcode)
    meta[i,1]  = "T_cell"
  if (meta[i,2] %in% B_Barcode)
    meta[i,1]  = "B_cell"  
  if (meta[i,2] %in% MP_Barcode)
    meta[i,1]  = "Macrophages"
  if (meta[i,2] %in% Glial)
    meta[i,1]  = "Glial"
  if (meta[i,2] %in% Pericyte)
    meta[i,1]  = "Pericyte"
  if (meta[i,2] %in% Epi)
    meta[i,1]  = "Epicardial"
}
meta = data.frame(meta)
proj5$subtypename  = meta$Subtype
table(proj5$subtypename)

nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI-Epi")
 B_cell Cardiomyocytes    Endothelial     Epicardial     Fibroblast 
           334           2521          28055            405          20831 
         Glial    granulocyte    Macrophages       Pericyte   Smoothmuscle 
           241            214           3367            144           4109 
        T_cell 
           199 
"#350b40","#1f441e","#121013"
p1 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", 
	name = "subtypename", 
	#col=c("#350b40","#ffd3b4","#0B70AC","#ffcc29","#D78427","#da7f8f","#1f441e","#1F7E4D","#3d84b8","#7876AF","#121013"),
	embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggAlignPlots(p1,p4,p3, type = "h")
plotPDF(p1,p4, p3,name = "allcelltype-rescue-UMAP.pdf", ArchRProj = nproj, 	addDOC = FALSE, width = 5, height = 5)


pdf("allcelltype-rescuecells-UMAP.pdf")
celltype = proj5$subtypename
df = data.frame(proj5@embeddings$UMAP$df, proj5$subtypename)
colnames(df)<-c("UMAP1","UMAP2","celltype")
###sample
ggplot(df, aes(x = UMAP1, y = UMAP2,color = celltype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#94d0cc","#1f441e","#0B70AC","#ffcc29","#D78427","#da7f8f","#ffd3b4","#1F7E4D","#3d84b8","#7876AF","#81b214")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("celltype")
dev.off();

saveArchRProject(ArchRProj = proj5, overwrite = FALSE,outputDirectory = "Save-allcelltype-rescue", load = T)
