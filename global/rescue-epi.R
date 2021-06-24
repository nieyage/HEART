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

Epimarker<-c("Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8","Msln",
  "Kcnj8","Adgre1","Itgal",####Pericyte
  "Plp1")
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = Epimarker, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)
plotPDF(plotList = p, 
    name = "Epimarker-Genes-umap2.pdf", 
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

saveArchRProject(ArchRProj = proj5, overwrite = FALSE,outputDirectory = "Save-allcelltype-rescue2", load = T)

####3 UMAP plot
pdf("allcelltype-3UMAP.pdf")
> levels(UMAP$celltype)
 [1] "B_cell"         "Cardiomyocytes" "Endothelial"    "Epicardial"    
 [5] "Fibroblast"     "Glial"          "granulocyte"    "Macrophages"   
 [9] "Pericyte"       "Smoothmuscle"   "T_cell" 
 colorvalue = c("#48453d","#1f441e","#0072B5","#fbaba9",  "#E18727","#b92419","#66E1E6","#20854E",  "#8edcaf","#7876B1","#0F3057")

#20854E    绿 MP
#0072B5    蓝 EC
#E18727    黄 FB
#7876B1    紫 SMC
#FCF5B3 米白  EPI 改成  #fbaba9
#1f441e 暗绿  Per 改成  #8edcaf
#66E1E6 浅蓝  GUA
#0F3057 藏蓝  T 
#ffcbcb 粉    B   改成  #48453d
#21DB00 嫩绿  GLI 改成  #b92419
########UMAP1 and UMAP2
df<-read.table("allcelltypeUMAP.txt")

ggplot(df, aes(x = UMAP1, y = UMAP2,color = celltype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = colorvalue) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("UMAP1 and UMAP2")

ggplot(df, aes(x = UMAP1, y = UMAP3,color = celltype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = colorvalue) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("UMAP1 and UMAP3")

ggplot(df, aes(x = UMAP2, y = UMAP3,color = celltype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = colorvalue) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("UMAP2 and UMAP3")

######3D UMAP######
 colorvalue = c("#48453d","#1f441e","#0072B5","#fbaba9",  "#E18727","#b92419","#66E1E6","#20854E",  "#8edcaf","#7876B1","#0F3057")
library(scatterplot3d)
pdf("3D-UMAP.pdf")
 
colors <- colorvalue[as.factor(UMAP$celltype)]

 s3d <- scatterplot3d(UMAP[,c("UMAP1","UMAP2","UMAP3")],
 color = colors,angle = 60, cex.symbols = 0.3,pch=16)

legend("bottom", legend = c("B","CM","EC","EPI","FB","Glial","Gran","MP","Per","SMC","T"),
col =colorvalue, pch = 15,cex=1,
inset = -0.2,xpd = TRUE, horiz = TRUE)



 s3d <- scatterplot3d(UMAP[,c("UMAP1","UMAP2","UMAP3")],
 color = colors,angle = 90, cex.symbols = 0.3,pch=16)

legend("bottom", legend = c("B","CM","EC","EPI","FB","Glial","Gran","MP","Per","SMC","T"),
col =colorvalue, pch = 15,cex=1,
inset = -0.2,xpd = TRUE, horiz = TRUE)

MV411_diff_gene <-subset(MV411res, padj < 0.01 & abs(log2FoldChange)> 1) 
MV411_diff<-rownames(diff_gene)
