
library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 
addArchRGenome("mm10")
library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 
addArchRGenome("mm10")

inputFiles <- c("GSE153479_fragments.tsv.gz")
names(inputFiles)<-c("sample")

inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,   ####fragment files or bam files 
  sampleNames = names(inputFiles),
  filterTSS = 6, #Dont set this too high because you can always increase later
  filterFrags = 2800, 
  addTileMat = TRUE,
  force=TRUE,
  addGeneScoreMat = TRUE
)



doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "downstream",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
#########QC plot##########
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "cellreport-TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1,p2, name = "cellreport-QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

getAvailableMatrices(proj)

Name            Sample        Sample Index    Cell Barcode Suffix
P1D3MI          P1+3 dpi        SI-NA-G1              -1
P1D3Sham        P1+3 dps        SI-NA-H1              -2
P8D3MI          P8+3 dpi        SI-NA-A2              -3
P8D3Sham        P8+3 dps        SI-NA-B2              -4

P1D3MI<-grep(".*-1",rownames(proj@cellColData))
P1D3MI<-rownames(proj@cellColData)[P1D3MI]
P1D3Sham<-grep(".*-2",rownames(proj@cellColData))
P1D3Sham<-rownames(proj@cellColData)[P1D3Sham]

P8D3MI<-grep(".*-3",rownames(proj@cellColData))
P8D3MI<-rownames(proj@cellColData)[P8D3MI]
P8D3Sham<-grep(".*-4",rownames(proj@cellColData))
P8D3Sham<-rownames(proj@cellColData)[P8D3Sham]

sample<-proj@cellColData$Sample
sample<-as.character(sample)

i=1
for (i in 1:37737){
  if(rownames(proj@cellColData)[i]%in%P1D3MI ){sample[i]="P1D3MI"}
  if(rownames(proj@cellColData)[i]%in%P1D3Sham ){sample[i]="P1D3Sham"}
  if(rownames(proj@cellColData)[i]%in%P8D3MI ){sample[i]="P8D3MI"}
  if(rownames(proj@cellColData)[i]%in%P8D3Sham ){sample[i]="P8D3Sham"}
  i=i+1
}
  P1D3MI P1D3Sham   P8D3MI P8D3Sham 
    8590     9504     8024    11619 
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-cellreport-data", load = FALSE)
proj<-loadArchRProject(path = "Save-cellreport-data")


proj <- filterDoublets(ArchRProj = proj)






proj2 <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10), 
    varFeatures = 25000, 
    dimsToUse = 1:8
)

####cluster###
proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

table(proj2$Clusters)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.1, 
  force = TRUE,
    metric = "cosine"
)
####minDist越小，分的越开####

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "CellReport-Plot-UMAP-Sample-Clusters-all.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

proj3 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 50, 
    minDist = 0.5, 
  force = TRUE,
    metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p3,p4, name = "CellReport-Plot-UMAP-Sample-Clusters-all.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-cellreport-all-filter", load = FALSE)
proj3<-loadArchRProject("Save-cellreport-all-filter")
#dev.off();
###downstream####
markersGS <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
"Itgam","Lgals3","Lyzl1","Mrc1",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1","Tie1","Fabp4","Esam","Kdr","Tek",###endothelial
"Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra",##fibroblast
"Abcc9","Rgs5","Ano1","Acta2",####smoothmuscle
"Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8"#####epicardial
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "CellReport-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)




idxSample <- BiocGenerics::which(proj3$Clusters %in% c("C1","C2","C3","C4","C5","C6",
  "C7","C9","C10","C11","C18","C12","C13","C14","C15","C16","C17"))
cellsSample <- proj3$cellNames[idxSample]
proj4=proj3[cellsSample, ]
proj4 <- addImputeWeights(proj4)

cluster<-proj4$Clusters
table(proj4$Clusters)

for (i in 1:33854){
if(cluster[i]=="C1"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C2"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C3"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C4"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C5"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C6"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C7"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C8"){cluster[i]="C8"}
if(cluster[i]=="C9"){cluster[i]="Epicardial"}
if(cluster[i]=="C10"){cluster[i]="Endothelial"}
if(cluster[i]=="C11"){cluster[i]="Endothelial"}
if(cluster[i]=="C12"){cluster[i]="Endothelial"}
if(cluster[i]=="C13"){cluster[i]="Endothelial"}
if(cluster[i]=="C14"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C15"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C16"){cluster[i]="Fibroblast"}
if(cluster[i]=="C17"){cluster[i]="Fibroblast"}
if(cluster[i]=="C18"){cluster[i]="Fibroblast"}
}
table(cluster)
proj4$celltype<-cluster
p3 <- plotEmbedding(ArchRProj = proj4, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj4, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj4, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
ggAlignPlots(p3, p4,p5, type = "h")
plotPDF(p3,p4,p5 name = "CellReport-addcelltype-UMAP-all.pdf", ArchRProj = proj4, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj4, overwrite = F,outputDirectory = "Save-cellreport-all-filter-addcelltype", load = FALSE)

##########extract fibroblast##########
idxSample <- BiocGenerics::which(proj4$celltype %in% c("Fibroblast"))
cellsSample <- proj4$cellNames[idxSample]
proj5=proj4[cellsSample, ]
proj5 <- addImputeWeights(proj5)
saveArchRProject(ArchRProj = proj5, overwrite = F,outputDirectory = "Save-cellreport-fibroblast", load = FALSE)

nproj <- addIterativeLSI(
    ArchRProj = proj5,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-FB", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(1.2), 
        sampleCells = 6000, 
        n.start = 10
    ), force = TRUE,
    varFeatures = 50000, 
    dimsToUse = 1:15
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-FB",
    method = "Seurat",
    #force = TRUE,
    name = "ClustersFB2",
  force = TRUE,
    resolution = 0.5)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-FB", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.05, 
  force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "sample", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "ClustersFB2", embedding = "UMAP")

ggAlignPlots(p1, p2,p3, type = "h")
plotPDF(p1,p2,p3, name = "CellReport-fibroblast-UMAP-2.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

######rm C1 C2 #######

idxSample <- BiocGenerics::which(nproj$ClustersFB %in% c("C1","C2"))
cellsSample <- nproj$cellNames[-idxSample]
proj5=nproj[cellsSample, ]
#proj5 <- addImputeWeights(proj5)
saveArchRProject(ArchRProj = proj5, overwrite = F,outputDirectory = "Save-cellreport-fibroblast", load = FALSE)

######our top marker gene 在heatmap上的分布
ourFBproj<-loadArchRProject(path = "../Save-fibroblast2")

markersGS <- getMarkerFeatures(
    ArchRProj = ourFBproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
FB1<-markerList$FB1$name
FB2<-markerList$FB2$name
FB3<-markerList$FB3$name
FB4<-markerList$FB4$name
FBgene<-c(FB1,FB2,FB3,FB4)


markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "ClustersFB2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = FBgene,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "CellReport-FB-topmarker-GeneScores-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

cluster<-nproj$ClustersFB2
table(nproj$ClustersFB2)
for (i in 1:6913){
if(cluster[i]=="C1"){cluster[i]="F1"}
if(cluster[i]=="C2"){cluster[i]="F1"}
if(cluster[i]=="C3"){cluster[i]="F2"}
if(cluster[i]=="C4"){cluster[i]="F3"}
if(cluster[i]=="C5"){cluster[i]="F4"}
if(cluster[i]=="C6"){cluster[i]="F5"}
if(cluster[i]=="C7"){cluster[i]="F3"}
if(cluster[i]=="C8"){cluster[i]="F3"}
if(cluster[i]=="C9"){cluster[i]="F2"}
}
nproj$subtype<-cluster


markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = FBgene,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "CellReport-FB-topmarker-subtype", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

F1<-markerList$F1$name
F2<-markerList$F2$name
F3<-markerList$F3$name
F4<-markerList$F4$name
F5<-markerList$F5$name


#####GO analysis for F##########
library(clusterProfiler)
library(org.Mm.eg.db)

pdf("F-GO.pdf")
ego_geneID <- bitr(F5, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
###GO analysis####
ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"F5-BP.csv")
nproj<-loadArchRProject(path = "./publish_data/scATAC/Save-cellreport-all-filter-addcelltype")
FB<-loadArchRProject(path ="./publish_data/scATAC/Save-cellreport-fibroblast")
EC<-loadArchRProject(path ="./publish_data/scATAC/Save-cellreport-Endothelial")


marker_genes<-c("Brinp1","Pim1","Birc6","Il1rapl1","Cdh18","Grid2","Mdm4","Cebpd")######FB marker genes

####features 中不存在（其他几个成员均存在）

library(reshape2)
GM<-getMatrixFromProject(
  ArchRProj =FB,
  useMatrix = "GeneScoreMatrix"
)
data<-assays(GM)$GeneScoreMatrix
data<-as(as.matrix(data), 'sparseMatrix')
gene<-rowData(GM)$name
rownames(data)<-gene
data<-as.matrix(data)
#data<-as.data.frame(data)
vln.df=as.data.frame(data[marker_genes,])
vln.df$gene=rownames(vln.df)

vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=data.frame(CB = rownames(FB@cellColData),celltype=FB$sample)
#anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = marker_genes) #为了控制画图的基因顺序

pdf("scATAC-Vlnplot-FBmarkergene-sample.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+geom_boxplot(width=0.1,outlier.size=0,
    notch=TRUE,position=position_dodge(width=0.9))+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()

##########FB top gene violin#########
marker_genes<-c("Frmpd4"       ,"Hspa5"        ,"Aff2"         ,"Cxcl2"        ,
"Efhc2"        ,"Lrp1b"        ,"Fgf14"        ,"Ralyl"        ,
"Smoc1"        ,"4933402E13Rik","4930524N10Rik","Gm8787"       ,
"Lgals8"       ,"Ldhd"         ,"Ctsc"         ,"Il13ra2"      ,
"Gm5071"       ,"Rtl4"         ,"4930474H20Rik","Mageb18"      ,
"Mug2"         ,"Cntnap3"      ,"Ttr"          ,"Ppp1r2-ps9"   ,
"Mageb2"       ,"Npbwr1"       ,"Snca"         ,"Snhg14"       ,
"Psg18"        ,"Vmn2r104"     ,"Mir7669"      ,"Mageb1"       ,
"Hp"           ,"Gm15328"      ,"C130080G10Rik","Vmn2r120"     ,
"Gm5622"       ,"Serpina3h"    ,"Cpxcr1"       ,"Gm36633"      ,
"Zpld1"        ,"Tex13c1"      ,"Gm16445"      ,"Vmn1r221"     ,
"Olfr1341"     ,"Mir98"        ,"4930520P13Rik","Clec2f"       ,
"Rp1"          ,"Esp24"        ,"Dgkk"         ,"Mir6411"      ,
"Fhl5"         ,"4930515L19Rik","BC048507"     ,"Fabp12"       ,
"Mill1"        ,"Cypt15"       ,"Adad1"        ,"Olfr1323"     ,
"4930552N02Rik","Ssxb1"        ,"Gm5166"       ,"Serpina3g"    ,
"5430427O19Rik","1810014B01Rik","Vmn2r3"       ,"Samt2"        ,
"Slc39a4"      ,"Lce1d"        ,"4933433F19Rik","Pglyrp3"      ,
"Mir6901"      ,"Ros1"         ,"Ifi213"       ,"Olfr586"      ,
"Mir6996"      ,"Cypt14"       ,"Mir6954"      ,"Slitrk4"      ,
"1700031F05Rik","4930595M18Rik","Ctag2"        ,"Tubb2a-ps2"   ,
"H2afb1"       ,"BC061195"     ,"Mir7091"      ,"4932429P05Rik",
"Olfr589"      ,"Tex16"        ,"Usp17ld"      ,"Klra9"        ,
"Olfr1451"     ,"Slc15a3"      ,"Vmn2r68"      ,"Olfr1507"     ,
"Ssxb9"        ,"Tmem121b"     ,"Gm732"        ,"Lrrc14b"      ,
"A530058N18Rik","Tyr"          ,"Crisp1"       ,"Mageb4"       ,
"Tpsg1"        ,"Nup37"        ,"Plin4"        ,"Yipf7"        ,
"Dmrtc1a"      ,"Klre1"        ,"Olfm3"        ,"Mup5"         ,
"Olfr346"      ,"Cd209f")   
vln.df=as.data.frame(data[marker_genes,])
vln.df<-na.omit(vln.df)
vln<-vln.df
pdf("scATAC-Vlnplot-FBtopgene-mergeday13.pdf")

for(i in 1:19){
  j=i*6
  k=j-5
  vln.df=vln[k:j,]
  vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
anno=data.frame(CB = rownames(FB@cellColData),celltype=FB$sample)
#anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
vln.df=inner_join(vln.df,anno,by="CB")
#vln.df$gene=factor(vln.df$gene,levels = marker_genes) #为了控制画图的基因顺序

p<-vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+geom_boxplot(width=0.1,outlier.size=0,
    notch=TRUE,position=position_dodge(width=0.9))+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
print(p)
}
dev.off()





nega_marker_genes<-c("Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Tcf4","Col4a2")######EC3 marker genes
pro_marker_genes<-c("Plcg1","Acvrl","Flt1","Itgb3","Itgb1bp1")######EC3 marker genes

GM<-getMatrixFromProject(
  ArchRProj =EC,
  useMatrix = "GeneScoreMatrix"
)
data<-assays(GM)$GeneScoreMatrix
data<-as(as.matrix(data), 'sparseMatrix')
gene<-rowData(GM)$name
rownames(data)<-gene
data<-as.matrix(data)


vln.df=as.data.frame(data[nega_marker_genes,])
vln.df$gene=rownames(vln.df)
df<-vln.df
vln.df<-df[1:4,]
vln.df<-df[5:8,]

vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=data.frame(CB = rownames(EC@cellColData),celltype=EC$sample)
vln.df=inner_join(vln.df,anno,by="CB")
pdf("scATAC-Vlnplot-EC-nega_marker_genes-sample.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+geom_boxplot(width=0.1,outlier.size=0,
    notch=TRUE,position=position_dodge(width=0.9))+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()
pro_marker_genes<-gene[which(gene%in%pro_marker_genes)]
vln.df=as.data.frame(data[pro_marker_genes,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
vln.df=inner_join(vln.df,anno,by="CB")
pdf("Vlnplot-EC-pro_marker_genes-sample-each.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+geom_boxplot(width=0.1,outlier.size=0,
    notch=TRUE,position=position_dodge(width=0.9))+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()




GM<-getMatrixFromProject(
  ArchRProj =nproj,
  useMatrix = "GeneScoreMatrix"
)
data<-assays(GM)$GeneScoreMatrix
data<-as(as.matrix(data), 'sparseMatrix')
gene<-rowData(GM)$name
rownames(data)<-gene
data<-as.matrix(data)
#data<-as.data.frame(data)

PFgene<-c("Trap1","Cct2","Dnajb5","Dnajb11","Ppia","Ppil1","Stub1","Cdc37","Mkks","Grpel1","Fkbp8","Calr","P3h1","Ptges3l","Dnajb1","Hsp90ab1","Hsp90aa1","Hspd1","Pdcl3","Nudc","Dnaja1","Pfdn6","Mesd","Hspa1b","Hspe1","Hspa8","Hspa1a")

vln.df=as.data.frame(data[PFgene,])
vln.df$gene=rownames(vln.df)

vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")


library(RColorBrewer)
colourCount = length(unique(PFgene))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
library(ggridges)
library(ggplot2)

pdf("scATAC-Ridgeplot-PF-markergene-mergeday13.pdf")
ggplot(vln.df, aes(x = exp, y = gene, fill = gene)) +
  geom_density_ridges() +
  theme_ridges() + xlim(-0.5, 5)+
  theme(legend.position = "none")+
  scale_fill_manual(values = getPalette(colourCount))
 dev.off();
library(pheatmap)

vln.df=as.data.frame(data[PFgene,])
#vln.df<-apply(vln.df,2,as.numeric)
vln.df<-na.omit(vln.df)
count=t(scale(t(vln.df),scale = T,center = T))


pdf("Allcell_PFgene_heatmap-zscore.pdf")
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy","white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=F)

dev.off()








pdf("Allcell_PFgene_heatmap-zscore-bk3.pdf")
bk <- c(seq(-1,-0.1,by=0.001),seq(0,1,by=0.001))


pheatmap(count,
         scale = "none",
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         #legend_breaks=seq(-8,8,2),
         show_rownames=T,show_colnames=F,
         breaks=bk)


##########用gene promoter处accessibility 代替genescore############
library(RColorBrewer)
#########FB signal###########
FBmarker_genes<-c("Brinp1","Pim1","Birc6","Il1rapl1","Cdh18","Grid2","Mdm4","Cebpd")######FB marker genes
marker_genes<-c("Frmpd4"       ,"Hspa5"        ,"Aff2"         ,"Cxcl2"        ,
"Efhc2"        ,"Lrp1b"        ,"Fgf14"        ,"Ralyl"        ,
"Smoc1"        ,"4933402E13Rik","4930524N10Rik","Gm8787"       ,
"Lgals8"       ,"Ldhd"         ,"Ctsc"         ,"Il13ra2"      ,
"Gm5071"       ,"Rtl4"         ,"4930474H20Rik","Mageb18"      ,
"Mug2"         ,"Cntnap3"      ,"Ttr"          ,"Ppp1r2-ps9"   ,
"Mageb2"       ,"Npbwr1"       ,"Snca"         ,"Snhg14"       ,
"Psg18"        ,"Vmn2r104"     ,"Mir7669"      ,"Mageb1"       ,
"Hp"           ,"Gm15328"      ,"C130080G10Rik","Vmn2r120"     ,
"Gm5622"       ,"Serpina3h"    ,"Cpxcr1"       ,"Gm36633"      ,
"Zpld1"        ,"Tex13c1"      ,"Gm16445"      ,"Vmn1r221"     ,
"Olfr1341"     ,"Mir98"        ,"4930520P13Rik","Clec2f"       ,
"Rp1"          ,"Esp24"        ,"Dgkk"         ,"Mir6411"      ,
"Fhl5"         ,"4930515L19Rik","BC048507"     ,"Fabp12"       ,
"Mill1"        ,"Cypt15"       ,"Adad1"        ,"Olfr1323"     ,
"4930552N02Rik","Ssxb1"        ,"Gm5166"       ,"Serpina3g"    ,
"5430427O19Rik","1810014B01Rik","Vmn2r3"       ,"Samt2"        ,
"Slc39a4"      ,"Lce1d"        ,"4933433F19Rik","Pglyrp3"      ,
"Mir6901"      ,"Ros1"         ,"Ifi213"       ,"Olfr586"      ,
"Mir6996"      ,"Cypt14"       ,"Mir6954"      ,"Slitrk4"      ,
"1700031F05Rik","4930595M18Rik","Ctag2"        ,"Tubb2a-ps2"   ,
"H2afb1"       ,"BC061195"     ,"Mir7091"      ,"4932429P05Rik",
"Olfr589"      ,"Tex16"        ,"Usp17ld"      ,"Klra9"        ,
"Olfr1451"     ,"Slc15a3"      ,"Vmn2r68"      ,"Olfr1507"     ,
"Ssxb9"        ,"Tmem121b"     ,"Gm732"        ,"Lrrc14b"      ,
"A530058N18Rik","Tyr"          ,"Crisp1"       ,"Mageb4"       ,
"Tpsg1"        ,"Nup37"        ,"Plin4"        ,"Yipf7"        ,
"Dmrtc1a"      ,"Klre1"        ,"Olfm3"        ,"Mup5"         ,
"Olfr346"      ,"Cd209f")   
FBmarker_genes<-c(FBmarker_genes,marker_genes)

FB <- addGroupCoverages(ArchRProj = FB, groupBy = "sample")
pathToMacs2 <- findMacs2()
FB <- addReproduciblePeakSet(
    ArchRProj = FB, 
    groupBy = "sample", 
    pathToMacs2 = pathToMacs2
)
FBpeak<-getPeakSet(FB)

FBpeak<-FBpeak[which(FBpeak$nearestGene%in% FBmarker_genes),]
FBpeak<-FBpeak[which(FBpeak$peakType=="Promoter"),]

PM<-getMatrixFromProject(
  ArchRProj =FB,
  useMatrix = "PeakMatrix"
)
data<-assays(PM)$PeakMatrix
data<-as(as.matrix(data), 'sparseMatrix')
peakrange<-PM@ rowRanges

####peak matrix中第N行属于我们需要的peak#####
FBgene<-FBpeak$nearestGene
FBgene<-unique(FBgene)
P1MI<-rownames(FB@cellColData)[which(FB$sample=="P1D3MI")]
P1Sham<-rownames(FB@cellColData)[which(FB$sample=="P1D3Sham")]
P8MI<-rownames(FB@cellColData)[which(FB$sample=="P8D3MI")]
P8Sham<-rownames(FB@cellColData)[which(FB$sample=="P8D3Sham")]
pdf("scATAC-promoter-FB2topgene-ECDF.pdf",width=5,height=3)

for(i in 1:15){
  gene<-FBgene[i];
  peak<-FBpeak[which(FBpeak$nearestGene%in% gene),]
  o<-findOverlaps(peak,peakrange)@to
  print(o);
  vln.df=data[o,];
  if(length(o)>=2){
    vln<-colSums(vln.df)
  }else{
    vln<-vln.df
  };
  P1MIcount<-vln[P1MI]
  P1Shamcount<-vln[P1Sham]
  P8MIcount<-vln[P8MI]
  P8Shamcount<-vln[P8Sham]
  mydf = data.frame(
     value = c(P1MIcount,P8MIcount,P1Shamcount,P8Shamcount),
     sample =c(rep("P1MIcount",length(P1MIcount)),rep("P8MIcount",length(P8MIcount)),
      rep("P1Shamcount",length(P1Shamcount)),rep("P8Shamcount",length(P8Shamcount)))
  )
  
  p<-ggplot(mydf,aes(x = value)) + stat_ecdf(aes(colour = sample))+
      scale_color_manual(values =brewer.pal(8,"Set2")[1:4])+ggtitle(gene)+
      #scale_x_discrete("")+scale_y_continuous("")+
      theme_bw() +
      theme(panel.grid=element_blank())
  print(p)
}

dev.off();




#########EC signal###########

nega_marker_genes<-c("Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Tcf4","Col4a2")######EC3 marker genes
pro_marker_genes<-c("Plcg1","Acvrl","Flt1","Itgb3","Itgb1bp1")######EC3 marker genes
ECmarker_genes<-c(nega_marker_genes,pro_marker_genes)

EC <- addGroupCoverages(ArchRProj = EC, groupBy = "sample")
pathToMacs2 <- findMacs2()
EC <- addReproduciblePeakSet(
    ArchRProj = EC, 
    groupBy = "sample", 
    pathToMacs2 = pathToMacs2
)
ECpeak<-getPeakSet(EC)

ECpeak<-ECpeak[which(ECpeak$nearestGene%in% ECmarker_genes),]
ECpeak<-ECpeak[which(ECpeak$peakType=="Promoter"),]
EC <- addPeakMatrix(EC)
PM<-getMatrixFromProject(
  ArchRProj =EC,
  useMatrix = "PeakMatrix"
)
data<-assays(PM)$PeakMatrix
data<-as(as.matrix(data), 'sparseMatrix')
peakrange<-PM@ rowRanges

####peak matrix中第N行属于我们需要的peak#####
ECgene<-ECpeak$nearestGene
ECgene<-unique(ECgene)
P1MI<-rownames(EC@cellColData)[which(EC$sample=="P1D3MI")]
P1Sham<-rownames(EC@cellColData)[which(EC$sample=="P1D3Sham")]
P8MI<-rownames(EC@cellColData)[which(EC$sample=="P8D3MI")]
P8Sham<-rownames(EC@cellColData)[which(EC$sample=="P8D3Sham")]

pdf("scATAC-promoter-EC3topgene-ECDF.pdf",width=5,height=3)
for(i in 1:12){
  gene<-ECgene[i];
  peak<-ECpeak[which(ECpeak$nearestGene%in% gene),]
  o<-findOverlaps(peak,peakrange)@to
  print(o);
  vln.df=data[o,];
  if(length(o)>=2){
    vln<-colSums(vln.df)
  }else{
    vln<-vln.df
  };
  P1MIcount<-vln[P1MI]
  P1Shamcount<-vln[P1Sham]
  P8MIcount<-vln[P8MI]
  P8Shamcount<-vln[P8Sham]
  mydf = data.frame(
     value = c(P1MIcount,P8MIcount,P1Shamcount,P8Shamcount),
     sample =c(rep("P1MIcount",length(P1MIcount)),rep("P8MIcount",length(P8MIcount)),
      rep("P1Shamcount",length(P1Shamcount)),rep("P8Shamcount",length(P8Shamcount)))
  )
  
  p<-ggplot(mydf,aes(x = value)) + stat_ecdf(aes(colour = sample))+
      scale_color_manual(values =brewer.pal(8,"Set2")[1:4])+ggtitle(gene)+
      #scale_x_discrete("")+scale_y_continuous("")+
      theme_bw() +
      theme(panel.grid=element_blank())
  print(p) 
}

dev.off();
