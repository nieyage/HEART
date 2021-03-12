#######给全部cluster添加注释名，并根据细胞类型分类 #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 1) 
addArchRGenome("mm10")

####RNA#########
library(Seurat)
scRNA <- readRDS("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/integrative_analysis/cmc_sct.rds")
selRNA <- scRNA[,which(scRNA@meta.data$orig.ident == "A3" |
                         scRNA@meta.data$orig.ident == "A7" | scRNA@meta.data$orig.ident == "A14" |
                         scRNA@meta.data$orig.ident == "P7" | scRNA@meta.data$orig.ident == "P14" )]

selRNA<-selRNA[,which(selRNA@meta.data$cell_types == "Endotheliocytes" )]
DefaultAssay(selRNA) <- "RNA" # Create dummy new assay to demo switching default assays
selRNA <- NormalizeData(selRNA) ###Normalizing the data
#Identification of highly variable features (feature selection)
selRNA <- FindVariableFeatures(selRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(selRNA)
selRNA <- ScaleData(selRNA, features = all.genes)
#Perform linear dimensional reduction
selRNA <- RunPCA(selRNA, features = VariableFeatures(object = selRNA))

selRNA.integrated <- JackStraw(selRNA, num.replicate = 100, dims =50)
selRNA.integrated <- ScoreJackStraw(selRNA.integrated, dims = 1:50)

pdf("/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype/RNAplot/selscRNA.pdf")
JackStrawPlot(selRNA.integrated, dims = 1:50)
ElbowPlot(selRNA.integrated,ndims=50)
dev.off()

#Cluster the cells
selRNA.cluster.integrated <- FindNeighbors(selRNA, dims = 1:25)
selRNA.cluster.integrated <- FindClusters(selRNA.cluster.integrated, resolution = 0.5)

selRNA.cluster.integrated<-BuildClusterTree(selRNA.cluster.integrated)
Tool(object = selRNA.cluster.integrated, slot = 'BuildClusterTree')
pdf("/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype/RNAplot/PlotClusterTreeselscRNA.pdf")
PlotClusterTree(selRNA.cluster.integrated)
dev.off();
#Run non-linear dimensional reduction (UMAP/tSNE)
selRNA.cluster.integrated <- RunUMAP(object = selRNA.cluster.integrated, dims = 1:25)
pdf("/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype/RNAplot/EC_cluster_Umap.pdf")
DimPlot(object = selRNA.cluster.integrated, group.by = "orig.ident",label = T)
table(selRNA.cluster.integrated@meta.data$orig.ident,selRNA.cluster.integrated$seurat_clusters)
dev.off();
pdf("/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype/RNAplot/selscRNA_test_Umap.pdf")
DimPlot(object = selRNA.cluster.integrated, reduction = "umap",label = T)
dev.off()

##整合MP#######
#MP3proj<-loadArchRProject(path = "Save-Proj2-MP3-4Subtype")
#
#MP3<- addGeneIntegrationMatrix(
#  ArchRProj = MP3proj,
#  useMatrix = "GeneScoreMatrix",
#  matrixName = "GeneIntegrationMatrix",
#  seRNA = selRNA.cluster.integrated,
#  addToArrow = TRUE,
#  groupRNA = "cell_types",
#  force = TRUE,
#  nameCell = "predictedCell_Un",
#  nameGroup = "predictedGroup_Un",
#  nameScore = "predictedScore_Un"
#)
#getAvailableMatrices(MP3)

EC4proj<-loadArchRProject(path = "ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
EC4<- addGeneIntegrationMatrix(
  ArchRProj = EC4proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  seRNA = selRNA.cluster.integrated,
  addToArrow = TRUE,
  groupRNA = "cell_types",
  nameCell = "predictedCell_Un",
  force = TRUE,
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
getAvailableMatrices(EC4)

#SMC1<- addGeneIntegrationMatrix(
#  ArchRProj = SMC1proj,
#  useMatrix = "GeneScoreMatrix",
#  matrixName = "GeneIntegrationMatrix",
#  seRNA = selRNA.cluster.integrated,
#  addToArrow = TRUE,
#   reducedDims = "IterativeLSI-dim15-SMC-1.2-80000-4",
#  groupRNA = "cell_types",
#  nameCell = "predictedCell_Un",
#  nameGroup = "predictedGroup_Un",
#  force = TRUE,
#  nameScore = "predictedScore_Un"
#)
#getAvailableMatrices(SMC1)
#
#
#all <- addGeneIntegrationMatrix(
#  ArchRProj = all,
#  useMatrix = "GeneScoreMatrix",
#  matrixName = "GeneIntegrationMatrix",
#  reducedDims = "IterativeLSI-dim8",
#  seRNA = selRNA.cluster.integrated,
#  addToArrow = TRUE,
#  groupRNA = "cell_types",
#  nameCell = "predictedCell_Un",
#   useImputation = FALSE,
#   force = TRUE,
#  nameGroup = "predictedGroup_Un",
#  nameScore = "predictedScore_Un"
#)

#######EC3中target gene的表达情况#########
######negative regulation of vasculature development######
Adamts1/Ccn6/Cd36/Col4a2/Foxj2/Itgb1bp1/Meox2/Nfatc1/Rgcc/Tcf4/Zeb1


negative_gene<-"Adamts1/Ccn6/Cd36/Col4a2/Foxj2/Meox2/Tcf4/Xdh/Zeb1"
Xdh:negative regulation of endothelial cell proliferation
Zeb1:negative regulation of endothelial cell differentiation
Rgcc:negative regulation of angiogenesis
Nfatc1:negative regulation of VSMC formation 
Acvrl1/Flt1/Itgb1bp1/Itgb3/Plcg1/Rgcc/Xdh
Cd36:negative regulation of angiogenesis

ECproliferation_gene<-"Acvrl1/Flt1/Itgb1bp1/Itgb3/Plcg1/Acvrl1"
negative_gene<-strsplit(negative_gene,'/')
ECproliferation_gene<-strsplit(ECproliferation_gene,'/')
negative_gene<-negative_gene[[1]]
ECproliferation_gene<-ECproliferation_gene[[1]]

EC3target<-c("1600002D24Rik","1700026D11Rik","1700049E15Rik","1700061I17Rik","2310061N02Rik","4921539E11Rik","4930459C07Rik","4930486F22Rik","4930503H13Rik",
	"4930524C18Rik","8430437L04Rik","A230028O05Rik","Abr","Actbl2","Acvrl1","Adam28","Adamts1","Adamts5","Adarb1","Afap1l2","Agpat4","AI427809","Aldh8a1",
	"Ank3","Ankrd55","Ap1s3","Apbb1ip","Apold1","Arf2","Arhgap18","Arhgap28","Arhgef2","Arhgef7","Asap1","Atp10d","Casz1","Ccdc12","Ccn6","Cd36","Cenpi",
	"Cfap69","Clrn1","Col4a1","Col4a2","Crhbp","Crp","Cubn","Cyyr1","Dgki","Dip2c","Dock4","Dtl","E230013L22Rik","Ebf1","Efcc1","Ell","Emb","Eml6","Ercc4",
	"Ets1","Ets2","Eva1c","Exoc6","Fam129b","Fam187b","Fam76b","Fbxo33","Fig4","Flt1","Foxa2","Foxj2","Fut2","Galnt2","Gimap5","Gm1821","Gm20740","Gm29536",
	"Gm32865","Gm39653","Gm6116","Gm6559","Gramd4","H2afy3","Hand2os1","Heatr5a","I730028E13Rik","Igfbp3","Inpp5f","Irx3","Itgb1bp1","Itgb3","Krtap13",
	"Lama4","Lmnb1","Ltbp1","Ly75","Lyrm1","Man2a2","Mbnl1","Mctp1","Mctp2","Meox2","Mir101b","Mir135b","Mir1951","Mir3967","Mir467h","Mir5114","Mro",
	"Myo5c","N4bp2l1","Nab1","Nav1","Nckap5","Ncoa2","Nde1","Nebl","Nfatc1","Nid1","Nipsnap3b","Ntn5","Papss1","Phf2","Plcb1","Plcg1","Plekha7","Plet1",
	"Plpp1","Plppr4","Plxna2","Pnrc1","Pter","Ptpn3","Ptprg","Rasgef1b","Reln","Rgcc","Rnf125","Rnf128","Rnf19a","Rrp1b","Rtl5","S1pr1","Sema3a","Sh3bp5",
	"Shroom3","Slc17a6","Slc36a4","Slc45a1","Slc46a3","Sorbs2","St3gal1","St3gal6","St8sia1","Stk24","Tank","Tapt1","Tcf4","Tlr13","Tmem161b","Tnfrsf9",
	"Tsga13","Tshz1","Ttc30a2","Ttpal","Tubb2a","Tulp1","Tyk2","Xdh","Xiap","Xylt1","Zeb1","Zfp326","Zfp608","Zfp648","Zfp697")

p<-getGroupSE(
  ArchRProj = EC4proj,
  useMatrix = "GeneIntegrationMatrix",
  groupBy = "EC4proj$Subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = EC4proj,
  useMatrix = "GeneIntegrationMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
count<-as.data.frame(count)
count=as.data.frame(lapply(count,as.numeric))
rownames(count)<-gene

#ECcount=t(scale(t(count),scale = T,center = T))


allcount<-count[which(rownames(count)%in% EC3target),]

#######protein_folding genes heatmap in FB MP EC ####
library(pheatmap)
ECcount=t(scale(t(allcount),scale = T,center = T))
ECcount<-na.omit(ECcount)
pdf("EC3target-heatmap.pdf")
pheatmap(ECcount,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 20, #cellheight = 4,
         show_rownames=F,show_colnames=T)

pheatmap(ECcount,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 20, #cellheight = 4,
         show_rownames=F,show_colnames=T)
dev.off();

#####针对负调控基因#######

negacount<-ECcount[which(rownames(ECcount)%in% negative_gene),]
pdf("negative_gene_heatmap.pdf")
pheatmap(negacount,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 20, cellheight = 20,
         show_rownames=T,show_colnames=T)
dev.off()
ECproliferation_genecount<-ECcount[which(rownames(ECcount)%in% ECproliferation_gene),]
pdf("ECproliferation_gene_heatmap.pdf")
pheatmap(ECproliferation_genecount,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 20, cellheight = 20,
         show_rownames=T,show_colnames=T)
dev.off()

######violin plot ###########
EC4<-addImputeWeights(EC4proj)
p<- plotGroups(
    ArchRProj = EC4, 
    groupBy = "EC4$Subtype", 
    colorBy = "GeneIntegrationMatrix", 
	name = negative_gene,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    pal =c("#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F"),#######EC12345 color 
    addBoxPlot = TRUE,
    size=0.2
   )
"#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F"

#Rearrange for grid plotting

p2 <- lapply(p, function(x){x })
pdf("negative_gene_violin_9fig.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))
dev.off();
plotPDF(plotList = p,     name = "negative_gene_violin.pdf", ArchRProj = EC4, 
    addDOC = FALSE, width = 5, height = 5)

p<- plotGroups(
    ArchRProj = EC4, 
    groupBy = "EC4$Subtype", 
    colorBy = "GeneScoreMatrix", 
	name = negative_gene,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    pal =c("#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F"),#######EC12345 color 
    addBoxPlot = TRUE,
    size=0.2
   )

#Rearrange for grid plotting
p2 <- lapply(p, function(x){x })
pdf("negative_gene_violin_accessbility_9fig.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))
dev.off();


######this 9 genes's track#######
p <- plotBrowserTrack(
    ArchRProj = EC4, 
    groupBy = "Subtype", 
    geneSymbol = negative_gene, 
    upstream = 50000,
    downstream = 50000,
    pal=c("#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F")
)

plotPDF(plotList = p, 
    name = "negative_gene_track.pdf", 
    ArchRProj = EC4, 
    addDOC = FALSE, width = 5, height = 5)

##########overlap genes experssion in 4 subtypes #######
all <-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/Save-allcelltype" )
all<-addImputeWeights(all,reducedDims="IterativeLSI-dim8")
overlap<-c("Trap1","Cct2","Dnajb5","Dnajb11","Ppia","Ppil1","Stub1","Cdc37","Mkks","Grpel1","Fkbp8","Calr","P3h1",
  "Ptges3l","Dnajb1","Hsp90ab1","Hsp90aa1","Hspd1","Pdcl3","Nudc","Dnaja1","Pfdn6","Mesd","Hspa1b","Hspe1","Hspa8","Hspa1a")
getAvailableMatrices(all)

#########4个subtype在所有细胞umap上的分布#########
FB3proj<-loadArchRProject(path = "Save-fibroblast2")
EC4proj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
MP3proj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")
SMC1proj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")

######提取这四个cluster的barcode#######
cellNames <- EC4proj$cellNames
m<-EC4proj$Subtype
EC4proj <- addCellColData(ArchRProj = EC4proj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)

cellNames <- MP3proj$cellNames
m<-MP3proj$Subtype
MP3proj <- addCellColData(ArchRProj = MP3proj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)

FB3barcode<-rownames(FB3proj@cellColData[which(FB3proj@cellColData$subtype=="FB3"),])
EC4barcode<-rownames(EC4proj@cellColData[which(EC4proj@cellColData$subtype=="EC4"),])
MP3barcode<-rownames(MP3proj@cellColData[which(MP3proj@cellColData$subtype=="MP3"),])
SMC1barcode<-rownames(SMC1proj@cellColData[which(SMC1proj@cellColData$subtype=="SMC1"),])

cellsSample <- c(FB3barcode,EC4barcode,MP3barcode,SMC1barcode)
nproj=proj5[cellsSample, ]
nproj<-addImputeWeights(nproj,reducedDims="IterativeLSI-dim8")
subtype<-c(rep("other",9278))
m<-rownames(nproj@cellColData)
for (i in 1:9278){
if(m[i] %in% FB3barcode)
    subtype[i]="FB3"
if(m[i] %in% MP3barcode)
    subtype[i]="MP3"
if(m[i] %in% EC4barcode)
    subtype[i]="EC4"
if(m[i] %in% SMC1barcode)
    subtype[i]="SMC1"    
    }
table(subtype)

cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(subtype),
    cells=cellNames, name = "subtype",force = TRUE)
 EC4  FB3  MP3 SMC1 
4908 2796  267 1307 
nproj <- addGeneIntegrationMatrix(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI-dim8",
  seRNA = selRNA.cluster.integrated,
  addToArrow = TRUE,
  groupRNA = "cell_types",
  nameCell = "predictedCell_Un",
   useImputation = FALSE,
   force = TRUE,
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
getAvailableMatrices(nproj)

p<- plotGroups(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    colorBy = "GeneIntegrationMatrix", 
	name = overlap,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    pal =c("#0072B5FF", "#E18727FF","#20854EFF","#7876B1FF"),#######EC12345 color 
    addBoxPlot = TRUE,
    size=0.2
   )
"#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F"

#Rearrange for grid plotting
p2 <- lapply(p, function(x){x })
pdf("overlap_gene_violin_9fig.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))
dev.off();
plotPDF(plotList = p,name="/md01/nieyg/scATAC-ArchR/overlap_gene_violin.pdf", ArchRProj = EC4, addDOC = FALSE, width = 5, height = 5)



#######all subtypes########
FBproj<-loadArchRProject(path = "Save-fibroblast2")
ECproj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
MPproj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")
SMCproj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")

######提取这四个cluster的barcode#######
cellNames <- ECproj$cellNames
m<-ECproj$Subtype
ECproj <- addCellColData(ArchRProj = ECproj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)

cellNames <- MPproj$cellNames
m<-MPproj$Subtype
MPproj <- addCellColData(ArchRProj = MPproj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)


FB1barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB1"),])
FB2barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB2"),])
FB3barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB3"),])
FB4barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB4"),])

EC1barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC1"),])
EC2barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC2"),])
EC3barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC3"),])
EC4barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC4"),])
EC5barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC5"),])

MP1barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP1"),])
MP2barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP2"),])
MP3barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP3"),])
MP4barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP4"),])

SMC1barcode<-rownames(SMCproj@cellColData[which(SMCproj@cellColData$subtype=="SMC1"),])
SMC2barcode<-rownames(SMCproj@cellColData[which(SMCproj@cellColData$subtype=="SMC2"),])
SMC3barcode<-rownames(SMCproj@cellColData[which(SMCproj@cellColData$subtype=="SMC3"),])



########No CMs####
cellsSample <- c(FB1barcode,FB2barcode,FB3barcode,FB4barcode,
	EC1barcode,EC2barcode,EC3barcode,EC4barcode,EC5barcode,
	MP1barcode,MP2barcode,MP3barcode,MP4barcode,
	SMC1barcode,SMC2barcode,SMC3barcode)
nproj=proj5[cellsSample, ]
nproj<-addImputeWeights(nproj,reducedDims="IterativeLSI-dim8")

subtype<-c(rep("other",56362))
m<-rownames(nproj@cellColData)
for (i in 1:56362){
if(m[i] %in% FB3barcode)
    subtype[i]="FB3"
if(m[i] %in% MP3barcode)
    subtype[i]="MP3"
if(m[i] %in% EC4barcode)
    subtype[i]="EC4"
if(m[i] %in% SMC1barcode)
    subtype[i]="SMC1"    
if(m[i] %in% FB1barcode)
    subtype[i]="FB1"
if(m[i] %in% MP1barcode)
    subtype[i]="MP1"
if(m[i] %in% EC1barcode)
    subtype[i]="EC1"
if(m[i] %in% FB2barcode)
    subtype[i]="FB2"
if(m[i] %in% MP2barcode)
    subtype[i]="MP2"
if(m[i] %in% EC2barcode)
    subtype[i]="EC2"
if(m[i] %in% SMC2barcode)
    subtype[i]="SMC2" 
if(m[i] %in% EC3barcode)
    subtype[i]="EC3"
if(m[i] %in% SMC3barcode)
    subtype[i]="SMC3"    
if(m[i] %in% FB4barcode)
    subtype[i]="FB4"
if(m[i] %in% MP4barcode)
    subtype[i]="MP4"
if(m[i] %in% EC5barcode)
    subtype[i]="EC5"
    }
table(subtype)

cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(subtype),
    cells=cellNames, name = "subtype",force = TRUE)

##########add RNA data###########
library(Seurat)
scRNA <- readRDS("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/integrative_analysis/cmc_sct.rds")
selRNA <- scRNA[,which(scRNA@meta.data$orig.ident == "A3" |
                         scRNA@meta.data$orig.ident == "A7" | scRNA@meta.data$orig.ident == "A14" |
                         scRNA@meta.data$orig.ident == "P7" | scRNA@meta.data$orig.ident == "P14" )]
selRNA<-selRNA[,which(selRNA@meta.data$cell_types != "Cardiomyocytes" )]

DefaultAssay(selRNA) <- "RNA" # Create dummy new assay to demo switching default assays
selRNA <- NormalizeData(selRNA) ###Normalizing the data
#Identification of highly variable features (feature selection)
selRNA <- FindVariableFeatures(selRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(selRNA)
selRNA <- ScaleData(selRNA, features = all.genes)
#Perform linear dimensional reduction
selRNA <- RunPCA(selRNA, features = VariableFeatures(object = selRNA))

selRNA.integrated <- JackStraw(selRNA, num.replicate = 100, dims =50)
selRNA.integrated <- ScoreJackStraw(selRNA.integrated, dims = 1:50)
#Cluster the cells
selRNA.cluster.integrated <- FindNeighbors(selRNA, dims = 1:25)
selRNA.cluster.integrated <- FindClusters(selRNA.cluster.integrated, resolution = 0.5)
#Run non-linear dimensional reduction (UMAP/tSNE)
selRNA.cluster.integrated <- RunUMAP(object = selRNA.cluster.integrated, dims = 1:25)

nproj<- addGeneIntegrationMatrix(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI-dim8",
  seRNA = selRNA.cluster.integrated,
  addToArrow = TRUE,
  groupRNA = "cell_types",
  nameCell = "predictedCell_Un",
   useImputation = FALSE,
   force = TRUE,
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

#########plot violin plot #######

p<- plotGroups(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    colorBy = "GeneIntegrationMatrix", 
	name = overlap,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    pal =c("#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F",
           "#41668F","#3383B1","#65A4D1","#B4D1E1",
           "#3D5B4A","#458455","#6A9D62","#C0D2A9",
           "#AAA1CC","#9C63A8","#C39FCB"),
    addBoxPlot = TRUE,
    size=0.2
   )

"#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F" ######EC
"#41668F","#3383B1","#65A4D1","#B4D1E1" ######FB
"#3D5B4A","#458455","#6A9D62","#C0D2A9" #####MP
"#AAA1CC","#9C63A8","#C39FCB" ####SMC

#Rearrange for grid plotting
p2 <- lapply(p, function(x){x })
pdf("overlap_gene_violin_9fig_allsubtype.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2[10:18]))
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2[19:26]))
dev.off();
plotPDF(plotList = p, name = "overlap_gene_violin_allsubtype.pdf", ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)


common_TF <- c("Jun","Fosl1","Batf","Fosl2","Atf3","Bach1","Bach2","Junb")
p<- plotGroups(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    colorBy = "GeneIntegrationMatrix", 
	name = common_TF,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    pal =c("#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F",
           "#41668F","#3383B1","#65A4D1","#B4D1E1",
           "#3D5B4A","#458455","#6A9D62","#C0D2A9",
           "#AAA1CC","#9C63A8","#C39FCB"),
    addBoxPlot = TRUE,
    size=0.2
   )
p2 <- lapply(p, function(x){x })
pdf("common_TF_violin_9fig_allsubtype.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))
dev.off();
plotPDF(plotList = p, name = "common_TF_violin_allsubtype.pdf", ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)




#########plot violin plot #######

p<- plotGroups(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    colorBy = "GeneScoreMatrix", 
  name = overlap,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    pal =c("#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F",
           "#41668F","#3383B1","#65A4D1","#B4D1E1",
           "#3D5B4A","#458455","#6A9D62","#C0D2A9",
           "#AAA1CC","#9C63A8","#C39FCB"),
    addBoxPlot = TRUE,
    size=0.2
   )

"#A92B14","#EC431A","#FF7701","#F1C149","#F4E48F" ######EC
"#41668F","#3383B1","#65A4D1","#B4D1E1" ######FB
"#3D5B4A","#458455","#6A9D62","#C0D2A9" #####MP
"#AAA1CC","#9C63A8","#C39FCB" ####SMC

#Rearrange for grid plotting
p2 <- lapply(p, function(x){x })
pdf("overlap_gene_violin_9fig_allsubtype_access.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2[10:18]))
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2[19:26]))
dev.off();
plotPDF(plotList = p, name = "overlap_gene_violin_allsubtype.pdf", ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

