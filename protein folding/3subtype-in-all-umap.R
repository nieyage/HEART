#######获取各个亚型的protein folding gene的FC情况######
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")
proj5 <- addImputeWeights(proj5)
#######Hspa8在所有细胞类型上的分布状态#########
p<-plotEmbedding(
    ArchRProj = proj5, 
    colorBy = "GeneScoreMatrix", 
    name = "Hspa8", 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj5)
)
plotPDF(plotList = p, 
    name = "Hspa8-umap.pdf", 
    ArchRProj = proj5, 
    addDOC = FALSE, width = 5, height = 5)
#####所有PP基因的分布########
PP<-c("Map3k5","Atf4","Atf6","Bak1","Bcap31","Bax","Bcl2","Capn1","Capn2","Casp12","Stub1","Ddit3","Canx",
	"Calr","Calr3","Calr4","Cul1","Derl2","Derl3","Derl1","Plaa","Ubqln1","Ubqln2","Ubqln3","Ubqln4","Edem1",
	"Edem2","Edem3","Eif2ak1","Eif2ak2","Eif2ak4","Lman1","Lman1l","Man1a","Man1a2","Man1b1","Man1c1",
	"Pdia3","Ppp1r15a","Amfr","Herpud1","Syvn1","Dnaja1","Dnaja2","Dnajb1","Dnajb11","Dnajb12",
	"Dnajb14","Dnajb2","Dnajc1","Dnajc10","Dnajc3","Dnajc5","Dnajc5b","Dnajc5g","Hspa1a","Hspa1b","Hspa1l",
	"Hspa2","Hspa4l","Hspa5","Hspa8","Hsp90aa1","Hsp90ab1","Hsp90b1","Ern1","Mapk8","Map2k7","Nploc4","Nfe2l2",
	"Os9","Rrbp1","Vcp","Prkn","Pdia4","Pdia6","Eif2ak3","Ngly1","Rad23a","Rad23b","Rbx1","Rnf5",
	"Mbtps1","Mbtps2","Sar1a","Sar1b","Preb","Sec13","Sec23a","Sec23b","Sec24a","Sec24b","Sec24c","Sec24d",
	"Sec31a","Sec31b","Sec61a1","Sec61a2","Sec61b","Sec61g","Sec62","Sec63","Sel1l","Sel1l2","Svip",
	"Traf2","Tram1","Ssr1","Ssr2","Ssr3","Ssr4","Ube2j1","Ube2j2","Ube2g1","Ube2g2","Ube2d1","Ube2d2a",
	"Ube2d2b","Ube2d3","Ube2dnl1","Ube2dnl2","Ubxn1","Ubxn2a","Ubxn4","Ubxn6","Ubxn8","Ufd1","Ube4b",
	"Selenos","Lman2","Wfs1","Xbp1","Erlec1")
p<-plotEmbedding(
    ArchRProj = proj5, 
    colorBy = "GeneScoreMatrix", 
    name = PP, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj5)
)
plotPDF(plotList = p, 
    name = "PP-umap.pdf", 
    ArchRProj = proj5, 
    addDOC = FALSE, width = 5, height = 5)

#########3个subtype在所有细胞umap上的分布#########
FB3proj<-loadArchRProject(path = "Save-fibroblast2")
EC4proj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
MP3proj<-readRDS("/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/MP3.rds")
######提取这三个cluster的barcode#######
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

######free umap df#######
#####re draw umap 

df = data.frame(proj5@embeddings$UMAP$df, subtype=c(rep("other",60420)))
df$subtype<-as.character(df$subtype)
subtype<-df$subtype
m<-rownames(df)
for (i in 1:60420){
if(m[i] %in% FB3barcode)
    df[i,3]="FB3"
if(m[i] %in% MP3barcode)
    df[i,3]="MP3"
if(m[i] %in% EC4barcode)
    df[i,3]="EC4"
    }

###timepoint
pdf("threesub_allumap.pdf")
ggplot(df, aes(x = IterativeLSI.dim8.UMAP_Dimension_1, y = IterativeLSI.dim8.UMAP_Dimension_2, color = subtype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c('#BD956A', '#8C549C', '#585658',"grey")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Timepoint")

p<-getGroupSE(
  ArchRProj = FB3proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "time"
)
GM<-getMatrixFromProject(
  ArchRProj = FB3proj,
  useMatrix = "GeneScoreMatrix"
)

gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
count<-as.data.frame(count)
data=as.data.frame(lapply(count,as.numeric))
write.table(data,"PF-Use-matrix-FB.txt")

p<-getGroupSE(
  ArchRProj = EC4proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "genotype"
)
GM<-getMatrixFromProject(
  ArchRProj = EC4proj,
  useMatrix = "GeneScoreMatrix"
)

gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
count<-as.data.frame(count)
data=as.data.frame(lapply(count,as.numeric))
write.table(data,"PF-Use-matrix-EC.txt")


p<-getGroupSE(
  ArchRProj = MP3proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "genotype"
)
GM<-getMatrixFromProject(
  ArchRProj = MP3proj,
  useMatrix = "GeneScoreMatrix"
)

gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
count<-as.data.frame(count)
data=as.data.frame(lapply(count,as.numeric))
write.table(data,"PF-Use-matrix-MP.txt")

uprgene<-uprgene[order(uprgene$progress),]
FBcount<-FB[which(rownames(FB) %in% uprgene$gene),]
rnamat=t(scale(t(FBcount),scale = T,center = T))
head(rnamat)
library(pheatmap)
# 添加注释信息
annotation_col = data.frame(
  Sample =factor( c(rep("FB1",3),rep("FB2",3),rep("FB3",3),rep("FB4",3) )), 
  time = rep(c("Day14","Day3","Day7"),4)
    )
rownames(annotation_col) = colnames(rnamat)
head(annotation_col)
annotation_row = data.frame(progress=uprgene$progress)
rownames(annotation_row) = uprgene$gene
ann_colors = list(
  Sample = c(FB1="#A20056B2",FB2="#3B4992B2",FB3="#F48639",FB4="#D51F26"),
  genotype = c(AR="#E64B35B2",NAR="#00468BB2",P="#42B540B2"),
  progress =c(ATF6="#E64B35B2", IRE1="#4DBBD5B2" ,PERK="#00A087B2", "PERK/ATF6"="#3C5488B2")
)

ann_colors = list(
  Sample = c(FB1="#A20056B2",FB2="#3B4992B2",FB3="#F48639",FB4="#D51F26"),
  time = c(Day14="#8491B4B2",Day3="#91D1C2B2",Day7="#DC0000B2"),
  progress =c(ATF6="#E64B35B2", IRE1="#4DBBD5B2" ,PERK="#00A087B2", "PERK/ATF6"="#3C5488B2")
)

head(ann_colors)
rnamat<-rnamat[rownames(annotation_row),]

pdf("UPR-FB-time-heatmap.pdf")
pheatmap(rnamat,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_row = annotation_row,
         show_rownames=T,show_colnames=T)
dev.off();


ERADgene<-trans[which(trans$progress=="ERAD"),]
FBcount<-FB[which(rownames(FB) %in% ERADgene$gene),]
rnamat=t(scale(t(FBcount),scale = T,center = T))
head(rnamat)
library(pheatmap)

pdf("ERAD-FB-heatmap.pdf")
pheatmap(rnamat,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         #annotation_row = annotation_row,
         show_rownames=F,show_colnames=T)
dev.off();


markersGS <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

getFeatures(
  ArchRProj = proj5,
  useMatrix = "GeneScoreMatrix"
)

m<-getMarkers(markersGS,cutOff = "FDR <= 0.01 & Log2FC >= 2")
m<-rownames(FBcount)
ego_geneID <- bitr(m, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
p <- pathview(gene.data = y[, 1], pathway.id = "04141", species = "mmu", out.suffix = "PF_pathway")


p<-getGroupSE(
  ArchRProj = MP3proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "MP3proj$Subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = MP3proj,
  useMatrix = "GeneScoreMatrix"
)


gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
count<-as.data.frame(count)
data=as.data.frame(lapply(count,as.numeric))

MPcount<-data[which(rownames(data) %in% trans$gene),]
rnamat=t(scale(t(MPcount),scale = T,center = T))
head(rnamat)
ego_geneID <- bitr(rownames(rnamat), fromType = "SYMBOL", toType =  "ENTREZID", OrgDb = org.Mm.eg.db,drop = T)
rnamat<-rnamat[ego_geneID$SYMBOL,]
rnamat<-cbind(ego_geneID,rnamat)