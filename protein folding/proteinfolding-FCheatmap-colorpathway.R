#######获取各个亚型的protein folding gene的FC情况######
#####################################################
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")

#####FB3 protein folding gene#####
g1<-"Hspa8/Hspa1l/Hspa1a/Hspa1b/Hsph1/Dnaja1/Hspb6/Calr/H2-DMa/Hspd1/Grn/Emc6/Hspe1/Dnajb1/Hsp90ab1/Bag1/Calr3"
g2<-"Hspa8/Hspa1l/Hspa1a/Ddit3/Hspa1b/Hsph1/Tbl2/Calr/Hspd1/Atf3/Mfn2/Herpud1/Calr3"
g3<-"Hspa8/Hspa1l/Hspa1a/Hspa1b/Hsph1/H2-DMa/Hspd1/Hspe1/Dnajb1/Bag1"
g4<-"Hspa8/Hspa1l/Hspa1a/Hspa1b/Hsph1/Hspb6/H2-DMa/Hspe1/Dnajb1/Bag1"
g5<-"Hspa8/Hspa1l/Hspa1a/Hspa1b/Hsph1/H2-DMa/Hspe1/Dnajb1/Bag1"
g1<-strsplit(g1,'/')
g2<-strsplit(g2,'/')
g3<-strsplit(g3,'/')
g4<-strsplit(g4,'/')
g5<-strsplit(g5,'/')
g1<-g1[[1]]
g2<-g2[[1]]
g3<-g3[[1]]
g4<-g4[[1]]
g5<-g5[[1]]
FB3<-union(g1,union(g2,union(g3,union(g4,g5))))
#####EC4 gene#########
g1<-"Hspa1a/Hspa8/Dnajb14/Hspa1l"
g1<-strsplit(g1,'/')
EC4<-g1[[1]]
#####MP3#####
g1<-"Hspa1b/Dnajb1/Hspa8/Emc6/Bag5"
g1<-strsplit(g1,'/')
MP3<-g1[[1]]
####all protein folding gene#######
all<-union(FB3,union(MP3,EC4))
######构建protein folding matrix##########
#######FB matrix#######
nproj<-loadArchRProject(path = "Save-fibroblast2")
markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
heatmap<-plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & (Log2FC >= 0.5|Log2FC <= -0.5)",
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  labelMarkers = NULL,
  returnMat = TRUE,
  transpose = FALSE
)
FBmat<-heatmap[which(rownames(heatmap)%in% all),]
########EC4 matrix######
nproj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "nproj$Subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
ECheatmap<-plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & (Log2FC >= 0.5|Log2FC <= -0.5)",
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  labelMarkers = NULL,
  returnMat = TRUE,
  transpose = FALSE
)
ECmat<-ECheatmap[which(rownames(ECheatmap)%in% all),]

########MP3 matrix########
/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/MP3.rds
MP3proj<-readRDS("/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/MP3.rds")
MPmarkersGS <- getMarkerFeatures(
    ArchRProj = MP3proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "MP3proj$Subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
MPheatmap<-plotMarkerHeatmap(
  seMarker = MPmarkersGS,
  cutOff = "FDR <= 0.01 & (Log2FC >= 0.5|Log2FC <= -0.5)",
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  labelMarkers = NULL,
  returnMat = TRUE,
  transpose = FALSE
)
MPmat<-MPheatmap[which(rownames(MPheatmap)%in% all),]
#####all matrix####
FB3<-FBmat[,3]
MP3<-MPmat[,3]
EC4<-ECmat[,4]
EC4<-as.data.frame(EC4)
EC4$name<-rownames(EC4)
MP3<-as.data.frame(MP3)
MP3$name<-rownames(MP3)
FB3<-as.data.frame(FB3)
FB3$name<-rownames(FB3)
c<-merge(EC4,FB3,by="name",all=T)
c[is.na(c)] <- 0
pdf("protein_folding.pdf")
pheatmap(c,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(colors = c("white","red"))(100),
         cellwidth = 12, cellheight = 12,
         show_rownames=T,show_colnames=T)
dev.off()

library(clusterProfiler)
library(org.Mm.eg.db)
ego_geneID <- bitr(all, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)

ECmat<-ECheatmap[which(rownames(ECheatmap)%in% all),]

#mmu    highlight
15482	red
15505	red
15510	red
15528	red
11910	red
13198	red
66048	red
14824	red
14998	red
15511	red
193740	red
15516	red
15502	red
12017	red
170731	red
27368	red
243912	red
73316	red
81489	red
12317	red
64209	red
15481	red
70604	red
70369	red


PP<-c("Hspa1l","Hsph1","Hspd1","Hspe1","Atf3","Ddit3","Emc6","Grn","H2-DMa","Hspa1b","Hspa1a","Hsp90ab1","Dnaja1","Bag1","Mfn2","Tbl2","Hspb6",
	"Calr3","Dnajb1","Calr","Herpud1","Hspa8","Dnajb14","Bag5","Sec61a1","Sec61a2","Sec61b","Sec61g","Gm15266","Gm30534","Sec62","Sec63","Rpn1",
	"Rpn2","Dad1","Tusc3","Ddost","Stt3a","Stt3b","Mogs","Ckap4","Rrbp1","Sil1","Hyou1","Hspa5","Dnajb11","Dnajc1","Dnajc3","Dnajc10","Hsp90b1",
	"Ganab","Prkcsh","Canx","Pdia3","Calr","Calr4","Man1a","Man1a2","Man1c1","Man1b1","Lman2","Lman1","Lman1l","Preb","Sar1a","Sar1b","Sec13",
	"Sec31a","Sec31b","Sec23b","Sec23a","Sec24b","Sec24a","Sec24c","Sec24d","Uggt2","Uggt1","Edem1","Edem2","Edem3","P4hb","Pdia4","Pdia6","Erp29",
	"Txndc5","Ero1a","Ero1b","Os9","Erlec1","Ssr1","Ssr2","Ssr3","Ssr4","Bcap31","Tram1","Tram1l1","Derl1","Derl2","Derl3","Ubxn4","Ubxn1","Ubxn2a",
	"Ubxn6","Ubxn8","Nsfl1c","Svip","Vcp","Nploc4","Ufd1","Hspa8","Hspa1l","Hspa1a","Hspa2","Hspa1b","Dnaja1","Dnaja2","Dnajb1","Dnajb2","Dnajb12",
	"Dnajc5","Dnajc5b","Dnajc5g","Hsp90aa1","Hsp90ab1","Hsph1","Hspa4l","Bag1","Bag2","Hspbp1","Cryaa","Cryab","Yod1","Plaa","Rad23b","Rad23a",
	"Ubqln1","Ubqln3","Ubqln2","Ubqln4","Ngly1","Atxn3","Ube4b","Eif2ak1","Eif2ak2","Eif2ak3","Eif2ak4","Eif2s1","Nfe2l2","Atf4","Ppp1r15a",
	"Ddit3","Bcl2","Atf6","Atf6b","Wfs1","Mbtps1","Mbtps2","Xbp1","Ern1","Traf2","Map3k5","Map2k7","Mapk8","Mapk9","Mapk10","Bax","Bak1","Capn1",
	"Capn2","Casp12","Marchf6","Ube2j1","Ube2j2","Ube2g1","Ube2g2","Syvn1","Rnf5","Rnf185","Selenos","Sel1l","Sel1l2","Herpud1","Amfr","Stub1",
	"Ube2d1","Ube2d2a","Ube2d3","Ube2d2b","Ube2dnl1","Ube2dnl2","Prkn","Rbx1","Cul1","Skp1","Fbxo2","Fbxo6")

pdf("PP.pdf")
pheatmap(c,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         cellwidth = 8, cellheight = 7,fontsize=7,
         show_rownames=T,show_colnames=T)
dev.off()

ego_geneID <- bitr(rownames(c), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)


