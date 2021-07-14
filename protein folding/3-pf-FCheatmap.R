library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")
proj5 <- addImputeWeights(proj5,reducedDims="IterativeLSI-dim8")

FB3proj<-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/Save-fibroblast2")
EC4proj<-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
MP3proj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")
SMC1proj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")
p<-getGroupSE(
  ArchRProj = FB3proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = FB3proj,
  useMatrix = "GeneScoreMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
count<-as.data.frame(count)
count=as.data.frame(lapply(count,as.numeric))
rownames(count)<-gene
FBcount=t(scale(t(count),scale = T,center = T))

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
count<-as.data.frame(count)
count=as.data.frame(lapply(count,as.numeric))
rownames(count)<-gene
MPcount=t(scale(t(count),scale = T,center = T))
#write.table(count,"PF-FCmatrix-MP.txt")

p<-getGroupSE(
  ArchRProj = EC4proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "EC4proj$Subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = EC4proj,
  useMatrix = "GeneScoreMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
count<-as.data.frame(count)
count=as.data.frame(lapply(count,as.numeric))
rownames(count)<-gene
ECcount=t(scale(t(count),scale = T,center = T))
#write.table(count,"PF-FCmatrix-EC.txt")

#####SMC#######
p<-getGroupSE(
  ArchRProj = SMC1proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "subtype"
)
GM<-getMatrixFromProject(
  ArchRProj = SMC1proj,
  useMatrix = "GeneScoreMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
count<-as.data.frame(count)
count=as.data.frame(lapply(count,as.numeric))
rownames(count)<-gene
SMCcount=t(scale(t(count),scale = T,center = T))
#write.table(count,"PF-FCmatrix-SMC.txt")

PP_pathway<-c("Map3k5","Atf4","Atf6","Bak1","Bcap31","Bax","Bcl2","Capn1","Capn2","Casp12","Stub1","Ddit3","Canx",
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
protein_folding<-c("Dnajb1","Dnaja2","P3h1","Cct4","Ppid","Hsp90b1","Qsox1","Cryaa","Cct3","Ppih","Pfdn2","Vbp1","Dnajb11",
	"Fkbp6","Pfdn1","Clpx","Trap1","Ppib","Dnaja4","Pfdn4","Dnajb4","Cct6a","Pfdn5","P4hb","Zmynd10","Ppig","Ero1l","Pdrg1",
	"Tbcd","Cryab","Tcp1","Clgn","Ranbp2","Cwc27","Erp27","Grpel2","Hspe1","Hspd1","Cdc37l1","Ptges3l","Pdcl","Mesd",
	"Hsp90aa1","Ppil3","Dnajc25","Hsp90ab1","Dnajb13","Hspa8","Qsox2","Calr4","Fkbp8","Cdc37","Ppic","Nktr","Ppif","Dnajb5",
	"Fkbp9","Pdilt","Cct8l1","Mfsd13b","Dnajb6","Hspa9","Grpel1","Ppie","Pfdn6","Calr3","Tbce","Mkks","Dnaja3","Nudcd2",
	"Ppia","Cct7","Dnaja1","Nudcd3","Dnajc1","Ube4b","Cct8","Dnlz","Hspa4l","Ptges3","Tbca","Ric3","Cct5","Erp44","Ppil2",
	"Nppa","Hspa1b","Hspa1a","Pdia2","Pdia3","Nppc","Cct2","Pdcl3","Cd74","Ahsa1","Calr","Ppil1","Cct6b","Hspe1-rs1","Canx",
	"Stub1","Nudc","Gm7879","Hspa4l")

allcount<-cbind(FBcount,MPcount)
allcount<-cbind(allcount,ECcount)
allcount<-cbind(allcount,SMCcount)
pfcount<-allcount[which(rownames(allcount)%in% protein_folding),]
#######protein_folding genes heatmap in FB MP EC ####
library(pheatmap)
# 添加注释信息
annotation_col = data.frame(
  celltype =factor( c(rep("FB",4),rep("MP",4),rep("EC",5),rep("SMC",3)  ))
    )
rownames(annotation_col) = colnames(pfcount)
head(annotation_col)

ann_colors = list(
  celltype = c(FB="#E58965",MP="#519A6D",EC="#6B58B3",SMC="#945E9C")
 )
head(ann_colors)
pdf("protein_folding-heatmap_for4celltypes.pdf")
pheatmap(pfcount,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 4,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         #annotation_row = annotation_row,
         show_rownames=F,show_colnames=T)
dev.off();

#######ORDER HEATMAP####
fb<-pfcount[which(pfcount[,3]>0.5),]
fb<-fb[order(fb[,3]),]
ec<-pfcount[which(pfcount[,12]>0.5),]
ec<-ec[order(ec[,12]),]
mp<-pfcount[which(pfcount[,7]>0.5),]
mp<-mp[order(mp[,7]),]
smc<-pfcount[which(pfcount[,14]>0.5),]
smc<-smc[order(smc[,14]),]

fb<-rownames(fb)
ec<-rownames(ec)
mp<-rownames(mp)
smc<-rownames(smc)
######high in EC4 FB3 MP3######
overlap<-intersect(ec,fb)
overlap<-intersect(overlap,mp)
overlap<-intersect(overlap,smc)
pdf("protein_folding-heatmap_for_overlap.pdf")
overlap
pheatmap(pfcount[overlap,],cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         #annotation_row = annotation_row,
         show_rownames=T,show_colnames=T)
dev.off();
