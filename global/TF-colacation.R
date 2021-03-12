#######获取各个细胞类型中peak在基因组上的位置信息######
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
FB3proj<-loadArchRProject(path = "Save-fibroblast2")
EC4proj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
MP3proj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")
SMC1proj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")

#######FB######
markersPeaks <- getMarkerFeatures(ArchRProj = FB3proj, useMatrix = "PeakMatrix", 
                                  groupBy = "subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList
markerList$FB1[,8]<-1:length(rownames(markerList$FB1))
markerList$FB2[,8]<-1:length(rownames(markerList$FB2))
markerList$FB3[,8]<-1:length(rownames(markerList$FB3))
markerList$FB4[,8]<-1:length(rownames(markerList$FB4))

write.table(markerList$FB1[,c(8,1,3,4)],"FB1-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$FB2[,c(8,1,3,4)],"FB2-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$FB3[,c(8,1,3,4)],"FB3-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$FB4[,c(8,1,3,4)],"FB4-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)

#####SMC###########

markersPeaks <- getMarkerFeatures(ArchRProj = SMC1proj, useMatrix = "PeakMatrix", 
                                  groupBy = "subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList
markerList$SMC1[,8]<-1:length(rownames(markerList$SMC1))
markerList$SMC2[,8]<-1:length(rownames(markerList$SMC2))
markerList$SMC3[,8]<-1:length(rownames(markerList$SMC3))
write.table(markerList$SMC1[,c(8,1,3,4)],"SMC1-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$SMC2[,c(8,1,3,4)],"SMC2-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$SMC3[,c(8,1,3,4)],"SMC3-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)

##############MP#############
markersPeaks <- getMarkerFeatures(ArchRProj = MP3proj, useMatrix = "PeakMatrix", 
                                  groupBy = "MP3proj$Subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList
markerList$MP1[,8]<-1:length(rownames(markerList$MP1))
markerList$MP2[,8]<-1:length(rownames(markerList$MP2))
markerList$MP3[,8]<-1:length(rownames(markerList$MP3))
markerList$MP4[,8]<-1:length(rownames(markerList$MP4))

write.table(markerList$MP1[,c(8,1,3,4)],"MP1-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$MP2[,c(8,1,3,4)],"MP2-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$MP3[,c(8,1,3,4)],"MP3-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$MP4[,c(8,1,3,4)],"MP4-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)

markerList$EC2[,8]<-1:length(rownames(markerList$EC2))
markerList$EC4[,8]<-1:length(rownames(markerList$EC4))
markerList$EC5[,8]<-1:length(rownames(markerList$EC5))
write.table(markerList$EC2[,c(8,1,3,4)],"EC2-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$EC4[,c(8,1,3,4)],"EC4-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)




sed -i 's/"//g' *
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' FB2-FDR03-peak.txt>FB2-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' FB3-FDR03-peak.txt>FB3-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' FB4-FDR03-peak.txt>FB4-FDR03-peak.bed

awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' SMC2-FDR03-peak.txt>SMC2-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' SMC3-FDR03-peak.txt>SMC3-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' SMC1-FDR03-peak.txt>SMC1-FDR03-peak.bed

awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' MP4-FDR03-peak.txt>MP4-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' MP1-FDR03-peak.txt>MP1-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' MP2-FDR03-peak.txt>MP2-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' MP3-FDR03-peak.txt>MP3-FDR03-peak.bed

awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' EC2-FDR03-peak.txt>EC2-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' EC4-FDR03-peak.txt>EC4-FDR03-peak.bed
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' EC5-FDR03-peak.txt>EC5-FDR03-peak.bed

annotatePeaks.pl SMC1-FDR03-peak.bed mm10 > SMC1_peak.annotation.xls
annotatePeaks.pl SMC3-FDR03-peak.bed mm10 > SMC3_peak.annotation.xls
annotatePeaks.pl SMC2-FDR03-peak.bed mm10 > SMC2_peak.annotation.xls

annotatePeaks.pl MP1-FDR03-peak.bed mm10 > MP1_peak.annotation.xls
annotatePeaks.pl MP3-FDR03-peak.bed mm10 > MP3_peak.annotation.xls
annotatePeaks.pl MP2-FDR03-peak.bed mm10 > MP2_peak.annotation.xls
annotatePeaks.pl MP4-FDR03-peak.bed mm10 > MP4_peak.annotation.xls

annotatePeaks.pl FB3-FDR03-peak.bed mm10 > FB3_peak.annotation.xls
annotatePeaks.pl FB2-FDR03-peak.bed mm10 > FB2_peak.annotation.xls
annotatePeaks.pl FB4-FDR03-peak.bed mm10 > FB4_peak.annotation.xls

annotatePeaks.pl EC5-FDR03-peak.bed mm10 > EC5_peak.annotation.xls
annotatePeaks.pl EC2-FDR03-peak.bed mm10 > EC2_peak.annotation.xls
annotatePeaks.pl EC4-FDR03-peak.bed mm10 > EC4_peak.annotation.xls


######EC3中5个TF家族中各个成员的RNA表达情况##########

getAvailableMatrices(EC4proj)

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

#########
AP1<-c("Fosl1","Fosl2","Junb","Jund","Jun","Fosb","Fos")
Bach<-c("Bach1","Bach2")
Nfe<-c("Nfe2l1","Nfe2","Nfe2l2","Nfe2l3")
Batf<-c("Batf","Batf2","Batf3")

             EC1         EC2         EC3        EC4        EC5
Batf3 0.03829594 0.021542535 0.026357862 0.08265418 0.13356197
Batf  0.01212321 0.005288383 0.053219333 0.01998048 0.12970402
Batf2 0.01054446 0.005876381 0.009947658 0.01406209 0.02645608


getAvailableMatrices(FB3proj)


library(Seurat)
scRNA <- readRDS("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/integrative_analysis/cmc_sct.rds")

selRNA <- scRNA[,which(scRNA@meta.data$orig.ident == "A3" |
                         scRNA@meta.data$orig.ident == "A7" | scRNA@meta.data$orig.ident == "A14" |
                         scRNA@meta.data$orig.ident == "P7" | scRNA@meta.data$orig.ident == "P14" )]
table(selRNA@meta.data$orig.ident)

selRNA<-selRNA[,which(selRNA@meta.data$cell_types == "Fibroblasts" )]

FB3<- addGeneIntegrationMatrix(
  ArchRProj = FB3proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  seRNA = selRNA,
  addToArrow = TRUE,
  groupRNA = "cell_types",
  nameCell = "predictedCell_Un",
  force = TRUE,
  reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
p<-getGroupSE(
  ArchRProj = FB3proj,
  useMatrix = "GeneIntegrationMatrix",
  groupBy = "subtype"
)
GM<-getMatrixFromProject(
  ArchRProj =FB3proj,
  useMatrix = "GeneIntegrationMatrix"
)
gene<-rowData(GM)$name
count<-assay(p)
rownames(count)<-gene
Cebp<-c("Cebpa","Cebpb","Cebpd","Cebpe","Cebpz","Cebpg")

               FB1          FB2          FB3          FB4
Cebpe 1.702905e-75 3.959651e-81 6.977537e-63 4.438249e-62
Cebpz 1.456573e-45 2.970417e-41 5.041934e-51 5.536910e-51
Cebpb 5.302227e-44 4.642958e-39 4.093816e-48 2.356864e-48
Cebpg 1.368868e-45 8.573249e-41 7.622777e-51 8.454093e-51
Cebpa 1.519211e-47 2.896135e-42 5.288550e-51 2.507314e-51


#############TF 共定位信息#############

markersPeaks <- getMarkerFeatures(ArchRProj = EC4proj, useMatrix = "PeakMatrix", 
                                  groupBy = "EC4proj$Subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.3 & Log2FC >= 0.1")
markerList
EC3peak<-markerList$EC3
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = EC4proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.3 & Log2FC >= 0.1"
  )
AP1<-c("Fosl1","Fosl2","Junb","Jund","Jun","Fosb","Fos")
Bach<-c("Bach1","Bach2")
Nfe<-c("Nfe2l1","Nfe2","Nfe2l2","Nfe2l3")
Batf<-c("Batf","Batf3")

TFname<-c(AP1,Bach,Nfe,Batf)
TFid<-which(getPeakAnnotation(EC4proj)$motifSummary$name %in% TFname)
position<-readRDS("/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype/Annotations/Motif-Positions-In-Peaks.rds")
library(GenomicRanges)
peak<-GRanges(
	seqnames=EC3peak$seqnames,
	ranges=IRanges(start=EC3peak$start,end=EC3peak$end,names = 1:546	))
m<-matrix(nrow=546,ncol=15)
rownames(m)<-names(peak)
colnames(m)<-getPeakAnnotation(EC4proj)$motifSummary$name[TFid]

#########构建共定位矩阵##############
p=0;

for  (i in TFid){
  p<- p+1;
  print(p)
  tf<-position[[i]];
  hits<-findOverlaps(peak,tf)@ from
  for (j in 1:546){
  	if (j %in% hits )
  	{m[j,p]=1;}
  	else
  	{m[j,p]=0;}
  }
 }
o<-m[which(rowSums(m)>0),]

library(pheatmap)
pdf("TF-colocation-1.pdf")

out<-pheatmap(o,cluster_cols = T,cluster_rows = T,
         color = c("white", "red"),cutree_rows = 7,
         show_rownames=FALSE,show_colnames=TRUE)
rownames(o[out$tree_row[["order"]],])
row_cluster <- cutree(out$tree_row,k=7)
Cluster2<-peak[which(row_cluster==2)]

#######peak 注释到基因########
library(org.Mm.eg.db)
peakAnnoList <- lapply(Cluster2, TxDb=txdb,tssRegion=c(-1000, 1000), 
  verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Mm.eg.db")


annotatePeak, 
aCR<-assignChromosomeRegion(sampleA, nucleotideLevel=FALSE,  precedence=c("Promoters", "immediateDownstream","fiveUTRs", "threeUTRs", "Exons", "Introns"),
  TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)barplot(aCR$percentage, las=3)






#########EC1 target gene heatmap ###########



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
rownames(count)<-gene
Bach1<-c("4930542M03Rik","4930572O13Rik","AI847159","Ano6","Apoa1","Arf1","Arhgef12","Atp2a2","Babam2","Bcl2l14","Capg","Cdip1","Chd7","Cmip","Col17a1",
	"Col18a1","Col4a2","Crybg3","Ctla2a","Ctnnbip1","Cyp11a1","Dst","Fam220a","Flt1","Fosl1","Fscn1","Gm11978","Gm15713","Gm20544","Gm35206","Gm8615",
	"Herc2","Il20rb","Irf2bpl","Jdp2","Kdm6b","Lpp","Map4k4","Micall2","Mir12182","Mir23a","Ndrg1","Nfib","Orai1","Oser1","Pak6","Pald1","Papss2","Parvb",
	"Peak1","Pip5k1c","Plk2","Plpp1","Plxna1","Prickle2","Prr5l","Ptk2","Ptprg","Scd2","Sdcbp2","Sema7a","Sipa1l3","Snx11","Spag9","Sparc","Spata24","Synpo",
	"Tbc1d1","Tbce","Tet2","Tgfbr1","Tmbim1","Tmem44","Tnfaip2","Traf1","Trio","Trp53cor1","Wdr1","Wdr60","Zfand3")
heat<-count[which(rownames(count)%in% Nfe2),]
heat=t(scale(t(heat),scale = T,center = T))
heat1<-na.omit(heat)
pdf("EC1-Nfe2—targetgene-heatmap.pdf")
pheatmap(heat1,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #ellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()


Batf<-c("0610040F04Rik","0610040F04Rik","1110002J07Rik","1110036E04Rik","1600010M07Rik","1600020E01Rik","1700001C19Rik","1700006H21Rik","1700007L15Rik","1700016G22Rik","1700020M21Rik","1700020N01Rik","1700022A22Rik","1700034J05Rik","1700096J18Rik","1700099I09Rik","1700099I09Rik","1700099I09Rik","1700112H15Rik","1700113A16Rik","2510009E07Rik","2610507B11Rik","2810025M15Rik","4930431P22Rik","4930445N18Rik","4930459C07Rik","4930511A08Rik","4930523C07Rik","4930542M03Rik","4930554I06Rik","4930572O13Rik","4933430N04Rik","4933433G15Rik","5031425E22Rik","5830416I19Rik","6030407O03Rik","6820426E19Rik","9330111N05Rik","A330009N23Rik","A530006G24Rik","A530013C23Rik","A730017L22Rik","A930003A15Rik","Aak1","Abca1","Abca1","Abca4","Abca4","Ablim1","Ablim1","Ablim1","Ablim3","Actb","Acvr1","Adamts2","Adamts9","Adarb1","Adgrg1","Adgrg3","Adgrl2","Adgrl2","Adipor2","Aebp1","Afap1l1","Agfg2","Ahctf1","Ahnak","Ahnak","AI847159","Akap2","Akirin2","Akt3","Amotl1","Ampd3","Amz1","Angpt2","Angptl4","Ankrd11","Ankrd12","Ankrd2","Ankrd33b","Ankrd55","Anks1","Ano6","Ano6","Apbb2","Apoa1","App","Aqp1","Aqp7","Aqp7","Arap3","Arap3","Arf1","Arf6","Arhgap28","Arhgap31","Arhgap31","Arhgap35","Arhgef10l","Arhgef12","Arhgef2","Arhgef3","Arhgef5","Arhgef7","Arid5a","Arid5b","Armc2","Arntl","Arpc3","Arsa","Asah2","Asap2","Asap2","Asap2","Aspa","Atp2a2","Atp2b4","B230217C12Rik","B230319C09Rik","B3gnt2","Baalc","Babam2","Bach2","Bach2","Bahcc1","Bak1","Bbs12","Bbs9","Bcar3","Bcas1os2","Bcl2l11","Bcl2l14","Bcor","Bin3","C130050O18Rik","C130071C03Rik","C1qtnf9","Cacng3","Camk2a","Capg","Card10","Card10","Carmn","Casz1","Cavin1","Ccdc12","Ccdc50","Ccdc68","Ccdc81","Ccdc81","Ccdc97","Cd300lg","Cd44","Cd44","Cd68","Cd83","Cd83","Cd93","Cdc42bpg","Cdc42bpg","Cdc42ep2","Cdc42ep3","Cdca7l","Cdh15","Cdh15","Cdh5","Cdip1","Cdk14","Cdkl4","Cds2","Ceacam16","Ceacam16","Cel","Cemip2","Cflar","Cflar","Cflar","Chd7","Chd7","Chd9","Chst1","Chst15","Chst2","Chst2","Cirbp","Cldn5","Cldn5","Cltb","Cmip","Cmip","Cnksr3","Cnn3","Col15a1","Col17a1","Col18a1","Col18a1","Col4a1","Col4a1","Col4a2","Col4a2","Col4a2","Coro6","Cpsf4l","Creb3l2","Crybg3","Crybg3","Crybg3","Csgalnact1","Csnk1d","Cspg4","Csrnp1","Csrp1","Cst6","Ctla2a","Ctnnbip1","Cx3cl1","Cx3cl1","Cxcl17","Cxxc5","Cyp11a1","Cyp11a1","Cyp2d22","D330050G23Rik","D830013O20Rik","Daam1","Dab2","Dab2ip","Dclk2","Ddi1","Dgke","Dhrs3","Dhrs3","Dhrs3","Dhrs3","Dip2c","Dip2c","Dipk2a","Disc1","Dlc1","Dock1","Dock4","Dock4","Dock4","Dock4","Dpysl3","Dst","Dusp1","Dusp1","Dusp1","Dusp1","Dusp5","Dyrk1a","E030018B13Rik","E230013L22Rik","E230013L22Rik","Eaf1","Efna1","Egr1","Egr3","Egr3","Ehd1","Eif1","Eif2ak3","Elmo1","Elmo1","Elmsan1","Eml1","Eng","Entpd1","Ercc1","Erg","Ergic1","Esrrb","Ets1","Ets1","Etv6","Etv6","Exoc5","Ext1","Ext1","F2r","Fa2h","Faap20","Fads3","Fam110a","Fam129a","Fam220a","Fam43a","Fat4","Fbn1","Fbxo28","Fez2","Fez2","Fgd5","Fgd5","Filip1","Filip1","Fkbp1a","Fkbp1a","Flnb","Flnb","Flt1","Flt1","Flt1","Flt1","Fmo9","Fos","Fosb","Fosl1","Fosl1","Foxc1","Foxo1","Foxp2","Frmd4a","Frmd4a","Frmd4a","Frmd4a","Frmd4a","Fscn1","Fscn1","Furin","Gadd45b","Gadd45g","Galnt15","Gas7","Gcnt2","Gde1","Gdpd5","Gfpt2","Ghr","Ghsr","Git2","Gli2","Glrx","Gm10389","Gm10415","Gm10532","Gm10532","Gm10532","Gm10556","Gm10941","Gm10941","Gm10941","Gm11978","Gm12216","Gm13872","Gm14634","Gm15569","Gm15713","Gm16701","Gm16701","Gm16793","Gm17619","Gm18409","Gm18409","Gm19276","Gm19510","Gm19557","Gm20544","Gm20750","Gm29682","Gm29682","Gm29684","Gm30108","Gm31508","Gm32168","Gm35206","Gm36221","Gm45924","Gm4832","Gm5069","Gm5129","Gm5820","Gm6116","Gm6117","Gm7444","Gm8615","Gm8817","Gnai2","Gns","Gosr2","Gpr15","Gprc5b","Gpx4","Grap","Grrp1","Gse1","Gse1","Gsg1","Gxylt1","H1foo","H1foo","H2afy3","Halr1","Herc2","Hexa","Hhex","Hif1a","Hip1r","Hk1","Hk1","Hk1","Hmbox1","Hmga1","Hmox1","Hnrnpl","Hook2","Hotairm1","Hoxb1","Hsd3b7","Hsf2bp","Hsp90aa1","Hspg2","Hspg2","Ibtk","Id1","Id3","Ier5","Ifitm10","Ifitm5","Ift122","Il20rb","Il34","Inhbb","Inpp5d","Insyn2a","Irf2","Irf2","Irf2bpl","Irx3","Itga2","Itga3","Itgav","Itgav","Itgav","Itgb1","Itpkb","Itprid2","Itpripl2","Jak1","Jcad","Jdp2","Kansl1l","Kdm1a","Kdm6b","Kdm6b","Kdr","Kdr","Kdr","Kdr","Kdr","Kidins220","Kitl","Kitl","Klf10","Klf6","Klhl21","Knop1","Ksr1","Lama4","Lamc1","Laptm5","Lfng","Lhx6","Lifr","Lipc","Lipc","Lmna","Lmnb1","Lmnb1","Lmod3","Lnx1","Lnx1","LOC102634389","Lpar6","Lpp","Lrrc32","Lrrc32","Lrrc8c","Lrrfip2","Lrrtm2","Lysmd4","Mab21l3","Maff","Mafk","Man1c1","Manba","Map2k6","Map3k1","Map4k4","Map7d1","Mapkap1","Marcks","Mark2","Mbnl2","Mcl1","Mdfi","Me3","Me3","Med13","Med13l","Med13l","Mef2c","Mertk","Mex3c","Mfn2","Mfng","Mfsd2a","Micall2","Midn","Mir1190","Mir12182","Mir12196","Mir12198","Mir139","Mir139","Mir181b-1","Mir1901","Mir1904","Mir1943","Mir1954","Mir1956","Mir1962","Mir1962","Mir21a","Mir21c","Mir23a","Mir3074-2","Mir3092","Mir3109","Mir342","Mir3962","Mir3965","Mir5098","Mir574","Mir6368","Mir6377","Mir6391","Mir6481","Mir6897","Mir6899","Mir6951","Mir6951","Mir7008","Mir7008","Mir7018","Mir7018","Mir7032","Mir705","Mir7077","Mir7089","Mir7090","Mir762","Mir7659","Mir7661","Mir7681","Misp","Mknk1","Mlxip","Mmd","Mmp15","Mnt","Mocs2","Mrm3","Mrpl42","Msh5","Msn","Myo10","Myo1b","Myo1d","Myo1e","Myo1e","Myo5b","Myo9b","Nadk","Nav2","Nav3","Ncor2","Ncor2","Ndrg1","Ndrg1","Nedd4","Nedd9","Nek6","Neu3","Nfib","Nfib","Nfkb1","Nfyb","Nhlrc3","Nhsl1","Nhsl1","Nmt2","Noct","Nr1d2","Nrp2","Nrp2","Nt5c3","Ntn5","Nudt16l1","Oacyl","Oaz2","Olfm2","Omg","Orai1","Oser1","P2ry2","Pak6","Pald1","Pan3","Papss2","Parp14","Parvaos","Parvb","Pcdh1","Pcdh1",
	"Pcdh1","Pcdh12","Pde4b","Pde4b","Pdgfa","Pdgfb","Pdlim7","Pdzd2","Peak1","Peak1","Pecam1","Pecam1","Pecam1","Pfdn1","Pfkfb3","Pfn1","Pgf","Pgk1","Pgpep1","Pgs1","Phf2","Phldb1","Piezo1","Pik3r5","Pip5k1c","Pitpnc1","Pitpnc1","Pkm","Pla2g5","Plaat5","Plaur","Plbd2","Plcg1","Plekhg1","Plekhg1","Plekhg1","Plk2","Pllp","Plp2","Plpp1","Plpp1","Plpp1","Plpp6","Pls3","Plxna1","Plxna1","Plxna1","Plxna1","Plxnd1","Pop5","Por","Ppp1r9b","Ppp2r3c","Ppp3ca","Pqlc1","Prelp","Prex1","Prg2","Prickle2","Prkch","Prnp","Prpf39","Prr13","Prr5l","Prr5l","Prss52","Ptk2","Ptpa","Ptprb","Ptpre","Ptprg","Ptprm","Ptprm","Ptprm","Pvr","Pwwp2b","Rab8b","Ralb","Ralgds","Rapgef3","Rapgef4","Rara","Rasa4","Rassf3","Rassf8","Rassf9","Rbbp6","Rbfox1","Rbp7","Resf1","Rftn1","Rfx5","Rfx8","Rhoc","Rhoq","Rin2","Rin3","Rmi2","Rnf125","Rnf130","Rnf157","Rnf19a","Rnf216","Rnf24","Rpl37rt","Rplp0","Rplp0","Rptoros","Rptoros","Rptoros","Rreb1","Rreb1","Rtkn2","S1pr1","Sash1","Saysd1","Saysd1","Scarb1","Scarletltr","Scd2","Scd4","Scgb1a1","Schip1","Scimp","Scnn1a","Sdc3","Sdc3","Sdcbp2","Sec14l1","Sec24d","Sema3f","Sema6d","Sema7a","Sema7a","Serpind1","Sgce","Sgk1","Sgk1","Sgk1","Sgpl1","Sh3glb1","Sh3pxd2b","Sh3tc1","Shank3","Shroom3","Sik1","Sipa1l3","Skil","Slc16a11","Slc16a2","Slc22a16","Slc25a13","Slc25a24","Slc43a3","Slc44a1","Slc44a1","Smad6","Smad7","Smchd1","Smco4","Smpdl3a","Smtn","Smurf2","Snora3","Snord83b","Snrk","Snx11","Socs5","Soga1","Sos2","Sox11","Sox7","Spaca6","Spag9","Spag9","Sparc","Sparc","Spata13","Specc1l","Spop","Spop","Spsb1","Sqstm1","Srgap1","Srsf3","Stk10","Stk10","Stk32c","Strap","Stum","Sumo3","Svil","Sybu","Syne1","Synj2","Synpo","Synpo","Tbc1d1","Tbc1d1","Tbc1d2","Tbc1d9b","Tbcc","Tbce","Tcf4","Tcf4","Tcof1","Tdp1","Tead1","Tek","Tekt1","Tes","Tet2","Tet2","Tet3","Tex43","Tgfbr1","Tgfbr1","Tgfbr1","Tgif1","Thbd","Ticam1","Timp4","Tinagl1","Tm4sf1","Tm4sf1","Tmbim1","Tmbim1","Tmem17","Tmem204","Tmem241","Tmem243","Tmem28","Tmem44","Tmem51os1","Tmpo","Tnfaip2","Tnfaip2","Tnfrsf10b","Tnfrsf1b","Tnfrsf23","Tnip1","Tnk1","Tnk2","Tnk2os","Tnni1","Tnnt3","Tns1","Tor4a","Tpcn1","Tpm3","Tpst1","Tpst2","Traf1","Traf1","Trib1","Trib2","Trim24","Trim6","Trim8","Trim8","Trio","Trp53cor1","Tsen2","Tshz2","Tspan14","Tspan18","Ttc28","Ttll7","Ttpal","Ttpal","Txndc11","Uaca","Ubl3","Ubtf","Usp6nl","Utrn","Vav1","Vcl","Vps25","Vps45","Wdr1","Wdr60","Wdr73","Wdr91","Wnt7b","Wsb2","Xaf1","Xdh","Xpo6","Ybx3","Ypel2","Ythdf3","Ywhag","Zbtb2","Zcchc14","Zeb1","Zeb1","Zeb1","Zfand3","Zfp1","Zfp219","Zfp469","Zfp541","Zfp572","Zfp608","Zfp608","Zfp827","Zfp957","Zfyve1","Znrf1","Znrf1")

Fosl2<-c("0610040F04Rik","1110002J07Rik","1600010M07Rik","1600020E01Rik","1700001C19Rik","1700016G22Rik","1700020M21Rik","1700020N01Rik","1700022A22Rik","1700034J05Rik","1700099I09Rik","1700099I09Rik","1700112H15Rik","1700113A16Rik","2510009E07Rik","2610507B11Rik","4930445N18Rik","4930459C07Rik","4930511A08Rik","4930523C07Rik","4930542M03Rik","4930554I06Rik","4930572O13Rik","4933433G15Rik","6030407O03Rik","6820426E19Rik","9330111N05Rik","A530013C23Rik","A930003A15Rik","Abca4","Abca4","Ablim1","Ablim1","Ablim3","Actb","Adamts2","Adamts9","Adarb1","Adgrg1","Adgrg3","Aebp1","Afap1l1","Agfg2","Ahctf1","Ahnak","Ahnak","AI847159","Akt3","Angpt2","Angptl4","Ankrd12","Ankrd2","Ankrd33b","Ankrd55","Ano6","Apoa1","App","Aqp1","Aqp7","Arap3","Arap3","Arf1","Arf6","Arhgap28","Arhgap31","Arhgap31","Arhgef10l","Arhgef12","Arhgef2","Arhgef3","Arhgef5","Arid5a","Armc2","Arntl","Arsa","Asap2","Asap2","Aspa","Atp2a2","Atp2b4","B230319C09Rik","Baalc","Babam2","Bach2","Bach2","Bahcc1","Bak1","Bbs12","Bbs9","Bcar3","Bcl2l14","Bcor","Bin3","C130050O18Rik","C130071C03Rik","Cacng3","Camk2a","Capg","Card10","Card10","Carmn","Cavin1","Ccdc12","Ccdc50","Ccdc68","Ccdc81","Ccdc97","Cd300lg","Cd68","Cd83","Cd83","Cdc42bpg","Cdc42bpg","Cdc42ep2","Cdc42ep3","Cdca7l","Cdh15","Cdh5","Cdip1","Cdkl4","Ceacam16","Cel","Cemip2","Cflar","Cflar","Chd7","Chd9","Chst2","Cirbp","Cldn5","Cltb","Cmip","Cmip","Col15a1","Col17a1","Col18a1","Col18a1","Col4a1","Col4a1","Col4a2","Col4a2","Col4a2","Cpsf4l","Crybg3","Csgalnact1","Csnk1d","Cspg4","Csrnp1","Cst6","Ctla2a","Ctnnbip1","Cx3cl1","Cx3cl1","Cxcl17","Cyp11a1","Cyp11a1","Cyp11a1","D330050G23Rik","D830013O20Rik","Daam1","Dclk2","Dhrs3","Dhrs3","Dhrs3","Dip2c","Disc1","Dlc1","Dock4","Dock4","Dpysl3","Dst","Dusp1","Dusp1","Dyrk1a","E030018B13Rik","E230013L22Rik","E230013L22Rik","Efna1","Egr1","Egr3","Elmo1","Elmsan1","Eml1","Entpd1","Ercc1","Erg","Ets1","Etv6","Etv6","Ext1","F2r","Faap20","Fads3","Fam110a","Fam129a","Fam220a","Fam43a","Fat4","Fez2","Fez2","Fgd5","Filip1","Flnb","Flt1","Flt1","Flt1","Fmo9","Fos","Fosb","Fosl1","Foxc1","Foxo1","Foxp2","Frmd4a","Frmd4a","Frmd4a","Fscn1","Fscn1","Furin","Gadd45b","Gadd45g","Galnt15","Gcnt2","Gde1","Gdpd5","Gfpt2","Ghr","Git2","Glrx","Gm10532","Gm10556","Gm10941","Gm10941","Gm11978","Gm13872","Gm14634","Gm15569","Gm15713","Gm16793","Gm17619","Gm18409","Gm19276","Gm19557","Gm20544","Gm20750","Gm29682","Gm30108","Gm32168","Gm45924","Gm4832","Gm5069","Gm5129","Gm5820","Gm6116","Gm7444","Gm8615","Gm8817","Gnai2","Gprc5b","Gpx4","Grap","Grrp1","Gse1","Gxylt1","H1foo","H1foo","Herc2","Hexa","Hip1r","Hk1","Hk1","Hmga1","Hnrnpl","Hook2","Hotairm1","Hsd3b7","Hsf2bp","Hsp90aa1","Hspg2","Ibtk","Ifitm5","Il20rb","Il34","Inhbb","Irf2","Irf2","Irf2bpl","Irx3","Itga2","Itga3","Itgav","Itgb1","Itpkb","Itprid2","Jak1","Jdp2","Kansl1l","Kdm1a","Kdm6b","Kdm6b","Kdr","Kdr","Kdr","Kidins220","Kitl","Klf10","Knop1","Ksr1","Lama4","Lamc1","Laptm5","Lfng","Lhx6","Lifr","Lipc","Lipc","Lmna","Lnx1","LOC102634389","Lpp","Lrrc32","Lrrc32","Lrrfip2","Lysmd4","Maff","Man1c1","Manba","Map3k1","Map4k4","Marcks","Mark2","Mbnl2","Mcl1","Mdfi","Me3","Me3","Med13","Med13l","Med13l","Mertk","Mex3c","Mfn2","Mfsd2a","Micall2","Midn","Mir1190","Mir12182","Mir12196","Mir12198","Mir139","Mir139","Mir181b-1","Mir1904","Mir1943","Mir1954","Mir1956","Mir1962","Mir1962","Mir21c","Mir23a","Mir3092","Mir342","Mir3962","Mir5098","Mir574","Mir6377","Mir6481","Mir6899","Mir6951","Mir6951","Mir7008","Mir7008","Mir7018","Mir7018","Mir7032","Mir705","Mir7077","Mir762","Mir7659","Mir7661","Misp","Mknk1","Mlxip","Mmd","Mmp15","Mrm3","Mrpl42","Myo10","Myo1b","Myo1d","Myo1e","Myo9b","Nav3","Ncor2","Nedd4","Nek6","Neu3","Nfib","Nfkb1","Nfyb","Nhlrc3","Nhsl1","Nhsl1","Noct","Nrp2","Nt5c3","Ntn5","Nudt16l1","Olfm2","Omg","Orai1","Oser1","Pak6","Pald1","Pan3","Papss2","Parvaos","Parvb","Pcdh1")

	"Pcdh12","Pde4b","Pdgfa","Pdgfb","Pdlim7","Peak1","Peak1","Pecam1","Pecam1","Pfdn1","Pfkfb3","Pfn1","Pgk1","Pgs1","Phf2","Phldb1","Pik3r5","Pip5k1c","Pitpnc1","Pkm","Plaat5","Plbd2","Plcg1","Plekhg1","Plekhg1","Plp2","Plpp1","Plpp1","Plpp6","Pls3","Plxna1","Plxna1","Plxna1","Plxna1","Plxnd1","Pop5","Por","Ppp1r9b","Ppp3ca","Prelp","Prex1","Prg2","Prickle2","Prnp","Prpf39","Prr13","Prr5l","Prr5l","Prss52","Ptk2","Ptprb","Ptpre","Ptprg","Ptprm","Ptprm","Pvr","Pwwp2b","Ralgds","Rasa4","Rassf3","Rbfox1","Rftn1","Rin2","Rin3","Rmi2","Rnf125","Rnf157","Rnf19a","Rnf24","Rpl37rt","Rplp0","Rptoros","Rptoros","Rptoros","S1pr1","Sash1","Saysd1","Scarb1","Scarletltr","Scd2","Scd4","Scgb1a1","Scimp","Scnn1a","Sdc3","Sdc3","Sdcbp2","Sec24d","Sema6d","Sema7a","Sgce","Sgk1","Sgk1","Sgpl1","Sh3glb1","Sh3pxd2b","Sh3tc1","Shroom3","Sipa1l3","Skil","Slc16a2","Slc25a13","Slc25a24","Slc43a3","Slc44a1","Smad6","Smad7","Smco4","Smpdl3a","Smtn","Snora3","Snord83b","Snrk","Snx11","Soga1","Sorbs2","Sos2","Sox11","Sox7","Spaca6","Spag9","Sparc","Sparc","Spata13","Specc1l","Spsb1","Sqstm1","Srsf3","Stk10","Stk10","Stum","Svil","Sybu","Syne1","Synj2","Synpo","Synpo","Tbc1d1","Tbcc","Tbce","Tcf4","Tcf4","Tcof1","Tet2","Tet2","Tgfbr1","Tgfbr1","Tgif1","Ticam1","Timp4","Tm4sf1","Tmbim1","Tmbim1","Tmem17","Tmem241","Tmem44","Tmem51os1","Tmpo","Tnfaip2",
	"Tnfrsf1b","Tnfrsf23","Tnip1","Tnk1","Tnk2","Tnk2os","Tnni1","Tnnt3","Tns1","Tor4a","Tpcn1","Tpm3","Tpst2","Traf1","Traf1","Trib2","Trim24","Trim6","Trim8","Trio","Trp53cor1","Tsen2","Tspan14","Ttll7","Txndc11","Ubtf","Usp6nl","Utrn","Vav1","Vcl","Vps25","Wdr1","Wdr60","Wdr73","Wdr91","Wsb2","Xaf1","Ybx3","Ypel2","Ywhag","Zcchc14","Zeb1","Zeb1","Zfand3","Zfp1","Zfp219","Zfp469","Zfp541","Zfp572","Zfp608","Znrf1")


Junb<-c("5830454E08Rik","A630001O12Rik","Aak1","Adgrf5","Afap1","AI847159","Arap1","Atg16l2","Capzb","Carhsp1","Ccdc81","Col4a2","Dusp16","Dync2li1","Dyrk1a","Ednrb","Efnb1","Egr3","Elovl5","Entpd1","Esm1",
	"Gm15713","Gm5083","Gmpr","Gnl3l","Hip1r","Hs6st1","Ifitm10","Irak2","Jak1","Klf7","Lhx8","Lrrc8d","Mapk14","Mef2c","Mir1962","Mir6516","Mir7661","Mob2","Ncor2","Nfkbid","Nhsl2","Oaz2","Prg3","Rftn1","Rin3","Rptoros","Sdc3","Sema7a","Smad6","Smchd1","Snx11","Sox18","Srgap1","St5","Stum","Tet2","Tm2d1","Tnni1","Wwtr1","Zfp958")

Nfe2<-c("1700022A22Rik","4930542M03Rik","Aebp1","AI847159","Ano6","Apoa1","Arhgef12","Atp2a2","Babam2","Bach2","Capg","Cdc42bpg","Cdip1","Cds2","Cflar","Cflar","Chd7","Col17a1","Col18a1","Col4a1","Col4a2","Col4a2","Ctla2a","Ctnnbip1","Cyp11a1","Daam1","Dipk2a","Dst","Egr3","Fam220a","Fosl1","Frmd4a","Gfpt2","Gm11978","Gm15713","Gm35206","Gm8615","Herc2","Irf2bpl","Jdp2","Klf7","Lpp","Man1c1","Map4k4","Marcks","Mir12198","Ndrg1","Nek6","Nfib","Ntn5","Olfm2","Orai1","Oser1","Pak6","Pald1","Papss2","Parvb","Peak1","Pecam1","Phf21b","Plpp1","Plxna1","Ptk2","Ptprg","Scd2","Sema7a","Sipa1l3","Slc43a3","Smad7","Snx11","Spag9","Tac4","Tbc1d1","Tbce","Tgfbr1","Tmbim1","Tmem44","Tnfaip2","Tpm3","Trim6","Trio","Trp53cor1","Ttpal","Wdr60","Zfand3","Zmiz1")


