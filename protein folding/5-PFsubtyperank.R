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

protein_folding<-c("Dnajb1","Dnaja2","P3h1","Cct4","Ppid","Hsp90b1","Qsox1","Cryaa","Cct3","Ppih","Pfdn2","Vbp1","Dnajb11",
  "Fkbp6","Pfdn1","Clpx","Trap1","Ppib","Dnaja4","Pfdn4","Dnajb4","Cct6a","Pfdn5","P4hb","Zmynd10","Ppig","Ero1l","Pdrg1",
  "Tbcd","Cryab","Tcp1","Clgn","Ranbp2","Cwc27","Erp27","Grpel2","Hspe1","Hspd1","Cdc37l1","Ptges3l","Pdcl","Mesd",
  "Hsp90aa1","Ppil3","Dnajc25","Hsp90ab1","Dnajb13","Hspa8","Qsox2","Calr4","Fkbp8","Cdc37","Ppic","Nktr","Ppif","Dnajb5",
  "Fkbp9","Pdilt","Cct8l1","Mfsd13b","Dnajb6","Hspa9","Grpel1","Ppie","Pfdn6","Calr3","Tbce","Mkks","Dnaja3","Nudcd2",
  "Ppia","Cct7","Dnaja1","Nudcd3","Dnajc1","Ube4b","Cct8","Dnlz","Hspa4l","Ptges3","Tbca","Ric3","Cct5","Erp44","Ppil2",
  "Nppa","Hspa1b","Hspa1a","Pdia2","Pdia3","Nppc","Cct2","Pdcl3","Cd74","Ahsa1","Calr","Ppil1","Cct6b","Hspe1-rs1","Canx",
  "Stub1","Nudc","Gm7879","Hspa4l")


markersGS <- getMarkerFeatures(
    ArchRProj = FB3proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
FBpval<-assays(markersGS)$Pval
rownames(FBpval)<-rowData(markersGS)$name

markersGS <- getMarkerFeatures(
    ArchRProj = MP3proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "MP3proj$Subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
MPpval<-assays(markersGS)$Pval
rownames(MPpval)<-rowData(markersGS)$name

markersGS <- getMarkerFeatures(
    ArchRProj = EC4proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "EC4proj$Subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
ECpval<-assays(markersGS)$Pval
rownames(ECpval)<-rowData(markersGS)$name

markersGS <- getMarkerFeatures(
    ArchRProj = SMC1proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
SMCpval<-assays(markersGS)$Pval
rownames(SMCpval)<-rowData(markersGS)$name

pval<-cbind(FBpval,ECpval)
pval<-cbind(pval,MPpval)
pval<-cbind(pval,SMCpval)

#####only use 27 overlap protein folding genes######
overlap<-c("Trap1","Cct2","Dnajb5","Dnajb11","Ppia","Ppil1","Stub1","Cdc37","Mkks","Grpel1","Fkbp8","Calr","P3h1",
 "Ptges3l","Dnajb1","Hsp90ab1","Hsp90aa1","Hspd1","Pdcl3","Nudc","Dnaja1","Pfdn6","Mesd","Hspa1b","Hspe1","Hspa8","Hspa1a")
pval<-pval[which(rownames(pval)%in% protein_folding),]
p<--log2(colMeans(pval))
      FB1       FB2       FB3       FB4       EC1       EC2       EC3       EC4 
1.7788596 1.2367612 3.1759865 1.0761329 1.1676111 1.3945948 2.0981198 2.0134154 
      EC5       MP1       MP2       MP3       MP4      SMC1      SMC2      SMC3 
1.4479617 1.2893713 0.9134264 2.1849630 1.2323930 3.8401142 1.9767158 1.6833626 



FBmarker<-getMarkers(markersGS,cutOff = "FDR <= 0.05 & Log2FC >= 0.1",returnGR = F)
FBmarker<-as.data.frame(FBmarker)
FB<-FBmarker[which(FBmarker$name %in% protein_folding),]

ECmarker<-getMarkers(markersGS,cutOff = "FDR <= 0.05 & Log2FC >= 0.1",returnGR = F)
ECmarker<-as.data.frame(ECmarker)
EC<-ECmarker[which(ECmarker$name %in% protein_folding),]

SMCmarker<-getMarkers(markersGS,cutOff = "FDR <= 0.05 & Log2FC >= 0.1",returnGR = F)
SMCmarker<-as.data.frame(SMCmarker)
SMC<-SMCmarker[which(SMCmarker$name %in% protein_folding),]


MPmarker<-getMarkers(markersGS,cutOff = "FDR <= 0.05 & Log2FC >= 0.1",returnGR = F)
MPmarker<-as.data.frame(MPmarker)
MP<-MPmarker[which(MPmarker$name %in% protein_folding),]

data<-data.frame()


#######27 overlap protein folding gene#######
             p       ratio RANK
MP3  2.1849630 0.020689655    1
EC4  2.0134154 0.018867925    2
SMC1 3.8401142 0.008201203    3
FB3  3.1759865 0.006564551    4
FB1  1.7788596 0.000000000    5
FB2  1.2367612 0.000000000    6
FB4  1.0761329 0.000000000    7
EC1  1.1676111 0.000000000    8
EC2  1.3945948 0.000000000    9
EC3  2.0981198 0.000000000   10
MP1  1.2893713 0.000000000   11
MP2  0.9134264 0.000000000   12
MP4  1.2323930 0.000000000   13
SMC2 1.9767158 0.000000000   14
SMC3 1.6833626 0.000000000   15

ggplot(data,aes(x=RANK,y=ratio))+
  geom_point(aes(size=p,color=color),alpha=0.6)+
  #scale_size(range=c(1,12))+
  theme_bw()+
  theme( )+
  geom_text_repel(
    data = data,
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE )+theme(panel.grid=element_blank())

######all protein folding gene ######
  p(mlog2pval)        ratio RANK   ID color
FB1  1.285042 0.0000000000   13  FB1 black
FB2  1.240302 0.0025445293    6  FB2 black
FB3  1.665372 0.0098468271    4  FB3   red
FB4  1.076701 0.0000000000   14  FB4 black
EC1  1.327418 0.0025641026    5  EC1 black
EC2  1.122320 0.0014771049   10  EC2 black
EC3  1.676403 0.0008474576   12  EC3 black
EC4  1.346810 0.0226415094    1  EC4   red
MP1  1.197658 0.0020408163    7  MP1 black
MP2  1.085216 0.0013642565   11  MP2 black
MP3  1.296905 0.0206896552    2  MP3   red
MP4  1.304237 0.0016992353    9  MP4 black
SMC1 2.084002 0.0103881903    3 SMC1   red
SMC2 1.407202 0.0019821606    8 SMC2 black
SMC3 1.219885 0.0000000000   15 SMC3 black






#######计算Fold Enrichment:GeneRatio/BgRatio########

inPF<-c(0,1,18,0,2,1,1,6,1,1,6,2,19,4,0)
marker<-c(181,393,1828,455,780,677,1180,265,490,733,290,1177,1829,2018,1185)
noPF<-23210-104
PF<-104
data<-data.frame(inPF=inPF,marker=marker,noPF=noPF,PF=PF)
data$PFm1<-data$inPF-1

data$p<-1-phyper(data$PFm1,104,23106,data$marker)
data$mlog2pval<- -log2(data$p)
data$geneRatio<-data$inPF/data$marker
data$BgRatio<-104/23106
data$FoldEnrichment<-data$geneRatio/data$BgRatio


> data
   inPF marker  noPF  PF PFm1            p    mlog2pval    geneRatio      BgRatio FoldEnrichment  
1     0    181 23106 104   -1 1.0000000000  0.000000000 0.0000000000  0.004500995      0.0000000  
2     1    393 23106 104    0 0.8313670497  0.266442525 0.0025445293  0.004500995      0.5653259  
3    18   1828 23106 104   17 0.0012111857  9.689364237 0.0098468271  0.004500995      2.1876999  
4     0    455 23106 104   -1 1.0000000000  0.000000000 0.0000000000  0.004500995      0.0000000  
5     2    780 23106 104    1 0.8686525342  0.203148888 0.0025641026  0.004500995      0.5696746  
6     1    677 23106 104    0 0.9542971848  0.067489478 0.0014771049  0.004500995      0.3281729  
7     1   1180 23106 104    0 0.9956558122  0.006280991 0.0008474576  0.004500995      0.1882823  
8     6    265 23106 104    5 0.0012485866  9.645488445 0.0226415094  0.004500995      5.0303338  
9     1    490 23106 104    0 0.8918358165  0.165149955 0.0020408163  0.004500995      0.4534144  
10    1    733 23106 104    0 0.9647389554  0.051789472 0.0013642565  0.004500995      0.3031011  
11    6    290 23106 104    5 0.0019706445  8.987116708 0.0206896552  0.004500995      4.5966844  
12    2   1177 23106 104    1 0.9709988414  0.042458521 0.0016992353  0.004500995      0.3775243  
13   19   1829 23106 104   18 0.0004535135 11.106566887 0.0103881903  0.004500995      2.3079762  
14    4   2018 23106 104    3 0.9832866878  0.024315984 0.0019821606  0.004500995      0.4403827  
15    0   1185 23106 104   -1 1.0000000000  0.000000000 0.0000000000  0.004500995      0.0000000  

data$ID<-c("FB1","FB2","FB3","FB4","EC1","EC2","EC3","EC4","MP1","MP2","MP3","MP4","SMC1","SMC2","SMC3"    )

data<-data[order(data$FoldEnrichment,decreasing=TRUE),]
data$RANK<-1:15
library(ggrepel)

pdf("PF-Foldenrichment-p.pdf")
ggplot(data,aes(x=RANK,y=FoldEnrichment))+
  geom_point(aes(size=mlog10Pvalue,color="b"),alpha=0.6)+
  #scale_size(range=c(1,12))+
  theme_bw()+
  theme( )+
  geom_text_repel(
    data = data,
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE )+theme(panel.grid=element_blank())
dev.off();


"#48453d","#20854EFF" ,"#0072B5FF" ,"#E18727FF", 
   "#7876B1FF" ,"#fbaba9","#8edcaf", "#66E1E6","#0F3057", "#F06C5D", "#b92419","#EB31FA",
   "#525c19","#059E8D","#FA6200","#9F7717","#f38370","#e2ad15",
   "#afd961","#7e6db2" ,"#6E9E6D","#6D749E"

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175')
pdf("celltype-cluster-Plot-UMAP2.pdf")
p1 <- plotEmbedding(ArchRProj =  EC4proj,size = 0.5, labelSize=2,legendSize = 0,colorBy = "cellColData", name = "Clusters_3", embedding = "UMAP",
    pal=my36colors[1:13])
p2 <- plotEmbedding(ArchRProj =  FB3proj,size = 0.5, labelSize=2,legendSize = 0,colorBy = "cellColData", name = "Clusters3", embedding = "UMAP",
    pal=my36colors[14:24])
p3 <- plotEmbedding(ArchRProj = SMC1proj,size = 0.5,labelSize=2,legendSize = 0, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP",
    pal=my36colors[25:36])
p1 
p2 
p3 
ggAlignPlots(p1, p2,p3, type = "h")

dev.off()



SMC<-SMCmarker[SMCmarker$group_name=="SMC1",]$name
EC<-ECmarker[ECmarker$group_name=="EC4",]$name
MP<-MPmarker[MPmarker$group_name=="MP3",]$name
FB<-FBmarker[FBmarker$group_name=="FB3",]$name
vennplot<-venn.diagram(
  x = list(SMC,FB,EC,MP,protein_folding),
  category.names = c("SMC" , "FB" , "EC","MP","protein_folding"),
  filename =NULL,
  fill = c('#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3'),
  alpha = 0.50,
  output=TRUE
)
pdf("markergene_PFgene_vennplot.pdf")
grid.draw(vennplot)
dev.off()
overlap<-intersect(FB,SMC)
overlap<-intersect(overlap,MP)
overlap<-intersect(overlap,EC)
overlap<-intersect(overlap,protein_folding)


FBPF<-intersect(FB,protein_folding)
MPPF<-intersect(MP,protein_folding)
ECPF<-intersect(EC,protein_folding)
SMCPF<-intersect(SMC,protein_folding)

vennplot<-venn.diagram(
  x = list(SMCPF,FBPF,ECPF,MPPF),
  category.names = c("SMC" , "FB" , "EC","MP"),
  filename =NULL,
  fill = c( '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3'),
  alpha = 0.50,
  output=TRUE
)
pdf("markergene_PFgene_overlap_vennplot.pdf")
grid.draw(vennplot)
dev.off()

vennplot<-venn.diagram(
  x = list(SMC,FB,EC,MP),
  category.names = c("SMC" , "FB" , "EC","MP"),
  filename =NULL,
  fill = c('#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3'),
  alpha = 0.50,
  output=TRUE
)
pdf("markergene_overlap_vennplot.pdf")
grid.draw(vennplot)
dev.off()

overlap<-intersect(FB,SMC)
overlap<-intersect(overlap,MP)
overlap<-intersect(overlap,EC)

