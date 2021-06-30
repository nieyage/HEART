
library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 
addArchRGenome("mm10")

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
