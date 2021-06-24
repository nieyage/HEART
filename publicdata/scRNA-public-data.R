########scRNA public data#########
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
hisall<-readRDS("addpubanno-scRNA.rds")
all<-readRDS("scRNAseq.rds")
DefaultAssay(all) <- "RNA" 

umap<-read.table("heart.rna.umap.coord.txt")
anno<-read.table("heart.rna.anno.txt")
#ourcell<-rownames(all@meta.data)
#ourcell<-gsub("_","",ourcell)
#ourcell<-gsub("-1","",ourcell)

hiscell<-rownames(anno)
hiscell<-sub("_","",hiscell)
hiscell<-paste(hiscell,seq="-1")
hiscell<-sub(" ","",hiscell)
####ourcell:18862 hiscell:17320####

hisall<-subset(all,cells=hiscell)
#order<-colnames(hisall)
#rownames(anno)<-hiscell

hisall@meta.data$celltype2<-anno$FineID
#######提取FB##########
hisall@active.ident<-hisall$celltype2
FB<-subset(hisall,idents=c("FB1","FB2","FB3","FB4","Pro.FB"))

marker_genes<-c("Brinp1","Pim1","Birc6","Il1rapl1","Cdh18","Grid2","Mdm4")######FB marker genes

"Cebpd"####features 中不存在（其他几个成员均存在）
CEBP<-c("Cebpa","Cebpb","Cebpg","Cebpe","Cebpz")
library(reshape2)
data<-as.data.frame(FB[["RNA"]]@data)
vln.df=as.data.frame(data[CEBP,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=data.frame(CB = colnames(FB),celltype=FB$orig.ident)
#anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = CEBP) #为了控制画图的基因顺序

pdf("Vlnplot-FBmarkergene-sample.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()

meta<-matrix(data=NA,nrow=5639,ncol=3);
colnames(meta)=c("group","cellbarcode","subtypename");
meta[,1] = FB$orig.ident
meta[,2] = rownames(FB@meta.data)
#cl <- meta[,1]
for ( i in 1:nrow(meta))
{
  if (meta[i,1] =="P11MI"|  meta[i,1] =="P13MI")
    meta[i,3] = "P1MI"
  if (meta[i,1] =="P11Sham"|  meta[i,1] =="P13Sham")
    meta[i,3] = "P1Sham"
  if (meta[i,1] =="P81MI"|  meta[i,1] =="P83MI")
    meta[i,3] = "P8MI"
  if (meta[i,1] =="P81Sham"|  meta[i,1] =="P83Sham")
    meta[i,3] = "P8Sham"
}
meta = data.frame(meta)
FB$subtypename  = meta$subtypename
table(FB$subtypename)


anno=data.frame(CB = colnames(FB),celltype=FB$subtypename)
#anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = CEBP) #为了控制画图的基因顺序

pdf("Vlnplot-CEBPgene-mergeday13.pdf")
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


marker_genes<-c("Brinp1","Pim1","Birc6","Il1rapl1","Cdh18","Grid2","Mdm4")######FB marker genes
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
data<-as.data.frame(FB[["RNA"]]@data)
vln.df=as.data.frame(data[marker_genes,])
vln.df<-na.omit(vln.df)
vln<-vln.df
pdf("Vlnplot-FBtopgene-mergeday13.pdf")

for(i in 1:8){
  j=i*6
  k=j-5
  vln.df=vln[k:j,]
  vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
anno=data.frame(CB = colnames(FB),celltype=FB$subtypename)
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


########EC
EC<-subset(hisall,idents=c("Endo","VEC1","VEC2","VEC3","Pro.EC"))

nega_marker_genes<-c("Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Tcf4","Col4a2")######EC3 marker genes
pro_marker_genes<-c("Plcg1","Acvrl","Flt1","Itgb3","Itgb1bp1")######EC3 marker genes

library(reshape2)
data<-as.data.frame(EC[["RNA"]]@data)
vln.df=as.data.frame(data[nega_marker_genes,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=data.frame(CB = colnames(EC),celltype=EC$orig.ident)
#anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = nega_marker_genes) #为了控制画图的基因顺序

pdf("Vlnplot-EC-nega-markergene-sample.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()


vln.df=as.data.frame(data[pro_marker_genes,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=data.frame(CB = colnames(EC),celltype=EC$orig.ident)
#anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = pro_marker_genes) #为了控制画图的基因顺序

pdf("Vlnplot-EC-pro-markergene-sample.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()


meta<-matrix(data=NA,nrow=5337,ncol=3);
colnames(meta)=c("group","cellbarcode","subtypename");
meta[,1] = EC$orig.ident
meta[,2] = rownames(EC@meta.data)

for ( i in 1:nrow(meta))
{
  if (meta[i,1] =="P11MI"|  meta[i,1] =="P13MI")
    meta[i,3] = "P1MI"
  if (meta[i,1] =="P11Sham"|  meta[i,1] =="P13Sham")
    meta[i,3] = "P1Sham"
  if (meta[i,1] =="P81MI"|  meta[i,1] =="P83MI")
    meta[i,3] = "P8MI"
  if (meta[i,1] =="P81Sham"|  meta[i,1] =="P83Sham")
    meta[i,3] = "P8Sham"
}
meta = data.frame(meta)
EC$subtypename  = meta$subtypename
table(EC$subtypename)



data<-as.data.frame(EC[["RNA"]]@data)
vln.df=as.data.frame(data[nega_marker_genes,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=data.frame(CB = colnames(EC),celltype=EC$subtypename)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = nega_marker_genes) #为了控制画图的基因顺序

pdf("Vlnplot-EC-nega-markergene-mergeday13.pdf")
vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+geom_boxplot(width=0.1,outlier.size=0,
    notch=TRUE,position=position_dodge(width=0.9))+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 0,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()


vln.df=as.data.frame(data[pro_marker_genes,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
#anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = pro_marker_genes) #为了控制画图的基因顺序

pdf("Vlnplot-EC-pro-markergene-mergeday13.pdf")
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



# library
library(ggridges)
library(ggplot2)
# Diamonds dataset is provided by R natively
#head(diamonds)
 PFgene<-c("Trap1","Cct2","Dnajb5","Dnajb11","Ppia","Ppil1","Stub1","Cdc37","Mkks","Grpel1","Fkbp8","Calr","P3h1","Ptges3l","Dnajb1","Hsp90ab1","Hsp90aa1","Hspd1","Pdcl3","Nudc","Dnaja1","Pfdn6","Mesd","Hspa1b","Hspe1","Hspa8","Hspa1a")
 data<-as.data.frame(hisall[["RNA"]]@data)

vln.df=as.data.frame(data[PFgene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
##anno=data.frame(CB=rownames(anno),celltype=anno$celltype)
#vln.df=inner_join(vln.df,anno,by="CB")
#vln.df$gene=factor(vln.df$gene,levels = PF_genes) #为了控制画图的基因顺序

# basic example
library(RColorBrewer)
colourCount = length(unique(PFgene))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
 


pdf("Ridgeplot-EC-pro-markergene-mergeday13.pdf")
ggplot(vln.df, aes(x = exp, y = gene, fill = gene)) +
  geom_density_ridges() +
  theme_ridges() +  
  theme(legend.position = "none")+
  scale_fill_manual(values = getPalette(colourCount))
 dev.off();

#######是否是同一群细胞表达###########

pdf("Allcell_PFgene_umap.pdf")
FeaturePlot(hisall, features = PFgene[1:9])
FeaturePlot(hisall, features = PFgene[10:18])
FeaturePlot(hisall, features = PFgene[19:27])
dev.off();

#####用heatmap观察 是否是同一群细胞表达PFgene######

PFgene<-c("Trap1","Cct2","Dnajb5","Dnajb11","Ppia","Ppil1","Stub1","Cdc37","Mkks","Grpel1","Fkbp8","Calr","P3h1","Ptges3l","Dnajb1","Hsp90ab1","Hsp90aa1","Hspd1","Pdcl3","Nudc","Dnaja1","Pfdn6","Mesd","Hspa1b","Hspe1","Hspa8","Hspa1a")
data<-as.data.frame(hisall[["RNA"]]@data)
vln.df=as.data.frame(data[PFgene,])
count<-t(vln.df)
count<-count[,-23]
for(i in 1:26){
	for(j in 1:17320){
		if (count[j,i]== 0){count[j,i]=0}
		else{count[j,i]=1}
	}
}

library(pheatmap)
count=t(scale(t(vln.df),scale = T,center = T))

pdf("Allcell_PFgene_heatmap-zscore.pdf")
count<-na.omit(count)
#breaks

pdf("Allcell_PFgene_heatmap-zscore-bk3.pdf")
vln.df<-na.omit(vln.df)
count=t(scale(t(vln.df),scale = T,center = T))
bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))

pheatmap(count,
         scale = "none",
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         #legend_breaks=seq(-8,8,2),
         show_rownames=T,show_colnames=F,
         breaks=bk)



pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy","white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=F)

dev.off()

saveRDS(hisall,"addpubanno-scRNA.rds")


#######提取FB中的P1MI 是否表达FB2基因的细胞是同一群#######

meta<-matrix(data=NA,nrow=5639,ncol=3);
colnames(meta)=c("group","cellbarcode","subtypename");
meta[,1] = FB$orig.ident
meta[,2] = rownames(FB@meta.data)
#cl <- meta[,1]
for ( i in 1:nrow(meta))
{
  if (meta[i,1] =="P11MI"|  meta[i,1] =="P13MI")
    meta[i,3] = "P1MI"
  if (meta[i,1] =="P11Sham"|  meta[i,1] =="P13Sham")
    meta[i,3] = "P1Sham"
  if (meta[i,1] =="P81MI"|  meta[i,1] =="P83MI")
    meta[i,3] = "P8MI"
  if (meta[i,1] =="P81Sham"|  meta[i,1] =="P83Sham")
    meta[i,3] = "P8Sham"
}
meta = data.frame(meta)
FB$subtypename  = meta$subtypename
table(FB$subtypename)
FB@active.ident<-FB$subtypename
P1MIFB<-subset(FB,idents=c("P1MI"))
marker_genes<-c("Brinp1","Pim1","Birc6","Il1rapl1","Cdh18","Grid2","Mdm4")######FB marker genes

CEBP<-c("Cebpa","Cebpb","Cebpg","Cebpe","Cebpz")
library(reshape2)
data<-as.data.frame(FB[["RNA"]]@data)
vln.df=as.data.frame(data[marker_genes,])
vln.df<-na.omit(vln.df)
count=t(scale(t(vln.df),scale = T,center = T))
count<-na.omit(count)

pdf("P1MIFB_gene_heatmap-zscore-bk3.pdf")
bk <- c(seq(-3,-0.1,by=0.001),seq(0,3,by=0.001))
library(pheatmap)
pheatmap(count,
         scale = "none",
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         #legend_breaks=seq(-8,8,2),
         show_rownames=T,show_colnames=F,
         breaks=bk)
dev.off()





pdf("P1MIFB_topgene_heatmap-zscore-bk5.pdf")
bk <- c(seq(-5,-0.1,by=0.001),seq(0,5,by=0.001))
library(pheatmap)
pheatmap(count,
         scale = "none",       
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         #legend_breaks=seq(-8,8,2),
         show_rownames=T,show_colnames=F,
         breaks=bk)
dev.off()


###########分不同sample 不同cell type进行细胞的PF 的heatmap绘制

PFgene<-c("Trap1","Cct2","Dnajb5","Dnajb11","Ppia","Ppil1","Stub1","Cdc37","Mkks","Grpel1","Fkbp8","Calr","P3h1","Ptges3l","Dnajb1","Hsp90ab1","Hsp90aa1","Hspd1","Pdcl3","Nudc","Dnaja1","Pfdn6","Hspa1b","Hspe1","Hspa8","Hspa1a")
hisall <- ScaleData(hisall, features =rownames(hisall))

pdf("PFgene_groupby_celltype_heatmap.pdf")
DoHeatmap(hisall, features = PFgene,group.by = "celltype2",combine = F) 
dev.off()

meta<-matrix(data=NA,nrow=17320,ncol=3);
colnames(meta)=c("group","cellbarcode","subtypename");
meta[,1] = hisall$orig.ident
meta[,2] = rownames(hisall@meta.data)
for ( i in 1:nrow(meta))
{
  if (meta[i,1] =="P11MI"|  meta[i,1] =="P13MI")
    meta[i,3] = "P1MI"
  if (meta[i,1] =="P11Sham"|  meta[i,1] =="P13Sham")
    meta[i,3] = "P1Sham"
  if (meta[i,1] =="P81MI"|  meta[i,1] =="P83MI")
    meta[i,3] = "P8MI"
  if (meta[i,1] =="P81Sham"|  meta[i,1] =="P83Sham")
    meta[i,3] = "P8Sham"
}
meta = data.frame(meta)
hisall$sample  = meta$subtypename
pdf("PFgene_groupby_sample_heatmap.pdf")
DoHeatmap(hisall, features = PFgene,group.by = "subtypename",combine = F) 
dev.off()

#####merge subtype#######
    Art.EC    B cells         CM    DC-like       Endo        EPI        FB1 
      1114        107        278         64        276        420       2280 
       FB2        FB3        FB4      glial        Gra Macrophage   Monocyte 
      1336        706        333         31        412       1299        354 
  Pericyte     Pro.EC     Pro.FB        SMC    T cells       VEC1       VEC2 
      1036       1049        984       1112        117       2615       1229 
      VEC3 
       168 


meta<-matrix(data=NA,nrow=17320,ncol=3);
colnames(meta)=c("group","cellbarcode","subtypename");
meta[,1] = as.character(hisall$celltype2)
meta[,2] = rownames(hisall@meta.data)
meta[,3] = "other"
for ( i in 1:nrow(meta))
{
  if (meta[i,1] %in% c("Art.EC","Endo","VEC1","VEC2","VEC3","Pro.EC"))
    meta[i,3] = "EC"
  if (meta[i,1] %in% c("FB1","FB2","FB3","FB4","Pro.FB"))
    meta[i,3] = "FB"
  if (meta[i,1] %in% c("B cells","DC-like","glial","Gra","Macrophage","Monocyte","T cells"))
    meta[i,3] = "IC"
  if (meta[i,1] %in% c("CM"))
    meta[i,3] = "CM"
  if (meta[i,1] %in% c("EPI"))
    meta[i,3] = "EPI"
  if (meta[i,1] %in% c("Pericyte"))
    meta[i,3] = "Pericyte"
  if (meta[i,1] %in% c("SMC"))
    meta[i,3] = "SMC"
}
meta = data.frame(meta)
hisall$subtypename  = meta$subtypename
table(hisall$subtypename)
hisall@active.ident<-hisall$subtypename
hisall <- ScaleData(hisall, features =rownames(hisall))

pdf("PFgene_groupby_mergecelltype_heatmap.pdf")
DoHeatmap(hisall, features = PFgene,group.by = "subtypename",combine = F) 
dev.off()





#########split by sample and merged celltype

meta<-matrix(data=NA,nrow=17320,ncol=3);
colnames(meta)=c("group","cellbarcode","subtypename");
meta[,1] = as.character(hisall$subtypename)
meta[,2] = as.character(hisall$sample)
meta[,3] = "other"
for ( i in 1:nrow(meta))
{
  if (meta[i,1]=="EC" & meta[i,2]=="P1MI")
    meta[i,3] = "EC_P1MI"
  if (meta[i,1]=="EC" & meta[i,2]=="P8MI")
    meta[i,3] = "EC_P8MI"
  if (meta[i,1]=="EC" & meta[i,2]=="P1Sham")
    meta[i,3] = "EC_P1Sham"
  if (meta[i,1]=="EC" & meta[i,2]=="P8Sham")
    meta[i,3] = "EC_P8Sham"
  if (meta[i,1]=="FB" & meta[i,2]=="P1MI")
    meta[i,3] = "FB_P1MI"
  if (meta[i,1]=="FB" & meta[i,2]=="P8MI")
    meta[i,3] = "FB_P8MI"
  if (meta[i,1]=="FB" & meta[i,2]=="P1Sham")
    meta[i,3] = "FB_P1Sham"
  if (meta[i,1]=="FB" & meta[i,2]=="P8Sham")
    meta[i,3] = "FB_P8Sham"
  if (meta[i,1]=="SMC" & meta[i,2]=="P1MI")
    meta[i,3] = "SMC_P1MI"
  if (meta[i,1]=="SMC" & meta[i,2]=="P8MI")
    meta[i,3] = "SMC_P8MI"
  if (meta[i,1]=="SMC" & meta[i,2]=="P1Sham")
    meta[i,3] = "SMC_P1Sham"
  if (meta[i,1]=="SMC" & meta[i,2]=="P8Sham")
    meta[i,3] = "SMC_P8Sham"
  if (meta[i,1]=="IC" & meta[i,2]=="P1MI")
    meta[i,3] = "IC_P1MI"
  if (meta[i,1]=="IC" & meta[i,2]=="P8MI")
    meta[i,3] = "IC_P8MI"
  if (meta[i,1]=="IC" & meta[i,2]=="P1Sham")
    meta[i,3] = "IC_P1Sham"
  if (meta[i,1]=="IC" & meta[i,2]=="P8Sham")
    meta[i,3] = "IC_P8Sham"
  if (meta[i,1]=="EPI" & meta[i,2]=="P1MI")
    meta[i,3] = "EPI_P1MI"
  if (meta[i,1]=="EPI" & meta[i,2]=="P8MI")
    meta[i,3] = "EPI_P8MI"
  if (meta[i,1]=="EPI" & meta[i,2]=="P1Sham")
    meta[i,3] = "EPI_P1Sham"
  if (meta[i,1]=="EPI" & meta[i,2]=="P8Sham")
    meta[i,3] = "EPI_P8Sham"
  if (meta[i,1]=="CM" & meta[i,2]=="P1MI")
    meta[i,3] = "CM_P1MI"
  if (meta[i,1]=="CM" & meta[i,2]=="P8MI")
    meta[i,3] = "CM_P8MI"
  if (meta[i,1]=="CM" & meta[i,2]=="P1Sham")
    meta[i,3] = "CM_P1Sham"
  if (meta[i,1]=="CM" & meta[i,2]=="P8Sham")
    meta[i,3] = "CM_P8Sham"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P1MI")
    meta[i,3] = "Pericyte_P1MI"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P8MI")
    meta[i,3] = "Pericyte_P8MI"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P1Sham")
    meta[i,3] = "Pericyte_P1Sham"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P8Sham")
    meta[i,3] = "Pericyte_P8Sham"

}


meta = data.frame(meta)
hisall$subtypename  = meta$subtypename
table(hisall$subtypename)
hisall@active.ident<-hisall$subtypename
hisall <- ScaleData(hisall, features =rownames(hisall))

pdf("PFgene_groupby_sample_mergecelltype_heatmap.pdf")
DoHeatmap(hisall, features = PFgene,group.by = "subtypename",combine = F) 
dev.off()

#####分开cell type######
saveRDS(hisall,"hisall.rds")
CM<-subset(hisall, idents = c( "CM"))
EC<-subset(hisall, idents = c( "EC"))
SMC<-subset(hisall, idents = c( "SMC"))
IC<-subset(hisall, idents = c( "IC"))
Pericyte<-subset(hisall, idents = c( "Pericyte"))
EPI<-subset(hisall, idents = c( "EPI"))
FB<-subset(hisall, idents = c( "FB"))

cols=brewer.pal(8,"Set2")[1:4]
pdf("PFgene_splitcelltype_heatmap.pdf")
DoHeatmap(CM,       features = PFgene,group.by = "sample",group.colors =cols,disp.min=-2,disp.max=2)+scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(EC,       features = PFgene,group.by = "sample",group.colors =cols,disp.min=-2,disp.max=2)+scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(IC,       features = PFgene,group.by = "sample",group.colors =cols,disp.min=-2,disp.max=2)+scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(FB,       features = PFgene,group.by = "sample",group.colors =cols,disp.min=-2,disp.max=2)+scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(SMC,      features = PFgene,group.by = "sample",group.colors =cols,disp.min=-2,disp.max=2)+scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(EPI,      features = PFgene,group.by = "sample",group.colors =cols,disp.min=-2,disp.max=2)+scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
DoHeatmap(Pericyte, features = PFgene,group.by = "sample",group.colors =cols,disp.min=-2,disp.max=2)+scale_fill_gradientn(colors = c("#4858A7","#788FC8","#D6DAE1","#F49B7C","#B51F29"))
dev.off()
