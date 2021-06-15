#########ECDF 分布 in FB each sample#########
#######提取FB##########
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
hisall<-readRDS("addpubanno-scRNA.rds")
FB<-subset(hisall,idents=c("FB1","FB2","FB3","FB4","Pro.FB"))
marker_genes<-c("Brinp1","Pim1","Birc6","Il1rapl1","Cdh18","Grid2","Mdm4",
  "Frmpd4"       ,"Hspa5"        ,"Aff2"         ,"Cxcl2"        ,
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
library(reshape2)
data<-as.data.frame(FB[["RNA"]]@data)
vln.df=as.data.frame(data[marker_genes,])
vln.df=na.omit(vln.df)

#vln.df$gene=rownames(vln.df)
allcell<-colnames(vln.df)
P1MI<-allcell[grep("^P11MI.*|^P13MI.*",allcell)]
P8MI<-allcell[grep("^P81MI.*|^P83MI.*",allcell)]
P1Sham<-allcell[grep("^P11Sham.*|^P13Sham.*",allcell)]
P8Sham<-allcell[grep("^P81Sham.*|^P83Sham.*",allcell)]
#allgene<-rownames(data)
marker_genes<-rownames(vln.df)
pdf("FB2topgene-ECDF.pdf",width=5,height=3)
for( i in 1:55){
  gene=marker_genes[i];
  data<-vln.df[which(rownames(vln.df)%in%gene),]
  P1MIcount<-data[which(names(data)%in%P1MI)]
  P1MIcount<-as.numeric(P1MIcount)
  P8MIcount<-data[which(names(data)%in%P8MI)]
  P8MIcount<-as.numeric(P8MIcount)
  P1Shamcount<-data[which(names(data)%in%P1Sham)]
  P1Shamcount<-as.numeric(P1Shamcount)
  P8Shamcount<-data[which(names(data)%in%P8Sham)]
  P8Shamcount<-as.numeric(P8Shamcount)
  mydf = data.frame(
     value = c(P1MIcount,P8MIcount,P1Shamcount,P8Shamcount),
     sample =c(rep("P1MIcount",length(P1MIcount)),rep("P8MIcount",length(P8MIcount)),
      rep("P1Shamcount",length(P1Shamcount)),rep("P8Shamcount",length(P8Shamcount)))
  )
  # print(mydf)
  p<-ggplot(mydf,aes(x = value)) + stat_ecdf(aes(colour = sample))+
    scale_color_manual(values =brewer.pal(8,"Set2")[1:4])+
    #scale_x_discrete("")+scale_y_continuous("")+
    theme_bw()+ggtitle(gene) +
    theme(
     panel.grid.major = element_blank(),panel.grid.minor = element_blank()
    )
    print(p);
}


data<-vln.df[which(rownames(vln.df)=="Pim1"),]
P1MIcount<-data[which(names(data)%in%P1MI)]
P1MIcount<-as.numeric(P1MIcount)
P8MIcount<-data[which(names(data)%in%P8MI)]
P8MIcount<-as.numeric(P8MIcount)
P1Shamcount<-data[which(names(data)%in%P1Sham)]
P1Shamcount<-as.numeric(P1Shamcount)
P8Shamcount<-data[which(names(data)%in%P8Sham)]
P8Shamcount<-as.numeric(P8Shamcount)
mydf = data.frame(
   value = c(P1MIcount,P8MIcount,P1Shamcount,P8Shamcount),
   sample =c(rep("P1MIcount",length(P1MIcount)),rep("P8MIcount",length(P8MIcount)),
    rep("P1Shamcount",length(P1Shamcount)),rep("P8Shamcount",length(P8Shamcount)))
)

pdf("FB2topgene-ECDF.pdf",width=5,height=3)
ggplot(mydf,aes(x = value)) + stat_ecdf(aes(colour = sample))+
    scale_color_manual(values =brewer.pal(8,"Set2")[1:4])+ggtitle(gene)+
    #scale_x_discrete("")+scale_y_continuous("")+
    theme_bw() +
    theme(panel.grid=element_blank())
dev.off()



#####EC marker gene's ECDF distribution#######$#

EC<-subset(hisall,idents=c("Endo","VEC1","VEC2","VEC3","Pro.EC"))

nega_marker_genes<-c("Adamts1","Xdh","Foxj2","Zeb1","Meox2","Cd36","Tcf4","Col4a2")######EC3 marker genes
pro_marker_genes<-c("Plcg1","Acvrl","Flt1","Itgb3","Itgb1bp1")######EC3 marker genes

data<-as.data.frame(EC[["RNA"]]@data)
vln.df=as.data.frame(data[pro_marker_genes,])
vln.df=na.omit(vln.df)
allcell<-colnames(vln.df)
P1MI<-allcell[grep("^P11MI.*|^P13MI.*",allcell)]
P8MI<-allcell[grep("^P81MI.*|^P83MI.*",allcell)]
P1Sham<-allcell[grep("^P11Sham.*|^P13Sham.*",allcell)]
P8Sham<-allcell[grep("^P81Sham.*|^P83Sham.*",allcell)]
pro_marker_genes<-rownames(vln.df)
library(RColorBrewer)
pdf("EC-pro_marker_genes-ECDF.pdf",width=5,height=3)
for( i in 1:5){
  gene=pro_marker_genes[i];
  data<-vln.df[which(rownames(vln.df)%in%gene),]
  P1MIcount<-data[which(names(data)%in%P1MI)]
  P1MIcount<-as.numeric(P1MIcount)
  P8MIcount<-data[which(names(data)%in%P8MI)]
  P8MIcount<-as.numeric(P8MIcount)
  P1Shamcount<-data[which(names(data)%in%P1Sham)]
  P1Shamcount<-as.numeric(P1Shamcount)
  P8Shamcount<-data[which(names(data)%in%P8Sham)]
  P8Shamcount<-as.numeric(P8Shamcount)
  mydf = data.frame(
     value = c(P1MIcount,P8MIcount,P1Shamcount,P8Shamcount),
     sample =c(rep("P1MIcount",length(P1MIcount)),rep("P8MIcount",length(P8MIcount)),
      rep("P1Shamcount",length(P1Shamcount)),rep("P8Shamcount",length(P8Shamcount)))
  )
  # print(mydf)
  p<-ggplot(mydf,aes(x = value)) + stat_ecdf(aes(colour = sample))+
    scale_color_manual(values =brewer.pal(8,"Set2")[1:4])+
    #scale_x_discrete("")+scale_y_continuous("")+
    theme_bw()+ggtitle(gene) +
    theme(
     panel.grid.major = element_blank(),panel.grid.minor = element_blank()
    )
    print(p);
}
