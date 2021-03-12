#######3 celltype
#######de novo' protein folding########
FB3<-"Hspa8/Hspa1l/Hspa1a/Hspa1b/Hsph1/H2-DMa/Hspd1/Hspe1/Dnajb1/Bag1"
EC4<-"Hspa1a/Hspa8/Dnajb14/Hspa1l"
MP3<-"Hspa1b/Dnajb1/Hspa8"
FB3<-strsplit(FB3,'/')
EC4<-strsplit(EC4,'/')
MP3<-strsplit(MP3,'/')
FB3<-FB3[[1]]
EC4<-EC4[[1]]
MP3<-MP3[[1]]
alltermgene<-union(FB3,EC4)
alltermgene<-union(alltermgene,MP3)
#####构建upsetR  matrix#########
meta<-matrix(data=NA,nrow=length(alltermgene),ncol=3);
colnames(meta)<-c("FB3","EC4","MP3");
rownames(meta)<-alltermgene;

for ( i in 1:length(alltermgene))
{ 
  if (alltermgene[i] %in% FB3)
    meta[i,1] = 1
    else 
    meta[i,1] = 0
  if (alltermgene[i] %in% EC4)
    meta[i,2] = 1
    else 
    meta[i,2] = 0
  if (alltermgene[i] %in% MP3)
    meta[i,3] = 1
    else 
    meta[i,3] = 0
}
######plot heatmap####
pdf("de_novo'_protein_folding.pdf")
pheatmap(meta,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(colors = c("white","red"))(100),
         cellwidth = 12, cellheight = 12,
         show_rownames=T,show_colnames=T)

#meta = data.frame(meta)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
library(UpSetR)
m<-as.data.frame(meta)
c<-cbind(alltermgene,m)
upset(c, mb.ratio = c(0.6, 0.4), order.by = "freq", 
      nsets = 7, number.angles = 0, point.size = 6, line.size = 1, mainbar.y.label = "Number of genes",
      sets.x.label = "Total number of genes in each cell type", text.scale = c(2, 2, 0.8, 2, 2, 3))


#######2 celltype
#######chaperone-mediated protein folding########
#FB3<-"Hspa8/Hspa1l/Hspa1a/Hspa1b/Hsph1/Hspb6/H2-DMa/Hspe1/Dnajb1/Bag1"
MP3<-"Hspa1b/Dnajb1/Hspa8"
#FB3<-strsplit(FB3,'/')
MP3<-strsplit(MP3,'/')
#FB3<-FB3[[1]]
MP3<-MP3[[1]]
EC4<-"Hspa1a/Hspa8/Dnajb14/Hspa1l"
EC4<-strsplit(EC4,'/')
EC4<-EC4[[1]]

alltermgene<-union(MP3,EC4)
#####构建upsetR  matrix#########
meta<-matrix(data=NA,nrow=length(alltermgene),ncol=2);
colnames(meta)<-c("MP3","EC4");
rownames(meta)<-alltermgene;

for ( i in 1:length(alltermgene))
{ 
  if (alltermgene[i] %in% MP3)
    meta[i,1] = 1
    else 
    meta[i,1] = 0
  if (alltermgene[i] %in% EC4)
    meta[i,2] = 1
    else 
    meta[i,2] = 0
}
######plot heatmap####
pdf("chaperone_cofactor-dependent_protein_refolding.pdf")
pheatmap(meta,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(colors = c("white","red"))(100),
         cellwidth = 12, cellheight = 12,
         show_rownames=T,show_colnames=T)

#meta = data.frame(meta)
require(ggplot2); require(plyr); require(gridExtra); require(grid);
library(UpSetR)
m<-as.data.frame(meta)
c<-cbind(alltermgene,m)
upset(c, mb.ratio = c(0.6, 0.4), order.by = "freq", 
      nsets = 7, number.angles = 0, point.size = 6, line.size = 1, mainbar.y.label = "Number of genes",
      sets.x.label = "Total number of genes in each cell type", text.scale = c(2, 2, 0.8, 2, 2, 3))

