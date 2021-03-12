#####find special ECM gene for subtype ####
####use archR to find overlap ####

library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.5 & Log2FC >= 0.3")
mart<-read.csv("mart_export.csv")
ECMgene<-mart[which(mart$GO.term.name=="extracellular matrix"),]
ECMgene<-ECMgene$Gene.name

FB1<-markerList$FB1$name
FB1ECM<-intersect(FB1,ECMgene)
FB1ECM

FB2<-markerList$FB2$name
FB2ECM<-intersect(FB2,ECMgene)
FB2ECM

FB3<-markerList$FB3$name
FB3ECM<-intersect(FB3,ECMgene)
FB3ECM

FB4<-markerList$FB4$name
FB4ECM<-intersect(FB4,ECMgene)
FB4ECM
g<-c(FB1ECM,FB2ECM,FB3ECM,FB4ECM)
count1<-subcount[which(rownames(subcount)%in%g),]
count1=t(scale(t(count1),scale = T,center = T))
pdf("ECM-overlap-ArchR-1.pdf",width=6,height=12)
list<-pheatmap(c, cluster_rows=F, 
               clustering_method="ward.D2",
               show_rownames=T,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)
dev.off()











pdf("FB234ECM-overlap-ArchR.pdf",width=3,height=4)
count1<-subcount[which(rownames(subcount)%in%FB4ECM),]
count1=t(scale(t(count1),scale = T,center = T))
pheatmap(count1, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               show_rownames=T,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)
dev.off()

####根据proline含量绘制ECDF图#####
FB1<-c(4.18,5.80,9.97 ,11.58 ,5.33 ,2.85 ,5.41 ,12.08 ,7.64 ,19.17 ,14.90 ,3.79 ,4.75 ,5.67)
FB2<-c(19.47 ,5.98 ,5.60 ,4.46 ,16.71 ,7.36 ,18.29 ,4.39 ,6.82 ,10.63 ,25.65 ,5.23 ,11.47 ,7.55 ,7.65 ,4.38 ,5.85 ,8.27 ,4.87)
FB3<-c(10.28 ,11.88 ,19.00 ,8.13 ,9.69 ,7.83 ,18.17 ,10.99 ,21.64 ,8.93 ,10.23 ,13.75 ,7.46 ,5.37 ,6.52 ,7.12 ,7.01 ,20.75 ,8.58 ,7.51 )
FB4<-c(11.14 ,7.96 ,4.41 ,6.94 ,16.62 ,13.76 ,6.10 ,7.40 ,7.99 ,8.43 ,6.84 ,2.62 ,5.79 ,6.67 ,8.00 ,22.45 ,8.30 ,5.59 ,5.09 ,7.41 ,9.92 ,4.77 ,15.40 ,7.69 ,6.58)


###提取FB1234中的AR，NAR，P组的barcode信息，根据细胞barcodes从snapATAC中提取cebp家族的表达量信息做violin图########
sample<-nproj$Sample
sub<-nproj$subtype
for (i in 1:21621){
if(sample[i]=="AR3"&&sub[i]=="FB1"){sample[i]="FB1-AR"}
if(sample[i]=="AR3"&&sub[i]=="FB2"){sample[i]="FB2-AR"}
if(sample[i]=="AR3"&&sub[i]=="FB3"){sample[i]="FB3-AR"}
if(sample[i]=="AR3"&&sub[i]=="FB4"){sample[i]="FB4-AR"}
if(sample[i]=="NAR3"&&sub[i]=="FB1"){sample[i]="FB1-NAR"}
if(sample[i]=="NAR3"&&sub[i]=="FB2"){sample[i]="FB2-NAR"}
if(sample[i]=="NAR3"&&sub[i]=="FB3"){sample[i]="FB3-NAR"}
if(sample[i]=="NAR3"&&sub[i]=="FB4"){sample[i]="FB4-NAR"}
if(sample[i]=="P3"&&sub[i]=="FB1"){sample[i]="FB1-P"}
if(sample[i]=="P3"&&sub[i]=="FB2"){sample[i]="FB2-P"}
if(sample[i]=="P3"&&sub[i]=="FB3"){sample[i]="FB3-P"}
if(sample[i]=="P3"&&sub[i]=="FB4"){sample[i]="FB4-P"}
if(sample[i]=="AR7"&&sub[i]=="FB1"){sample[i]="FB1-AR"}
if(sample[i]=="AR7"&&sub[i]=="FB2"){sample[i]="FB2-AR"}
if(sample[i]=="AR7"&&sub[i]=="FB3"){sample[i]="FB3-AR"}
if(sample[i]=="AR7"&&sub[i]=="FB4"){sample[i]="FB4-AR"}
if(sample[i]=="NAR7"&&sub[i]=="FB1"){sample[i]="FB1-NAR"}
if(sample[i]=="NAR7"&&sub[i]=="FB2"){sample[i]="FB2-NAR"}
if(sample[i]=="NAR7"&&sub[i]=="FB3"){sample[i]="FB3-NAR"}
if(sample[i]=="NAR7"&&sub[i]=="FB4"){sample[i]="FB4-NAR"}
if(sample[i]=="P7"&&sub[i]=="FB1"){sample[i]="FB1-P"}
if(sample[i]=="P7"&&sub[i]=="FB2"){sample[i]="FB2-P"}
if(sample[i]=="P7"&&sub[i]=="FB3"){sample[i]="FB3-P"}
if(sample[i]=="P7"&&sub[i]=="FB4"){sample[i]="FB4-P"}
if(sample[i]=="AR14"&&sub[i]=="FB1"){sample[i]="FB1-AR"}
if(sample[i]=="AR14"&&sub[i]=="FB2"){sample[i]="FB2-AR"}
if(sample[i]=="AR14"&&sub[i]=="FB3"){sample[i]="FB3-AR"}
if(sample[i]=="AR14"&&sub[i]=="FB4"){sample[i]="FB4-AR"}
if(sample[i]=="NAR14"&&sub[i]=="FB1"){sample[i]="FB1-NAR"}
if(sample[i]=="NAR14"&&sub[i]=="FB2"){sample[i]="FB2-NAR"}
if(sample[i]=="NAR14"&&sub[i]=="FB3"){sample[i]="FB3-NAR"}
if(sample[i]=="NAR14"&&sub[i]=="FB4"){sample[i]="FB4-NAR"}
if(sample[i]=="P14"&&sub[i]=="FB1"){sample[i]="FB1-P"}
if(sample[i]=="P14"&&sub[i]=="FB2"){sample[i]="FB2-P"}
if(sample[i]=="P14"&&sub[i]=="FB3"){sample[i]="FB3-P"}
if(sample[i]=="P14"&&sub[i]=="FB4"){sample[i]="FB4-P"}
}
nproj$genotype<-sample
cell<-nproj@cellColData
A1<-cell[cell$genotype=="FB1-AR",]
A2<-cell[cell$genotype=="FB2-AR",]
A3<-cell[cell$genotype=="FB3-AR",]
A4<-cell[cell$genotype=="FB4-AR",]
N1<-cell[cell$genotype=="FB1-NAR",]
N2<-cell[cell$genotype=="FB2-NAR",]
N3<-cell[cell$genotype=="FB3-NAR",]
N4<-cell[cell$genotype=="FB4-NAR",]
P1<-cell[cell$genotype=="FB1-P",]
P2<-cell[cell$genotype=="FB2-P",]
P3<-cell[cell$genotype=="FB3-P",]
P4<-cell[cell$genotype=="FB4-P",]

A1<-rownames(A1)
A2<-rownames(A2)
A3<-rownames(A3)
A4<-rownames(A4)
N1<-rownames(N1)
N2<-rownames(N2)
N3<-rownames(N3)
N4<-rownames(N4)
P1<-rownames(P1)
P2<-rownames(P2)
P3<-rownames(P3)
P4<-rownames(P4)

A1<-gsub("^.*#","",A1)
A2<-gsub("^.*#","",A2)
A3<-gsub("^.*#","",A3)
A4<-gsub("^.*#","",A4)
N1<-gsub("^.*#","",N1)
N2<-gsub("^.*#","",N2)
N3<-gsub("^.*#","",N3)
N4<-gsub("^.*#","",N4)
P1<-gsub("^.*#","",P1)
P2<-gsub("^.*#","",P2)
P3<-gsub("^.*#","",P3)
P4<-gsub("^.*#","",P4)


###snap###

library(SnapATAC);
library(GenomicRanges);
library(viridisLite);
x1.sp = readRDS("/public/home/nieyg/scATAC/denovo/fibroblast/fibroblast.rds");

######scRNA-seq based annotation

library(Seurat);
pbmc.rna = readRDS("/public/home/nieyg/scATAC/integrate/cmc_sct.rds");  
pbmc.rna$tech = "rna";
variable.genes = VariableFeatures(object = pbmc.rna);
genes.df = read.table("/public/home/nieyg/scATAC/findmarker/gencode.vM16.gene.bed");
genes.gr = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4]);
genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)];

## reload the bmat
x1.sp = addBmatToSnap(x1.sp);
x1.sp = createGmatFromMat(
    obj=x1.sp, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=10
  );
#####convert the snap object to Seurat object
pbmc.atac <- snapToSeurat(
    obj=x1.sp, 
    eigs.dims=1:6, 
    norm=TRUE,
    scale=TRUE
  );
  
transfer.anchors <- FindTransferAnchors(
    reference = pbmc.rna, 
    query = pbmc.atac, 
    features = variable.genes, 
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca"
  );
  
#####   Create psudo multiomics cells
refdata <- GetAssayData(
    object = pbmc.rna, 
    assay = "RNA", 
    slot = "data"
  );
imputation <- TransferData(
    anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = pbmc.atac[["SnapATAC"]], 
    dims = 1:6
  );  
x1.sp@gmat = t(imputation@data);  
rm(imputation); 
rm(refdata);    
rm(pbmc.rna);   
rm(pbmc.atac);  
genemat<-x1.sp@gmat
rownames(genemat)<-gsub("-.*$","-1",rownames(genemat))
cebp<-c("Cebpa","Cebpb","Cebpd","Cebpe","Cebpg","Cebpz")
genemat<-genemat[,which(colnames(genemat)%in%cebp)]
genemat<-as.data.frame(genemat)
rownames(genemat)<-gsub(".1$","-1",rownames(genemat))

A1mat<-genemat[which(rownames(genemat)%in%A1),]
A2mat<-genemat[which(rownames(genemat)%in%A2),]
A3mat<-genemat[which(rownames(genemat)%in%A3),]
A4mat<-genemat[which(rownames(genemat)%in%A4),]
N1mat<-genemat[which(rownames(genemat)%in%N1),]
N2mat<-genemat[which(rownames(genemat)%in%N2),]
N3mat<-genemat[which(rownames(genemat)%in%N3),]
N4mat<-genemat[which(rownames(genemat)%in%N4),]
P1mat<-genemat[which(rownames(genemat)%in%P1),]
P2mat<-genemat[which(rownames(genemat)%in%P2),]
P3mat<-genemat[which(rownames(genemat)%in%P3),]
P4mat<-genemat[which(rownames(genemat)%in%P4),]


####创建绘制violin图的矩阵####
Group<-c(rep("FB1",6240),rep("FB2",62),rep("FB3",1114),rep("FB4",1271))
Group<-factor(Group)
Attribute<-c(rep("AR",1324),rep("NAR",2575),rep("P",2341),rep("AR",0),rep("NAR",42),rep("P",20),
	rep("AR",0),rep("NAR",564),rep("P",550),rep("AR",0),rep("NAR",667),rep("P",604))
Attribute<-factor(Attribute)
cA1<-A1mat[,5]
cA2<-A2mat[,5]
cA3<-A3mat[,5]
cA4<-A4mat[,5]
cN1<-N1mat[,5]
cN2<-N2mat[,5]
cN3<-N3mat[,5]
cN4<-N4mat[,5]
cP1<-P1mat[,5]
cP2<-P2mat[,5]
cP3<-P3mat[,5]
cP4<-P4mat[,5]
value<-c(cA1,cN1,cP1,cA2,cN2,cP2,cA3,cN3,cP3,cA4,cN4,cP4)
data<-data.frame(Group=Group,Attribute=Attribute,value=value)
###创建summary函数##
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
###############
Data_summary <- summarySE(data, measurevar="value", groupvars=c("Group","Attribute"))
P2<- ggplot(data, aes(x=Group, y=value,fill=Attribute)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_point(data = Data_summary,aes(x=Group, y=value),pch=19,position=position_dodge(0.9),size=1.5)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) + 
  scale_fill_manual(values = c("#FEE500", "#8A9FD1","#C06CAB"))+ #设置填充的颜色,AR,NAR,P
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=0,hjust = 0,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("Cebpz")
pdf("Cebpz-exp.pdf")
P2
dev.off();






#####血管生成####
positive regulation of angiogenesis

negative regulation of angiogenesis

mart<-read.csv("mart_export.csv",header=TRUE)
negative<-mart[which(mart$GO.term.name=="negative regulation of angiogenesis"),]
negativegene<-negative$Gene.name
negativecount<-subcount[which(rownames(subcount)%in%negativegene),]

negativecount<-apply(negativecount,2,as.numeric)
pdf("negative-angiogenesis.pdf",width=5, height=8)

count=t(scale(t(negativecount),scale = T,center = T))
pheatmap(count, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)

dev.off()
positive<-mart[which(mart$GO.term.name=="positive regulation of angiogenesis"),]
positivegene<-positive$Gene.name
positivecount<-subcount[which(rownames(subcount)%in%positivegene),]

positivecount<-apply(positivecount,2,as.numeric)
pdf("positive-angiogenesis.pdf",width=5, height=8)

count=t(scale(t(positivecount),scale = T,center = T))
pheatmap(count, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)

dev.off()
#####提取accessbility矩阵
#####调取各个亚群的基因分数矩阵######
GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)

genename<-rowData(GM)$name ####共24333个基因
count<-GM@assays@data$GeneScoreMatrix 
A1<-cell[cell$genotype=="FB1-AR",]
A2<-cell[cell$genotype=="FB2-AR",]
A3<-cell[cell$genotype=="FB3-AR",]
A4<-cell[cell$genotype=="FB4-AR",]
N1<-cell[cell$genotype=="FB1-NAR",]
N2<-cell[cell$genotype=="FB2-NAR",]
N3<-cell[cell$genotype=="FB3-NAR",]
N4<-cell[cell$genotype=="FB4-NAR",]
P1<-cell[cell$genotype=="FB1-P",]
P2<-cell[cell$genotype=="FB2-P",]
P3<-cell[cell$genotype=="FB3-P",]
P4<-cell[cell$genotype=="FB4-P",]

A1<-rownames(A1)
A2<-rownames(A2)
A3<-rownames(A3)
A4<-rownames(A4)
N1<-rownames(N1)
N2<-rownames(N2)
N3<-rownames(N3)
N4<-rownames(N4)
P1<-rownames(P1)
P2<-rownames(P2)
P3<-rownames(P3)
P4<-rownames(P4)


A1count<-count[,which(colnames(count)%in%A1)]
A2count<-count[,which(colnames(count)%in%A2)]
A3count<-count[,which(colnames(count)%in%A3)]
A4count<-count[,which(colnames(count)%in%A4)]
N1count<-count[,which(colnames(count)%in%N1)]
N2count<-count[,which(colnames(count)%in%N2)]
N3count<-count[,which(colnames(count)%in%N3)]
N4count<-count[,which(colnames(count)%in%N4)]
P1count<-count[,which(colnames(count)%in%P1)]
P2count<-count[,which(colnames(count)%in%P2)]
P3count<-count[,which(colnames(count)%in%P3)]
P4count<-count[,which(colnames(count)%in%P4)]

A1count<-as.matrix(A1count)
A2count<-as.matrix(A2count)
A3count<-as.matrix(A3count)
A4count<-as.matrix(A4count)
N1count<-as.matrix(N1count)
N2count<-as.matrix(N2count)
N3count<-as.matrix(N3count)
N4count<-as.matrix(N4count)
P1count<-as.matrix(P1count)
P2count<-as.matrix(P2count)
P3count<-as.matrix(P3count)
P4count<-as.matrix(P4count)

A1count<-rowMeans(A1count)
A2count<-rowMeans(A2count)
A3count<-rowMeans(A3count)
A4count<-rowMeans(A4count)
N1count<-rowMeans(N1count)
N2count<-rowMeans(N2count)
N3count<-rowMeans(N3count)
N4count<-rowMeans(N4count)
P1count<-rowMeans(P1count)
P2count<-rowMeans(P2count)
P3count<-rowMeans(P3count)
P4count<-rowMeans(P4count)

subcount<-cbind(A1count,N1count)
subcount<-cbind(subcount,P1count)
subcount<-cbind(subcount,A2count)
subcount<-cbind(subcount,N2count)
subcount<-cbind(subcount,P2count)
subcount<-cbind(subcount,A3count)
subcount<-cbind(subcount,N3count)
subcount<-cbind(subcount,P3count)
subcount<-cbind(subcount,A4count)
subcount<-cbind(subcount,N4count)
subcount<-cbind(subcount,P4count)

subcount<-cbind(genename,subcount)
rownames(subcount)<-subcount[,1]
subcount<-subcount[,-1]
subcount<-apply(subcount,2,as.numeric)
n<-c("FB1-AR","FB1-NAR","FB1-P","FB2-AR","FB2-NAR","FB2-P","FB3-AR","FB3-NAR","FB3-P","FB4-AR","FB4-NAR","FB4-P")
colnames(subcount)<-n
rownames(subcount)<-genename

negativecount<-subcount[which(rownames(subcount)%in%negativegene),]

pdf("negative-angiogenesis-archR.pdf",width=5, height=8)
count=t(scale(t(negativecount),scale = T,center = T))
pheatmap(count, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)
dev.off()

positivecount<-subcount[which(rownames(subcount)%in%positivegene),]
pdf("positive-angiogenesis-archR.pdf",width=5, height=8)
count=t(scale(t(positivecount),scale = T,center = T))
pheatmap(count, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)
dev.off()





#############更改umap图片的颜色
Subtype = nproj@cellColData$subtype
###timepoint
pdf("FB-4-popu-changecolor.pdf",width=5,height=5)
Subtype = nproj@cellColData$subtype
df = data.frame(nproj@embeddings$UMAP$df,nproj@cellColData$subtype)
colnames(df)<-c("UMAP.Dimension.1","UMAP.Dimension.2","subtype")
ggplot(df, aes(x = UMAP.Dimension.1, y = UMAP.Dimension.2, color = Subtype)) +
  geom_point(size=0.6,alpha=1) + 
  theme_bw() + 
  scale_color_manual(values = c("#A20056B2","#3B4992B2","#FA8639","#D51F26")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")+ggtitle("FB")
dev.off()
pdf("SMC-5-popu-changecolor.pdf",width=5,height=5)
Subtype = proj@cellColData$subtype
df = data.frame(proj@embeddings$UMAP$df,proj@cellColData$subtype)
colnames(df)<-c("UMAP.Dimension.1","UMAP.Dimension.2","subtype")
ggplot(df, aes(x = UMAP.Dimension.1, y = UMAP.Dimension.2, color = Subtype)) +
  geom_point(size=0.6,alpha=1) + 
  theme_bw() + 
  scale_color_manual(values = c("#A6D587","#5C89CB","#EC9A9D","#AB67B0")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")+ggtitle("FB")
dev.off();


#################################################
##########FB2's genes overlap with terms ########
################################################# 


markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
pdf("FB234term-overlap-ArchR.pdf")

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.1")
mart<-read.csv("mart_export.csv")
FB2<-markerList$FB2
FB2<-FB2$name
ECMgene<-mart[which(mart$GO.term.name=="cardiac muscle tissue regeneration"),]

ECMgene<-ECMgene$Gene.name
gg<-intersect(FB2,ECMgene)
count1<-subcount[which(rownames(subcount)%in%gg),]

count1=t(scale(t(count1),scale = T,center = T))

pheatmap(count1, cluster_rows=TRUE, 
               clustering_method="ward.D2",
               cellwidth = 50, cellheight =25,
               show_rownames=T,border_color = NA,
               colorRampPalette(colors = c("blue","white","red"))(100),
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)








dev.off()












####创建绘制violin图的矩阵####
Group<-c(rep("FB1",10889),rep("FB2",1370),rep("FB3",2796),rep("FB4",5776))
Group<-factor(Group)
Attribute<-c(rep("AR",2440),rep("NAR",4221),rep("P",4228),rep("AR",1179),rep("NAR",142),rep("P",49),
	rep("AR",905),rep("NAR",963),rep("P",928),rep("AR",2265),rep("NAR",1384),rep("P",2127))
Attribute<-factor(Attribute)
cA1<-A1mat[,1]
cA2<-A2mat[,1]
cA3<-A3mat[,1]
cA4<-A4mat[,1]
cN1<-N1mat[,1]
cN2<-N2mat[,1]
cN3<-N3mat[,1]
cN4<-N4mat[,1]
cP1<-P1mat[,1]
cP2<-P2mat[,1]
cP3<-P3mat[,1]
cP4<-P4mat[,1]
value<-c(cA1,cN1,cP1,cA2,cN2,cP2,cA3,cN3,cP3,cA4,cN4,cP4)
data<-data.frame(Group=Group,Attribute=Attribute,value=value)
d<-data[data$value>5,]
c<-data[-which(rownames(data)%in% rownames(d)),]
Data_summary <- summarySE(c, measurevar="value", groupvars=c("Group","Attribute"))
P2<- ggplot(c, aes(x=Group, y=value,fill=Attribute)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_point(data = Data_summary,aes(x=Group, y=value),pch=19,position=position_dodge(0.9),size=1.5)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) + 
  expand_limits(y=c(0,2))+
  scale_fill_manual(values = c("#FEE500", "#8A9FD1","#C06CAB"))+ #设置填充的颜色,AR,NAR,P
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=0,hjust = 0,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("Cebpe")
pdf("Cebpe-acc.pdf")
P2 
dev.off();


write.table(A1mat,"A1mat.txt")
write.table(A2mat,"A2mat.txt")
write.table(A3mat,"A3mat.txt")
write.table(A4mat,"A4mat.txt")
write.table(N1mat,"N1mat.txt")
write.table(N2mat,"N2mat.txt")
write.table(N3mat,"N3mat.txt")
write.table(N4mat,"N4mat.txt")
write.table(P1mat,"P1mat.txt")
write.table(P2mat,"P2mat.txt")
write.table(P3mat,"P3mat.txt")
write.table(P4mat,"P4mat.txt")



A1mat<-read.table("A1mat.txt")
A2mat<-read.table("A2mat.txt")
A3mat<-read.table("A3mat.txt")
A4mat<-read.table("A4mat.txt")
N1mat<-read.table("N1mat.txt")
N2mat<-read.table("N2mat.txt")
N3mat<-read.table("N3mat.txt")
N4mat<-read.table("N4mat.txt")
P1mat<-read.table("P1mat.txt")
P2mat<-read.table("P2mat.txt")
P3mat<-read.table("P3mat.txt")
P4mat<-read.table("P4mat.txt")
1