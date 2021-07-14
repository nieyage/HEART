library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")
proj5 <- addImputeWeights(proj5,reducedDims ="IterativeLSI-dim8")



########overlap的基因list#########
overlap<-c("Trap1","Cct2","Dnajb5","Dnajb11","Ppia","Ppil1","Stub1","Cdc37","Mkks","Grpel1","Fkbp8","Calr","P3h1",
	"Ptges3l","Dnajb1","Hsp90ab1","Hsp90aa1","Hspd1","Pdcl3","Nudc","Dnaja1","Pfdn6","Mesd","Hspa1b","Hspe1","Hspa8","Hspa1a")


##########提取violin plot 需求的4个subtype的信息#######
GM<-getMatrixFromProject(
  ArchRProj = proj5,
  useMatrix = "GeneScoreMatrix"
)
count<-assays(GM)$GeneScoreMatrix
gene<-rowData(GM)$name
rownames(count)<-gene
df = data.frame(proj5@embeddings$UMAP$df, subtype=c(rep("other",60420)),treatment=c(rep("other",60420)))
df$subtype<-as.character(df$subtype)
subtype<-df$subtype
m<-rownames(df)
#######subtype#######
for (i in 1:60420){
if(m[i] %in% FB3barcode)
    df[i,3]="FB3"
if(m[i] %in% MP3barcode)
    df[i,3]="MP3"
if(m[i] %in% EC4barcode)
    df[i,3]="EC4"
if(m[i] %in% SMC1barcode)
    df[i,3]="SMC1"    
    }
table(df$subtype)

#######treatment####
allcell<-rownames(df)
Pcell<-grep("^P.*",allcell,value=T)
ARcell<-grep("^AR.*",allcell,value=T)
NARcell<-grep("^NAR.*",allcell,value=T)

for (i in 1:60420){
if(m[i] %in% ARcell)
    df[i,4]="2AR"
if(m[i] %in% NARcell)
    df[i,4]="3NAR"
if(m[i] %in% Pcell)
    df[i,4]="1P"
    }
table(df$treatment)
colnames(df)<-c("UMAP_Dimension_1","UMAP_Dimension_2","subtype","treatment")

df<-df[which(df$subtype!="other"),]
count<-count[,which(colnames(count)%in% rownames(df))]
df$treatment<-as.factor(df$treatment)
######构建violin plot函数####
#violin plots represent gene expression for each subpopulation. The color of each violin represents the mean gene expression after log2 transformation.
#gene: Gene name of interest as string. DATAuse: Gene expression matrix with rownames containing gene names. tsne.popus = dbscan output, axis= if FALSE no axis is printed. legend_position= default "none" indicates where legend is plotted. gene_name = if FALSE gene name will not be included in the plot.
plot.violin2 <- function(gene, DATAuse, umap.popus, axis=FALSE, legend_position="none", gene_name=FALSE){

 testframe<-data.frame(expression=as.numeric(DATAuse[paste(gene),]), Population=umap.popus$subtype,Treatment=umap.popus$treatment)
 testframe$Population <- as.factor(testframe$Population)
 testframe$Treatment <- testframe$Treatment
 colnames(testframe)<-c("expression", "Population","Treatment")

 #testframe$expression<-log2(testframe$expression+1)

 p <- ggplot(testframe, aes(x=Population, y=expression, fill= Treatment))+
 geom_violin(scale="width") +geom_boxplot(width=0.2,position=position_dodge(0.9))+
 scale_fill_manual(values = c("#E64B35B2","#00468BB2","#42B540B2"))+
 labs(title=paste(gene))+
  ylab("log2(GeneScore)")+ #设置x轴和y轴的标题
 theme_classic() 
 #theme(axis.title.y = element_blank())+
 #theme(axis.ticks= element_blank())+
 #theme(axis.line.y = element_blank())
 #theme(axis.text.y = element_blank())+
 #theme(axis.title.x = element_blank())
p
}


#################画图###########
pdf("PF-overlaplist_violin.pdf")


p1<-plot.violin2(gene = overlap[25], DATAuse = count, umap.popus = df)
p2<-plot.violin2(gene = overlap[26], DATAuse = count, umap.popus = df)
p3<-plot.violin2(gene = overlap[27], DATAuse = count, umap.popus = df)

plot_grid(p1,p2,p3,ncol = 1, nrow = 3)
