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

A2<-rownames(A2)
N2<-rownames(N2)
P2<-rownames(P2)

genename<-rowData(GM)$name ####共24333个基因
count<-GM@assays@data$GeneScoreMatrix   #####

A2count<-count[,which(colnames(count)%in%A2)]
N2count<-count[,which(colnames(count)%in%N2)]
P2count<-count[,which(colnames(count)%in%P2)]
A2count<-as.matrix(A2count)
N2count<-as.matrix(N2count)
P2count<-as.matrix(P2count)
A2count<-rowMeans(A2count)
N2count<-rowMeans(N2count)
P2count<-rowMeans(P2count)
subcount<-cbind(A2count,N2count)
subcount<-cbind(subcount,P2count)

mart<-read.csv("mart_export.csv",header=TRUE)
negative<-mart[which(mart$GO.term.name=="negative regulation of angiogenesis"),]
negativegene<-negative$Gene.name
negativecount<-subcount[which(rownames(subcount)%in%negativegene),]
pdf("FB2-negative.pdf")
negativecount=t(scale(t(negativecount),scale = T,center = T))

list<-pheatmap(negativecount, cluster_rows=TRUE, 
        clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)
list
dev.off();


count3=t(scale(t(A2count),scale = T,center = T))
count3<-na.omit(count3)
pdf("FB2-AR-positive.pdf",width=24, height=10)
pheatmap(count3, cluster_rows=TRUE, 
        clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=TRUE,cutree_rows =1,
               show_colnames = F)
dev.off();

rowlist<-rownames(b[list$tree_row[["order"]],])

write.table(A2count,"A2count.txt")
write.table(N2count,"N2count.txt")
write.table(P2count,"P2count.txt")

######根据AR组的行顺序排序，不对行进行cluster########
pdf("FB2-NAR-1.pdf",,width=3, height=10)
count3=t(scale(t(N2count),scale = T,center = T))
b<-na.omit(count3)
#b<-b[rowlist,]
#pheatmap(b, cluster_rows=F, 
#        clustering_method="ward.D2",
#               show_rownames=F,border_color = NA,
#               cluster_cols=T,cutree_rows =1,
#               show_colnames = F)
#dev.off()
pheatmap(b, cluster_rows=T, 
        clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=T,cutree_rows =1,
               show_colnames = F)

pdf("FB2-P-1.pdf",,width=24, height=10)
count1=t(scale(t(P2count),scale = T,center = T))
b<-count1[-is.na(count1),]
b<-b[rowlist,]
pheatmap(b, cluster_rows=T, 
        clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=T,cutree_rows =1,
               show_colnames = F)
dev.off()




#####re draw umap 
df = data.frame(nproj@embeddings$UMAP$df, nproj$treatment)
df = data.frame(nproj@embeddings$UMAP$df, nproj$timepoint)
df = data.frame(nproj@embeddings$UMAP$df, nproj$subtype)
Sample = nproj$treatment
Timepoint = nproj$timepoint
Subtype = nproj$subtype

###timepoint
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = Timepoint)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#FFA500","#0073C2FF", "#20854EFF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Timepoint")

##subtype umap
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = Subtype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#BC3C29FF","#E18727FF","#0072B5FF","#7876B1FF","#20854EFF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Subtype")
  
#"#A20056B2","#3B4992B2", "#F48639","#D51F26" ####FB
#"#A6D587","#5C89CB", "#EC9A9D","#AB67B0" ####SMC
Cardiomyocytes    Endothelial     Fibroblast    ImmuneCells   Smoothmuscle

"#BC3C29FF","#E18727FF","#0072B5FF","#7876B1FF","#20854EFF" IC 
###sample
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2,color = Sample)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#E64B35B2","#00468BB2","#42B540B2" )) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Genotype")


####用Arch R 中的violin看cebp的gene score######
pdf("cebp-archr-violin.pdf")
aa<-c("Cebpa","Cebpb","Cebpd","Cebpe","Cebpg","Cebpz")

nproj <-addImputeWeights(
       ArchRProj = nproj,
       reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
       dimsToUse = 1:15
     )
	 
  p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = aa,######genename
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )




for(i in 1:144){
t[i]=sd(positivecount[i,])
}
df<-data.frame(gene=rep(rownames(d),3),genotype=c(rep("FB2-AR",20),rep("FB2-NAR",20),rep("FB2-P",20)),
	count=c(rep(1179,20),rep(142,20),rep(49,20)),genescore=c(d[,1],d[,2],d[,3])
)
df<-data.frame(gene=rep(rownames(f),3),genotype=c(rep("FB2-AR",30),rep("FB2-NAR",30),rep("FB2-P",30)),
	count=c(rep(1179,30),rep(142,30),rep(49,30)),genescore=c(f[,1],f[,2],f[,3])
)
library(ggplot2)
pdf("FB2-negative-bubbleplot.pdf")
ggplot(df, aes(x = genotype,y = gene),title="Gene score") +
  scale_color_gradient(low="blue",high = "red")+
  geom_point(aes(size=count,color=genescore))   +
  theme_bw();
dev.off()






markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.3")

markerGenes<-c("Cd36","Col3a1","Cav1","Ptn","Dcn","Egfl7",###early
	"Pdgfra","Ddr2","Gfpt2","Mfap5","Tcf21","Sema3c",####later
	"Il1rapl1","Cdh18","Grid2","Tenm2","Cntn5","Nlgn1","Cdh8",######cell-cell adhesion ###FB2
	"Egr1","Pdgfrb","Wt1",####cell proliferation  ###FB2
     "Kcnq5", "Tnp1","Nckap5","Tbx19", "Kif26b","Stxbp6", "Gm7550", "Esrrb",#####FB1
     "Sulf1","Zfp451", "Ackr3", "Cops8", "Pam", "Gpr39",####FB4
     "Hspa8","Hspa1l","Hspa1a","Hspa1b","Hsph1","Dnaja1","Hspb6","Calr","H2-DMa","Hspd1","Grn","Emc6","Bag1","Calr3"#####protein folding
	)
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.3", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
pdf("fibroblast-GeneScores-labelMarker-Heatmap.pdf",width = 10, height = 4)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "fibroblast-GeneScores-labelMarker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
dev.off()




nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI-dim15-fibroblast-1.2-40000")

markerGenes  <- c(
    "Il1rapl1","Cdh18","Grid2","Tenm2","Cntn5","Nlgn1","Cdh8",######cell-cell adhesion ###FB2
	"Egr1","Pdgfrb","Wt1",####cell proliferation
	"Trp63","Brinp1","Pim1","Birc6"####cell cycle

  )

p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
    name = "FB2-functiongene-UMAP-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)


p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = markerGenes,######genename
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )

ridgeScale




plotPDF(plotList = p2, 
    name = "FB2-functiongene-violin.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 10, height = 3)

cebp<-c("Cebpa","Cebpb","Cebpd","Cebpe","Cebpg","Cebpz")


p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = cebp,######genename
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 5,
    addBoxPlot = TRUE,
    size=0.2
   )

ridgeScale




plotPDF(plotList = p2, 
    name = "cebp-violin.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 10, height = 3)


 