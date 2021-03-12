


library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")


nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI-dim15-fibroblast-1.2-40000")

markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
nproj <- addMotifAnnotations(ArchRProj = nproj, motifSet = "cisbp", name = "Motif")

####To prepare this data for plotting with ggplot we can create a simplified data.frame object containing the motif names,
#### the corrected p-values, and the significance rank.b -
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

motiflist<-c("Klf4_143","Klf1_214","Klf3_804","Klf2_819","Klf8_194","Klf12_236",
           "Sp3_802","Smad5_883","Sp2_150","Sp5_238","Sp6_191","Sp4_167","Klf7_171","Klf5_145","E2f4_257","Klf14_237","Klf11_798","Klf13_817","Bach2_119",
		   "Tcf3_31","Tcf12_59","Ascl2_23","Klf15_806","Klf16_874","Sp7_222","Sp8_207","Sp9_231","Zfp219_816","Gata5_385","Gata4_386","Gata1_387","Mesp2_57",
		   "Gata2_383","Sox9_725","Egr2_188","Egr3_183","Zfp161_209","Zic5_195","Tcfap2c_3","Gata3_384","Zic3_229","Zbtb7a_808","Zbtb7c_812","Elk1_274",
		   "Zfp281_193","Zfp740_204","Plagl1_152","Tcfap2a_1","Zic2_225","LINE3878_878","Egr4_234","Zic4_187","Tcfl5_71","Klf10_810","Smad4_860",
		   "E2f6_264","Zic1_181","Fosb_98","Nfe2l2_101","Jund_135","Batf_790","Fosl1_107",
		   "Nfe2l1_117","Nfe2l3_871","Jun_126","Cebpb_130","Cebpg_129","Zfp202_169","Tbpl2_773","Nfatc1_703","Nfatc3_702","Nfata4_857")
motif<-heatmapEM@ column_names_param$ labels 
motiflist<-c(motif,motiflist)

idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% motiflist)
cellsSample <- enrichMotifs@NAMES[-idxSample]
enrichMotifs1=enrichMotifs[cellsSample, ]
heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 50, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "FB-Motifs-FDR01",width = 12, height = 4,ArchRProj = nproj, addDOC = FALSE)


#####正向selectTF困难，直接用反向选择#########
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 1 & Log2FC >= 0.5"
  )

df <- data.frame(TF = rownames(enrichMotifs), mlog10Padj = assay(enrichMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

FB4<-c("Smarcc1_843","Fos_104","Bach1_108","Batf3_789","Junb_127")
FB3<-c("Sp2_150","Wt1_872","Klf7_171","Tcfap2d_5","Smad1_863")
FB2<-c("Cebpd_97","LINE10076_380","Tead4_868")
FB1<-c("Tcf21_79","Nfatc2_700","Irf1_631","Atoh7_64","Tbp_722","Bhlha15_87")
fb<-c(FB1,FB2,FB3,FB4)
idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% fb)
cellsSample <- enrichMotifs@NAMES[idxSample]
enrichMotifs1=enrichMotifs[cellsSample, ]
heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 100, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "FB-Motifs-FDR1-2",width = 4, height = 8,ArchRProj = nproj, addDOC = FALSE)


                          FB1        FB2 FB3       FB4
Tcf21_79 (35)      1.00000000 0.00000000   0 0.0000000
Nfatc2_700 (30)    1.00000000 0.03011820   0 0.0000000
Irf1_631 (27)      1.00000000 0.61014972   0 0.2411853
Atoh7_64 (23)      1.00000000 0.00414133   0 0.0000000
Bhlha15_87 (23)    1.00000000 0.02433855   0 0.0000000
Cebpd_97 (35)      0.00000000 1.00000000   0 0.1043674
LINE10076_380 (24) 0.50252680 1.00000000   0 0.6171055
Sp2_150 (113)      0.00000000 0.00000000   1 0.0000000
Wt1_872 (111)      0.00000000 0.00000000   1 0.0000000
Tcfap2d_5 (98)     0.00000000 0.00000000   1 0.0000000
Smad1_863 (90)     0.00000000 0.00000000   1 0.0000000
Klf7_171 (72)      0.00000000 0.00000000   1 0.0000000
Smarcc1_843 (243)  0.01293673 0.00000000   0 1.0000000
Fos_104 (220)      0.03829238 0.00000000   0 1.0000000
Bach1_108 (149)    0.03737454 0.00000000   0 1.0000000
Junb_127 (89)      0.11385678 0.00000000   0 1.0000000
Batf3_789 (84)     0.07229696 0.00000000   0 1.0000000


df <- data.frame(TF = rownames(enrichMotifs), assay(enrichMotifs)[,c(1,2,3,4)])
 
                      TF        FB1         FB2       FB3        FB4
Tcf21_79           Tcf21_79 34.7058575  0.00000000   0.00000   0.000000
Nfatc2_700       Nfatc2_700 29.6185575  0.89205749   0.00000   0.000000
Irf1_631           Irf1_631 26.5586575 16.20475749   0.00000   6.405557
Atoh7_64           Atoh7_64 23.4604575  0.09715749   0.00000   0.000000
Bhlha15_87       Bhlha15_87 23.0727575  0.56155749   0.00000   0.000000
Cebpd_97           Cebpd_97  0.0000000 34.59565749   0.00000   3.610657
LINE10076_380 LINE10076_380 12.0080575 23.89535749   0.00000  14.745957
Tead4_868         Tead4_868  0.4047575  9.85465749   0.00000   0.000000
Sp2_150             Sp2_150  0.0000000  0.00000000 113.38016   0.000000
Wt1_872             Wt1_872  0.0000000  0.00000000 110.76536   0.000000
Klf7_171           Klf7_171  0.0000000  0.00000000  71.55516   0.000000
Tcfap2d_5         Tcfap2d_5  0.0000000  0.00000000  98.31636   0.000000
Smad1_863         Smad1_863  0.0000000  0.00000000  89.52786   0.000000
Smarcc1_843     Smarcc1_843  3.1455575  0.00000000   0.00000 243.149257
Fos_104             Fos_104  8.4425575  0.00000000   0.00000 220.476157
Bach1_108         Bach1_108  5.5754575  0.00000000   0.00000 149.177957
Batf3_789         Batf3_789  6.0873575  0.00000000   0.00000  84.199357
Junb_127           Junb_127 10.1225575  0.00000000   0.00000  88.906057

##针对原始mlog10adj矩阵直接进行heatmap

pdf("tf-df-mlog10adj-archR.pdf",width=4, height=6)
df=t(scale(t(df),scale = T,center = F))
pheatmap(df, cluster_rows=FALSE, 
               clustering_method="ward.D2",
               show_rownames=TRUE,border_color = "black",
               cluster_cols=FALSE,cutree_rows =1,
               show_colnames = TRUE)
dev.off()

findMotifsGenome.pl FB2-2-2.bed mm10 ./ -find motif1.motif > ./motif1.txt;
annotatePeaks.pl motif1-1.bed mm10 > motif1_peak.annotation.xls

targetgene<-c("Pim1","Rhob","Mdm4","Dram1","Dusp6","D030024E09Rik","Plac8","1700016H13Rik","Atp6v1b2","Phactr1","4930549C15Rik","Celf2",
              "Cdkl4","Fabp12","4930563F08Rik","Auh","Uhrf1bp1l","Mir21a","Pag1","Slc44a3","Gm5079","Nfkbia","Ociad1","Hspa5")
#####"Gm5079"找不到##########
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = targetgene, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p,     name = "cebpd-targetgene-UMAP-Imputation.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

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
p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = targetgene,
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )
plotPDF(plotList = p2, name = "cebpd-targetgene-violinplot.pdf",   ArchRProj = nproj,    addDOC = FALSE, width = 10, height = 3)
#####大部分gene在FB2处富集#####
GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)

genename<-rowData(GM)$name ####共24333个基因
count<-GM@assays@data$GeneScoreMatrix
rownames(count)<-genename
targetcount<-count[which(rownames(count)%in%targetgene),]
#####求平均####
targetcount<-colMeans(targetcount)
barcodeorder<-rownames(nproj@embeddings$UMAP$df)
targetcount<-targetcount[barcodeorder]
m<-log10(targetcount+1)
df = data.frame(nproj@embeddings$UMAP$df, m)
colnames(df)<-c("IterativeLSI.UMAP_Dimension_1","IterativeLSI.UMAP_Dimension_2","logtargetcount")
pdf("cebpd-targetgene-avg-umap.pdf")
#m<-log10(targetcount+1)
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = df[,3])) +
  geom_point(size=0.8,alpha=0.6) + 
  theme_bw() + 
  scale_color_gradient(low="blue",high="yellow") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Cebpd Target Genes")
dev.off()
ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = targetcount) +
  geom_point(size=0.8,alpha=0.6) + 
  theme_bw() + 
  scale_color_gradient(low="grey",high="red") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("D030024E09Rik")

######use ECDF plot#####
FB1<-rownames(nproj@cellColData[which(nproj@cellColData$subtype=="FB1"),])
FB1<-targetcount[which(names(targetcount)%in%FB1)]
FB2<-rownames(nproj@cellColData[which(nproj@cellColData$subtype=="FB2"),])
FB2<-targetcount[which(names(targetcount)%in%FB2)]
FB3<-rownames(nproj@cellColData[which(nproj@cellColData$subtype=="FB3"),])
FB3<-targetcount[which(names(targetcount)%in%FB3)]
FB4<-rownames(nproj@cellColData[which(nproj@cellColData$subtype=="FB4"),])
FB4<-targetcount[which(names(targetcount)%in%FB4)]
FB1<-as.numeric(FB1)
FB2<-as.numeric(FB2)
FB3<-as.numeric(FB3)
FB4<-as.numeric(FB4)

FB1<-ecdf(FB1)
FB2<-ecdf(FB2)
FB3<-ecdf(FB3)
FB4<-ecdf(FB4)
pdf("cebpd-targetgene-ECDF.pdf")
plot(FB1,col="#A20056B2")
lines(FB2, col = "#3B4992B2")
lines(FB3, col = "#F48639")
lines(FB4, col = "#D51F26")

ks.test(FB2,FB1)
ks.test(FB2,FB3)
ks.test(FB2,FB4)





TM<-getMatrixFromProject(ArchRProj = nproj,useMatrix = "TileMatrix")
