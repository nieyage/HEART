########2020.08.13##########
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")
#####change color for treatment umap####
sample<-nproj$treatment
sample[1]=""
sample[2]="  "
sample[3]=" "
sample[4]="   "
sample[5]="    "
table(nproj$sample)
nproj$treatment2<-sample
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "treatment2", embedding = "UMAP")
ggAlignPlots(p3, type = "h")
plotPDF(p3, name = "FB-UMAP-treatment-changecolor.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

#######add early and later gene for marker gene heatmap############


markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes  <- c( "Mdk","Mfap2","Col1a1","Ptn","Col3a1","Fabp4","Gpihbp1","Cav1","Egfl7","Cd36",#####early
	"Mfap5","Pi16","Sema3c","Dcn","Gfpt2","Tcf21","Pdgfra","Ddr2"####Muture
)
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 0.1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "FB-early-muture-Marker-Heatmap", width = 6, height = 4, ArchRProj = nproj, addDOC = FALSE)

########New marker peak########
markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0", returnGR = TRUE)
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "FB-Peak-Heatmap-zongxiang", width = 4, height = 6, ArchRProj = nproj, addDOC = FALSE)

#######单个基因的小提琴图(区分不同的genotype)########


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
VEGFa，FGF，PDGF
例如TGFβ[转化生长因子β]，AngII [血管紧张素II]）
nproj <-addImputeWeights(
       ArchRProj = nproj,
       reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
       dimsToUse = 1:15
     )

aa<-c("Shh","Lox","Acan","Fn1","Lama1","Lama5","Lamb1","Lamc1","Lamc2","Postn",######再生相关ECM基因
"Hic1","Runx1","Rara","Rarb",####RUNA
"VEGFa","Pdgfra",#####血管再生相关
"Mmp3","Mmp10","Mmp1b","Mmp16","Mmp2","Mmp9","Mmp7",##ECM降解相关（MMP家族）
"Tnc","Erbb2","Acta2",
"Col5a1","Col5a2","Col5a3")

aa<-{"Pdgfa","Pdgfb"}
p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = aa,######genename
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

pdf("violinplot.pdf")
p2;
dev.off();
plotPDF(p2, name = "FB-violin-xueguan", width = 6, height =4, ArchRProj = nproj, addDOC = FALSE)

#######gene Umap ######


	 
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = aa, 
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
    name = "FB--UMAP-Marker-Genes-xueguan.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 4, height = 4)

MMPGenes  <- c( "Mmp1","Mmp2","Mmp3","Mmp4","Mmp5","Mmp6","Mmp7","Mmp8","Mmp9","Mmp10",
	            "Mmp11","Mmp12","Mmp13","Mmp14","Mmp15","Mmp16","Mmp17","Mmp18","Mmp19","Mmp20",
	            "Mmp21","Mmp22","Mmp23","Mmp24","Mmp25","Mmp26","Mmp28","Mmp23b","Mmp27"
	)

###find marker peak###
markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", returnGR = TRUE)
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "FB-sub-Peak-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
nproj <- addMotifAnnotations(ArchRProj = nproj, motifSet = "cisbp", name = "Motif",force = TRUE)




####marker TF#####

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 50, transpose = F,returnMatrix = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "FB-sub-Motifs", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
  

####去掉同一个基因家族的TF，仅保留一个，添加尽可能多的基因家族进来

  idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% c(
           "Sp3_802","Smad5_883","Sp2_150","Sp5_238","Sp6_191","Sp4_167","Bach2_119",
		   "Tcf3_31","Tcf12_59","Ascl2_23","Sp7_222","Sp8_207","Sp9_231","Zfp219_816","Gata5_385","Gata4_386","Gata1_387","Mesp2_57",
		   "Gata2_383","Sox9_725","Egr2_188","Egr3_183","Zfp161_209","Zic5_195","Tcfap2c_3","Gata3_384",
		   "Zfp281_193","Zfp740_204","Plagl1_152","Tcfap2a_1","Zic2_225","LINE3878_878",
		   "E2f6_264","Zic1_181","Fosb_98","Nfe2l2_101","Jund_135","Batf_790","Fosl1_107",
		   "Nfe2l1_117","Nfe2l3_871","Jun_126","Zfp202_169","Tbpl2_773","Nfatc1_703","Nfatc3_702","Nfata4_857"))



  idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% c("Zfp281_193","Zfp740_204","Plagl1_152","Tcfap2a_1","Zic2_225","LINE3878_878",
		   "E2f6_264","Zic1_181","Fosb_98","Nfe2l2_101","Jund_135","Batf_790","Fosl1_107",
		   "Nfe2l1_117","Nfe2l3_871","Jun_126","Cebpb_130","Cebpg_129"))


cellsSample <- enrichMotifs@NAMES[-idxSample]
enrichMotifs1=enrichMotifs[cellsSample, ]


heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 50, transpose = F,returnMatrix = TRUE)


pheatmap(heatmapEM,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("#E6E7E8", "#3A97FF", "#8816A7","black"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)



  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "FB-Motifs-rmrep2",width = 8, height = 4,ArchRProj = nproj, addDOC = FALSE)

#########Hic1 etc.TFS#

markerGenes  <- c(
    "Hic1","Runx1","Rara","Rarb"
  )

p <- plotBrowserTrack(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    geneSymbol = markerGenes, 
    upstream = 1000,
    downstream = 50000,
    ylim=c(0,1.2)
)
plotPDF(plotList = p, 
    name = "FB-Tracks-1.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)



count$FB1count<-as.numeric(count$FB1count)
count$FB2count<-as.numeric(count$FB2count)
count$FB3count<-as.numeric(count$FB3count)
count$FB4count<-as.numeric(count$FB4count)



count<-t(scale(t(count),scale = T,center = T))

pdf("ecm-type.pdf",width=6,height=12)
pheatmap(count,cluster_cols = F,show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,show_rownames=T,show_colnames=T)


proj2 <- addGeneIntegrationMatrix(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = RNA,
    addToArrow = TRUE,
    nameCell = "predictedCell_Un",
    nameScore = "predictedScore_Un"
)
nproj<- addGeneIntegrationMatrix(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = RNA,
    addToArrow = FALSE,
    groupRNA = "seurat_clusters",##seurat_clusters
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)


RNA <- FindClusters(
  object = RNA,
  reduction.type = "umap",
  dims.use = 1:9,
  resolution = 0.6,
  print.output = 0,
  save.SNN = TRUE)




n<-nproj@cellColData
FB1<-n[which(n$subtype=="FB1"),]
FB1<-rownames(FB1)
FB1<-gsub("^.*#","",FB1)
FB1<-gsub("-1$","",FB1)
FB1<-exp[which(rownames(exp)%in%FB1),]

FB2<-as.data.frame(FB2)
FB1<-colMeans(FB1)



