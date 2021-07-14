########2020.08.13##########
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")

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
nproj <-addImputeWeights(
       ArchRProj = nproj,
       reducedDims = "IterativeLSI-dim15-fibroblast-1.2-40000",
       dimsToUse = 1:15
     )

#######gene Umap ######
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
    name = "FB--UMAP-Marker-Genes.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 4, height = 4)

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
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "FB-Motifs-rmrep",width = 8, height = 4,ArchRProj = nproj, addDOC = FALSE)

