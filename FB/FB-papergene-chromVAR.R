
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-fibroblast2")


nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI-dim15-fibroblast-1.2-40000")

markerGenes  <- c(
    "Scn7a","Bmper","Acsm1",###paper-FB1
    "Cfh","ID4","Kcnt2",###paper-FB2
    "Ptx3","Osmr","Il6st",###FB3，细胞因子受体 OSMR,ILST6
    "Postn","Tnc","Fap",###FB4,genes responsive to TGF-signaling
    "Fbln2","Pcolce2",####FB5,ECM产生,重塑和降解.
    "Cd36","Egflam","Ftl1"####FB6
      )

p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "FB-papergene-UMAP-Imputation.pdf", 
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
plotPDF(plotList = p2, 
    name = "FB2-papergene-violin.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 10, height = 3)
#Rearrange for grid plotting

pFB5<-c("PCOLCE2","MFAP5","FBN1","CFD","DCN","IGFBP6","S100A6","UAP1","GSN","GFPT2","TIMP1","CREB5","ABI3BP",
	"SCARA5","TNXB","FSTL1","CCDC80","HTRA3","ITM2B","S100A4","CST3","ACKR3","B2M","S100A10","S100A11","FBLN2",
	"ADAMTS5","ACTB","LAPTM4A","SERPINF1","CD63","PLA2G2A","LTBP4","ANXA2","C1R","PLAC9","VIM","SEMA3C","H3F3B",
	"VCAN","SERPING1","FTL","EMP3","LINC01133","PSAP","PTGIS","TMSB10","PMP22","CD55","C1S","TIMP3","CD34","MMP2",
	"NTN4","AXL","TAGLN2","ITGA11","IFITM3","MFAP4","MT-CO1","FTH1","SH3BGRL3","DST","SDC2","EIF1","SDK1","HLA-C",
	"LUM","RPLP1","SOD3","CYP4B1","TIMP2","LDHA","KLF4","ZBTB7C","MT-ND4","RPS15","PRDX1","TMSB4X","MT-CYB","PTMA",
	"RPL11","ANXA1","FP700111.1","EEF1A1","RPL15","MT-CO3","RPL41","RPL13A","ECM1","RPL34","PI16","GSTP1","JUN","UBB",
	"GRK5","GAPDH","DBN1","IFITM2","RPL36")
markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.3")
> a
 [1] "Gfpt2" "Grk5"  "Dst"   "Uap1"  "Sod3"  "Ackr3" "Mfap5" "Fbn1"  "Ptgis"
[10] "Fstl1" "Htra3"

gf<-getFeatures(
       ArchRProj = nproj,
       useMatrix = "GeneScoreMatrix",
       select = NULL,
       ignoreCase = TRUE
     )
c<-intersect(pFB5,gf)
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = c, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "FB-papergene-pFB5-UMAP-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = c,######genename
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )
plotPDF(plotList = p2, 
    name = "FB2-papergene-pFB5-violin.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 10, height = 3)



######FB2-AR and paper-FB3#####
pFB3<-c("JUNB","NAMPT","MT2A","ZFP36","NR4A1","GADD45B","FOS","EGR1","THBS1","SAT1","ATP1B3","ELL2","ADAMTS4","GFPT2",
	"SOCS3","3-Mar","VMP1","STAT3","MT1M","CCNL1","C11orf96","CYR61","MT1X","CCL2","DUSP1","SPSB1","DDX21","IER2","UAP1",
	"SLC2A3","SLC39A14","NFIL3","CRISPLD2","PLA2G2A","PDE4D","YBX3","ACSL4","ATP13A3","ARID5B","OSMR","AKAP12","GALNT2",
	"KAZN","CD44","BTG2","LITAF","MT1A","FGF7","PNRC1","NNMT","CBLB","CFH","COL4A1","PVT1","NEAT1","CHD1","FOSL2","PGAP1",
	"MTHFD1L","STEAP4","VCAN","C1QTNF1","MCL1","IL1R1","MEDAG","PTPRG","ITPKC","GPRC5A","NABP1","TMEM165","H3F3B","HIVEP2",
	"PLPP3","CEBPD","IFI16","CHSY1","SLC4A7","SIK3","ZFP36L1","SLC2A13","SLC19A2","THBS2","JUN","IGFBP4","FSIP1","MT1E",
	"COL4A2","BHLHE40","ID2","IL4R","B4GALT1","RUNX1","PITPNC1","CDKN1A","TNFAIP6","ZNF331","RND3","PABPC1","CTSL","DENND5A")

pFB3<-tolower(pFB3)
pFB3<-capitalize(pFB3)
c<-intersect(pFB3,gf)
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = c, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "FB-papergene-pFB3-UMAP-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name = c,######genename
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )
plotPDF(plotList = p2, 
    name = "FB2-papergene-pFB3-violin.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 10, height = 3)


markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

PM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "PeakMatrix",
  verbose = FALSE
)
pmat<-PM@assays@data$PeakMatrix
peaklo<-rowData(markersPeaks)
n<-getPeakSet(nproj)
v<-n[which(n@ranges@NAMES=="FB2"),]

FB2peakbed<-v %>% {paste0(seqnames(.), "    ", start(.),"    ",end(.),idx(.))}

######ArchR chromVAR#######



if("Motif" %ni% names(nproj@peakAnnotation)){
    nproj <- addMotifAnnotations(ArchRProj = nproj, motifSet = "cisbp", name = "Motif")
}
nproj <- addBgdPeaks(nproj)
nproj <- addDeviationsMatrix(
  ArchRProj = nproj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(nproj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "FB-Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = nproj, addDOC = FALSE)
motifs <- c("TCF21", "NFATC2", "ATOH7","TBP","MEF2C", ####FB1
	"CEBPD","LINE10076",###FB2
	"WT1","ZFP148","SMAD1","ZIC4","E2F1","TCF4","SP1",###FB3
	"SMARCC1", "FOS", "BACH1","NFE2","JUNB","BATF3")####FB4

motifs <- c("TCF21", "NFATC2", "IRF1","ATOH7","TBP","BHLHA15","AHCTF1","MSC","NFATC4","PHF21A","MEF2C")####FB1



markerMotifs <- getFeatures(nproj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs
p <- plotGroups(ArchRProj = nproj, 
  groupBy = "genotype", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(nproj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
plotPDF(p2, name = "FB-chromVAR-Deviations-density", width = 5, height = 5, ArchRProj = nproj, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)
plotPDF(plotList = p, 
    name = "FB-chromVAR-devscore-umap-1.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 10, height = 3)




perl ~/softwares/miniconda3/bin/findMotifsGenome.pl encc-enhancer-atac.promt.Homer.bed hg38 promt.motif -p 10 -size 200 -find jaspar.motifs > promt.jaspar.txt

###find motif localtion in DEP
for ((i=1 ;i<=9;i++));
do
findMotifsGenome.pl newCM1.bed mm10 motif/CM1/homerResults/ -find motif/CM1/homerResults/motif1.motif > ./Motif_position/CM1/homer_motif1.txt;
done

for ((i=1 ;i<=30;i++));
do
findMotifsGenome.pl newCM1.bed mm10 motif/CM1/knownResults/ -find motif/CM1/knownResults/known1.motif > ./Motif_position/CM1/motif1.txt;
done




b<-c("D030024E09Rik","Bnipl")
p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = b, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "cebptargetgene-UMAP-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

p2 <- plotGroups(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    colorBy = "GeneScoreMatrix", 
	name =b,######genename
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 4,
    addBoxPlot = TRUE,
    size=0.2
   )
plotPDF(plotList = p2, 
    name = "cebptargetgene-violin.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
    ArchRProj = nproj, 
    groupBy = "genotype", 
    geneSymbol = c("D030024E09Rik"), 
    upstream =  300000,
    downstream = -10000)
plotPDF(plotList = p, 
    name = "cebptargetgene-trackD030024E09Rik.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

p<-ggplot(df, aes(x = d1, y = d2, color = bnipl)) +
  geom_point(size=0.8,alpha=0.5) + 
  theme_bw()  +  scale_color_gradient(low="white",high="red")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("D030024E09Rik")





plotPDF(plotList = p, 
    name = "D0rik.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)





colorRampPalette(colors = c("white","blue"))(100)


df = data.frame(nproj@embeddings$UMAP$df, pmat1[3240,],pmat1[15305,],pmat1[24209,],pmat1[35996,],pmat1[43014,],pmat1[50530,],pmat1[56862,],pmat1[64190,],pmat1[69895,],
	pmat1[77144,],pmat1[83249,],pmat1[88152,],pmat1[100942,],pmat1[109570,],pmat1[120337,],pmat1[130811,],pmat1[139743,],pmat1[149755,],pmat1[158630,])
p = pmat1[3240,]

ggplot(df, aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = df[,9])) +
  geom_point(size=0.8,alpha=0.6) + 
  theme_bw() + 
  scale_color_gradient(low="grey",high="red") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("D030024E09Rik")


3240  15305  24209  35996  43014  50530  56862  64190  69895  77144
83249  88152 100942 109570 120337 130811 139743 149755 158630




