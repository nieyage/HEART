library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")

idxSample <- BiocGenerics::which(proj5$Clusters %in% c("C14","C15","C17","C16","C18","C19","C20","C21"))
cellsSample <- proj5$cellNames[idxSample]
nproj=proj5[cellsSample, ]
saveArchRProject(ArchRProj = nproj, outputDirectory = "Save-EC", load = T)
nproj<-loadArchRProject(path = "Save-EC")

nproj
sample<-nproj$Sample
table(nproj$Sample)
for (i in 1:28055){
if(sample[i]=="AR3"){sample[i]="Day3"}
if(sample[i]=="AR7"){sample[i]="Day7"}
if(sample[i]=="AR14"){sample[i]="Day14"}
if(sample[i]=="NAR3"){sample[i]="Day3"}
if(sample[i]=="NAR7"){sample[i]="Day7"}
if(sample[i]=="NAR14"){sample[i]="Day14"}
if(sample[i]=="P3"){sample[i]="Day3"}
if(sample[i]=="P7"){sample[i]="Day7"}
if(sample[i]=="P14"){sample[i]="Day14"}
}

cellNames <- nproj$cellNames
nproj$timepoint<-sample2
nproj <- addCellColData(ArchRProj = nproj, data = paste0(sample),
    cells=cellNames, name = "timepoint",force = TRUE)

sample2<-nproj$Sample
table(nproj$Sample)
for (i in 1:28055){
if(sample2[i]=="AR3"){sample2[i]="AR"}
if(sample2[i]=="AR7"){sample2[i]="AR"}
if(sample2[i]=="AR14"){sample2[i]="AR"}
if(sample2[i]=="NAR3"){sample2[i]="NAR"}
if(sample2[i]=="NAR7"){sample2[i]="NAR"}
if(sample2[i]=="NAR14"){sample2[i]="NAR"}
if(sample2[i]=="P3"){sample2[i]="P"}
if(sample2[i]=="P7"){sample2[i]="P"}
if(sample2[i]=="P14"){sample2[i]="P"}
}
nproj$treatment<-sample2
nproj <- addCellColData(ArchRProj = nproj, data = paste0(sample2),
    cells=cellNames, name = "treatment",force = TRUE)
	table(nproj$treatment)

nproj <- addIterativeLSI(
    ArchRProj = nproj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI-dim15-fibroblast-1.2-50000", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(1.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 50000, 
    dimsToUse = 1:15
)

####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-dim15-fibroblast-1.2-50000",
    method = "Seurat",
    name = "Clusters2",
	force = TRUE,
    resolution = 0.8
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-dim15-fibroblast-1.2-50000", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.05, 
	force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "treatment", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")

ggAlignPlots(p1, p2,p3,p4, type = "h")
plotPDF(p1,p2,p3,p4, name = "Endothelial-Plot-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("Cd34","Flt1","Flk1","Kdr","Prom1","Ptprc","Eng","Vcam1","Kit","Kdr","Mcam","Pecam1","Cdh5","Nos3","Vwf","Plvap")

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)


GM<-getMatrixFromProject(
  ArchRProj = nproj,
  useMatrix = "GeneScoreMatrix"
)
head(GM@assays@data$GeneScoreMatrix)[,1:10]
coldata<-nroj@cellColData
gene<-GM@assays@data$GeneScoreMatrix

C1<-coldata[which(coldata$Clusters2=="C1"),]
C1gene<-gene[,which(colnames(gene)%in% rownames(C1))]
C1<-rowMeans(C1gene)
C2<-coldata[which(coldata$Clusters2=="C2"),]
C2gene<-gene[,which(colnames(gene)%in% rownames(C2))]
C2<-rowMeans(C2gene)
C3<-coldata[which(coldata$Clusters2=="C3"),]
C3gene<-gene[,which(colnames(gene)%in% rownames(C3))]
C3<-rowMeans(C3gene)
C4<-coldata[which(coldata$Clusters2=="C4"),]
C4gene<-gene[,which(colnames(gene)%in% rownames(C4))]
C4<-rowMeans(C4gene)

C5<-coldata[which(coldata$Clusters2=="C5"),]
C5gene<-gene[,which(colnames(gene)%in% rownames(C5))]
C5<-rowMeans(C5gene)

C6<-coldata[which(coldata$Clusters2=="C6"),]
C6gene<-gene[,which(colnames(gene)%in% rownames(C6))]
C6<-rowMeans(C6gene)

C7<-coldata[which(coldata$Clusters2=="C7"),]
C7gene<-gene[,which(colnames(gene)%in% rownames(C7))]
C7<-rowMeans(C7gene)
C8<-coldata[which(coldata$Clusters2=="C8"),]
C8gene<-gene[,which(colnames(gene)%in% rownames(C8))]
C8<-rowMeans(C8gene)
C9<-coldata[which(coldata$Clusters2=="C9"),]
C9gene<-gene[,which(colnames(gene)%in% rownames(C9))]
C9<-rowMeans(C9gene)

C10<-coldata[which(coldata$Clusters2=="C10"),]
C10gene<-gene[,which(colnames(gene)%in% rownames(C10))]
C10<-rowMeans(C10gene)

C11<-coldata[which(coldata$Clusters2=="C11"),]
C11gene<-gene[,which(colnames(gene)%in% rownames(C11))]
C11<-rowMeans(C11gene)

C12<-coldata[which(coldata$Clusters2=="C12"),]
C12gene<-gene[,which(colnames(gene)%in% rownames(C12))]
C12<-rowMeans(C12gene)

C13<-coldata[which(coldata$Clusters2=="C13"),]
C13gene<-gene[,which(colnames(gene)%in% rownames(C13))]
C13<-rowMeans(C13gene)

all<-cbind(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13)

tall<-t(all)
tdist=dist(tall,method="euclidean")
hc=hclust(tdist,method="complete")

plclust(hc)
plot(hc, hang=-1, xlab="");

######用EPC marker 注释newEC######


markerEPC<-c("Cd34","Flt1","Kdr","Prom1","Ptprc","Eng","Vcam1","Kit")
nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI-dim15-fibroblast-1.2-50000")


p <- plotEmbedding(
    ArchRProj = nproj, 
    colorBy = "GeneScoreMatrix", 
    name = markerEPC, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(nproj)
)

plotPDF(plotList = p, 
    name = "EC-EPCmarker-UMAP-Imputation.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerEPC,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-EPCMarker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)



artery<-c("Rbp7","Ly6a","8430408G22Rik","Id1","Stmn2","Fbln5","Glul","Cxcl12","Sox17","Hey1","Mgll","Gadd45g","Hspa1a","Dusp1","Alpl","Btg2","Ier2",
	"Mecom","Slc6a6","Klf4","Plk2","Fos","Tsc22d1","Crip1","Vegfc","Cdkn1a","Jun","Cyr61","Junb","Rgs5","Igfbp3","Aqp7","Amd1","Klf6","Tinagl1",
	"Azin1","Acvrl1","Jund","Fbln2","Sema3g","Adamts1","Sat1","Zfp36","Hist1h2bc","Ebf1","Nebl","Nrarp","Dll4","Gja5","Ptprr")

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = artery,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-arteryMarker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

capillary<-c("Cd36","Fabp4","Aqp1","Rgcc","Gpihbp1","Aplnr","Lpl","C1qtnf9","Sparcl1","Car4","Sparc","Tcf15","AW112010","Tspan13","Timp4","Sgk1",
"Kdr","Trp53i11","Slc28a2","Cd300lg","Cav1","Cyp4b1","Meox2","Cd200","Fscn1","Hspb1","Rsad2","Vwa1","Kcna5","Jam3","Dok4","Sox4",
"Dhrs3","Pdzd2","Tspan7","Gm12002","Fmo1","BC028528","Ifi203","Cav2","Ramp3","Ccdc85a","Plpp1","Myadm","Tuba1a","Arhgap18","Thrsp",
"Klk8","Cxcl9","Ssh2")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = capillary,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-capillaryMarker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)


vein<-c("Mgp","Cfh","Apoe","Cpe","Cytl1","Bgn","Dcn","Plvap","Ctsh","Vwf","Rbp1","Fabp5","Npr3","H19","Vcam1","Tmem108","Tm4sf1","Id2",
	"Tmem176b","Hmcn1","Ptgs1","Ccdc80","F2r","Igfbp4","Il6st","Ctla2a","Emcn","Cd59a","Clu","Ahsg","Spint2","Smoc1","Vim","Comp","Tmem176a",
	"Mgst1","Plagl1","Lum","Igf2","Gpr182","Colec11","Cxcl16","Tmem2","Tmem100","Cdh11","Ablim1","Bmx","Ramp2","Pdlim3","Ednrb")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = vein,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-veinMarker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)





lymphatic<-c("Ccl21a","Mmrn1","Fgl2","Prss23","Fxyd6","Thy1","Igfbp5","Fth1","Lcn2","Lyve1","Prelp","Nts","S100a6","Ifi27l2a","Cd9",
	"Pard6g","Lrg1","Flt4","Timp2","Lbp","Cp","Reln","Cd63","Tmsb10","Pdpn","Serpine2","Bcr","Nsg1","Maf","Scn1b","Marcks","Clca3a1",
	"Marcksl1","Lmo2","Timp3","Gpm6a","Stab1","Ntn1","Meox1","Arl4a","Anxa1","Cd24a","Abi3bp","Pglyrp1","Slco3a1","Rnase4","Nr2f2",
	"Tgfbr2","Ptn","Dapk2")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = lymphatic,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-lymphaticMarker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)


