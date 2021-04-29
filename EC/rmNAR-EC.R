

#######提取Endothelial-cluster #####

library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")

proj2<-loadArchRProject(path = "Save-rmNAR")

idxSample <- BiocGenerics::which(proj2$celltype %in% c("Endothelial"))
cellsSample <- proj2$cellNames[idxSample]
nproj=proj2[cellsSample, ]
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-rmNAR-Endothelial", load = T)
nproj<-loadArchRProject(path = "Save-rmNAR-Endothelial")
table(nproj$Sample)

nproj <- addIterativeLSI(
    ArchRProj = nproj,
    name = "IterativeLSI-rmNAR-Endothelial",
    useMatrix = "TileMatrix",
    iterations = 2, 
    clusterParams = list(resolution = c(1.2),
    sampleCells = 10000, n.start = 10), 
    varFeatures = 50000, 
    force = TRUE,
    dimsToUse = 1:15
)


####cluster###
nproj <- addClusters(
    input = nproj,
    reducedDims = "IterativeLSI-rmNAR-Endothelial",
    name = "Clusters2", force = TRUE, resolution = 1, dimsToUse = 1:15
)

nproj <- addUMAP(
    ArchRProj = nproj, 
    reducedDims = "IterativeLSI-rmNAR-Endothelial", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.05, 
    force = TRUE,
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")

ggAlignPlots(p1, p4, type = "h")
plotPDF(p1, p4, name = "rmNAR-Endothelial-Plot-UMAP.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "rmNAR-Endothelial-GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)



cluster<-nproj$Clusters2
table(nproj$Clusters2)

for (i in 1:24437){
if(cluster[i]=="C1"){cluster[i]="EC5"}
if(cluster[i]=="C2"){cluster[i]="EC4"}
if(cluster[i]=="C3"){cluster[i]="EC5"}
if(cluster[i]=="C4"){cluster[i]="EC5"}
if(cluster[i]=="C5"){cluster[i]="EC4"}
if(cluster[i]=="C6"){cluster[i]="EC1"}
if(cluster[i]=="C7"){cluster[i]="EC1"}
if(cluster[i]=="C8"){cluster[i]="EC1"}
if(cluster[i]=="C9"){cluster[i]="EC3"}
if(cluster[i]=="C10"){cluster[i]="EC1"}
if(cluster[i]=="C11"){cluster[i]="EC2"}
if(cluster[i]=="C12"){cluster[i]="EC2"}
if(cluster[i]=="C13"){cluster[i]="EC1"}
if(cluster[i]=="C14"){cluster[i]="EC1"}
if(cluster[i]=="C15"){cluster[i]="EC1"}
}
table(cluster)
nproj$subtype<-cluster
p3 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = nproj, colorBy = "cellColData", name = "subtype", embedding = "UMAP")
ggAlignPlots(p3, p4,p5, type = "h")
plotPDF(p3,p4,p5, name = "rmNAR-addsubtype-UMAP-EC.pdf", ArchRProj = nproj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = nproj, overwrite = F,outputDirectory = "Save-rmNAR-Endothelial", load = FALSE)
nproj<-loadArchRProject("Save-rmNAR-Endothelial")
table(nproj$subtype,nproj$Sample)
markersGS <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerGenes<-c("Cd36","Fabp4","Aqp1","Rgcc","Gpihbp1","Aplnr","Lpl","C1qtnf9","Sparcl1","Car4","Sparc","Tcf15","AW112010","Tspan13","Timp4","Sgk1",
"Kdr","Trp53i11","Slc28a2","Cd300lg","Cav1","Cyp4b1","Meox2","Cd200","Fscn1","Hspb1","Rsad2","Vwa1","Kcna5","Jam3","Dok4","Sox4",
"Dhrs3","Pdzd2","Tspan7","Gm12002","Fmo1","BC028528","Ifi203","Cav2","Ramp3","Ccdc85a","Plpp1","Myadm","Tuba1a","Arhgap18","Thrsp",
"Klk8","Cxcl9","Ssh2",#######C
"Mgp","Cfh","Apoe","Cpe","Cytl1","Bgn","Dcn","Plvap","Ctsh","Vwf","Rbp1","Fabp5","Npr3","H19","Vcam1","Tmem108","Tm4sf1","Id2",
  "Tmem176b","Hmcn1","Ptgs1","Ccdc80","F2r","Igfbp4","Il6st","Ctla2a","Emcn","Cd59a","Clu","Ahsg","Spint2","Smoc1","Vim","Comp","Tmem176a",
  "Mgst1","Plagl1","Lum","Igf2","Gpr182","Colec11","Cxcl16","Tmem2","Tmem100","Cdh11","Ablim1","Bmx","Ramp2","Pdlim3","Ednrb",######Vein
  "Rbp7","Ly6a","Id1","Stmn2","Fbln5","Glul","Cxcl12","Sox17","Hey1","Mgll","Gadd45g","Hspa1a","Dusp1","Alpl","Btg2","Ier2",
  "Mecom","Slc6a6","Klf4","Plk2","Fos","Tsc22d1","Crip1","Vegfc","Cdkn1a","Jun","Cyr61","Junb","Rgs5","Igfbp3","Aqp7","Amd1","Klf6","Tinagl1",
  "Azin1","Acvrl1","Jund","Fbln2","Sema3g","Adamts1","Sat1","Zfp36","Hist1h2bc","Ebf1","Nebl","Nrarp","Dll4","Gja5","Ptprr"#####A
  )

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-ACVMarker-Heatmap", width = 8, height = 6, ArchRProj =nproj, addDOC = FALSE)

markerGenes<-c("Mgp","Cfh","Apoe","Cpe","Cytl1","Bgn","Dcn","Plvap","Ctsh","Vwf","Rbp1","Fabp5","Npr3","H19","Vcam1","Tmem108","Tm4sf1","Id2",
  "Tmem176b","Hmcn1","Ptgs1","Ccdc80","F2r","Igfbp4","Il6st","Ctla2a","Emcn","Cd59a","Clu","Ahsg","Spint2","Smoc1","Vim","Comp","Tmem176a",
  "Mgst1","Plagl1","Lum","Igf2","Gpr182","Colec11","Cxcl16","Tmem2","Tmem100","Cdh11","Ablim1","Bmx","Ramp2","Pdlim3","Ednrb",######Vein EC5
  "Smad6","Litaf","Nfkbia","Noct","Abca1","Trib1","Nod2",#####EC1
  "Cd36","Dach1","Apln","Mc4r","Uchl1","Htr1a","Pyy",#####EC2
  "Tenm1","Il1rapl1","Pcdh15","Ptprd","Grid2","Cntn5",###EC3
  "Hspa1a","Hspa8","Dnajb14","Hspa1l","Hspa1b"##EC4
  )

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-ECMarker-Heatmap", width = 8, height = 6, ArchRProj =nproj, addDOC = FALSE)


markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
EC1<-markerList$EC1$name
EC2<-markerList$EC2$name
EC3<-markerList$EC3$name
EC4<-markerList$EC4$name
EC5<-markerList$EC5$name




library(clusterProfiler)
library(org.Mm.eg.db)

pdf("/md01/nieyg/scATAC-ArchR/Save-rmNAR-Endothelial/Plots/EC-GO.pdf")
ego_geneID <- bitr(EC4, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
###GO analysis####
ego <- enrichGO(ego_geneID$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"/md01/nieyg/scATAC-ArchR/Save-rmNAR-Endothelial/Plots/EC5-BP.csv")

      AR14  AR3  AR7  P14   P3   P7
  EC1 6005  203  752 4565  439  645
  EC2   62   67 2749   46 1255 1835
  EC3    1 1342   38   11   53    4
  EC4  373   43  815  400  214  429
  EC5  155  557  300  161  482  436


#######
pathToMacs2 <- findMacs2()
nproj <- addGroupCoverages(ArchRProj = nproj, groupBy = "subtype")

nproj <- addReproduciblePeakSet(
    ArchRProj = nproj, 
    groupBy = "subtype", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(nproj)
nproj <- addPeakMatrix(nproj)


markersPeaks <- getMarkerFeatures(
    ArchRProj = nproj, 
    useMatrix = "PeakMatrix", 
    groupBy = "subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "rmNAR-ENdothelial-Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)

nproj <- addMotifAnnotations(ArchRProj = nproj, force = TRUE,motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = nproj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "EC-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = nproj, addDOC = FALSE)
  
####去掉同一个基因家族的TF，仅保留一个，添加尽可能多的基因家族进来

  idxSample <- BiocGenerics::which(enrichMotifs@NAMES %in% c(
           "Sp3_802","Smad5_883","Sp2_150","Sp5_238","Sp6_191","Sp4_167","Bach2_119",
       "Tcf3_31","Tcf12_59","Ascl2_23","Sp7_222","Sp8_207","Sp9_231","Zfp219_816","Gata5_385","Gata4_386","Gata1_387","Mesp2_57",
       "Gata2_383","Sox9_725","Egr2_188","Egr3_183","Zfp161_209","Zic5_195","Tcfap2c_3","Gata3_384",
       "Zfp281_193","Zfp740_204","Plagl1_152","Tcfap2a_1","Zic2_225","LINE3878_878",
       "E2f6_264","Zic1_181","Fosb_98","Nfe2l2_101","Jund_135","Batf_790","Fosl1_107",
       "Nfe2l1_117","Nfe2l3_871","Jun_126","Zfp202_169","Tbpl2_773","Nfatc1_703","Nfatc3_702","Nfata4_857"))
cellsSample <- enrichMotifs@NAMES[-idxSample]
enrichMotifs1=enrichMotifs[cellsSample, ]


heatmapEM <- plotEnrichHeatmap(enrichMotifs1, n = 10, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "rmNAR-fibroblast-Motifs",width = 12, height = 6,ArchRProj = nproj, addDOC = FALSE)
