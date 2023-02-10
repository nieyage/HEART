library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
nproj<-loadArchRProject(path = "Save-Proj2-EC5-subtype")
nproj <- addImputeWeights(nproj,reducedDims="IterativeLSI")
markerEPC<-c("Cd34","Flt1","Kdr","Prom1","Ptprc","Eng","Vcam1","Kit")

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
    ArchRProj = ArchR.EC, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "ArchR.EC$Subtype",
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

ArchR.EC <- loadArchRProject(path = "Save-Proj2-EC5-subtype")

artery<-c("Rbp7","Ly6a","Id1","Stmn2","Fbln5","Glul","Cxcl12","Sox17","Hey1","Mgll","Gadd45g","Hspa1a","Dusp1","Alpl","Btg2","Ier2",
	"Mecom","Slc6a6","Klf4","Plk2","Fos","Tsc22d1","Crip1","Vegfc","Cdkn1a","Jun","Cyr61","Junb","Rgs5","Igfbp3","Aqp7","Amd1","Klf6","Tinagl1",
	"Azin1","Acvrl1","Jund","Fbln2","Sema3g","Adamts1","Sat1","Zfp36","Hist1h2bc","Ebf1","Nebl","Nrarp","Dll4","Gja5","Ptprr")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = artery,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-arteryMarker-Heatmap-1", width = 8, height = 6, ArchRProj =ArchR.EC, addDOC = FALSE)

 ArchR.EC <- addImputeWeights( ArchR.EC)

p <- plotEmbedding(
    ArchRProj = ArchR.EC, 
    colorBy = "GeneScoreMatrix", 
    name = artery, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(ArchR.EC)
)

plotPDF(plotList = p, 
    name = "EC-artery-UMAP-Imputation.pdf", 
    ArchRProj = ArchR.EC, 
    addDOC = FALSE, width = 5, height = 5)
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
plotPDF(heatmapGS, name = "Endothelial-capillaryMarker-Heatmap-1", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)
p <- plotEmbedding(
    ArchRProj = ArchR.EC, 
    colorBy = "GeneScoreMatrix", 
    name = capillary, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(ArchR.EC)
)

plotPDF(plotList = p, 
    name = "EC-capillary-UMAP-Imputation.pdf", 
    ArchRProj = ArchR.EC, 
    addDOC = FALSE, width = 5, height = 5)



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
plotPDF(heatmapGS, name = "Endothelial-veinMarker-Heatmap-1", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)


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
plotPDF(heatmapGS, name = "Endothelial-lymphaticMarker-Heatmap-1", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)
Angiogenic<-c("Tmsb10","Ift122","Sparc","Cldn5","Col4a2","Apln","Vim","Igfbp7","Vwa1","Mest","Sparcl1","Adm","Aplnr","Meox1","Trp53i11","Prnp",
"Nrp2","Mycn","Tnfaip8l1","Ccnd1","Gnb4","Sox4","Dapk2","Bcl6b","Pgf","Kcna5","Col15a1","Bst2","Nudt4","Nid2","Actb","Ets2","Cyp4b1","Pcdh17",
"Akap13","Rgs16","Prcp","Stmn1","Knop1","Jam3","Lgals1","Itga6","Hmgb2","Tubb6","Plxnd1","Arl4c","Noct","Lamc1","Kit","Hcls1")


heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = Angiogenic,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-AngiogenicMarker-Heatmap-1", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)

Interferon<-c("Rgcc","Lpl","Isg15","Car4","Ifit3","Rtp4","Ifi203","Rsad2","Tcf15","Ifit1","Ly6a","Mndal","Fabp4","Iigp1","Aqp1","Ifit3b",
	"AW112010","Timp4","Gpihbp1","Gbp7","Cdkn1c","Sgk1","Xdh","Irf7","Ifit2","Cmpk2","Slc28a2","Hspb1","Tgtp2","Gbp2","Stat1","Oasl2",
	"Lpar6","Gm12002","Pdzd2","Irgm1","Gbp4","Efnb1","Usp18","Igtp","Tgtp1","Parp14","Ifi44","Cmtm8","Trim12a","Ifi47","Dok4","Thrsp","Irf9","Cd274")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = Interferon,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-InterferonMarker-Heatmap-1", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)

Capillaryvenous<-c("Vcam1","Ier3","Pltp","Pi16","Fmo2","Fmo1","Calcrl","Emp1","AU021092","Bsg","Nfkbia","Eln","Bhlhe40","Icam1","Klk8","Socs3",
	"Slfn2","2200002D01Rik","Rnd1","Emp2","Arl4d","Rasa4","Ier5","Atf3","Ramp3","Hs3st1","Actn1","Rhob","Lhx6","Lmo1","Ankrd37","Col13a1","Csf1",
	"Dll1","Hist1h1c","Lfng","Egr1","Kcnb1","Nuak1","Prpf40b","Bmp4","Socs2","Gm12216","Entpd1","Tob1","Icosl","Fam71a","P4ha2","Mgat4a","Cyp26b1")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = Capillaryvenous,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "Endothelial-CapillaryvenousMarker-Heatmap-1", width = 8, height = 6, ArchRProj = ArchR.EC, addDOC = FALSE)



