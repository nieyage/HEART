#######给全部cluster添加注释名，并根据细胞类型分类 #####
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")
nproj<-loadArchRProject(path = "Save-nProj")

######添加celltype信息###
cluster<-proj5$Clusters
table(proj5$Clusters)

for (i in 1:60420){
if(cluster[i]=="C1"){cluster[i]="Cardiomyocytes"}
if(cluster[i]=="C2"){cluster[i]="Fibroblast"}
if(cluster[i]=="C3"){cluster[i]="Fibroblast"}
if(cluster[i]=="C4"){cluster[i]="Fibroblast"}
if(cluster[i]=="C5"){cluster[i]="Fibroblast"}
if(cluster[i]=="C6"){cluster[i]="Fibroblast"}
if(cluster[i]=="C7"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C8"){cluster[i]="Smoothmuscle"}
if(cluster[i]=="C9"){cluster[i]="Fibroblast"}
if(cluster[i]=="C10"){cluster[i]="Fibroblast"}
if(cluster[i]=="C11"){cluster[i]="Fibroblast"}
if(cluster[i]=="C12"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C13"){cluster[i]="ImmuneCells"}
if(cluster[i]=="C14"){cluster[i]="Endothelial"}
if(cluster[i]=="C15"){cluster[i]="Endothelial"}
if(cluster[i]=="C16"){cluster[i]="Endothelial"}
if(cluster[i]=="C17"){cluster[i]="Endothelial"}
if(cluster[i]=="C18"){cluster[i]="Endothelial"}
if(cluster[i]=="C19"){cluster[i]="Endothelial"}
if(cluster[i]=="C20"){cluster[i]="Endothelial"}
if(cluster[i]=="C21"){cluster[i]="Endothelial"}
}
table(cluster)

cluster
           C10            C11            C14            C15             C9 
           813           2261           1523           2945            897 
Cardiomyocytes    Endothelial     Fibroblast    ImmuneCells   Smoothmuscle 
          2521          23587          17479           4285           4109 
cluster
Cardiomyocytes    Endothelial     Fibroblast    ImmuneCells   Smoothmuscle 
          2521          28055          21621           4114           4109
		  

proj5$celltype<-cluster
cellNames <- proj5$cellNames
proj5 <- addCellColData(ArchRProj = proj5, data = paste0(cluster),
    cells=cellNames, name = "celltype",force = TRUE)


#########add timepoint/treatment#####
sample<-proj5$Sample
table(proj5$Sample)
for (i in 1:60420){
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
proj5$timepoint<-sample
proj5 <- addCellColData(ArchRProj = proj5, data = paste0(sample),
    cells=cellNames, name = "timepoint",force = TRUE)

sample2<-proj5$Sample
table(proj5$Sample)
for (i in 1:60420){
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
proj5$treatment<-sample2
proj5 <- addCellColData(ArchRProj = proj5, data = paste0(sample2),
    cells=cellNames, name = "treatment",force = TRUE)
	table(proj5$treatment)

p1 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "celltype2", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj5, colorBy = "cellColData", name = "treatment", embedding = "UMAP")

ggAlignPlots(p1, p2, p3, type = "h")
plotPDF(p1,p2,p3, name = "New-UMAP-celltype2-timepoint-treatment.pdf", ArchRProj = proj5, addDOC = FALSE, width = 5, height = 5)

#####add marker gene######
####have other genes###
markerGenes  <- c(
"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1",####xinji
"Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
"Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk",###Macrophages
"Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat",###Tcell
"Cd79a","Cd79b","Mzb1","Ly6d",####Bcell
"Cd74","Cd83","Cd86","Flt3","Cd209a",####DC
"Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
"Tie1","Fabp4","Esam","Kdr","Tek","Eng",###endothelial
"Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
"Dclk1","Fkbp5","Akt3",##fibroblast
"Rgs5","Ano1","Acta2",####smoothmuscle
"Prr15","Slc9a3r1","Lrrn4","Slc39a8","Krt8","Anxa8",#####epicardial
"Pdgfrb","Abcc9",#Pericyte
"Pde3b","Ghr","Sik2","Pparg","Pde8b","Mgst1","Mapk10","Mme","Lama4",
"Fasn","Igf1","Dhrs3","Adra1a","Adh1b","Rtn4","Ptger3",#####Adipocyte
"Mbnl1","Fn1","Pde4d",####VSMC
"Scn7a","Bcl2",####Neuronal
"Ccnd3","Ptprc"#####Lymphocyte
)
C9:
	Fcer1a, Sertad4, Gm15867, Aifm2, Rufy2, Gm3285, Odf3l2, Efna2, Gamt, Mir6913, Mir1191b, Gli1, Pmel, Mir7239, Mrpl55, Rnf227, Ctdnep1, Mir451b, Hils1, Rundc3a
C10:
	Gsta3, Lyg2, Unc80, Olfr12, Cdh19, Ascl5, Mroh3, Gpr25, Rxrg, Trappc3l, Slc35f1, Lrrc75b, Slc6a15, Lgr5, Erbb3, Mir6917, Adcy1, Gm11985, Vstm2a, Slit3
C11:
	Fn1, Cyp27a1, Resp18, Myoc, 9030612E09Rik, Apc2, C1qtnf2, 1810065E05Rik, Smtnl2, Proca1, Mycbpap, Gm11627, 2810433D01Rik, Sstr2, Tpo, Mir673, Mir493, Mir3544, Mir434, Mir432
C14:
	Higd1b, Nrip2, Sox17, Mrpl15, Sgk3, Ppp1r42, Prex2, Gm17644, Sulf1, Slco5a1, Msc, Trpa1, Jph1, Pi15, Tfap2d, Lincmd1, Mir133b, Gsta3, Mir30c-2, Gm29669
C15:
	Cfh, Ptgs2os, Sele, Selp, Tbx19, Sft2d2, 1700007P06Rik, Ptpn14, Slc35d3, Pawr, 4921513I03Rik, Zpbp, Cldn7, Agr2, Slc25a21, Lrrc9, Smoc1, 1700018A04Rik, Foxf2, Serpinb9d
Cardiomyocytes:
	Jph1, Tfap2d, Mir30c-2, Bag2, Mir6896, Ankrd23, Mir6349, Fhl2, 1500015O10Rik, Ndufb3, Acadl, Erbb4, Tuba4a, Inha, Slc4a3, Irs1, Slc19a3, Ramp1, Klhl30, Ppip5k2
Endothelial:
	Sox17, Mrpl15, Prex2, Neurl3, Sema4c, Mgat4a, 4930594C11Rik, Mfsd6, Casp8, 9530026F06Rik, Nrp2, Klf7, Mir6899, Fam124b, Sag, St8sia4, 4930598F16Rik, Rnf152, Ralb, Platr22
Fibroblast:
	Ppp1r42, Gm17644, Sulf1, Slco5a1, Trpa1, Lincmd1, Mir133b, Gm29669, 1700001G17Rik, Aff3, Il1rl2, 4930521E06Rik, Col3a1, Col5a2, Aox1, Aox2, Nif3l1, Crygb, Gm8883, Tnp1
ImmuneCells:
	Sgk3, Bend6, Ptpn18, Prss39, Cfc1, Cnot11, Rnf149, Creg2, Gm16894, Il18r1, Il18rap, Mir7681, Cflar, Cryga, Slc11a1, 1700016L21Rik, Dock10, Pid1, Mir6353, Dner
Smoothmuscle:
	Igfbp5, Ptprn, Asic4, Krtap28-13, Gm19461, Lmod1, Rgs5, Rgs4, Mark1, 5033404E19Rik, 4930567K20Rik, Taar5, Taar6, Hey2, Mical1, Stk11, Jsrp1, Usp44, Mir331, Mir7211
	


markersGS <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "New-GeneScores-Heatmap-haveother", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
####save marker list to GO ananlysis####

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
write.csv(markerList,"allsample-markerList.csv")

####have 5 kinds genes###
markerGenes  <- c(
"Pdgfrb","Abcc9","Rgss","Fam65c","Lhfp","Scn3a","Tbx2","Pdzd2",#Pericyte
"Pde3b","Ghr","Sik2","Pparg","Pde8b","Mgst1","Mapk10","Mme","Lama4",
"Fasn","Igf1","Dhrs3","Adra1a","Adh1b","Rtn4","Ptger3","Slc1a3",#####Adipocyte
"Mbnl1","Fn1","Pde4d",####VSMC
"Scn7a","Bcl2","Nrxn1","Nrxn3","Negr1","Lgi4","Nkain3","Ncam2",####Neuronal
"Ccnd3","Ptprc","Aoah","Rhoh","Camk4","Sytl3","Skap1","Aim1","Ikzf3"#####Lymphocyte
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Heatmap-5-last", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)

###find marker peak###
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "PeakMatrix", 
    groupBy = "celltype2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "celltype2-Peak-Heatmap", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif")

####marker TF#####

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 8, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapEM, name = "celltype2-Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)
  
  if("Motif" %ni% names(proj5@peakAnnotation)){
    proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif")
}
proj5 <- addBgdPeaks(proj5)

proj5 <- addDeviationsMatrix(
  ArchRProj = proj5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj5, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotPDF(plotVarDev, name = "celltype2-Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj5, addDOC = FALSE)


