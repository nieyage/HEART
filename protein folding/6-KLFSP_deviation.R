library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-allcelltype-rescue2")
proj5 <- addImputeWeights(proj5,reducedDims="IterativeLSI-dim8")
####Remove CMs##########
 B_cell Cardiomyocytes    Endothelial     Epicardial     Fibroblast 
           334           2521          28055            405          20831 
         Glial    granulocyte    Macrophages       Pericyte   Smoothmuscle 
           241            214           3367            144           4109 
        T_cell 
           199 


proj6<-proj5[which(proj5$subtypename %in%c( "B_cell","T_cell","Glial","Endothelial","Epicardial",
	"Fibroblast","granulocyte",
"Macrophages","Pericyte","Smoothmuscle" )),]
    proj6 <- addMotifAnnotations(ArchRProj = proj6, motifSet = "cisbp", name = "Motif")

proj6<-addPeakMatrix(proj6, force=TRUE)
proj6 <- addBgdPeaks(proj6,force = TRUE)
proj6 <- addDeviationsMatrix(
  ArchRProj = proj6, 
  peakAnnotation = "Motif",
  force = TRUE
)
getAvailableMatrices(proj6)
m<-getFeatures(proj6,"MotifMatrix")
TF<-c(m[grep("deviations:Klf.*",m)],m[grep("deviations:Sp.*",m)])

pdf("KLF-SP-deviations.pdf")
p <- plotEmbedding(
    ArchRProj = proj6, 
    colorBy = "MotifMatrix", 
    name = TF, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj6)
)

dev.off()








