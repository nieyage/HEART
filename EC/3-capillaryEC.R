library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
ArchR.EC<-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
ArchR.EC <- addImputeWeights(ArchR.EC,reducedDims="IterativeLSI")
##############除EC5之外都是capillary EC#########
capillary<-c("Cd36","Fabp4","Aqp1","Rgcc","Gpihbp1","Aplnr","Lpl","C1qtnf9","Sparcl1","Car4","Sparc","Tcf15","AW112010","Tspan13","Timp4","Sgk1",
"Kdr","Trp53i11","Slc28a2","Cd300lg","Cav1","Cyp4b1","Meox2","Cd200","Fscn1","Hspb1","Rsad2","Vwa1","Kcna5","Jam3","Dok4","Sox4",
"Dhrs3","Pdzd2","Tspan7","Fmo1","BC028528","Ifi203","Cav2","Ramp3","Ccdc85a","Plpp1","Myadm","Tuba1a","Arhgap18","Thrsp",
"Klk8","Cxcl9","Ssh2")

p <- plotEmbedding(
    ArchRProj = ArchR.EC, 
    colorBy = "GeneScoreMatrix", 
    name = capillary, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(ArchR.EC)
)
plotPDF(plotList = p, 
    name = "capillary-Genes-umap.pdf", 
    ArchRProj = ArchR.EC, 
    addDOC = FALSE, width = 5, height = 5)
