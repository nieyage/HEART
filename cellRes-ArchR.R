library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 
addArchRGenome("mm10")

inputFiles <- c("GSE153479_fragments.tsv.gz")
names(inputFiles)<-c("sample")

inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,   ####fragment files or bam files 
  sampleNames = names(inputFiles),
  filterTSS = 6, #Dont set this too high because you can always increase later
  filterFrags = 2800, 
  addTileMat = TRUE,
  force=TRUE,
  addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)




proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "downstream",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
getAvailableMatrices(proj)
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-cellres-data", load = FALSE)

