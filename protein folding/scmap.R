########scmap########
library(SingleCellExperiment)
library(scmap)
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
all <-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/Save-allcelltype" )
all<-addImputeWeights(all,reducedDims="IterativeLSI-dim8")

######extract cell-gene-matrix

#######all subtypes########
FBproj<-loadArchRProject(path = "Save-fibroblast2")
ECproj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
MPproj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")
SMCproj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")

######提取这四个cluster的barcode#######
cellNames <- ECproj$cellNames
m<-ECproj$Subtype
ECproj <- addCellColData(ArchRProj = ECproj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)

cellNames <- MPproj$cellNames
m<-MPproj$Subtype
MPproj <- addCellColData(ArchRProj = MPproj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)


FB1barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB1"),])
FB2barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB2"),])
FB3barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB3"),])
FB4barcode<-rownames(FBproj@cellColData[which(FBproj@cellColData$subtype=="FB4"),])

EC1barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC1"),])
EC2barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC2"),])
EC3barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC3"),])
EC4barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC4"),])
EC5barcode<-rownames(ECproj@cellColData[which(ECproj@cellColData$subtype=="EC5"),])

MP1barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP1"),])
MP2barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP2"),])
MP3barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP3"),])
MP4barcode<-rownames(MPproj@cellColData[which(MPproj@cellColData$subtype=="MP4"),])

SMC1barcode<-rownames(SMCproj@cellColData[which(SMCproj@cellColData$subtype=="SMC1"),])
SMC2barcode<-rownames(SMCproj@cellColData[which(SMCproj@cellColData$subtype=="SMC2"),])
SMC3barcode<-rownames(SMCproj@cellColData[which(SMCproj@cellColData$subtype=="SMC3"),])



########No CMs####
cellsSample <- c(FB1barcode,FB2barcode,FB3barcode,FB4barcode,
	EC1barcode,EC2barcode,EC3barcode,EC4barcode,EC5barcode,
	MP1barcode,MP2barcode,MP3barcode,MP4barcode,
	SMC1barcode,SMC2barcode,SMC3barcode)
nproj=all[cellsSample, ]
nproj<-addImputeWeights(nproj,reducedDims="IterativeLSI-dim8")

subtype<-c(rep("other",56362))
m<-rownames(nproj@cellColData)
for (i in 1:56362){
if(m[i] %in% FB3barcode)
    subtype[i]="FB3"
if(m[i] %in% MP3barcode)
    subtype[i]="MP3"
if(m[i] %in% EC4barcode)
    subtype[i]="EC4"
if(m[i] %in% SMC1barcode)
    subtype[i]="SMC1"    
if(m[i] %in% FB1barcode)
    subtype[i]="FB1"
if(m[i] %in% MP1barcode)
    subtype[i]="MP1"
if(m[i] %in% EC1barcode)
    subtype[i]="EC1"
if(m[i] %in% FB2barcode)
    subtype[i]="FB2"
if(m[i] %in% MP2barcode)
    subtype[i]="MP2"
if(m[i] %in% EC2barcode)
    subtype[i]="EC2"
if(m[i] %in% SMC2barcode)
    subtype[i]="SMC2" 
if(m[i] %in% EC3barcode)
    subtype[i]="EC3"
if(m[i] %in% SMC3barcode)
    subtype[i]="SMC3"    
if(m[i] %in% FB4barcode)
    subtype[i]="FB4"
if(m[i] %in% MP4barcode)
    subtype[i]="MP4"
if(m[i] %in% EC5barcode)
    subtype[i]="EC5"
    }
table(subtype)

cellNames <- nproj$cellNames
nproj <- addCellColData(ArchRProj = nproj, data = paste0(subtype),
    cells=cellNames, name = "subtype",force = TRUE)
GSM<-getMatrixFromProject(ArchRProj =nproj, useMatrix = "GeneScoreMatrix")
head(GSM@assays@data$GeneScoreMatrix)[,1:10]
genename<-rowData(GSM)$name ####共24333个基因
count<-GSM@assays@data$GeneScoreMatrix[,which(colnames(GSM) %in% rownames(nproj))]
count<-as.matrix(count)
rownames(count) = genename
anno<-nproj$subtype
anno<-as.data.frame(anno)
rownames(anno)<-rownames(nproj@cellColData)
count<-as.data.frame(count)
sce <- SingleCellExperiment(assays = list(counts = as.matrix(count)), colData = anno)
logcounts(sce) <- log2(counts(sce) + 1)
# use gene names as feature symbols
rowData(sce)$feature_symbol <- rownames(sce)
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
# remove features with duplicated names
sce <- sce[!duplicated(rownames(sce)), ]
sce
sce <- selectFeatures(sce,  n_features = 15000,suppress_plot = FALSE)

## Warning in linearModel(object, n_features): Your object does not contain
## counts() slot. Dropouts were calculated using logcounts() slot...

table(rowData(sce)$scmap_features)


sce <- indexCluster(sce,cluster_col="anno")
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))
head(metadata(sce)$scmap_cluster_index)
scmapCluster(projection = NULL, index_list = NULL, threshold = 0.7)
mydata <- yan[apply(yan,1,function(x) sum(x > 0) > 1),] 


scmapCluster_results <- scmapCluster(
  projection = sce, 
  index_list = list(
    mydata=mydata
  ),threshold = 0.5
)
getSankey(reference, clusters, plot_width = 400, plot_height = 600,
       colors = NULL)
reference: reference clustering labels

clusters: clustering labels under investigations


p<-plot(
  getSankey(
    colData(sce)$anno, 
    scmapCluster_results$scmap_cluster_labs[,'mydata'],
    plot_height = 400
  )
)




p1<-  p(
    colData(sce)$anno, 
    scmapCluster_results$scmap_cluster_labs[,'mydata'],
    plot_height = 400
  )

plot(
  getSankey(
    colData(sce)$cell_type1, 
    scmapCluster_results$scmap_cluster_labs[,'yan'],
    plot_height = 400
  )
)
pdf("test-sankey.pdf")
ggplot(as.data.frame(p1),
       aes(y = `# of cells`, axis1 = From, axis2 = To)) +
geom_alluvium(aes(fill = From), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") 


pdf("test-sankey.pdf")
  ggplot(data = p1,
       aes(axis1 = From, axis2 = To,
           weight = `# of cells`)) +
  scale_x_discrete(limits = c("From", "To")) +
  geom_alluvium(aes(fill = From)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  #theme_minimal() +
  #theme_tropical()
  ggtitle("Taxonomy abundance in each group")
