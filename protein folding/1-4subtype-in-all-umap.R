#######获取各个亚型的protein folding gene的FC情况######
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
proj5<-loadArchRProject(path = "Save-Proj5")
proj5 <- addImputeWeights(proj5,reducedDims ="IterativeLSI-dim8")
#######Hspa8等overlap gene在所有细胞类型上的分布状态#########
#########################################
p<-plotEmbedding(
    ArchRProj = proj5, 
    colorBy = "GeneScoreMatrix", 
    name = overlap, ######include Hspa8
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj5)
)

######九张放在一页图#######
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
pdf("/md01/nieyg/scATAC-ArchR/stressed-overlap-genelist-umap-9fig.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))

p2 <- lapply(p[10:18], function(x){
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
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))

p2 <- lapply(p[19:27], function(x){
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
do.call(cowplot::plot_grid, c(list(ncol = 3,nrow=3),p2))
dev.off();

#####所有PP基因的分布########
PP<-c("Map3k5","Atf4","Atf6","Bak1","Bcap31","Bax","Bcl2","Capn1","Capn2","Casp12","Stub1","Ddit3","Canx",
	"Calr","Calr3","Calr4","Cul1","Derl2","Derl3","Derl1","Plaa","Ubqln1","Ubqln2","Ubqln3","Ubqln4","Edem1",
	"Edem2","Edem3","Eif2ak1","Eif2ak2","Eif2ak4","Lman1","Lman1l","Man1a","Man1a2","Man1b1","Man1c1",
	"Pdia3","Ppp1r15a","Amfr","Herpud1","Syvn1","Dnaja1","Dnaja2","Dnajb1","Dnajb11","Dnajb12",
	"Dnajb14","Dnajb2","Dnajc1","Dnajc10","Dnajc3","Dnajc5","Dnajc5b","Dnajc5g","Hspa1a","Hspa1b","Hspa1l",
	"Hspa2","Hspa4l","Hspa5","Hspa8","Hsp90aa1","Hsp90ab1","Hsp90b1","Ern1","Mapk8","Map2k7","Nploc4","Nfe2l2",
	"Os9","Rrbp1","Vcp","Prkn","Pdia4","Pdia6","Eif2ak3","Ngly1","Rad23a","Rad23b","Rbx1","Rnf5",
	"Mbtps1","Mbtps2","Sar1a","Sar1b","Preb","Sec13","Sec23a","Sec23b","Sec24a","Sec24b","Sec24c","Sec24d",
	"Sec31a","Sec31b","Sec61a1","Sec61a2","Sec61b","Sec61g","Sec62","Sec63","Sel1l","Sel1l2","Svip",
	"Traf2","Tram1","Ssr1","Ssr2","Ssr3","Ssr4","Ube2j1","Ube2j2","Ube2g1","Ube2g2","Ube2d1","Ube2d2a",
	"Ube2d2b","Ube2d3","Ube2dnl1","Ube2dnl2","Ubxn1","Ubxn2a","Ubxn4","Ubxn6","Ubxn8","Ufd1","Ube4b",
	"Selenos","Lman2","Wfs1","Xbp1","Erlec1")
p<-plotEmbedding(
    ArchRProj = proj5, 
    colorBy = "GeneScoreMatrix", 
    name = PP, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj5)
)
plotPDF(plotList = p, 
    name = "PP-umap.pdf", 
    ArchRProj = proj5, 
    addDOC = FALSE, width = 5, height = 5)

#########4个subtype在所有细胞umap上的分布#########
FB3proj<-loadArchRProject(path = "Save-fibroblast2")
EC4proj<-loadArchRProject(path = "./ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
MP3proj<-loadArchRProject(path = "/md01/Chenzh275/ori/Data/scATAC/ArchR-Immune_Cell/Save-Proj2-MP3-4Subtype")
SMC1proj<-loadArchRProject(path="/md01/nieyg/scATAC-ArchR/Save-SMC-merge12")

######提取这四个cluster的barcode#######
cellNames <- EC4proj$cellNames
m<-EC4proj$Subtype
EC4proj <- addCellColData(ArchRProj = EC4proj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)

cellNames <- MP3proj$cellNames
m<-MP3proj$Subtype
MP3proj <- addCellColData(ArchRProj = MP3proj, data = paste0(m),
    cells=cellNames, name = "subtype",force = TRUE)

FB3barcode<-rownames(FB3proj@cellColData[which(FB3proj@cellColData$subtype=="FB3"),])
EC4barcode<-rownames(EC4proj@cellColData[which(EC4proj@cellColData$subtype=="EC4"),])
MP3barcode<-rownames(MP3proj@cellColData[which(MP3proj@cellColData$subtype=="MP3"),])
SMC1barcode<-rownames(SMC1proj@cellColData[which(SMC1proj@cellColData$subtype=="SMC1"),])


######free umap df#######
#####re draw umap 

df = data.frame(proj5@embeddings$UMAP$df, subtype=c(rep("other",60420)))
df$subtype<-as.character(df$subtype)
subtype<-df$subtype
m<-rownames(df)
for (i in 1:60420){
if(m[i] %in% FB3barcode)
    df[i,3]="FB3"
if(m[i] %in% MP3barcode)
    df[i,3]="MP3"
if(m[i] %in% EC4barcode)
    df[i,3]="EC4"
if(m[i] %in% SMC1barcode)
    df[i,3]="SMC1"    
    }
table(df$subtype)

pdf("four_sub_allumap.pdf")
colnames(df)<-c("UMAP_Dimension_1","UMAP_Dimension_2","subtype")
ggplot(df, aes(x = UMAP_Dimension_1, y = UMAP_Dimension_2, color = subtype)) +
  geom_point(size=0.8,alpha=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#0072B5FF", "#E18727FF", "#20854EFF","gray87","#7876B1FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Timepoint")
