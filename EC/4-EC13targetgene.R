library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome("mm10")
ArchR.EC<-loadArchRProject(path = "/md01/nieyg/scATAC-ArchR/ArchR-Endothelial_Cell/Save-Proj2-EC5-subtype")
############寻找EC3与EC1的差异peak#############
markersPeaks <- getMarkerFeatures(ArchRProj = ArchR.EC, useMatrix = "PeakMatrix", 
                                  groupBy = "ArchR.EC$Subtype", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.3 & Log2FC >= 0.1")
markerList
write.table(markerList$EC1[,c(8,1,3,4)],"EC1-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)
write.table(markerList$EC3[,c(8,1,3,4)],"EC3-FDR03-peak.txt",row.names=FALSE,col.names=FALSE)

sed -i 's/"//g' EC1-FDR05-peak.txt
awk '{print $1 "\t" $2 "\t" $3"\t" $4 "\t" "+"}' EC1-FDR03-peak.txt>EC1-FDR03-peak.bed

#####注释bed 文件##########
annotatePeaks.pl EC1-FDR03-peak.bed mm10 > EC1_peak.annotation.xls
annotatePeaks.pl EC3-FDR03-peak.bed mm10 > EC3_peak.annotation.xls
######寻找motif富集#####
findMotifsGenome.pl EC1-FDR03-peak.bed mm10 EC1-FDR03-peak_motif -len 6,8,10,12,13
findMotifsGenome.pl EC3-FDR03-peak.bed mm10 EC3-FDR03-peak_motif -len 6,8,10,12,13
#####寻找motif的结合位点##
findMotifsGenome.pl EC1-FDR03-peak.bed mm10 ./ -find ./EC1-FDR03-peak_motif/knownResults/known63.motif > ./EC1/Junb.txt;
findMotifsGenome.pl EC3-FDR03-peak.bed mm10 ./ -find ./EC3-FDR03-peak_motif/knownResults/known1.motif > ./EC3/Junb.txt;

findMotifsGenome.pl EC1-FDR03-peak.bed mm10 ./ -find ./EC1-FDR03-peak_motif/knownResults/known7.motif > ./EC1/Fosl2.txt;
findMotifsGenome.pl EC3-FDR03-peak.bed mm10 ./ -find ./EC3-FDR03-peak_motif/knownResults/known6.motif > ./EC3/Fosl2.txt;

findMotifsGenome.pl EC1-FDR03-peak.bed mm10 ./ -find ./EC1-FDR03-peak_motif/knownResults/known38.motif > ./EC1/Bach1.txt;
findMotifsGenome.pl EC3-FDR03-peak.bed mm10 ./ -find ./EC3-FDR03-peak_motif/knownResults/known43.motif > ./EC3/Bach1.txt;

findMotifsGenome.pl EC1-FDR03-peak.bed mm10 ./ -find ./EC1-FDR03-peak_motif/knownResults/known39.motif > ./EC1/Nfe2.txt;
findMotifsGenome.pl EC3-FDR03-peak.bed mm10 ./ -find ./EC3-FDR03-peak_motif/knownResults/known60.motif > ./EC3/Nfe2.txt;

findMotifsGenome.pl EC1-FDR03-peak.bed mm10 ./ -find ./EC1-FDR03-peak_motif/knownResults/known3.motif > ./EC1/Batf.txt;
findMotifsGenome.pl EC3-FDR03-peak.bed mm10 ./ -find ./EC3-FDR03-peak_motif/knownResults/known5.motif > ./EC3/Batf.txt;

#####注释结合位点peak，寻找靶基因#########
Junb
Bach1.txt
Nfe2.txt
Batf.txt

#!/bin/bash

awk '{print $1 }' Fosl2.txt |  sort -k 1 |uniq|grep -v "PositionID" | while read id;
do    read_id="^"$id$'\t';
    grep "$read_id" EC3-FDR03-peak.bed >> Fosl2-EC3.bed;done;



awk '{print $1 }' Batf.txt |  sort -k 1 |uniq|grep -v "PositionID" | while read id;
do    read_id="^"$id$'\t';
    grep "$read_id" EC1-FDR03-peak.bed >> Batf-EC1.bed;done;


annotatePeaks.pl Batf-EC1.bed mm10 > EC1_Batf.annotation.xls
annotatePeaks.pl Batf-EC3.bed mm10 > EC3_Batf.annotation.xls

annotatePeaks.pl Junb-EC1.bed mm10 > EC1_Junb.annotation.xls
annotatePeaks.pl Nfe2-EC1.bed mm10 > EC1_Nfe2.annotation.xls
annotatePeaks.pl Bach1-EC1.bed mm10 > EC1_Bach1.annotation.xls
annotatePeaks.pl Fosl2-EC1.bed mm10 > EC1_Fosl2.annotation.xls


annotatePeaks.pl Junb-EC3.bed mm10 > EC3_Junb.annotation.xls
annotatePeaks.pl Nfe2-EC3.bed mm10 > EC3_Nfe2.annotation.xls
annotatePeaks.pl Bach1-EC3.bed mm10 > EC3_Bach1.annotation.xls
annotatePeaks.pl Fosl2-EC3.bed mm10 > EC3_Fosl2.annotation.xls

#####注释靶基因的功能#########
library(clusterProfiler)
library(org.Mm.eg.db)
pdf("EC1-Batf-GO.pdf")
gene.df <- bitr(gene, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                readable = TRUE)

ego
barplot(ego, showCategory=20)

write.csv(ego,"EC3-Batf-GO-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Mm.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                readable = TRUE)

ego
barplot(ego, showCategory=20)


ggo <- groupGO(gene     = gene.df$ENTREZID,
               OrgDb    = org.Mm.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

head(ggo)


ggo<-as.data.frame(ggo)
ggo<-ggo[order(ggo$Count,decreasing=T),]
write.csv(ggo,"EC3-Nfe2-ggo.csv")

