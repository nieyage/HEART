ecm<-c("Tecta","Thbs2","Col20a1","Ccbe1","Lum","Wnt4","Hpse2","Col24a1","Thsd4","Col3a1","2300002M23Rik",
"Ndp","Lama1","Col6a4","Col4a1","Thbs1","Lrrtm4","Lingo2","Col4a2","Svep1","Col11a1","Itih1","Spon1","Ahsg",
"Amelx","Tnr","Col28a1","Oc90","Acan","Wnt3","Egflam","Mmp19","Shh","Ltbp3","Ltbp4","Col1a1","Efemp2",
"Adamtsl4","Fn1","Col5a1","Loxl1","Col10a1","Gfod2","Rell2","Enam","Emilin1","Mmp23","Bgn","Wnt6","Mmp12",
"Bcan","Lrrc24","Adamtsl2","Vasn","Ltbp1","Lamc1","Mfap5","Col5a3","Col15a1","Fbn1","Smoc2","Tgfbr3","Optc",
"Crispld2","Matn2","Ntn1","Ccbe1","Ccdc80","Col8a1","Adamts2","Lrrc32","Tnc","Cthrc1","Tnxb","Alpl","Abi3bp","Tgfb1","Col6a3")

b=t(scale(t(b),scale = T,center = T))

p<-pheatmap(b, cluster_rows=F, 
               clustering_method="ward.D2",
               show_rownames=T,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)

plotPDF(plotList = p, 
    name = "FB-ecm-overlap.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 4, height = 8)



c<-pheatmap(v4, cluster_rows=T, 
               clustering_method="ward.D2",
               show_rownames=T,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)

rowlist4<-rownames(v4[c$tree_row[["order"]],])


mat=t(scale(t(mat),scale = T,center = T))

p<-pheatmap(mat, cluster_rows=T, 
               clustering_method="ward.D2",
               show_rownames=F,border_color = NA,
               cluster_cols=F,cutree_rows =1,
               show_colnames = T)

plotPDF(plotList = p, 
    name = "FB2special-peak.pdf", 
    ArchRProj = nproj, 
    addDOC = FALSE, width = 4, height = 8)











