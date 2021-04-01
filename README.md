# HEART
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))#####小鼠：mmusculus_gene_ensembl
my_ensembl_gene_id<-counts$target_id
head(my_ensembl_gene_id)
m<-gsub("_.*","",my_ensembl_gene_id)
ENSEMBL <- gsub("\\.\\d*", "", m) 
head(ENSEMBL)

head(mart)
options(timeout = 4000000)
#my_ensembl_gene_id<-row.names(raw_count_filt)
mms_symbols<- getBM(attributes=c('ensembl_transcript_id','external_gene_name',"description"),
                    filters = 'ensembl_transcript_id', values= ENSEMBL,useCache=F, mart = mart)
head(listAttributes(mart))
#readcount<-read.csv(file="raw_count_filt.csv",header = TRUE)
head(mms_symbols)
head(counts)
counts$gene<-rownames(counts)
counts$target_id<-ENSEMBL
mms_symbols<-mms_symbols[counts$target_id,]
counts<-counts %>% distinct(target_id,.keep_all=T)
counts<-counts[-1,]
geneid<-counts$target_id
mms_symbols<-mms_symbols[geneid,]
last<-cbind(mms_symbols,counts)
write.table(last,"lastdata.txt")