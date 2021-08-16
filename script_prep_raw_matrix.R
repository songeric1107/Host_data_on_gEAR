library(GEOquery)
library(biomaRt)
library(tidyverse)
library("dplyr")
library(Matrix)

###prepare count file ##################
setwd("~/Downloads")

gse_num = 'GSE98969'
outfile_prefix = paste('ID', gse_num, sep = "")
data <- getGEO(gse_num, GSEMatrix=TRUE)



###load supplementary file
sfiles = getGEOSuppFiles(gse_num,fetch_files=T)
fnames = rownames(sfiles)
fnames
# Check the fnames to idenitfy the index of the file and replace the index
countsFileName = fnames[1]

###get the list of content
supple=untar(countsFileName,list=TRUE)


###unzip the zipped folder into gse_num folder
untar(countsFileName,exdir=gse_num )




setwd(gse_num)
path=getwd()

file.names <- dir(pattern ="tsv.gz")


###combine each individual dataset together
out.file<-""
for(i in 1:length(file.names)){
  
  file <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE,row.names=1)
  #names(file)[ncol(file)]=paste(file.names[i],names(file)[ncol(file)],sep="_")
  if (nrow(file)!=0) {
    file$gene=row.names(file)
   
    out.file <- cbind(out.file, file)}}



##get expression matrix
count = out.file[,!grepl("gene", colnames(out.file))]

###get gene annotation.
id= (out.file[,grepl("gene", colnames(out.file))])

names(id)=rep("gene", ncol(id))

id.uniq=data.frame(id[, !duplicated(colnames(id))])
names(id.uniq)="gene"

count.m=count[-1]
row.names(count.m)=id.uniq$gene


###prepare the ensembl annotation db

mart = useMart( 'ensembl' )
datasets <- listDatasets(mart)
mart = useDataset( 'mmusculus_gene_ensembl' , mart = mart )
ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') , mart=mart)
names(ensembl)[2] = "gene_symbol"


gene=ensembl[!duplicated(ensembl$gene_symbol),]



###merge to raw expression matrix

raw.ann=merge(gene,count.m,by.x="gene_symbol",by.y=0)

raw_count1=raw.ann[-c(1:2)]

row.names(raw_count1)=raw.ann$gene_symbol

sparse.matrix=as.sparse(raw_count1)


writeMM(obj = sparse.matrix, file="matrix.mtx")



barcode=data.frame(colnames(raw_count1))

genes=raw.ann[,c(2,1)]

write.table(barcode,"barcodes.tsv",sep="\t",row.names=F,col.names=F,quote=F)

write.table(genes,"genes.tsv",sep="\t",row.names=F,col.names=F,quote=F)

###zip the 3 files for uploading

system('COPYFILE_DISABLE=1 tar -cf upload1.tar *.tsv *.mtx')


##########################################################################




