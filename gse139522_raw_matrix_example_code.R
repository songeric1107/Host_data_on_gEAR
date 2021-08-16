#setwd("/Users/ysong/Downloads/GSE139522_RAW/")
#path="/Users/ysong/Downloads/GSE139522_RAW/"

setwd("~/Download")
path="~/Download"

library("ggplot2");library(GEOquery)
suppressMessages(library("monocle"))
suppressMessages(library("qlcMatrix"))

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("Matrix"))
suppressMessages(library("MAST"))
#sessionInfo()
#packageVersion("Seurat")


gse_num = 'GSE139522'
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


###unzip the zipped folder into GSM  folder
untar(countsFileName,exdir=gse_num )

setwd(gse_num)
path=getwd()
path1=list.dirs(path = path, full.names = TRUE, recursive = TRUE)[2]

out.file<-""
# list all files in the directory
list=as.matrix(list.files(pattern = ".gz", recursive = TRUE))
foo <- as.character(unique(data.frame(do.call('rbind', strsplit(as.character(list),'_',fixed=TRUE)))[1])$X1)


subfolder_names <- foo

for (j in seq_along(subfolder_names)){
  folder<-dir.create(paste0(subfolder_names[j]))}


#create directory and move files belong to each GSE number to individual folder.

library(filesstrings)
move_files(supple[1:3], foo[1], overwrite=TRUE)

move_files(supple[4:6], foo[2], overwrite=TRUE)
move_files(supple[7:9], foo[3], overwrite=TRUE)
move_files(supple[10:12], foo[4], overwrite=TRUE)

###load each dataset
pathm1=paste(path,gse_num,foo[1],sep="/")
pathm2=paste(path,gse_num,foo[2],sep="/")
pathm3=paste(path,gse_num,foo[3],sep="/")
pathm4=paste(path,gse_num,foo[4],sep="/")

#file.rename(list.files(path=pathm1,pattern = as.character(foo[1])), str_replace(list.files(path=pathm1),pattern = as.character(foo[1]), ""))
##check each folder,make sure there are 3 files named as barcodes.tsv.gz ,features.tsv.gz,matrix.mtx.gz,the files without the prefix name
##for filename in *; do echo mv \"$filename\" \"${filename//GSM4142870_Patient1_/}\"; done > mv.pbs


data1 <- Read10X(data.dir = pathm1)
data1 <- CreateSeuratObject(counts = data1 , project = "patient1",min.cells = 0, min.features = 0)

###load dataset2
data2 <- Read10X(data.dir = foo[2])
data2 <- CreateSeuratObject(counts = data2 , project = "patient2",min.cells = 0, min.features = 0)

###load dataset3
data3 <- Read10X(data.dir = foo[3])
data3 <- CreateSeuratObject(counts = data3 , project = "patient3",min.cells = 0, min.features = 0)

###load dataset4
data4 <- Read10X(data.dir = foo[4])
data4 <- CreateSeuratObject(counts = data4 , project = "patient4",min.cells = 0, min.features = 0)

####merge all the datasets together
data.combined <- merge(data1, y = c(data2, data3,data4), add.cell.ids = c("pateint1","patient2", "patient3", "patient4"),  project = "GSE139522")

library(Matrix)


#saveRDS(data.combined,file="combine.rds")



raw_count <- as.matrix(GetAssayData(data.combined, slot = "counts"))
mart = useMart( 'ensembl' )
datasets <- listDatasets(mart)
mart = useDataset( 'mmusculus_gene_ensembl' , mart = mart )
ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') , mart=mart)
names(ensembl)[2] = "gene_symbol"



gene=ensembl[!duplicated(ensembl$gene_symbol),]
raw.ann=merge(gene,count.m,by.x="gene_symbol",by.y=0)

raw_count1=raw.ann[-c(1:2)]

row.names(raw_count1)=raw.ann$gene_symbol

sparse.matrix=as.sparse(raw_count1)

library(Matrix)
                         

writeMM(obj = sparse.matrix, file="matrix.mtx")

barcode=data.frame(colnames(raw_count1))

genes=raw.ann[,c(2,1)]

write.table(barcode,"barcodes.tsv",sep="\t",row.names=F,col.names=F,quote=F)

write.table(genes,"genes.tsv",sep="\t",row.names=F,col.names=F,quote=F)
