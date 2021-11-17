getwd()
R.version

library(GEOquery)
library(biomaRt)
library(tidyverse)
library("dplyr")

# Retrieve the expression matrix from GEO
# Change the GSE number in red
gse_num = 'GSE102562'
outfile_prefix = paste('ID', gse_num, sep = "")
data <- getGEO(gse_num, GSEMatrix=TRUE)


###prepare count file ##################

###load supplementary file
sfiles = getGEOSuppFiles(gse_num,fetch_files=T)
fnames = rownames(sfiles)
fnames


# Check the fnames to idenitfy the index of the file and replace the index
countsFileName = fnames[1]


 ####load in count file##########
count = read.delim(sep="\t",countsFileName, header=T, check.names=F)
count
names(count)[1] = "gene"
rownames(count) = count$gene
head(count)
dim(count)

######generate the observation file##################

# Retrieve the column metadata
options(repr.matrix.max.rows=50, repr.matrix.max.cols=200)

meta = pData(phenoData(data[[1]]))
head(meta,10)
glimpse(meta)


#view the content of meta
meta1=data.frame((names(meta)))
#meta1



# check the column names of metadata file.
meta1 <- meta1 %>% mutate(id = row_number())
meta1

#Check if column has more than one variable 
#No spaces otherwise in quotes

#unique(meta$'title')


foo2 = data.frame(meta[,1])
foo2
colnames(foo2)=c("title")
foo2$x1=sub("SOD1:TREM2-","",foo2$title)
foo2



foo2$x1=gsub(" ","_",foo2$x1)
foo2
foo3 <- data.frame(do.call('rbind', strsplit(as.character(foo2[,2]),('_'))))
foo3

#print(class(meta[,1]))
#foo <- data.frame(do.call('rbind', strsplit(as.character(meta[,1]),('-|_| '))))
#foo
#delete space between sample and 1
#foo$X1=sub(" ","",foo$X1)
#foo

# select the columns you are interested in
obs =cbind(meta[1], foo3[c(1)], meta[41])
obs


#re-name the columns
names(obs) = c("observations","genotype", "gender")
obs



##another way to create the replicate column
#group_by(column names) should be what is in names(obs)
#keep replicate name the same
#pipe the dataframe to the group_by function and send output to mutate function to count row numbers and assign them to the variable replicate
obs <- obs %>% group_by(genotype, gender) %>% mutate(replicate = row_number())
obs = data.frame(obs)
obs


############To add replicate column##############################
##for one column
#obs$replicate <- ave( 1:nrow(obs), obs[,2], factor (obs[,2]), FUN=function(x) 1:length(x)) 
##below is for 2 columns
#obs$replicate <- ave( 1:nrow(obs), paste(obs[,2],obs[,3], sep="_"), factor( paste(obs[,2],obs[,3], sep="_")), FUN=function(x) 1:length(x) )
    
#obs$observations=foo$X1
#obs


write.table(obs, file=paste(outfile_prefix, "_COLmeta.tab", sep = ""), sep="\t", row.names=F, quote=F)


#names(count)

###check whether there is ensemble in the count file
temp=count[grepl("ENS",count[,1]),]
temp
names(count)
if (nrow(temp)>0) {
  print('ensembl ID are detected')
 

} else {
  
 print('Error: no ensemble detected')

}

temp1=count[grepl("gene|Gene|ID",names(count[,1])),]

#temp1=count[grepl("gene",count[,1])|grepl("Gene",count[,1]),]
temp1

if (nrow(temp1)>0) {
  print('gene symbols are detected')
 

} else {
  
 print('no gene symbols are detected')

}


###if ensemble ID detection and gene symbol detection are true, continue to block1

###if ensemble ID detection or gene sybmol dection is not true, continue to block2


dim(count)

count[1]


###block1: use the annotation from the count file directly

foo <- data.frame(do.call('rbind', strsplit(as.character(as.matrix(count[,1])),'.',fixed=TRUE)))


#foo <- data.frame(do.call('rbind', strsplit(as.character(count[，1]),'.',fixed=TRUE)))
foo
count[,1]=foo$X1
count


count_to_use = count[, -c(2:5)] ## delete duplicate
#colnames(count_to_use)[1] = 'ensembl_id' ## rename first column containing ensembl id's
names(count_to_use)[1]="ensembl_ID"

count_to_use$ensembl_ID=count_to_use$ensembl_ID

head(count_to_use)
write.table(count_to_use, file=paste(outfile_prefix, "_DataMTX.tab", sep = ""), sep="\t", quote=F, row.names=F)

genes_to_use=count[c(1,5)]

write.table(genes_to_use,file=paste(outfile_prefix, "_ROWmeta.tab", sep = ""), sep="\t", quote=F,row.names=F)


#count[1:3,1:3]

####block2: annotate the count file
# get ensemble annotation_info. Change the database to point to the right database name.
# Here we are using the mouse annotation database
mart = useMart( 'ensembl' )
datasets <- listDatasets(mart)
mart = useDataset( 'mmusculus_gene_ensembl' , mart = mart )
ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') , mart=mart)
names(ensembl)[2] = "gene_symbol"
head(ensembl)
####only keep the first ensemble ID for each gene
ensembl.m=ensembl[!duplicated(ensembl$gene_symbol),]

#merge by gene_symbol or ensembl ID # out the correct one below
#count.ann = merge(ensembl, count, by.x="gene_symbol", by.y="gene_symbol")

##merge by ensemble id           
##Y= column name in matrix
count.ann = merge(ensembl.m, count, by.x="gene_symbol", by.y="gene")
count.ann

#names(count.ann) = c("gene_symbol", "ensembl_ID")
names(count.ann)[c(2,1)] = c("ensembl_ID", "gene_symbol")

##remove the row with NA in ensembl or gene symbol

#count.ann=na.omit(count.ann)

##remove gen_symbol column for the final output
##if normalisation require comment out the write command below
write.table(count.ann[-1], file=paste(outfile_prefix, "_DataMTX.tab", sep = ""), sep="\t", quote=F, row.names=F)
dim(count.ann)
genes_to_use=count.ann[c(2,1)]
genes_to_use
write.table(genes_to_use,file=paste(outfile_prefix, "_ROWmeta.tab", sep = ""), sep="\t", quote=F,row.names=F)


count.ann#if the count file is raw count, basic normalization need to be completed before uploading
count_to_norm=count.ann[-c(1:2)]
count_to_norm
row.names(count_to_norm)=count.ann[,1]
count_to_norm #count=count_to_norm
cpm=sweep(count_to_norm,2,colSums(count_to_norm),FUN="/")*10^6
#cpm=cbind(count.ann[1:2],cpm)
cpm
#cpm=sweep(count[-2], 2, colSums(count[-1]), FUN='/')*10^6
#row.names(cpm)=count[1]


##quantile normalization of cpm

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
  
    nq=quantile_normalisation(data.frame(cpm))
nq
### you can use the script in block2 to annotate the "nq" files with ensembl_ID or gene_symbol
##row.names=T true to output the row names
write.table(nq, file=paste(outfile_prefix, "_DataMTX.tab", sep = ""), sep="\t", quote=F, row.names=T)



