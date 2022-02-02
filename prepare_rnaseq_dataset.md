 ####load in count file##########
 countsFileName="count.txt"
 
count = read.delim(sep="\t",countsFileName, header=T, check.names=F)
count
names(count)[1] = "gene"
rownames(count) = count$gene


######generate the observation file##################

# Retrieve the column metadata
options(repr.matrix.max.rows=50, repr.matrix.max.cols=200)

meta = read.delim("phenotype.txt")


#re-name the columns
names(obs) = c("observations","genotype", "gender")


#keep replicate name the same
#pipe the dataframe to the group_by function and send output to mutate function to count row numbers and assign them to the variable replicate
obs <- obs %>% group_by(genotype, gender) %>% mutate(replicate = row_number())
obs = data.frame(obs)



write.table(obs, file=paste(outfile_prefix, "observation.tab", sep = ""), sep="\t", row.names=F, quote=F)


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




###block1: use the annotation from the count file directly



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


write.table(count.ann[-1], file=paste(outfile_prefix, "_DataMTX.tab", sep = ""), sep="\t", quote=F, row.names=F)

genes_to_use=count.ann[c(2,1)]
genes_to_use
write.table(genes_to_use,file=paste(outfile_prefix, "genes.tab", sep = ""), sep="\t", quote=F,row.names=F)


count.ann#if the count file is raw count, basic normalization need to be completed before uploading
count_to_norm=count.ann[-c(1:2)]

row.names(count_to_norm)=count.ann[,1]

cpm=sweep(count_to_norm,2,colSums(count_to_norm),FUN="/")*10^6


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


write.table(nq, file=paste(outfile_prefix, "expression.tab", sep = ""), sep="\t", quote=F, row.names=T)

