Prepare 3 tab format for RNAseq dataset uploading


 1. load in count file
	 countsFileName="count.txt"
 
		count = read.delim(sep="\t",countsFileName, header=T, check.names=F)

		names(count)[1] = "gene"
		rownames(count) = count$gene


2. generate the observation file


		options(repr.matrix.max.rows=50, repr.matrix.max.cols=200)

		meta = read.delim("phenotype.txt")


			#re-name the columns
			names(obs) = c("observations","genotype", "gender")


#keep replicate name the same
#pipe the dataframe to the group_by function and send output to mutate function to count row numbers and assign them to the variable replicate

		obs <- obs %>% group_by(genotype, gender) %>% mutate(replicate = row_number())
		obs = data.frame(obs)
		write.table(obs, file= "observation.tab", sep="\t", row.names=F, quote=F)




3. annotaiton. check whether there is ensemble in the count file

			temp=count[grepl("ENS",count[,1]),]
			
			if (nrow(temp)>0) {
 			 print('ensembl ID are detected')

					} else {
  			print('Error: no ensemble detected')}


 Here we are using the mouse annotation database
 
	                library("ensembldb");library(biomaRt)
 
			mart = useMart( 'ensembl' )
			datasets <- listDatasets(mart)
			mart = useDataset( 'mmusculus_gene_ensembl' , mart = mart )
			ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') , mart=mart)
			names(ensembl)[2] = "gene_symbol"
			head(ensembl)

	#only keep the first ensemble ID for each gene

			ensembl.m=ensembl[!duplicated(ensembl$gene_symbol),]

	#merge by gene_symbol or ensembl ID 

			count.ann = merge(ensembl, count, by.x="gene_symbol", by.y="gene_symbol")

			names(count.ann)[c(2,1)] = c("ensembl_ID", "gene_symbol")


			write.table(count.ann[-1], file= "expression.tab", sep="\t", quote=F, row.names=F)

			genes_to_use=count.ann[c(2,1)]
	
			write.table(genes_to_use,file= "genes.tab", sep="\t", quote=F,row.names=F)


if the count file is raw count, basic normalization need to be completed before uploading

			count_to_norm=count.ann[-c(1:2)]

			row.names(count_to_norm)=count.ann[,1]

			cpm=sweep(count_to_norm,2,colSums(count_to_norm),FUN="/")*10^6


##quantile normalization of cpm

		quantile_normalisation <- function(df){
  		df_rank <- apply(df,2,rank,ties.method="min")
  		df_sorted <- data.frame(apply(df, 2, sort))
  		df_mean <- apply(df_sorted, 1, mean)
  		index_to_mean <- function(my_index, my_mean){
    	return(my_mean[my_index])}
  		df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  		rownames(df_final) <- rownames(df)
  		return(df_final)}
   	
		nq=quantile_normalisation(data.frame(cpm))
		write.table(nq, file=paste(outfile_prefix, "expression.tab", sep = ""), sep="\t", quote=F, row.names=T)


####compress files together for uploading
 
 	system( 'tar -czvf upload.tar.gz *.tab')

