suppressMessages(library("Seurat"));library(dplyr)



#Idents(data) <-"data.set"


metadata <- function(object,group){
  meta <- as.data.frame(cbind(row.names(object@meta.data), object@meta.data, Embeddings(object[["umap"]])))
  colnames(meta)[1] <- "observations"
 
  
  return(meta)
}


expression_matrix <- function(obj,ensembl){
  
  t1.m=data.frame(GetAssayData(object = obj))
  
  t1.m.ann=merge(ensembl,t1.m,by.y=0,by.x="gene_symbol")
  
  t3=t1.m.ann[grep("EN",t1.m.ann$ensembl_ID),]
  t4=t3[-c(1:2)]
  t4.s=t4[order(names(t4))]
  t5=cbind(t3[,2],t4.s)
  names(t5)[1]=c("ensembl_ID")
  
  return(t5)
  
}

genes_creation <- function(obj,ensembl){
  
  t1.m=data.frame(GetAssayData(object = obj))
  
  t1.m.ann=merge(ensembl,t1.m,by.y=0,by.x="gene_symbol")
  
  t3=t1.m.ann[grep("EN",t1.m.ann$ensembl_ID),]
  t4=t3[c(2,1)]
  names(t4)=c("ensembl_ID","gene_symbol")
  
  return(t4)
  
}

###zb fish not use
#mart = useMart( 'ensembl' )
#datasets <- listDatasets(mart)
#mart = useDataset( 'drerio_gene_ensembl' , mart = mart )
#mart = useDataset( 'mmusculus_gene_ensembl' , mart = mart )
#ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') , mart=mart)
#names(ensembl) = c("ensembl_ID","gene_symbol")
#head(ensembl)


#########################

#genes=read.delim(file="/local/projects-t3/XRHUM2017/gear_dataset/gear_dataset_upload/process/batch11/mouse_ensembl_v99.txt",sep="\t",header=T,fill=T)
mart = useMart( 'ensembl' )
datasets <- listDatasets(mart)
mart = useDataset( 'mmusculus_gene_ensembl' , mart = mart )
ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') , mart=mart)
names(ensembl)[2] = "gene_symbol"


names(ensembl) = c("ensembl_ID","gene_symbol")
ensembl.dedup=ensembl[!duplicated(ensembl$gene_symbol),]







path="/local/projects-t3/XRHUM2017/gear_dataset/gear_dataset_upload/process/batch11/SGNs/sgn_withoutunknown.RDS"
data <- readRDS(path)


##create metadata sheet using function above
obs <- metadata(data)

#add replicate column based on cell_type
obs <- obs %>% group_by(cell_type) %>% mutate(replicate = row_number())

##add cluster_label to able the primary analysis function
obs$cluster_label=obs$cell_type

##if you have the hexcode for each cell type to color, name that column as "cell_type_colors"



#########################create count matrix

count <- expression_matrix(data, ensembl.dedup)

#names(count)=gsub("X","",names(count))

colnames(count)=gsub("[.]","-",colnames(count))


###############create gene matrix
genes_to_use=genes_creation(data, ensembl.dedup)






 
 write.table(obs, "observations.tab", sep = "\t", quote =  FALSE, row.names = FALSE)
 write.table(genes_to_use, "genes.tab", sep = "\t", quote =  FALSE, row.names = FALSE)
 write.table(count, "expression.tab", sep = "\t", quote =  FALSE, row.names = FALSE, col.names = TRUE)
 
 
 
 ####zip file together for uploading
 
 system( 'tar -czvf upload.tar.gz *.tab')
 
 