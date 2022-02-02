

## Prepare the files using the raw matrix from 10x platform

1.	you can locate the files from the 10x cellranger output folder under raw_feature_bc_matrix folder, there are usually 3 compressed files, barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz

	There are 3 columns in the features.tsv.gz, please only keep the first 2 columns
```{bash}

zcat features.tsv.gz |cut -f1,2 > genes.tsv
```


	Decompress all the the zip files
```{bash}
gunzip *.gz
```

 renames features.tsv to genes.tsv:
```{bash}
mv features.tsv genes.tsv
```


## there should be 3 files for next step (MUST BE NAMED AS FOLLOWING): barcodes.tsv genes.tsv matrix.mtx

Compress all the files
```{bash}
tar -cf upload1.tar *.tsv *.mtx (windows)

COPYFILE_DISABLE=1 tar -cf upload1.tar *.tsv *.mtx (mac)

```


2.	Prepare the files based on seurat processed output (all the command below are R script)

```{r cars}
suppressMessages(library("Seurat"));suppressMessages(library(dplyr)); suppressMessages(library(biomaRt))

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


#path to your RDS file of Seurat analysis result

path="~/demo.RDS"
data <- readRDS(path)


##create metadata sheet using function above
obs <- metadata(data)

#add replicate column based on cell_type
obs <- obs %>% group_by(cell_type) %>% mutate(replicate = row_number())

##add cluster_label to able the primary analysis function
obs$cluster_label=obs$cell_type

##if you have the hexcode for each cell type to color, name that column as "cell_type_colors"

## Make sure the cluster names are not in numeric format
# In order to enable the “primary analysis” function of single cell workbench, one of the column names MUST be in this list. ['cluster', 'cell_type', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot']. As long as one of the tag is found in this list, the primary analysis will use that tag to display the cluster to compare with.

# Don’t put “,” or “  “ in the column name, otherwise, the multi-genes display function will have problems. 


########create count matrix

count <- expression_matrix(data, ensembl.dedup)

#names(count)=gsub("X","",names(count))

#if the “-“ is changed to “.” by R, using this script to replace “-“ to “.”
colnames(count)=gsub("[.]","-",colnames(count))



####prepare annotation

mart = useMart( 'ensembl' )
datasets <- listDatasets(mart)
mart = useDataset( 'mmusculus_gene_ensembl' , mart = mart )
ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') , mart=mart)
names(ensembl)[2] = "gene_symbol"
##human reference
##mart = useDataset( 'hsapiens_gene_ensembl' , mart = mart )


names(ensembl) = c("ensembl_ID","gene_symbol")
ensembl.dedup=ensembl[!duplicated(ensembl$gene_symbol),]



###############create gene matrix
genes_to_use=genes_creation(data, ensembl.dedup)

write.table(obs, "observations.tab", sep = "\t", quote =  FALSE, row.names = FALSE)
 write.table(genes_to_use, "genes.tab", sep = "\t", quote =  FALSE, row.names = FALSE)
 write.table(count, "expression.tab", sep = "\t", quote =  FALSE, row.names = FALSE, col.names = TRUE)
 
 ####zip file together for uploading
 
 system( 'tar -czvf upload.tar.gz *.tab')
 
```
3.	Prepare the files from SingleCellExperiment object 



```{r}

library(Rtsne);library(SingleCellExperiment);library(CellTrails);library(scuttle)
BP_Data <- readRDS("example.rds")
class(BP_Data)

#SCNorm, get the normalized count matrix

data <- scuttle::logNormCounts(BP_Data )

exp=data.frame(logcounts(data ))

ann=rowData(data)
ann$gene_symbol=row.names(ann)

exp.ann=merge(ann,exp,by.y=0,by.x=0)

names(exp.ann)[1]="gene_symbol"
write.table(exp.ann[-1],file="expression.tab",sep="\t",row.names=F,quote=F)

write.table(exp.ann[c(2,1)],file="genes.tab",sep="\t",row.names=F,quote=F)

raw.obs=data.frame(names(exp))
names(raw.obs)="observations"

obs=data.frame(BP_Data@colData)

obs1=cbind(data.frame(row.names(obs)),obs[,])
names(obs1)[1]="observations"
obs1$cell_type=obs1$CellTrails.state
reducedDimNames(BP_Data)

tsne.cor=data.frame(reducedDim(BP_Data, "CellTrails.tSNE")[,1:2])

colnames(tsne.cor)=c("tSNE_1","tSNE_2")

obs2=merge(obs1,tsne.cor,by.x=0,by.y=0)

#rearrange the columns
obs3=obs2[c(2,3,4,5,6,7)]
write.table(obs3=,file="observations.tab",sep="\t",row.names=F,quote=F)

####compress files together for uploading
 
 system( 'tar -czvf upload.tar.gz *.tab')


```


