####Prepare uploadint matrix based on cloupe files from cellranger pipeline.




	1.Identify the TSNE/UMAP coordinate from cellranger output (sample_ID/ANALYSIS/20220201_cellranger.6.1.2_6451/sampleID_ILxxxxxxxxx/outs/analysis/umap/), or you can export projection from Loupe Broswer version > 6.0

	2.Identify the generate cluster name from "graphclust" folder(sample_ID/ANALYSIS/20220201_cellranger.6.1.2_6451/sampleID_ILxxxxxxxxx/outs/analysis/graphclust) or  export the cluster name from cloupe by click the 3 dots on the right top corner of "Graph-Based" window,or 

	3.create the observation file with the 1st column named as "observations" and include umap/tsne coordinate. 

	4.also include a column named as "cluster_lable" for primary analysis.

	

####Prepare the expression matrix:




```{r}
	library(Seurat)
	
	pathm1="filtered_feature_bc_matrix/"

	data1 <- Read10X(data.dir = pathm1)

	data1 <- CreateSeuratObject(counts = data1 , project = "patient1",min.cells = 0, min.features = 0)
	data1 <- NormalizeData(object = data1)

	count=data.frame(GetAssayData(object = data1))

	ann=read.table("pathm1/features.tsv.gz",sep="\t",header=T)

	count.ann=merge(ann, count,by.y=0,by.x="V2")

	write.table(count.ann[-1],"expression.tab",sep="\t",quote=F,row.names=F)

	write.table(count.ann[c(2,1)],"genes.tab",sep="\t",quote=F,row.names=F)

```

