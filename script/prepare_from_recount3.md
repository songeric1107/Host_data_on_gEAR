#if (!requireNamespace("BiocManager", quietly = TRUE))

#BiocManager::install("recount3")




#########using recount###
library("recount3")

## Obtain all available projects
## Obtain all available projects
projects.all <- rbind(
  recount3::available_projects("human"),
  recount3::available_projects("mouse")
)




#proj_info1=subset(projects.all,organism=="human" & file_source == "gtex")

###subset nerve dataset

proj_info <- subset(
  projects.all,
  project == "NERVE" & project_type == "data_sources" )

# extract the dataset bsed on SRP ID
#proj_info <- subset(projects.all,project == "SRP009615" & project_type == "data_sources")




rse_gene_nerve <- create_rse(proj_info)

assay(rse_gene_nerve, "counts") <- transform_counts(rse_gene_nerve)

count=data.frame(assay(rse_gene_nerve, "counts") )

gene.raw=data.frame(rowData(rse_gene_nerve))

gene=gene.raw[grep("ENSEMBL",gene.raw$source),]



foo <- data.frame(do.call('rbind', strsplit(as.character( gene$gene_id),'.',fixed=TRUE)))

gene1=cbind(gene,foo)

gene2=gene1[,c(11,7)]

obs=data.frame(colData(rse_gene_nerve))

count.all=merge(gene2,count,by=0)
names(count.all)[2]="ensembl_id"


exp.agg=aggregate(count.all[-c(1:3)],by=list(count.all$ensembl_id,count.all$gene_name),FUN=mean)

names(exp.agg)[1:2]=c("ensembl_ID","gene_symbol")

write.table(exp.agg,"expression.tab",sep="\t",quote=F)
write.table(exp.agg[,c(1,2)],"expression.tab",sep="\t",quote=F)
write.table(obs,"observation.tab",sep="\t",quote=F)








