#use r-3.6.0, this script is an example to prepare the uploading profiles from SingleCellExperiemnt rds files
library(Rtsne)
library(SingleCellExperiment)
library(CellTrails);library(scuttle)



BP_Data <- readRDS("sce.rds")
class(BP_Data)

#SCNorm

sce <- scuttle::logNormCounts(BP_Data )



exp=data.frame(logcounts(sce ))

ann=rowData(sce)
ann$gene_symbol=row.names(ann)

exp.ann=merge(ann,exp,by.y=0,by.x=0)

names(exp.ann)[1]="gene_symbol"
write.table(exp.ann[-1],file="expression.tab",sep="\t",row.names=F,quote=F)

write.table(exp.ann[c(2,1)],file="genes.tab",sep="\t",row.names=F,quote=F)


raw.obs=data.frame(names(exp))
names(raw.obs)="observations"

obs=data.frame(BP_Data@colData)



obs1=cbind(data.frame(row.names(obs)),obs[c(38:39)])
names(obs1)[1]="observations"
obs1$cell_type=obs1$CellTrails.state



reducedDimNames(BP_Data)

tsne.cor=data.frame(reducedDim(BP_Data,"CellTrails.tSNE"))

tsne.cor=reducedDim(BP_Data, "CellTrails.tSNE")[,1:2]

colnames(tsne.cor)=c("tSNE_1","tSNE_2")

obs2=merge(obs1,tsne.cor,by.x=0,by.y=0)

obs3=obs2[c(2,3,4,5,6,7)]

write.table(obs3=,file="observations.tab",sep="\t",row.names=F,quote=F)


