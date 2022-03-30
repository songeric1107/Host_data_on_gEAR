#########install packages###################################


##https://www.bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("GEOquery")

#BiocManager::install("BiocGenerics")

#BiocManager::install("biomaRt")





library("BiocManager")
library("Biobase")

library("GEOquery")
library("biomaRt")




# load series and platform data from GEO
geoid="GSE55656"

gset <- getGEO(geoid, GSEMatrix =TRUE, AnnotGPL=T)


eData <- gset[[1]]

#metadata

meta <- data.frame(pData(eData))

obs=meta[c(0,8)]

#foo <- data.frame(do.call('rbind', strsplit(as.character(obs$title),',',fixed=TRUE)))
#obs1=cbind(obs,foo$X1)

names(obs)[1]="type"


#replace the redundant information

obs$type=gsub("mouse MVN cell type ","",obs$type)

###add replicate####
obs$replicate <- ave( 1:nrow(obs), obs$type, factor( obs$type), FUN=function(x) 1:length(x) )
obs$title=row.names(obs)
names(obs)[ncol(obs)]="observations"
obs1=obs[c(3,1:2)]



##the plateform for GSE69644 is GPL13667
platform_id=unique(meta$platform_id)
if (length(gset) > 1) idx <- grep(platform_id, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

expression=data.frame(exprs(gset))




####################get genes annotation###########################

gene=fData(gset)

gene.use= gene[c(0,3)]
#genes.use$ensembl_ID <- sapply(strsplit(as.character(genes.use$Ensembl),'///'), "[", 1)

gene.use$gene_symbol=sapply(strsplit(as.character(gene.use$`Gene symbol`),'///'), "[", 1)

gene.use$probe_id=row.names(gene.use)


###get gene annotaiton version#######

gene.version=unique(gene$`Annotation Date`)

foo <- data.frame(do.call('rbind', strsplit(as.character(gene.version[1]),',',fixed=TRUE)))
foo1 <- data.frame(do.call('rbind', strsplit(as.character(foo$X1),' ',fixed=TRUE)))

new.name=paste(foo1$X1,foo$X2,sep="")


version = data.frame(listEnsemblArchives())

gene.version=unique(gene$`Annotation Date`)
foo <- data.frame(do.call('rbind', strsplit(as.character(gene.version[1]),',',fixed=TRUE)))
foo1 <- data.frame(do.call('rbind', strsplit(as.character(foo$X1),' ',fixed=TRUE)))

new.name=paste(foo1$X1,foo$X2,sep="")


version = data.frame(listEnsemblArchives())

###access to the biomart database 

mart=useMart("ensembl",host=version[grep(new.name, version$date),3], dataset = "mmusculus_gene_ensembl")

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl", 
                   version = version[grep(new.name, version$date),4])

ensembl = getBM( attributes = c('ensembl_gene_id','external_gene_name') ,values=gene.use$`Gene Symbol`, mart=mart)

names(ensembl) = c("ensembl_ID","gene_symbol")
ensembl=ensembl[!duplicated(ensembl$gene_symbol),]

######annotate microarray expression matrix

genes.ann = merge(gene.use,ensembl,  by.y="gene_symbol", by.x="gene_symbol")
genes.ann=genes.ann[c(3:4,1)]


expression.ann=merge(genes.ann,expression,by.y=0,by.x="probe_id")


expr.f=expression.ann[which(expression.ann$ensembl_ID!=" "|expression.ann$gene_symbol!="---"|expression.ann$ensembl_ID!="---"),]



###aggregate the duplicate values for same genes

exp.agg=aggregate(expr.f[-c(1:3)],by=list(expr.f$ensembl_ID,expr.f$gene_symbol),FUN=mean)

exp.s=exp.agg[which(exp.agg$Group.2!="---"|exp.agg$Group.2!=""),]

exp.final=exp.s[grep("EN",exp.s$Group.1),]

names(exp.final)[1:2]=c("ensembl_id","gene_symbol")

paste(geoid, "_EXPmeta.json", sep = "")

write.table(obs1,file = paste("ID",geoid, "_COLmeta.tab", sep = ""),sep="\t",quote=F,row.names=F)

write.table(exp.final[-2],file=paste("ID",geoid, "_DataMTX.ta", sep = ""),sep="\t",quote=F,row.names=F)
write.table(exp.final[1:2],file = paste("ID",geoid, "_ROWmeta.tab", sep = ""),sep="\t",quote=F,row.names=F)
