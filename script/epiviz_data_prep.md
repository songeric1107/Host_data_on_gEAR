
title: "epiviz dataset_repare"

---
prepare bigwig files from narrow peak file

```{r setup, include=FALSE}
cut -f1,2 mm10.dna.toplevel.fa.fai > mm10.dna.toplevel.chromosome.sizes
cut -f1,2,3,5 *narrowpeak > myFile.bedgraph
sort -k1,1 -k2,2n myFile.bedgraph > myFile_sorted.bedgraph
/usr/local/packages/ucsc_utils/bin/bedGraphToBigWig myFile_sorted.bedgraph mm10.dna.toplevel.chromosome.sizes myBigWig.bw

```

prepare bigbed files from narrow peak file


```{r cars}

cut -f1,2 mm10.dna.toplevel.fa.fai > mm10.dna.toplevel.chromosome.sizes
cut -f1,2,3,5 *narrowpeak > myFile.bedgraph
cat file.bedGraph | awk '{print $1 "\t" $2 "\t" $3}' > file.bed
#If your bedGraph file contains a header line, you can get rid of it with an inverse grep:
grep -v track file.bedGraph | awk '{print $1 "\t" $2 "\t" $3}' > file.bed

bedToBigBed file.bed mm10.chrom.sizes regions.bb
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}





```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
