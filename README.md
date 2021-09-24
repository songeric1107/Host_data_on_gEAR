# Host_data_on_gEAR

https://umgear.org

gEAR has two functions for single-cell datasets analysis. gEAR could display your processed results as you published ; or re-analyze your raw matrix using our single-cell workbench, which is a scanpy-based platform.


If you want to re-analyze your matrix, 4 files are required for uploading: *.mtx(sparse matrix), barcode.tsv,genes.tsv, metadata.xls


If you want to display your processed results as you published, 4 different files are required: expression.tab, genes.tab,observations.tab, metadata.xls
you can also put 3 tab files as individual sheet of an Excel file.


This page provides a few examples on preparing the files from different platform(microarray, rnaseq,scRNAseq) for gEAR uploading using GEO dataset as example.

All the scripts are example scripts only, need to be slightly modified based on different datasets, especially for microarray.

We highly recommendated prepare 3 tab delimited files(MUST BE NAMED AS genes.tab, observations.tab,expression.tab) and compress them for gEAR submission, especially for scRNAseq. 

