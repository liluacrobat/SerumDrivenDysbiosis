## Clear workspace to avoid conflicts with existing objects
rm(list=ls())
## Load required packages
require(R.matlab)
library( "DESeq2" )
library(ggplot2)


## Read DESeq2 input data generated in MATLAB
## dtmtx: count matrix
## dtgroup: group labels corresponding to samples
D<-readMat("../data/GO_DESeq2_5vs50.mat");

## Extract raw count matrix
countData <- D$dtmtx

## Construct sample metadata
## Group information is converted to a factor as required by DESeq2

metaData <-DataFrame(factor(D$dtgroup))
colnames(metaData)<-'group'


## Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData=round(countData),
                              colData=metaData, 
                              design=~group)

## Estimate size factors and obtain normalized counts
ddsSF <- estimateSizeFactors(dds)
normalized_counts <- counts(ddsSF, normalized=TRUE)

## Run DESeq2 differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

## Extract commonly used result metrics
baseMean<-res$baseMean
log2FoldChange <- res$log2FoldChange
lfcSE <- res$lfcSE
stat <- res$stat
pvalue <- res$pvalue
padj <- res$padj
