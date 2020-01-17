#Differential expression analysis R Pipeline
#v0.1- @J0MS

#Clear variables
rm(list=ls())

#Load packages
library('limma')
library('edgeR')

#Read file of samples
targets <- readTargets()
targets

#Read file of counts
x <- read.delim("total_counts.tsv", sep=' ', row.names=1, stringsAsFactors=FALSE)
head(x)
tail(x)

#Put the counts and other info into a DGEList object (A list-based S4 class for storing read counts and associated information from digital gene expression)
y <- DGEList(counts=x[,1:4], group=targets$Treatment)
colnames(y) <- targets$Label
dim(y)

# Filtering: We filter out very lowly expressed tags, keeping genes that are expressed at a reasonable level in at least one treatment condition. We keep genes that achieve at least one count per million (cpm) in at least two samples:
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep,]
dim(y)

#Re-compute the library sizes:
y$samples$lib.size <- colSums(y$counts)


#Normalizing
#Compute effective library sizes using TMM normalization:
#Robinson, M.D., Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11, R25 (2010) doi:10.1186/gb-2010-11-3-r25
#https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25
y <- calcNormFactors(y)
y$samples


#Data Exploration
#MDS plots shows distances in terms of biological coefficient of variation (BCV), between samples:
pdf('MDSplot_control-treatment.pdf')
plotMDS(y)
dev.off()

#Estimating dispersion
#The common dispersion estimates the overall BCV of the dataset, averaged over all genes:
y <- estimateCommonDisp(y, verbose=TRUE)

#Now estimate gene-specific dispersions:
y <- estimateTagwiseDisp(y)

#Plot the estimated dispersions:
pdf('BCVplot_control-treatment.pdf')
plotBCV(y)
dev.off()

#Differential Expression
#Compute exact genewise tests for differential expression between tumor and non-tumor treatments,  test for differential expression between two groups of count libraries (Robinson and Smyth (2008)) for a difference in mean between two groups of negative binomial random variables.
et <- exactTest(y)
top <- topTags(et)
top

#Check the individual cpm values for the top genes:
cpm(y)[rownames(top), ]

#The total number of DE genes at 5% FDR is given by
summary(de <- decideTestsDGE(et))

#Plot the log-fold-changes, highlighting the DE genes:
pdf('logFCplot_control-treatment.pdf')
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()

#To obtain genelists for all genes and FDR, you could try
tab<-topTags(et, n=nrow(y))
write.table(tab, file="FDR_control-treatment.txt", sep=" ")
