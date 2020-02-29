#-------IMPORTING RAW MICROARRAY DATA INTO R---------
library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(affy)
library(oligo)
library(huex10sttranscriptcluster.db)

# getting the decompressed GEO series 
gse <- getGEO('GSE143150', GSEMatrix = FALSE)
# getting the name of the samples in the dataset 
names(GSMList(gse)) 
gsm <- GSMList(gse)[[1]]
gsm

#-------CREATING & ATTACHING A PHENOTDATA OBJECT-------
print(Meta(gsm)[['characteristics_ch1']])
# peek into an example dataset
for (gsm in GSMList(gse)) {
  print(Meta(gsm)[['characteristics_ch1']])
}

# getting only genotype and treatments for phenodata
genotype <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
}
treatment <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][3]
}

pd <- data.frame(genotype=as.factor(sapply(GSMList(gse),genotype)), 
                 treatment=as.factor(sapply(GSMList(gse),treatment)))
pd

# cleaning up the table for neater data
pd$genotype <- as.factor(pd$genotype)
levels(pd$genotype) <- c("KO","WT")
pd$treatment <- as.factor(pd$treatment)
levels(pd$treatment) <- c("DM","FM")
pd

# creating the phenodata
files <- list.files("./data", pattern = ".CEL.gz")
files
celfiles <- paste0('data/',files)
celfiles
data <- read.celfiles(celfiles,phenoData = new("AnnotatedDataFrame",pd))
data
phenoData(data)

#-------MICROARRAY DATA PROCESSING (RMA)-------
eset <- oligo::rma(data)
# shows gene expression levels, after rma
exprs(eset) 
# generate 12 MAplots from oligo package, each against pseudo-median reference
oligo::MAplot(exprs(eset))
# generate plot density from affy package
plotDensity(exprs(eset),xlab='log intensity',main="feature level densities after RMA",lwd=2)

#-------IDENTIFY DIFFERENTIALLY EXPRESSED GENES USING LINEAR MODELS-------
pData(eset)

f <- paste(eset$genotype, eset$treatment, sep="_")
f <- factor(f)
f
m <- model.matrix(~0+f)
colnames(m) <- levels(f)
colnames(m)
m
fit <- lmFit(eset,m)
fit <- eBayes(fit)

# show constructed table
topTable(fit)

contrasts <- makeContrasts(Contrast_1=KO_DM-KO_FM, 
                           Contrast_2=WT_DM-WT_FM,
                           Contrast_3=KO_DM-WT_DM,
                           Contrast_4=KO_FM-WT_FM,
                           interaction = (KO_DM-WT_DM) - (KO_FM-WT_FM),
                           genotype = (KO_DM-WT_DM) + (KO_FM-WT_FM),
                           media = (WT_DM-WT_FM) + (KO_DM-KO_FM),
                           levels=c("KO_DM","KO_FM","WT_DM","WT_FM"))
contrasts

# fit the contrasts to the table
fitted.contrast <- contrasts.fit(fit,contrasts)
fitted.ebayes <- eBayes(fitted.contrast)

#-----------FROM FEATURES TO ANNOTATED GENE LISTS---------

#show constructed table 
topTable(fitted.ebayes, number = 10, p.value = 0.05, lfc = 2)
probesets <- rownames(topTable(fitted.ebayes))
probesets

interesting_genes_up_all <- rownames(topTable(fitted.ebayes, number = Inf, p.value = 0.05, lfc = 2))
interesting_genes_up_all
detail_all <- AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_up_all,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
detail_all

# JMJD2B_KO = JMJD2B_KO_DM - JMJD2B_KO_FM
topTable(fitted.ebayes, coef=1, number = Inf, p.value = 0.05, lfc = 2)
probesets <- rownames(topTable(fitted.ebayes, coef=1))
probesets

interesting_genes_up_1 <- rownames(topTable(fitted.ebayes, coef=1, number = Inf, p.value = 0.05, lfc = 2))
interesting_genes_up_1
detail_1 <- AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_up_1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
detail_1

# wildtype = wildtype_DM - wildtype_FM
topTable(fitted.ebayes, coef=2, number = Inf, p.value = 0.05, lfc = 1)
probesets_2 <- rownames(topTable(fitted.ebayes, coef=2))
probesets_2

interesting_genes_up_2 <- rownames(topTable(fitted.ebayes, coef=2, number = Inf, p.value = 0.05, lfc = 2))
interesting_genes_up_2
detail_2 <- AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_up_2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
detail_2

# DM = JMJD2B_KO_DM - wildtype_DM
topTable(fitted.ebayes, coef=3)
probesets_3 <- rownames(topTable(fitted.ebayes, coef=3))
probesets_3

# FM = JMJD2B_KO_FM - wildtype_FM
topTable(fitted.ebayes, coef=4)
probesets_4 <- rownames(topTable(fitted.ebayes, coef=4))
probesets_4

# interaction = (JMJD2B_KO_DM-wildtype_DM) - (JMJD2B_KO_FM-wildtype_FM)
topTable(fitted.ebayes, coef=5)
probesets_5 <- rownames(topTable(fitted.ebayes, coef=5))
probesets_5

# genotype = (JMJD2B_KO_DM-wildtype_DM) + (JMJD2B_KO_FM - wildtype_FM)
topTable(fitted.ebayes, coef=6)
probesets_6 <- rownames(topTable(fitted.ebayes, coef=6))
probesets_6

# media = (wildtype_DM-wildtype_FM) + (JMJD2B_KO_DM-JMJD2B_KO_FM)
topTable(fitted.ebayes, coef=7,number = Inf, p.value = 0.05, lfc = 2)
probesets_7 <- rownames(topTable(fitted.ebayes, coef=7))
probesets_7

#----------PRIMARY ANALYSIS-------------

# retrieve genes that are differentially expressed as interesting genes, looking at coefficient 1
interesting_genes_1 <- topTable(fitted.ebayes, coef=1, number = Inf, p.value = 0.05, lfc = 2)
interesting_genes_1
# plot volcano plot for analysis
volcanoplot(fitted.ebayes, coef = 1, main=sprintf("%d features pass our cutoffs", nrow(interesting_genes_1)))
points(interesting_genes_1[['logFC']], -log10(interesting_genes_1[['P.Value']]), col='red')
# to extract symbols, gene identifiers and gene names from interesting genes
interesting_genes_up_1 <- rownames(interesting_genes_1[interesting_genes_1$logFC > 2,])
interesting_genes_down_1 <- rownames(interesting_genes_1[interesting_genes_1$logFC < -2,])
interesting_genes_up_1
interesting_genes_down_1
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_up_1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_down_1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

# retrieve genes that are differentially expressed as interesting genes, looking at coefficient 2
interesting_genes_2 <- topTable(fitted.ebayes, coef=2, number = Inf, p.value = 0.05, lfc = 2)
interesting_genes_2
# plot volcano plot for analysis
volcanoplot(fitted.ebayes, coef = 2, main=sprintf("%d features pass our cutoffs", nrow(interesting_genes_2)))
points(interesting_genes_2[['logFC']], -log10(interesting_genes_2[['P.Value']]), col='red')
# to extract symbols, gene identifiers and gene names from interesting genes
interesting_genes_up_2 <- rownames(interesting_genes_2[interesting_genes_2$logFC > 2,])
interesting_genes_down_2 <- rownames(interesting_genes_2[interesting_genes_2$logFC < -2,])
interesting_genes_up_2
interesting_genes_down_2
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_up_2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_down_2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

# retrieve genes that are differentially expressed as interesting genes, looking at coefficient 7
interesting_genes_7 <- topTable(fitted.ebayes, coef=7, number = Inf, p.value = 0.01, lfc = 3)
interesting_genes_7_a <- rownames(topTable(fitted.ebayes, coef=7, number = Inf, p.value = 0.01, lfc = 3))
# plot volcano plot for analysis
volcanoplot(fitted.ebayes, coef = 7, main=sprintf("%d features pass our cutoffs", nrow(interesting_genes_up_7)))
points(interesting_genes_7[['logFC']], -log10(interesting_genes_7[['P.Value']]), col='red')
# to extract symbols, gene identifiers and gene names from interesting genes
interesting_genes_up_7 <- (interesting_genes_7[interesting_genes_7$logFC >= 3,])
interesting_genes_down_7 <- (interesting_genes_7[interesting_genes_7$logFC <= -3,])
interesting_genes_up_7
interesting_genes_down_7
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_up_7,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_down_7,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

#----------SECONDARY ANALYSIS-------------

# imposing strict bounds on genes with high changes of expression
interesting_genes_up_7 <- rownames(interesting_genes_7[interesting_genes_7$logFC > 3,])
interesting_genes_down_7 <- rownames(interesting_genes_7[interesting_genes_7$logFC < -3,])
interesting_genes_up_7
interesting_genes_down_7
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_up_7,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_down_7,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")


#----------DOWNSTREAM ANALYSIS OF MICROARRAY DATA-------------

# make new (and smaller) eset with our interesing_genes
eset_of_interest <- eset[rownames(interesting_genes_7),]
heatmap(exprs(eset_of_interest))
# heatmaps are defaulted to euclidean distance
# using correlation instead, and removed gene labels and changed color
library(RColorBrewer)
interesting_genes_all_7 <- rownames(topTable(fitted.ebayes, coef=7, number = Inf, p.value = 0.01, lfc = 3))
interesting_genes_all_7
detail7 <- AnnotationDbi::select(huex10sttranscriptcluster.db,interesting_genes_all_7,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
heatmap(exprs(eset_of_interest),
        labCol=f,labRow=(detail7[2][, 1]),
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x)))
)
