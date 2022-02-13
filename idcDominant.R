library(GEOquery)
library(limma)
library("seqinr") 
library(dplyr)

# loading dataset from GEO
gset = getGEO("GSE55715")

# data cleaning to extract required genes
genes = gset$GSE55715_series_matrix.txt.gz@featureData@data
genesEx = gset$GSE55715_series_matrix.txt.gz@assayData$exprs

which(genes$ILMN_Gene == "TP53")
which(genes$ILMN_Gene == "ERBB2")

tp53 = genes[19326,]
tp53 = tp53 %>% select(c("ID", "Symbol","RefSeq_ID"))
tp53Ex = genesEx[19326,]

erbb2 = genes[12683,]
erbb2 = erbb2 %>% select(c("ID", "Symbol","RefSeq_ID"))
erbb2Ex = genesEx[12683,]

# i) checking on any outliers
boxplot(tp53Ex, main = "Boxplot for Expression of Gene TP53")
boxplot(erbb2Ex, main = "Boxplot for Expression of Gene ERBB2")

# ii) checking on the normality
plotDensities(tp53Ex, main = "Expression Value Distribution Plot for TP53")
plotDensities(erbb2Ex, main = "Expression Value Distribution Plot for ERBB2")

qqnorm(tp53Ex, main = "Normal Q-Q Plot for Expression of Gene TP53")
qqline(tp53Ex)
qqnorm(erbb2Ex, main = "Normal Q-Q Plot for Expression of Gene ERBB2")
qqline(erbb2Ex)

# iii) computing the mean and median
mean(tp53Ex)
median(tp53Ex)
mean(erbb2Ex)
median(erbb2Ex)

# iv) GC content
tp53Fasta <- read.fasta(file = 'TP53.fasta')
erbb2Fasta <- read.fasta(file = 'ERBB2.fasta')
tp53Seq <- tp53Fasta[[1]]
erbb2Seq <- erbb2Fasta[[1]]
GC(tp53Seq)
GC(erbb2Seq)
