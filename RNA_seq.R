setwd("C:/Users/songdongyuan/group02_final_project_packet")

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("limma")
biocLite("Glimma")
biocLite("org.Mm.eg.db")
biocLite("RColorBrewer")
biocLite("DESeq2")
biocLite("DEFormats")

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(DESeq2)
library(DEFormats)

filtered <- read.csv("filtered.tsv", sep = "\t", row.names = 1, header= TRUE, stringsAsFactors = F)
filtered <- filtered[, c(4,2,6, 3,1,5)]

group <- c("Co", "Co", "Co", "Mono", "Mono", "Mono")
names(filtered) <- c("MM_HS5", "RPMI_HS5", "KMS11_HS5", "MM", "RPMI", "KMS11")

filtered_counts <- DGEList(filtered, group = group)
filtered_counts$samples$lib.size
barplot(filtered_counts$samples$lib.size,names=colnames(filtered_counts),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# Get log2 counts per million
logcounts <- cpm(filtered_counts,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
# MA plot
par(mfrow = c(2,3))
maPlot(filtered_counts$counts[,1], filtered_counts$counts[,2], lowess = T)
title("MA plot (unnormalised) of #1 & #2")
maPlot(filtered_counts$counts[,1], filtered_counts$counts[,3], lowess = T)
title("MA plot (unnormalised) of #1 & #3")
maPlot(filtered_counts$counts[,2], filtered_counts$counts[,3], lowess = T)
title("MA plot (unnormalised) of #2 & #3")
maPlot(filtered_counts$counts[,4], filtered_counts$counts[,5], lowess = T)
title("MA plot (unnormalised) of #4 & #5")
maPlot(filtered_counts$counts[,4], filtered_counts$counts[,6], lowess = T)
title("MA plot (unnormalised) of #4 & #6")
maPlot(filtered_counts$counts[,5], filtered_counts$counts[,6], lowess = T)
title("MA plot (unnormalised) of #5 & #6")

# Apply normalisation to DGEList object
filtered_counts_n <- calcNormFactors(filtered_counts, method = "TMM")

par(mfrow = c(2,3))
maPlot(filtered_counts_n$counts[,1], filtered_counts_n$counts[,2], lowess = T)
title("MA plot (normalised with TMM) of #1 & #2")
maPlot(filtered_counts_n$counts[,1], filtered_counts_n$counts[,3], lowess = T)
title("MA plot (normalised with TMM) of #1 & #3")
maPlot(filtered_counts_n$counts[,2], filtered_counts_n$counts[,3], lowess = T)
title("MA plot (normalised with TMM) of #2 & #3")
maPlot(filtered_counts_n$counts[,4], filtered_counts_n$counts[,5], lowess = T)
title("MA plot (normalised with TMM) of #4 & #5")
maPlot(filtered_counts_n$counts[,4], filtered_counts_n$counts[,6], lowess = T)
title("MA plot (normalised with TMM) of #4 & #6")
maPlot(filtered_counts_n$counts[,5], filtered_counts_n$counts[,6], lowess = T)
title("MA plot (normalised with TMM) of #5 & #6")


par(mfrow = c(1, 2))
# Get log2 counts per million
logcounts <- cpm(filtered_counts,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
# Get log2 counts per million
logcounts <- cpm(filtered_counts_n,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (TMM normalised)")



# Using limma

# Create design matrix
design <- model.matrix(~ 0 + group)

# Fit limma
logCPM <- cpm(filtered_counts_n, log=TRUE, prior.count=1)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

# Using DESeq2

filtered_counts_n$counts <- round(filtered_counts_n$counts)

# You can easily convert data format between edgeR and DESeq2
dds <- as.DESeqDataSet(filtered_counts)
res <- DESeq(dds)
res <- results( res )
summary(res)

resSig <- res[ which(res$padj < 0.1 ), ]

DEgene_list <- rownames(resSig)
write.table(DEgene_list, file = "DEgene_list.tsv", row.names = FALSE, sep = '\t', col.names = F)