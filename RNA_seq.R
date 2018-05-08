setwd("C:/Users/songdongyuan/group02_final_project_packet")

suppressMessages(install.packages("pheatmap"))
#source("https://bioconductor.org/biocLite.R")

#R3.5.0 has issues when installing data.table
#install.packages("https://socialsciences.mcmaster.ca/jfox/.Pickup/data.table_1.10.4-3.zip",
                 #repos=NULL, type="win.binary")

#biocLite("edgeR")
#biocLite("limma")
#biocLite("Glimma")
#biocLite("org.Mm.eg.db")
#biocLite("RColorBrewer")
#biocLite("DESeq2")
#biocLite("DEFormats")
#install.packages("gplots")
rm(list=ls())

suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(Glimma))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(RColorBrewer))
suppressMessages(library(DESeq2))
suppressMessages(library(DEFormats))
suppressMessages(library(pheatmap))


filtered <- read.csv("filtered.tsv", sep = "\t", row.names = 1, header= TRUE, stringsAsFactors = F)
filtered <- filtered[, c(3,1,5, 4,2,6)]

group <- c( "Mono", "Mono", "Mono", "Co", "Co", "Co")
names(filtered) <- c("MM", "RPMI", "KMS11", "MM_HS5", "RPMI_HS5", "KMS11_HS5")

filtered_counts <- DGEList(filtered, group = group, samples = c("MM", "RPMI", "KMS11", "MM", "RPMI", "KMS11") )

png("Barplot of library sizes.png")
#filtered_counts$samples$lib.size
barplot(filtered_counts$samples$lib.size,names=colnames(filtered_counts),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()

# Get log2 counts per million
logcounts <- cpm(filtered_counts,log=TRUE)

png("Boxplots of logCPMs (unnormalised).png")
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

# MA plot
png("MA plot (unnormalised).png")
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
dev.off()

# Apply normalisation to DGEList object
filtered_counts_n <- calcNormFactors(filtered_counts, method = "TMM")

png("MA plot (normalised with TMM).png")
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
dev.off()

png("Boxplots of logCPMs.png")
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
dev.off()


# Using limma

# Create design matrix
design <- model.matrix(~ 0 + group)

# Fit limma
#logCPM <- cpm(filtered_counts_n, log=TRUE, prior.count=1)
#fit <- lmFit(logCPM, design)
#fit <- eBayes(fit, trend=TRUE)

# Using DESeq2

filtered_counts_n$counts <- round(filtered_counts_n$counts)

# You can easily convert data format between edgeR and DESeq2
dds <- as.DESeqDataSet(filtered_counts_n)
res <- DESeq(dds)
res <- results(res)
summary(res)

vsd <- vst(dds, blind=FALSE)

png("PCA.png")
plotPCA(vsd, intgroup=c("group", "samples"))
dev.off()

#sampleDists <- dist(t(assay(vsd)))
#sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$samples, vsd$group, sep="-")
#colnames(sampleDistMatrix) <- NULL

#png("heatmap.png")
#pheatmap(sampleDistMatrix,
         #clustering_distance_rows=sampleDists,
         #clustering_distance_cols=sampleDists)
#dev.off()



resSig <- res[ which(res$padj < 0.1 ), ]

DEgene_list <- rownames(resSig)
write.table(DEgene_list, file = "DEgene_list.tsv", row.names = FALSE, sep = '\t', col.names = F, quote=F)




####################
# Auto pre-processing
####################
raw <- read.csv("expressionFile_counts_MM.csv", row.names = 1, header= TRUE, stringsAsFactors = F)
colnames(raw) <- c("RPMI", "RPMI(+HS5)", "MM1S", "MM1S(+HS5)", "KMS11", "KMS11(+HS5)")
png("heatmap2.png")
pheatmap(cor(raw))
dev.off()

raw <- round(raw[, c(3,1,5, 4,2,6)])

group <- c( "Mono", "Mono", "Mono", "Co", "Co", "Co")
names(raw) <- c("MM", "RPMI", "KMS11", "MM_HS5", "RPMI_HS5", "KMS11_HS5")
raw_counts <- DGEList(raw, group = group, samples = c("MM", "RPMI", "KMS11", "MM", "RPMI", "KMS11") )

dds <- as.DESeqDataSet(raw_counts)
res <- DESeq(dds)
res <- results(res)
summary(res)


resSig <- res[ which(res$padj < 0.1 ), ]

DEgene_list_auto <- rownames(resSig)
write.table(DEgene_list_auto, file = "DEgene_list_auto.tsv", row.names = FALSE, sep = '\t', col.names = F, quote=F)

raw_count <- raw[row.names(raw) %in% DEgene_list_auto,]

png("heatmap3.png")
pheatmap(cor(raw_count))
dev.off()
