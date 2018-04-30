setwd("C:/Users/songdongyuan/BST281_HW")

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

#-------------------------------------------------------------------
  #Remeber! Four Steps!
  # 1. Filter none/low expression genes (*has been finished in Python*)
  # 2. Normalization
  # 3. Cal differential expressed genes
  # 4. Extract gene lists and do GO/enrichment analysis (you can use web tools, e.g. use David)
#-------------------------------------------------------------------
setwd("C:/Users/songdongyuan/BST281_HW")
#
filtered <- read.csv("filtered.tsv", sep = "\t", row.names = 1, header= TRUE, stringsAsFactors = F)
names(filtered)
filtered <- filtered[, c(2,4,6, 1,3,5)]

group <- c("Co", "Co", "Co", "Mono", "Mono", "Mono")
names(filtered) <- c("RPMI_HS5", "MM_HS5", "KMS11_HS5", "RPMI", "MM", "KMS11")

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
# Notice: please change the order so that it can be case group and control group

# Notice! Has not been normalized!


# MA plot
maPlot(filtered_counts$counts[,1], filtered_counts$counts[,2], lowess = T)
title("MA plot (unnormalised)")
# Very flat

#plotMDS(filtered_counts)


#--------------------------------------------------------------------
# Apply normalisation to DGEList object
filtered_counts1 <- calcNormFactors(filtered_counts, method = "TMM")

filtered_counts2 <- calcNormFactors(filtered_counts, method = "upperquartile")

maPlot(filtered_counts1$counts[,1], filtered_counts1$counts[,2], lowess = T)
title("MA plot (normalised with TMM)")

maPlot(filtered_counts2$counts[,1], filtered_counts2$counts[,2], lowess = T)
title("MA plot (normalised with upperquartile)")

# After Normalization
# Get log2 counts per million
logcounts1 <- cpm(filtered_counts1,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts1, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts1),col="blue")
title("Boxplots of logCPMs (TMM normalised)")

par(mfrow=c(1, 2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
boxplot(logcounts1, xlab="", ylab="Log2 counts per million",las=2)

# Very little difference!

# Using limma

# Create design matrix
design <- model.matrix(~ 0 + group)
design

#
logCPM <- cpm(filtered_counts1, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))
summary(fit)
#---------------------------------------------------------------
# Use DESeq! Please see URL http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts
# Using DESeq2
# One thing that you need to notice: if you use DESeq2, you may need use raw counts and use DESeq2 from start  
  


str(filtered_counts)
filtered_counts$counts <- round(filtered_counts$counts)
filtered_counts1$counts <- round(filtered_counts1$counts)
filtered_counts2$counts <- round(filtered_counts2$counts)

# You can easily convert data format between edgeR and DESeq2
dds <- as.DESeqDataSet(filtered_counts)
res <- DESeq(dds)
res <- results( res )
summary(res)

sum( res$padj < 0.1, na.rm=TRUE )
resSig <- res[ which(res$padj < 0.1 ), ]
head( resSig[ order( resSig$log2FoldChange ), ] )

DEgene_list <- rownames(resSig)
write.csv(DEgene_list, file = "Significant_gene_list.csv", row.names = FALSE)
