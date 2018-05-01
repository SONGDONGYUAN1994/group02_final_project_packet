source("https://bioconductor.org/biocLite.R")
library(DiffBind)

# Obtaining differentially bound sites
workdir <- "C:/Users/songdongyuan/group02_final_project_packet"
setwd(workdir)

design <- dba(sampleSheet='H3K4me3_design.csv')

# Heatmap for initial clustering of the samples
dba.plotHeatmap(design, main="H3K4me3: Correlation heatmap using occupancy data")

# Counting reads
design_count <- dba.count(design,summits = 250)
# Correlation heatmap based on the aﬃnity scores
jpeg(filename = "H3K4me3_heat_ocu.jpeg",width = 800, quality = 100)
plot(design_count, main="H3K4me3: Correlation heatmap using affinity (read count) data")
dev.off()
# Establishing a contrast
designContrast <- dba.contrast(design_count, minMembers = 2,categories=DBA_CONDITION)

# Performing the differential analysis 
designDB <- dba.analyze(designContrast) #Using only differential binding sites

# Correlation heatmap based on te differential binding sites
jpeg(filename = "H3K4me3_heat_DB.jpeg",width = 800, quality = 100)
plot(designDB, contrast=1,main="H3K4me3: Correlation heatmap, using only significantly differentially bound sites")
dev.off()
# Retrieving the differentially bound sites
desired_FDR <- 0.05
desired_fc <- 2
design.DB <- dba.report(designDB,th=0.001);design.DB#, fold=2)
design.DB$score <- -10*log10(design.DB$FDR)

# MA plots: Visualize the eﬀect of normalization on data
jpeg(filename = "H3K4me3_MA.jpeg",width = 800, quality = 100)
dba.plotMA(designDB, contrast=1)
dev.off()

# Affinity binding boxplot:  highlight signiﬁcantly diﬀerentially bound sites and show their fold changes
jpeg(filename = "H3K4me3_boxplot.jpeg",width = 800, quality = 100)
dba.plotBox(designDB)
dev.off()

# Affinity heatmap: view the patterns of binding aﬃnity directly in the diﬀerentially bound sites 
jpeg(filename = "H3K4me3_heat_affinity.jpeg",width = 800, quality = 100)
dba.plotHeatmap(designDB, correlations=FALSE)
dev.off()

# Valcano Plot: highlight signiﬁcantly diﬀerentially bound sites and show their fold changes
design.all <- dba.report(designDB, th=1, fold = 0) #??????????????
num.DB <- -log10(design.all$FDR) >= -log10(desired_FDR) & abs(design.all$Fold)>=desired_fc

jpeg(filename = "H3K4me3_vocano.jpeg",width = 800, quality = 100)
plot(design.all$Fold, -log10(design.all$FDR), pch=20, cex=0.3,
     col= ifelse(num.DB, "red", "grey"),
     xlab = "log2 fold change (co-culture - monoculture)",
     ylab = " -log10(FDR)",
     main = paste("Contrast: monoculture vs co-culture [", sum(num.DB) ,
                  paste0("FDR<",desired_FDR," & FoldChange>2")))
dev.off()

########## H3K27ac #############
design_2 <- dba(sampleSheet='H3K27ac_design.csv') 

jpeg(filename = "H3K27ac_heat_ocu.jpeg",width = 800, quality = 100)
plot(design_2, main="H3K27ac: Correlation heatmap using occupancy data")
dev.off()

design_count_2 <- dba.count(design_2,summits = 250)

jpeg(filename = "H3K27ac_heat_DB.jpeg",width = 800, quality = 100)
plot(design_count_2, main="H3K27ac: Correlation heatmap using affinity (read count) data")
dev.off()

designContrast_2 <- dba.contrast(design_count_2, minMembers = 2,categories=DBA_CONDITION)
designDB_2 <- dba.analyze(designContrast_2)

design.DB_2 <- dba.report(designDB_2,th=0.001);design.DB_2#, fold=2) #####GRanges object
design.DB_2$score <- -10*log10(design.DB_2$FDR)

#MA plots 
jpeg(filename = "H3K27ac_MA.jpeg",width = 800, quality = 100)
dba.plotMA(designDB_2, contrast=1)
dev.off()

# affinity binding boxplot
jpeg(filename = "H3K27ac_boxplot.jpeg",width = 800, quality = 100)
dba.plotBox(designDB_2)
dev.off()

# affinity heatmap
dba.plotHeatmap(designDB_2,contrast = 1, correlations=FALSE)
dba.plotHeatmap(designDB,contrast = 1, correlations=FALSE)
#Valcano Plot
design.all_2 <- dba.report(designDB_2, th=1, fold = 0) #??????????????
num.DB_2 <- -log10(design.all_2$FDR) >= -log10(desired_FDR) & abs(design.all_2$Fold)>=desired_fc

jpeg(filename = "H3K27ac_volcano.jpeg",width = 800, quality = 100)
plot(design.all_2$Fold, -log10(design.all_2$FDR), pch=20, cex=0.3,
     col= ifelse(num.DB_2, "red", "grey"),
     xlab = "log2 fold change (co-culture - monoculture)",
     ylab = " -log10(FDR)",
     main = paste("Contrast: monoculture vs co-culture [", sum(num.DB_2) ,
                  paste0("FDR<",desired_FDR," & FoldChange>2")))
dev.off()

########### For motif analysis ##############
df_HH <- data.frame(seqnames=seqnames(design.DB),
                    starts=start(design.DB)-1,
                    ends=end(design.DB),
                    names=names(design.DB),
                    scores=design.DB$score,
                    strands=strand(design.DB),
                    Fold=design.DB$Fold)
write.table(df_HH, file=workdir, 
            quote=F, sep="\t", row.names=F, col.names=F)

df_222 <- data.frame(seqnames=seqnames(design.DB_2),
                     starts=start(design.DB_2)-1,
                     ends=end(design.DB_2),
                     names=names(design.DB_2),
                     scores=design.DB_2$score,
                     strands=strand(design.DB_2),
                     Fold=design.DB_2$Fold)
write.table(df_222, file=workdir, 
            quote=F, sep="\t", row.names=F, col.names=F)


