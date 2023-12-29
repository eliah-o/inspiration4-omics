# All codes are written by Ezequiel Dantas (ecd4001@med.cornell.edu)

library(dplyr)
library(Biobase)
library(NOISeq)
library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) 
library("BisqueRNA")
library("SeuratObject")
library("scater")
library(loomR)
library(ggplot2)
library(ggforce)
library(DESeq2)
library(rstatix)
library(EnhancedVolcano)
library(ggpubr)
library(stringr)
library(ggrepel)
library(sva)
library(mgsub)
library("bladderbatch")
library("Harman")
library("RUVSeq")
library(qsvaR)
library("limma")
library("org.Mm.eg.db")
library("GeneTonic")
library(GeneStructureTools)
library(vsn)
library("RNAseqQC")
library(GeneStructureTools)
library(pheatmap)

##%######################################################%##
#                                                          #
####                        STAR                        ####
#                                                          #
##%######################################################%##


data <- read.csv("STAR_TA.csv", stringsAsFactors = T, header = T)
metadata <- read.csv("metanasa.csv", header = T, stringsAsFactors = T)

metastar<- with(metadata, metadata[order(Condition, Sample),])

rownames(data) <- data$gene_name
data$gene_name <- NULL

rownames(metastar) <- metastar$Sample
metastar$Sample <- NULL

all(rownames(metastar) == colnames(data)) # if it is true move to analysis otherwise following the seline below
match(colnames(data), rownames(metastar))
idxstar <- match(rownames(metastar), colnames(data))
reorder_data <- data[,idxstar]
#str(reorder_metadata)
all(rownames(metastar) == colnames(reorder_data))

# Filtering
keep <- rowSums(reorder_data) > 20
data_fil <- reorder_data[keep,]
dim(data_fil)

str(data_fil)
#Library sizes
librarySizes <- colSums(data_fil)
barplot(librarySizes, 
  names=names(librarySizes), 
  las=2, 
  main="Barplot of library sizes")
abline(h=2e7, lty=2)

dds <- DESeqDataSetFromMatrix(countData = round(data_fil),
  colData = metastar,
  design = ~ Condition)

dds <- DESeq(dds, betaPrior = FALSE)
vsd <- vst(dds, blind=FALSE)
?DESeq
res<- results(dds, contrast =c("Condition", "Flight","GC"), alpha = 0.05, lfcThreshold = 0.32)


resultsNames(dds)
normcounts<- counts(dds, normalized=TRUE)
normcounts <- as.data.frame(normcounts)
write.csv(normcounts, "star_norm.csv")

cytostar <- subset(normcounts, rownames(normcounts) %in% cytos)

cytostar$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cytostar), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(cytostar) <- cytostar$symbol
cytostar$symbol <- NULL

pheatmap(cytostar, cluster_rows=F, cluster_cols=F, annotation_col = metastar,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= c(8, 17,26),
  cellheight=20, cellwidth = 35)


boxplot(assay(vsd)["ENSMUSG00000025746",] ~ colData(vsd)[,"Condition"])


####QC
mean_sd_plot(vsd)
plot_total_counts(dds)
plot_library_complexity(dds)
plot_gene_detection(dds)
plotDispEsts(dds)

pcaData <- plotPCA(vsd, intgroup = c("Condition"), ntop = 500, returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 6, alpha = 0.9) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  geom_text_repel(label= rownames(metadata))+
  scale_color_manual(breaks = c("Flight", "GC", "Vivarium"),values = c("blue", "red", "steelblue3")) +
  labs(legend = "Dose (mg)")+
  
  #scale_shape_manual(breaks = c("Batch_1", "Batch_2"),
  # values = c(22, 24)) +
  # stat_ellipse(type = "t") +
  
  theme(panel.border = element_blank(), legend.position = "top", # borders white
    panel.background = element_rect(fill = NA), # background withe
    plot.margin = ggplot2::margin(2, 4, 2, 4, "cm"),
    panel.grid.major = element_line(colour= "gray95"), # grid empty
    panel.grid.minor = element_blank(),  # grid empty
    axis.line = element_line(colour = "black", size = 0.5),
    text = element_text(color = "black", size = 20),
    axis.title.y= element_text(colour= "black", size= 20),
    axis.text.x= element_text(size=20, color = "black"),
    axis.text.y= element_text(size=20, color = "black"),
    aspect.ratio = 1)

##%######################################################%##
#                                                          #
####                       SALMON                       ####
#                                                          #
##%######################################################%##


dataSAL <- read.csv("salmon_merged_gene_counts.csv", stringsAsFactors = T, header = T, check.names = F)
#datasal <- gsub("-", "_", colnames(dataSAL))
metasal <- read.csv("meta_salmon.csv", header = T, stringsAsFactors = T)
metasal1<- with(metasal, metasal[order(Condition, Sample),])

rownames(dataSAL) <- dataSAL$gene_id
dataSAL$gene_id <- NULL

rownames(metasal1) <- metasal1$Sample
metasal1$Sample <- NULL

all(rownames(metasal1) == colnames(dataSAL)) # if it is true move to analysis otherwise following the seline below
match(colnames(dataSAL), rownames(metasal1))
idx <- match(rownames(metasal1), colnames(dataSAL))
reorder_datasal <- dataSAL[,idx]
#str(reorder_metadata)
all(rownames(metasal1) == colnames(reorder_datasal))

# Filtering
keepsal <- rowSums(reorder_datasal) > 20
datasal_fil <- reorder_datasal[keepsal,]
dim(datasal_fil)

#Library sizes
librarySizes_sal <- colSums(datasal_fil)

barplot(librarySizes_sal, 
  names=names(librarySizes_sal), 
  las=2, 
  main="Barplot of library sizes")
abline(h=2e7, lty=2)

ddssal <- DESeqDataSetFromMatrix(countData = round(datasal_fil),
  colData = metasal1,
  design = ~ Condition)

ddssal <- DESeq(ddssal, betaPrior = FALSE)
vsdsal <- vst(ddssal, blind=FALSE)

normcounts_sal<- counts(ddssal, normalized=TRUE)
normcounts_sal <- as.data.frame(normcounts_sal)
write.csv(normcounts_sal, "sal_norm.csv")


keytypes(org.Mm.eg.db)

cytos <- c("ENSMUSG00000018930", "ENSMUSG00000034855", "ENSMUSG00000061068", "ENSMUSG00000025746", "ENSMUSG00000016529", "ENSMUSG00000024402", "ENSMUSG00000036117", "ENSMUSG00000034394", "ENSMUSG00000094686", "ENSMUSG00000095675", "ENSMUSG00000094065", "ENSMUSG00000024401", "ENSMUSG00000031778", "ENSMUSG00000003206", "ENSMUSG00000044701", "ENSMUSG00000040329", "ENSMUSG00000035385", "ENSMUSG00000073888", "ENSMUSG00000096826", "ENSMUSG00000029371", "ENSMUSG00000000869")
cytosal <- subset(normcounts_sal, rownames(normcounts_sal) %in% cytos)

cytosal$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cytosal), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(cytosal) <- cytosal$symbol
cytosal$symbol <- NULL

pheatmap(cytosal, cluster_rows=F, cluster_cols=F, annotation_col = metasal1,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= c(9, 18,27),
  cellheight=20, cellwidth = 35)

#IL-6 ENSMUSG00000025746 ENSMUSG00000025746
boxplot(assay(vsdsal)["ENSMUSG00000025746",] ~ colData(vsdsal)[,"Condition"])


####QC
mean_sd_plot(vsdsal)
plot_total_counts(ddssal)
plot_library_complexity(ddssal)
plot_gene_detection(ddssal)
plotDispEsts(ddssal)

pcaDatasal <- plotPCA(vsdsal, intgroup = c("Condition"), ntop = 500, returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaDatasal, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 6, alpha = 0.9) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  geom_text_repel(label= rownames(metasal))+
  scale_color_manual(breaks = c("Flight", "GC", "Vivarium"),values = c("blue", "red", "steelblue3")) +
  labs(legend = "Dose (mg)")+
  
  #scale_shape_manual(breaks = c("Batch_1", "Batch_2"),
  # values = c(22, 24)) +
  # stat_ellipse(type = "t") +
  
  theme(panel.border = element_blank(), legend.position = "top", # borders white
    panel.background = element_rect(fill = NA), # background withe
    plot.margin = ggplot2::margin(2, 4, 2, 4, "cm"),
    panel.grid.major = element_line(colour= "gray95"), # grid empty
    panel.grid.minor = element_blank(),  # grid empty
    axis.line = element_line(colour = "black", size = 0.5),
    text = element_text(color = "black", size = 20),
    axis.title.y= element_text(colour= "black", size= 20),
    axis.text.x= element_text(size=20, color = "black"),
    axis.text.y= element_text(size=20, color = "black"),
    aspect.ratio = 1)

salmonres<- results(ddssal, contrast =c("Condition", "Flight","GC"), alpha = 0.05, lfcThreshold = 0.32)
Salres <- as.data.frame(salmonres)
write.csv(Salres, "576_Salmon_res.csv")
cytos <- c("ENSMUSG00000034855", "ENSMUSG00000061068", "ENSMUSG00000025746", "ENSMUSG00000016529", "ENSMUSG00000036117", "ENSMUSG00000024401", "ENSMUSG00000031778", "ENSMUSG00000003206", "ENSMUSG00000044701", "ENSMUSG00000040329", "ENSMUSG00000035385", "ENSMUSG00000029371", "ENSMUSG00000000869", "ENSMUSG00000026981", "ENSMUSG00000037942", "ENSMUSG00000021538", "ENSMUSG00000040770", "ENSMUSG00000000982", "ENSMUSG00000009185")

cytores <- subset(salmonres, rownames(salmonres) %in% cytos)

cytores <- as.data.frame(cytores)

cytores$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cytores), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")


rownames(cytores) <- cytores$symbol
cytores$symbol <- NULL
write.csv(cytores, "salmon_576_heatmapres.csv" )


#6*5

EnhancedVolcano(toptable = cytores,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(cytores),
  #selectLab = c("Mcpt4","Ccl2", "Ccl4", "Cxcl5", "Il6"),
  #xlim = c(-5, +6),
  labSize = 5,
  #ylim = c(0,20),
  FCcutoff = 1,
  pCutoff =0.05,
  drawConnectors = TRUE,
  maxoverlapsConnectors = Inf,
  widthConnectors = 0.20,
  colConnectors = "black",
  arrowheads = FALSE,
  typeConnectors = "closed",
  lengthConnectors = unit(1, "npc"),
  boxedLabels = FALSE,
  max.overlaps = Inf)

salmonres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(salmonres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(salmonres) <- salmonres$symbol

EnhancedVolcano(toptable = salmonres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(salmonres),
  #selectLab = c("Mcpt4","Ccl2", "Ccl4", "Cxcl5", "Il6"),
  #xlim = c(-5, +6),
  labSize = 5,
  #ylim = c(0,20),
  FCcutoff = 1,
  pCutoff =0.01,
  drawConnectors = TRUE,
  maxoverlapsConnectors = Inf,
  widthConnectors = 0.20,
  colConnectors = "black",
  arrowheads = FALSE,
  typeConnectors = "closed",
  lengthConnectors = unit(1, "npc"),
  boxedLabels = FALSE,
  max.overlaps = Inf)

pheat <- subset(normcounts_sal, rownames(normcounts_sal) %in% cytos)
pheat$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(pheat), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")


rownames(pheat) <- pheat$symbol
pheat$symbol <- NULL

pheat1 <- pheat[,1:18]
pheatmet <- metasal1[1:18, ]
pheatmet <- as.data.frame(pheatmet)
rownames(pheatmet) <- colnames(pheat1)

pheatmap(pheat1, cluster_rows=F, cluster_cols=F, annotation_col = pheatmet,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= 9,
  cellheight=20, cellwidth = 35)



