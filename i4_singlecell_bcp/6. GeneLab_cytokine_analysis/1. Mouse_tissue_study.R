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

##https://doi.org/10.1038/s42003-020-01227-2

data <- read.csv("counts.csv", stringsAsFactors = T, header = T, check.names = FALSE)
rownames(data) <- data$ensemble_id
relevant <- c("900005Y", "900006Y", "900009Y", "900010Y", "900011Y", "900012Y", "900013Y", "900014Y", "900015Y", "900019Y", "900020Y", "900021Y", "900027Y", "900029Y", "900030Y", "900031Y", "900032Y", "900034Y", "900039Y", "900041Y", "900042Y", "900043Y", "900044Y", "900046Y", "900051Y", "900053Y", "900054Y", "900055Y", "900056Y", "900058Y", "900063Y", "900065Y", "900066Y", "900067Y", "900068Y", "900070Y", "900075Y", "900077Y", "900078Y", "900079Y", "900080Y", "900082Y", "900087Y", "900089Y", "900090Y", "900091Y", "900092Y", "900094Y", "900099Y", "900101Y", "900102Y", "900103Y", "900104Y", "900106Y", "900111Y", "900113Y", "900114Y", "900115Y", "900116Y", "900118Y", "900123Y", "900125Y", "900126Y", "900127Y", "900128Y", "900130Y", "900135Y", "900137Y", "900138Y", "900139Y", "900140Y", "900142Y", "900145Y", "900146Y", "900148Y", "900151Y", "900152Y", "900154Y", "900158Y", "900159Y", "900160Y", "900165Y", "900166Y", "900167Y")
data1<- data[, relevant]

metadata <- read.csv("metadata.csv", header = T, stringsAsFactors = T)
metadata$genotype <- NULL

metadataor<- with(metadata, metadata[order(Tissue, condition, sample),])
rownames(metadataor) <- metadataor$sample
metadataor$sample <- NULL

metadataor1<- read.csv("metadataor.csv", header = T, stringsAsFactors = T, check.names = FALSE)
rownames(metadataor1) <- metadataor1$sample
metadataor1$sample <- NULL

all(rownames(metadataor) == colnames(data1)) # if it is true move to analysis otherwise following the seline below
match(colnames(data1), rownames(metadataor))
idx <- match(rownames(metadataor), colnames(data1))
reorder_data <- data1[,idx]
#str(reorder_metadata)
all(rownames(metadataor1) == colnames(reorder_data))

# Filtering
rownames(reorder_data)<- removeVersion(rownames(reorder_data))
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
  colData = metadataor1,
  design = ~ merged)

dds <- DESeq(dds, betaPrior = FALSE)
vsd <- vst(dds, blind=FALSE)

#res<- results(dds, contrast =c("Condition", "Flight","GC"), alpha = 0.05, lfcThreshold = 0.32)

normcounts<- counts(dds, normalized=TRUE)
normcounts <- as.data.frame(normcounts)
write.csv(normcounts, "norm_japan.csv")

cytos <- c("ENSMUSG00000034855", "ENSMUSG00000061068", "ENSMUSG00000025746", "ENSMUSG00000016529", "ENSMUSG00000036117", "ENSMUSG00000024401", "ENSMUSG00000031778", "ENSMUSG00000003206", "ENSMUSG00000044701", "ENSMUSG00000040329", "ENSMUSG00000035385", "ENSMUSG00000029371", "ENSMUSG00000000869", "ENSMUSG00000026981", "ENSMUSG00000037942", "ENSMUSG00000021538", "ENSMUSG00000040770", "ENSMUSG00000000982", "ENSMUSG00000009185")

cytostar <- subset(normcounts, rownames(normcounts) %in% cytos)

cytostar$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cytostar), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(cytostar) <- cytostar$symbol
cytostar$symbol <- NULL

rownames(metadata) <- metadata$sample
metadata$sample <- NULL
pheatmap(cytostar, cluster_rows=F, cluster_cols=F, annotation_col = metadata,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= c(6, 12, 15, 18,24,30,33,36,42,48,54,60,67,72,78),
  cellheight=20, cellwidth = 35)

pheatmap(cytostar, cluster_rows=F, cluster_cols=F, annotation_col = metadata,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= c(12, 18,30,36,48,60,72),
  cellheight=9, cellwidth = 10)

write.csv(metadataor, "metadataor.csv")
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

#Contrasts

wat<- results(dds, contrast =c("merged", "White_adipose_tissue_Flight","White_adipose_tissue_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
watres <- subset(wat, rownames(wat) %in% cytos)


watres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(watres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

watres <- as.data.frame(watres)
rownames(watres) <- watres$symbol
watres$symbol <- NULL
str(watres)

write.csv(watres, "wat_res.csv")

wat$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(wat), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

wat <- as.data.frame(wat)
#6*5

EnhancedVolcano(toptable = watres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(watres),
  selectLab = c("Mcpt4","Ccl2", "Ccl4", "Cxcl5", "Il6", "Ccl3"),
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


BAT<- results(dds, contrast =c("merged", "Brown_adipose_tissue_Flight","Brown_adipose_tissue_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
battres <- subset(BAT, rownames(BAT) %in% cytos)

battres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(battres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

battres <- as.data.frame(battres)
rownames(battres) <- battres$symbol
battres$symbol <- NULL
str(battres)

write.csv(battres, "bat_res.csv")
#6*5

EnhancedVolcano(toptable = battres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(watres),
  selectLab = c("Tnf","Ccl8", "Ccl2", "Crp", "Cx3cl1"),
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

kidney<- results(dds, contrast =c("merged", "Kidney_Flight","Kidney_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
kidneyres <- subset(kidney, rownames(kidney) %in% cytos)

kidneyres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(kidneyres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

kidneyres <- as.data.frame(kidneyres)
rownames(kidneyres) <- kidneyres$symbol
kidneyres$symbol <- NULL

write.csv(kidneyres, "kidney_res.csv")


#6*5

EnhancedVolcano(toptable = kidneyres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(kidneyres),
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

liver<- results(dds, contrast =c("merged", "Liver_Flight","Liver_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
liverres <- subset(liver, rownames(liver) %in% cytos)

liverres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(liverres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

liverres <- as.data.frame(liverres)
rownames(liverres) <- liverres$symbol
liverres$symbol <- NULL

write.csv(liverres, "liver_res.csv")

#6*5

EnhancedVolcano(toptable = liverres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(liverres),
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

man<- results(dds, contrast =c("merged", "Mandibular_bone_Flight","Mandibular_bone_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
manres <- subset(man, rownames(man) %in% cytos)

manres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(manres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

manres <- as.data.frame(manres)
rownames(manres) <- manres$symbol
manres$symbol <- NULL

write.csv(manres, "man_res.csv")


#6*5

EnhancedVolcano(toptable = manres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(manres),
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

spleen<- results(dds, contrast =c("merged", "Spleen_Flight","Spleen_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
spleenres <- subset(spleen, rownames(spleen) %in% cytos)

spleenres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(spleenres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

spleenres <- as.data.frame(spleenres)
rownames(spleenres) <- spleenres$symbol
spleenres$symbol <- NULL

write.csv(spleen, "spleen_res.csv")

#6*5

EnhancedVolcano(toptable = spleenres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(spleenres),
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

spleen<- results(dds, contrast =c("merged", "Spleen_Flight","Spleen_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
spleenres <- subset(spleen, rownames(spleen) %in% cytos)

spleenres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(spleenres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

spleenres <- as.data.frame(spleenres)
rownames(spleenres) <- spleenres$symbol
spleenres$symbol <- NULL

write.csv(spleenres, "spleen_res.csv")

#6*5

EnhancedVolcano(toptable = spleenres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(spleenres),
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

temporal<- results(dds, contrast =c("merged", "Temporal_bone_Flight","Temporal_bone_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
temporalres <- subset(temporal, rownames(temporal) %in% cytos)

temporalres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(temporalres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

temporalres <- as.data.frame(temporalres)
rownames(temporalres) <- temporalres$symbol
temporalres$symbol <- NULL

write.csv(temporalres, "temporal_res.csv")
#6*5

EnhancedVolcano(toptable = temporalres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(temporalres),
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

thymus<- results(dds, contrast =c("merged", "Thymus_Flight","Thymus_Ground_control"), alpha = 0.05, lfcThreshold = 0.32)
thymusres <- subset(thymus, rownames(thymus) %in% cytos)

thymusres$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(thymusres), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

thymusres <- as.data.frame(thymusres)
rownames(thymusres) <- thymusres$symbol
thymusres$symbol <- NULL

write.csv(thymusres, "thymus_res.csv")

#6*5

EnhancedVolcano(toptable = thymusres,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = rownames(thymusres),
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
