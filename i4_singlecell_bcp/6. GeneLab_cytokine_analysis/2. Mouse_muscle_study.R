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

counts <- read.csv("salmon_merged_gene_counts.csv", header = T, stringsAsFactors = T, check.names = F)
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

metadata <- read.csv("Metadata_ed.csv", header = T, stringsAsFactors = T, check.names = F)
metadataor<- with(metadata, metadata[order(paper, Group, Sample),])

met_nat <- subset(metadataor, paper == "natbio" & Group %in% c("WT_FLT", "WT_GC"))
met_sci <- subset(metadataor, paper == "Scirep" & Group %in% c("MG", "GC"))
rownames(met_nat) <- met_nat$Sample
met_nat$Sample <- NULL
met_nat$paper <- NULL

rownames(met_sci) <- met_sci$Sample
met_sci$Sample <- NULL
met_sci$paper <- NULL

scicount <- counts[,rownames(met_sci)]
all(rownames(met_sci) == colnames(scicount)) # if it is true move to analysis otherwise following the seline below

# Filtering

natkeep <- rowSums(natbio) > 20
natkept <- natbio[natkeep,]
dim(natkept)

scikeep <- rowSums(scicount) > 20
scikept <- scicount[scikeep,]
dim(scikept)

#Library sizes
librarySizes <- colSums(natkept)
barplot(librarySizes, 
  names=names(librarySizes), 
  las=2, 
  main="Barplot of library sizes")
abline(h=2e7, lty=2)

ddsnat <- DESeqDataSetFromMatrix(countData = round(natkept),
  colData = met_nat,
  design = ~ Group)

ddssci <- DESeqDataSetFromMatrix(countData = round(scikept),
  colData = met_sci,
  design = ~ Group)

dim(met_nat)

ddsnat <- DESeq(ddsnat, betaPrior = FALSE)

ddssci <- DESeq(ddssci, betaPrior = FALSE)

vsdnat <- vst(ddsnat, blind=FALSE)
vsdsci <- vst(ddssci, blind=FALSE)

####QC
mean_sd_plot(vsdnat)
plot_total_counts(ddsnat)
plot_library_complexity(ddsnat)
plot_gene_detection(ddsnat)
plotDispEsts(ddsnat)

pcaDatanat <- plotPCA(vsdnat, intgroup = c("Group"), ntop = 500, returnData = TRUE)
pcaDatasci <- plotPCA(vsdsci, intgroup = c("Group"), ntop = 500, returnData = TRUE)

percentVarnat <- round(100 * attr(pcaDatanat, "percentVar"))
percentVarsci <- round(100 * attr(pcaDatasci, "percentVar"))


ggplot(percentVarsci, aes(PC1, PC2, color = Group)) +
  geom_point(size = 6, alpha = 0.9) +
  xlab(paste0("PC1: ",percentVarsci[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarsci[2],"% variance")) + 
  coord_fixed() +
  geom_text_repel(label= rownames(met_sci))+
  scale_color_manual(breaks = c("GC", "MG"),values = c("steelblue3", "red")) +
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

res_nat<- results(ddsnat, contrast =c("Group", "WT_FLT","WT_GC"), alpha = 0.05, lfcThreshold = 0.32)
write.csv(res_nat, "res_sol_nat.csv")
res_sci<- results(ddssci, contrast =c("Group", "MG","GC"), alpha = 0.05, lfcThreshold = 0.32)
res_sci<- as.data.frame(res_sci)
cytos <- c("ENSMUSG00000034855", "ENSMUSG00000061068", "ENSMUSG00000025746", "ENSMUSG00000016529", "ENSMUSG00000036117", "ENSMUSG00000024401", "ENSMUSG00000031778", "ENSMUSG00000003206", "ENSMUSG00000044701", "ENSMUSG00000040329", "ENSMUSG00000035385", "ENSMUSG00000029371", "ENSMUSG00000000869", "ENSMUSG00000026981", "ENSMUSG00000037942", "ENSMUSG00000021538", "ENSMUSG00000040770", "ENSMUSG00000000982", "ENSMUSG00000009185")

normcounts_nat<- counts(ddsnat, normalized=TRUE)
normcounts_nat<- as.data.frame(normcounts_nat)

normcounts_sci<- counts(ddssci, normalized=TRUE)
normcounts_sci<- as.data.frame(normcounts_sci)

write.csv(normcounts_nat, "norm_japan.csv")
cyto_nat <- subset(normcounts_nat, rownames(normcounts_nat) %in% cytos)
cyto_sci <- subset(normcounts_sci, rownames(normcounts_sci) %in% cytos)
cyto_sci <- as.data.frame(cyto_sci)

cyto_nat$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cyto_nat), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

cyto_sci$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cyto_sci), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(cyto_sci) <- cyto_sci$symbol
cyto_sci$symbol <- NULL

sci <- rownames(cyto_sci)
ressci <- subset(res_sci, symbol %in% rownames(cyto_sci))
write.csv(ressci, "res_sci_heat.csv")

nat <- cyto_nat$symbol
resnat <- subset(resnat1, symbol %in% nat)
write.csv(resnat, "res_nat_heat.csv")

pheatmap(cyto_sci, cluster_rows=F, cluster_cols=F, annotation_col = met_sci,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= 3,
  cellheight=20, cellwidth = 35)


res_nat$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(res_nat), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

dup <- c("Arfip1", "Bcl2l2", "Bfar", "Erdr1", "Pakap", "Pms2", "Taf9", "Txnl4a", "NA")

res_nat <- as.data.frame(res_nat)
resnat1 <- subset(res_nat, !symbol %in% dup)
rownames(resnat1) <- resnat1$symbol



res_sci$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(res_sci), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

write.csv(res_sci, "res_sol_sci.csv")
#6*5

EnhancedVolcano(toptable = res_sci,               #associated with low count genes is removed 
  x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
  y = "pvalue",                     
  lab = res_sci$symbol,
  #selectLab = c("Tubb6","Ccl2", "Ankrd2", "Cd68"),
  xlim = c(-10, 10),
  labSize = 6,
  #ylim = c(0,20),
  FCcutoff = 2.5,
  pCutoff =0.00001,
  drawConnectors = TRUE,
  maxoverlapsConnectors = 20,
  widthConnectors = 0.20,
  colConnectors = "black",
  arrowheads = FALSE,
  typeConnectors = "closed",
  lengthConnectors = unit(1, "npc"),
  boxedLabels = FALSE,
  max.overlaps = 20)

session()

sessionInfo()
?EnhancedVolcano



