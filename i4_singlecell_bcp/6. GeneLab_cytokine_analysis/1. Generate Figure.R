library(dplyr)
library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(scater)
library(ggplot2)
library(ggforce)
library(DESeq2)
library(rstatix)
library(EnhancedVolcano)
library(ggpubr)
library(pheatmap)
library(ggrepel)
library(org.Mm.eg.db)
library(GeneStructureTools)

##%######################################################%##
#                                                          #
####               Extended Data Figure 2               ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####                  Panel A Heatmap                   ####
#                                                          #
##%######################################################%##
#Analyzed from Suzuki et al., 2020 (https://doi.org/10.1038/s42003-020-01227-2)
#Load the data
data <- read.csv("counts_suzuki.csv", stringsAsFactors = T, header = T, check.names = FALSE)
rownames(data) <- data$ensemble_id
relevant <- c("900005Y", "900006Y", "900009Y", "900010Y", "900011Y", "900012Y", "900013Y", "900014Y", "900015Y", "900019Y", "900020Y", "900021Y", "900027Y", "900029Y", "900030Y", "900031Y", "900032Y", "900034Y", "900039Y", "900041Y", "900042Y", "900043Y", "900044Y", "900046Y", "900051Y", "900053Y", "900054Y", "900055Y", "900056Y", "900058Y", "900063Y", "900065Y", "900066Y", "900067Y", "900068Y", "900070Y", "900075Y", "900077Y", "900078Y", "900079Y", "900080Y", "900082Y", "900087Y", "900089Y", "900090Y", "900091Y", "900092Y", "900094Y", "900099Y", "900101Y", "900102Y", "900103Y", "900104Y", "900106Y", "900111Y", "900113Y", "900114Y", "900115Y", "900116Y", "900118Y", "900123Y", "900125Y", "900126Y", "900127Y", "900128Y", "900130Y", "900135Y", "900137Y", "900138Y", "900139Y", "900140Y", "900142Y", "900145Y", "900146Y", "900148Y", "900151Y", "900152Y", "900154Y", "900158Y", "900159Y", "900160Y", "900165Y", "900166Y", "900167Y")
data1<- data[, relevant]

metadata <- read.csv("metadata_suzuki.csv", header = T, stringsAsFactors = T)
metadata$genotype <- NULL
metadataor<- with(metadata, metadata[order(Tissue, condition, sample),])
rownames(metadataor) <- metadataor$sample
metadataor$sample <- NULL

metadataor <- metadataor %>% 
  mutate(merged = paste(Tissue, condition, sep = "_"))

match(colnames(data1), rownames(metadataor))
idx <- match(rownames(metadataor), colnames(data1))
reorder_data <- data1[,idx]
all(rownames(metadataor) == colnames(reorder_data))

# Filtering
rownames(reorder_data)<- removeVersion(rownames(reorder_data))
keep <- rowSums(reorder_data) > 20
data_fil <- reorder_data[keep,]

#For Deseq2
dds <- DESeqDataSetFromMatrix(countData = round(data_fil),
  colData = metadataor,
  design = ~ merged)
dds <- DESeq(dds, betaPrior = FALSE)

normcounts<- counts(dds, normalized=TRUE)
normcounts <- as.data.frame(normcounts)

#Cytokines relevant for i4
cytos <- c("ENSMUSG00000034855", "ENSMUSG00000061068", "ENSMUSG00000025746", "ENSMUSG00000016529", "ENSMUSG00000036117", 
  "ENSMUSG00000024401", "ENSMUSG00000031778", "ENSMUSG00000003206", "ENSMUSG00000044701", "ENSMUSG00000040329", 
  "ENSMUSG00000035385", "ENSMUSG00000029371", "ENSMUSG00000000869", "ENSMUSG00000026981", "ENSMUSG00000037942", 
  "ENSMUSG00000021538", "ENSMUSG00000040770", "ENSMUSG00000000982", "ENSMUSG00000009185")

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
  size = 13, gaps_col= c(12, 18,30,36,48,60,72),
  cellheight=9, cellwidth = 10)

##%######################################################%##
#                                                          #
####              Panel B Mandibular Bone               ####
#                                                          #
##%######################################################%##


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

EnhancedVolcano(toptable = manres,               
  x = "log2FoldChange",           
  y = "pvalue",                     
  lab = rownames(manres),
  labSize = 5,
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

##%######################################################%##
#                                                          #
####            Panel C Brown Adipose Tissue            ####
#                                                          #
##%######################################################%##

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

EnhancedVolcano(toptable = battres,               
  x = "log2FoldChange",          
  y = "pvalue",                     
  lab = rownames(watres),
  selectLab = c("Tnf","Ccl8", "Ccl2", "Crp", "Cx3cl1"),
  labSize = 5,
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

##%######################################################%##
#                                                          #
####            Panel D White Adipose Tissue            ####
#                                                          #
##%######################################################%##

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

wat$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(wat), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

wat <- as.data.frame(wat)

EnhancedVolcano(toptable = watres,               
  x = "log2FoldChange",           
  y = "pvalue",                     
  lab = rownames(watres),
  selectLab = c("Mcpt4","Ccl2", "Ccl4", "Cxcl5", "Il6", "Ccl3"),
  labSize = 5,
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

##%######################################################%##
#                                                          #
####                      Panel F                       ####
#                                                          #
##%######################################################%##

#Load Data from Okada et.al. (doi: 10.1038/s41598-021-88392-4)

counts <- read.csv("counts_okada_hayashi.csv", header = T, stringsAsFactors = T, check.names = F)
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

metadata <- read.csv("metadata_okada_hayashi.csv", header = T, stringsAsFactors = T, check.names = F)
metadataor<- with(metadata, metadata[order(paper, Group, Sample),])

met_sci <- subset(metadataor, paper == "Scirep" & Group %in% c("MG", "GC"))

rownames(met_sci) <- met_sci$Sample
met_sci$Sample <- NULL
met_sci$paper <- NULL

scicount <- counts[,rownames(met_sci)]
all(rownames(met_sci) == colnames(scicount)) 

# Filtering


scikeep <- rowSums(scicount) > 20
scikept <- scicount[scikeep,]

ddssci <- DESeqDataSetFromMatrix(countData = round(scikept),
  colData = met_sci,
  design = ~ Group)

ddssci <- DESeq(ddssci, betaPrior = FALSE)
normcounts_sci<- counts(ddssci, normalized=TRUE)
normcounts_sci<- as.data.frame(normcounts_sci)
cyto_sci <- subset(normcounts_sci, rownames(normcounts_sci) %in% cytos)
cyto_sci <- as.data.frame(cyto_sci)

cyto_sci$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cyto_sci), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(cyto_sci) <- cyto_sci$symbol
cyto_sci$symbol <- NULL

pheatmap(cyto_sci, cluster_rows=F, cluster_cols=F, annotation_col = metadataor_okada,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= 3,
  cellheight=20, cellwidth = 35)

##%######################################################%##
#                                                          #
####                      Panel H                       ####
#                                                          #
##%######################################################%##

res_sci<- results(ddssci, contrast =c("Group", "MG","GC"), alpha = 0.05, lfcThreshold = 0.32)
res_sci<- as.data.frame(res_sci)

res_sci$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(res_sci), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

EnhancedVolcano(toptable = res_sci,               
  x = "log2FoldChange",          
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

##%######################################################%##
#                                                          #
####                      Panel E                       ####
#                                                          #
##%######################################################%##


met_nat <- subset(metadataor, paper == "natbio" & Group %in% c("WT_FLT", "WT_GC"))
rownames(met_nat) <- met_nat$Sample
met_nat$Sample <- NULL
met_nat$paper <- NULL

natcount <- counts[,rownames(met_nat)]
all(rownames(met_nat) == colnames(natcount)) 

natkeep <- rowSums(natcount) > 20
natkept <- natcount[natkeep,]

ddsnat <- DESeqDataSetFromMatrix(countData = round(natkept),
  colData = met_nat,
  design = ~ Group)

ddsnat <- DESeq(ddsnat, betaPrior = FALSE)

normcounts_nat<- counts(ddsnat, normalized=TRUE)
normcounts_nat<- as.data.frame(normcounts_nat)
cyto_nat <- subset(normcounts_nat, rownames(normcounts_nat) %in% cytos)

cyto_nat$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cyto_nat), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(cyto_nat) <- cyto_nat$symbol
cyto_nat$symbol <- NULL

pheatmap(cyto_nat, cluster_rows=F, cluster_cols=F, annotation_col = met_nat,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= 3,
  cellheight=20, cellwidth = 35)

##%######################################################%##
#                                                          #
####                      Panel G                       ####
#                                                          #
##%######################################################%##

res_nat<- results(ddsnat, contrast =c("Group", "WT_FLT","WT_GC"), alpha = 0.05, lfcThreshold = 0.32)

res_nat$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(res_nat), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

EnhancedVolcano(toptable = res_nat,               
  x = "log2FoldChange",           
  y = "pvalue",                     
  lab = res_nat$symbol,
  xlim = c(-10, 10),
  labSize = 6,
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

##%######################################################%##
#                                                          #
####                      Panel I                       ####
#                                                          #
##%######################################################%##

counts_osd <- read.csv("OSD_576.csv", stringsAsFactors = T, header = T, check.names = F)
metadata_osd <- read.csv("metadata_OSD576.csv", header = T, stringsAsFactors = T)

meta_osd<- with(metadata_osd, metadata_osd[order(Condition, Sample),])

rownames(counts_osd) <- counts_osd$gene_id
counts_osd$gene_id <- NULL

rownames(meta_osd) <- meta_osd$Sample
meta_osd$Sample <- NULL

all(rownames(meta_osd) == colnames(counts_osd))
match(colnames(counts_osd), rownames(meta_osd))
idxosd <- match(rownames(meta_osd), colnames(counts_osd))
reorder_counts_osd <- counts_osd[,idxosd]
all(rownames(meta_osd) == colnames(reorder_counts_osd))

# Filtering
keep_osd <- rowSums(reorder_counts_osd) > 20
kept_osd <- reorder_counts_osd[keep_osd,]

dds_osd <- DESeqDataSetFromMatrix(countData = round(kept_osd),
  colData = meta_osd,
  design = ~ Condition)

dds_osd <- DESeq(dds_osd, betaPrior = FALSE)
res_osd<- results(dds_osd, contrast =c("Condition", "Flight","GC"), alpha = 0.05, lfcThreshold = 0.32)

normcounts_osd<- counts(dds_osd, normalized=TRUE)
normcounts_osd <- as.data.frame(normcounts_osd)

cyto_osd <- subset(normcounts_osd, rownames(normcounts_osd) %in% cytos)

cyto_osd$symbol <- mapIds(org.Mm.eg.db,
  keys=rownames(cyto_osd), 
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals= "first")

rownames(cyto_osd) <- cyto_osd$symbol
cyto_osd$symbol <- NULL

pheat_osd <- cyto_osd[,1:18]
met_pheat <- meta_osd[1:18, ]
met_pheat <- as.data.frame(met_pheat)
rownames(met_pheat) <- colnames(pheat_osd)

pheatmap(pheat_osd, cluster_rows=F, cluster_cols=F, annotation_col = met_pheat,
  scale = "row", 
  color = colorRampPalette(c("navy", "white", "#B33434"))(50), 
  size = 13, gaps_col= 9,
  cellheight=20, cellwidth = 35)



