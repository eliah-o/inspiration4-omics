library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)


library(patchwork)



library(TFBSTools)

library(JASPAR2020)
set.seed(1234)


temp.sobj <- readRDS("path_to_Seurat_object")

temp.sobj$predicted.id <- ifelse(temp.sobj$celltype.l2 == 'CD4 Naive', 'CD4_Naive', 
                              ifelse(temp.sobj$celltype.l2 == 'CD4 TCM', 'CD4_TCM', 
                                     ifelse(temp.sobj$celltype.l2 == 'CD4 TEM', 'CD4_TEM', 
                                            ifelse(temp.sobj$celltype.l2 == 'CD4 Proliferating', 'CD4_Proliferating', 
                                                   ifelse(temp.sobj$celltype.l2 == 'CD4 CTL', 'CD4_CTL',
                                                          ifelse(temp.sobj$celltype.l2 == 'CD8 Naive', 'CD8_Naive',
                                                                ifelse(temp.sobj$celltype.l2 == 'CD8 TCM', 'CD8_TCM',
                                                                       ifelse(temp.sobj$celltype.l2 == 'CD8 Proliferating', 'CD8_Proliferating',
                                                                             ifelse(temp.sobj$celltype.l2 == 'dnT', 'dnT',
                                                                                   ifelse(temp.sobj$celltype.l2 == 'MAIT', 'MAIT', 
                                                                                         ifelse(temp.sobj$celltype.l2 == 'gdT', 'gdT',
                                                                                               ifelse(temp.sobj$celltype.l2 == 'Treg', 'Treg',
                                                                                                     ifelse(temp.sobj$celltype.l2 == 'NK', 'NK',
                                                                                                           ifelse(temp.sobj$celltype.l2 == 'NK_CD56bright', 'NK_CD56bright',
                                                                                                                  ifelse(temp.sobj$celltype.l2 == 'NK Proliferating', 'NK_Proliferating',
                                                                                                                         ifelse(temp.sobj$celltype.l2 == 'B naive', 'B_Naive',
                                                                                                                               ifelse(temp.sobj$celltype.l2 == 'B intermediate', 'B_intermediate',
                                                                                                                                     ifelse(temp.sobj$celltype.l2 == 'B memory', 'B_memory',
                                                                                                                                           ifelse(temp.sobj$celltype.l2 == 'CD14 Mono', 'CD14_Mono',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'CD16 Mono', 'CD16_Mono',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'cDC1', 'cDC1',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'cDC2', 'cDC2',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'pDC', 'pDC',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'ASDC', 'ASDC',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'Plasmablast', 'Plasmablast',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'HSPC', 'HSPC',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'Eryth', 'Eryth',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'ASDC', 'ASDC',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'ILC', 'ILC',
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'Platelet', 'Platelet', 
                                                                                                                                                 ifelse(temp.sobj$celltype.l2 == 'CD8 TEM', 'CD8_TEM', 'Doublet')))))))))))))))))))))))))))))))





cell_list <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD4_Proliferating", "CD4_CTL", "CD8_Naive", "CD8_TCM",
               "CD8_TEM", "CD8_Proliferating", "dnT", "MAIT", "gdT", "Treg", "NK", "NK_CD56bright", "NK_Proliferating",
               "B_Naive", "B_intermediate", "B_memory", "CD14_Mono", "CD16_Mono", "cDC1", "cDC2", "pDC",
               "ASDC", "Plasmablast", "HSPC", "Eryth", "ILC", "Platelet", "Doublet")

Idents(temp.sobj) <- 'predicted.id'

for (i in cell_list){
    x <- subset(temp.sobj, ident = i)
    assign(paste0(i,".sobj"), x)
}





for (i in cell_list){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_November_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    
    write.csv(a1, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    
    
    
}

cell_list2 <- c("CD8_Naive", "CD8_TCM",
               "CD8_TEM", "CD8_Proliferating", "dnT", "MAIT", "gdT", "Treg", "NK", "NK_CD56bright", "NK_Proliferating",
               "B_Naive", "B_intermediate", "B_memory", "CD14_Mono", "CD16_Mono", "cDC1", "cDC2", "pDC",
               "ASDC", "Plasmablast", "HSPC", "Eryth", "ILC", "Platelet", "Doublet")

for (i in cell_list2){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_November_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    
    write.csv(a1, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    
    
    
}

cell_list3 <- c("dnT", "MAIT", "gdT", "Treg", "NK", "NK_CD56bright", "NK_Proliferating",
               "B_Naive", "B_intermediate", "B_memory", "CD14_Mono", "CD16_Mono", "cDC1", "cDC2", "pDC",
               "ASDC", "Plasmablast", "HSPC", "Eryth", "ILC", "Platelet", "Doublet")

for (i in cell_list3){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_November_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    
    write.csv(a1, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    
    
    
}



cell_list4 <- c("B_Naive", "B_intermediate", "B_memory", "CD14_Mono", "CD16_Mono", "cDC1", "cDC2", "pDC",
               "ASDC", "Plasmablast", "HSPC", "Eryth", "ILC", "Platelet", "Doublet")

for (i in cell_list4){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_November_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    
    write.csv(a1, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    
    
    
}



cell_list5 <- c("cDC2", "pDC",
               "ASDC", "Plasmablast", "HSPC", "Eryth", "ILC", "Platelet", "Doublet")

for (i in cell_list5){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_November_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    
    write.csv(a1, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    
    
    
}



cell_list6 <- c("HSPC", "Eryth", "ILC", "Platelet", "Doublet")

for (i in cell_list6){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_November_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    
    write.csv(a1, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    
    
    
}



cell_list7 <- c("Platelet", "Doublet")

for (i in cell_list7){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_November_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    
    write.csv(a1, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_15_celltype.l2/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    
    
    
}







cell_list.update <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD8_Naive", "CD8_TCM",
               "CD8_TEM", "dnT", "MAIT", "gdT", "Treg", "NK", "NK_CD56bright", 
               "B_Naive", "B_intermediate", "B_memory", "CD14_Mono", "CD16_Mono", "cDC2", "pDC",
               "HSPC", "ILC", "Platelet")

for (i in cell_list.update){
    x <- get(paste0(i, ".sobj"))
    y <- get(paste0(i, ".enriched.motifs.Sept_post"))
    
    options(repr.plot.width=16, repr.plot.height=8)
    
    a <- MotifPlot(
      object = x,
      motifs = y$motif[1:12]
)+ 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size=15),
       strip.text = element_text(size =30))
    
    
    ggsave(paste0("2023_08_15_celltype.l2/motif.image/", i, ".tiff"), width=16, height=8, dpi=300, compression='lzw')
    ggsave(paste0("2023_08_15_celltype.l2/motif.image/", i, ".png"), width=16, height=8, dpi=300)
    ggsave(paste0("2023_08_15_celltype.l2/motif.image/", i, ".svg"), width=16, height=8, dpi=300)
    ggsave(paste0("2023_08_15_celltype.l2/motif.image/", i, ".pdf"), width=16, height=8, dpi=300)
    
    print(a)
}







sessionInfo()

rm(list=ls())










