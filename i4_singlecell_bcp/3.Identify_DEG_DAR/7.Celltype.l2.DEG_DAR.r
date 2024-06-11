library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)

library(data.table)

temp.sobj <- readRDS("path_to_Seurat_object")

Idents(temp.sobj) <- "celltype.l2"


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

Idents(temp.sobj) <- "predicted.id"




for (i in cell_list){
    x <- subset(temp.sobj, idents = paste0(i))
    assign(paste0(i, '.sobj'), x)
}











cell_list <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD8_Naive", "CD8_TCM",
               "CD8_TEM", "dnT", "MAIT", "gdT", "Treg", "NK", "NK_CD56bright",
               "B_Naive", "B_intermediate", "B_memory", "CD14_Mono", "CD16_Mono", "cDC2", "pDC",
               "HSPC", "Platelet")

for (i in cell_list){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "timepoint"
    DefaultAssay(x) <- "SCT"
    
   
    y1 = FindMarkers(x, ident.1 = "September Post-launch", 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    y8 = FindMarkers(x, ident.1 = c("November 2021", "December 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)

    assign(paste0("fcs_", i, '_1'), y1)
    assign(paste0("fcs_", i, '_8'), y8)
    
    write.csv(y1, file = paste0("pval/update/celltype.l2/differential/all.list/DEGs_",i,"_September_Post_JAS.csv"))
    write.csv(y8, file = paste0("pval/update/celltype.l2/differential/all.list/DEGs_",i,"_Nov+Dec_JAS.csv"))
    

    
    
}



sessionInfo()

rm(list=ls())






