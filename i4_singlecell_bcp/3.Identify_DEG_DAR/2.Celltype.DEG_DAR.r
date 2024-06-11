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

Idents(temp.sobj) <- "celltype"


temp.sobj$predicted.id <- ifelse(temp.sobj$celltype == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype == 'CD14 Mono', 'CD14_Mono',
                                                                 ifelse(temp.sobj$celltype == 'CD16 Mono', 'CD16_Mono',
                                                                         ifelse(temp.sobj$celltype == 'DC', 'DC','other'))))))))





cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "CD14_Mono", "CD16_Mono", "DC", "other")

Idents(temp.sobj) <- "predicted.id"


cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "CD14_Mono", "CD16_Mono", "DC", "other")

for (i in cell_list){
    x <- subset(temp.sobj, idents = paste0(i))
    assign(paste0(i, '.sobj'), x)
}







all_list <- c('pbmc', cell_list)

pbmc.sobj <- temp.sobj

rm(list=c('temp.sobj'))



for (i in c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B', 'NK', 'CD14_Mono', 'CD16_Mono', 'DC', 'other')){
    
    
    y1 <- read.csv(file = paste0("pval/update/all.list/update/DEGs_",i,"_R+45_preflight.csv"))
    y8 <- read.csv(file = paste0("pval/update/all.list/update/DEGs_",i,"_R+82_preflight.csv"))
    
    
    write.csv(y1, paste0("pval/update/all.list/update/DEGs_",i,"_RP1.csv"))
    write.csv(y8, paste0("pval/update/all.list/update/DEGs_",i,"_RP2.csv"))
    
    
    
    
    
    
}







for (i in c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B', 'NK', 'CD14_Mono', 'CD16_Mono', 'DC', 'other')){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "timepoint"
    DefaultAssay(x) <- "SCT"
    
   
    y1 = FindMarkers(x, ident.1 = "November 2021", 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    y8 = FindMarkers(x, ident.1 = "December 2021", 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y9 = FindMarkers(x, ident.1 = c('November 2021', "December 2021"), 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)

    assign(paste0("fcs_", i, '_1'), y1)
    assign(paste0("fcs_", i, '_8'), y8)
    assign(paste0("fcs_", i, '_9'), y9)
    
    write.csv(y1, file = paste0("pval/update/all.list/update/DEGs_",i,"_RP1.csv"))
    write.csv(y8, file = paste0("pval/update/all.list/update/DEGs_",i,"_RP2.csv"))
    write.csv(y9, file = paste0("pval/update/all.list/update/DEGs_",i,"_RP3.csv"))
    
    
    
    
    
    
}



for (i in c('CD8_T', 'other_T', 'B', 'NK', 'CD14_Mono', 'CD16_Mono', 'DC', 'other')){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "timepoint"
    DefaultAssay(x) <- "SCT"
    
   
    y1 = FindMarkers(x, ident.1 = "November 2021", 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    y8 = FindMarkers(x, ident.1 = "December 2021", 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y9 = FindMarkers(x, ident.1 = c('November 2021', "December 2021"), 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)

    assign(paste0("fcs_", i, '_1'), y1)
    assign(paste0("fcs_", i, '_8'), y8)
    assign(paste0("fcs_", i, '_9'), y9)
    
    write.csv(y1, file = paste0("pval/update/all.list/update/DEGs_",i,"_RP1.csv"))
    write.csv(y8, file = paste0("pval/update/all.list/update/DEGs_",i,"_RP2.csv"))
    write.csv(y9, file = paste0("pval/update/all.list/update/DEGs_",i,"_RP3.csv"))
    
    
    
    
    
    
}









for (i in c('other', 'DC', 'CD16_Mono','CD14_Mono', 'NK', 'B', 'other_T', 'CD8_T', 'CD4_T', 'pbmc')){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "timepoint"
    DefaultAssay(x) <- "peaks"
    
   
    y1 = FindMarkers(x, ident.1 = "November 2021", 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    y8 = FindMarkers(x, ident.1 = "December 2021", 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    

    assign(paste0("fcs_", i, '_1.dar'), y1)
    assign(paste0("fcs_", i, '_8.dar'), y8)
    
    write.csv(y1, file = paste0("pval/update/all.list/update/DARs_",i,"_RP1.csv"))
    write.csv(y8, file = paste0("pval/update/all.list/update/DARs_",i,"_RP2.csv"))
    
    
    
    
    
    
    
}



for (i in c('other', 'DC', 'CD16_Mono','CD14_Mono', 'NK', 'B', 'other_T', 'CD8_T', 'CD4_T', 'pbmc')){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "timepoint"
    DefaultAssay(x) <- "peaks"
    
   
    
    
    y9 = FindMarkers(x, ident.1 = c('November 2021', "December 2021"), 
                     ident.2 = c("September Post-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)

    
    assign(paste0("fcs_", i, '_9.dar'), y9)
    
    write.csv(y9, file = paste0("pval/update/all.list/update/DARs_",i,"_RP3.csv"))
    
    
    
    
    
    
}







sessionInfo()

rm(list=ls())






