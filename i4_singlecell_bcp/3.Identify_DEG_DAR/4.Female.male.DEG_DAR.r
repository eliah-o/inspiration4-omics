library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)

library(data.table)

temp.sobj <- readRDS("path_to_seurat_object")

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



temp.sobj$gender <- ifelse(temp.sobj$ID == 'C001', 'M', 
                              ifelse(temp.sobj$ID == 'C004', 'M', 
                                     ifelse(temp.sobj$ID == 'C002', 'F', 'F')))



Idents(temp.sobj) <- "gender"


gender <- c("F", "M")

for (i in gender){
    x <- subset(temp.sobj, idents = paste0(i))
    assign(paste0(i, '.sobj'), x)
}





Idents(F.sobj) <- "predicted.id"

Idents(M.sobj) <- "predicted.id"


cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "CD14_Mono", "CD16_Mono", "DC", "other")

for (i in cell_list){
    
    F <- get("F.sobj")
    M <- get("M.sobj")
    
    f <- subset(F, idents = paste0(i))
    assign(paste0(i, '.F.sobj'), f)
    
    m <- subset(M, idents = paste0(i))
    assign(paste0(i, '.M.sobj'), m)
}





pbmc.F.sobj <- F.sobj
pbmc.M.sobj <- M.sobj



all_list <- c('pbmc', cell_list)

for (i in c('pbmc')){
    F <- get(paste0(i, '.F.sobj'))
    M <- get(paste0(i, '.M.sobj'))
    
    
    
    
    Idents(F) <- "timepoint"
    DefaultAssay(F) <- "SCT"
    
    Idents(M) <- "timepoint"
    DefaultAssay(M) <- "SCT"
    
    
    f = FindMarkers(F, 
                    ident.1 = c("September Post-launch"), 
                    ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(f, paste0("gender/2023_11_30/", i, ".F.DEG.all.csv"))
    
    
    m = FindMarkers(M, 
                    ident.1 = c("September Post-launch"), 
                    ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(m, paste0("gender/2023_11_30/", i, ".M.DEG.all.csv"))
}



for (i in cell_list){
    F <- get(paste0(i, '.F.sobj'))
    M <- get(paste0(i, '.M.sobj'))
    
    
    
    
    Idents(F) <- "timepoint"
    DefaultAssay(F) <- "SCT"
    
    Idents(M) <- "timepoint"
    DefaultAssay(M) <- "SCT"
    
    
    f = FindMarkers(F, 
                    ident.1 = c("September Post-launch"), 
                    ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(f, paste0("gender/2023_11_30/", i, ".F.DEG.all.csv"))
    
    
    m = FindMarkers(M, 
                    ident.1 = c("September Post-launch"), 
                    ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(m, paste0("gender/2023_11_30/", i, ".M.DEG.all.csv"))
}







for (i in all_list){
    F <- get(paste0(i, '.F.sobj'))
    M <- get(paste0(i, '.M.sobj'))
    
    
    
    
    Idents(F) <- "timepoint"
    DefaultAssay(F) <- "SCT"
    
    Idents(M) <- "timepoint"
    DefaultAssay(M) <- "SCT"
    
    
    f = FindMarkers(F, 
                    ident.1 = c("November 2021", "December 2021"), 
                    ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    write.csv(f, paste0("gender/2023_11_30/", i, ".F.DEG.longterm.all.csv"))
    
    
    m = FindMarkers(M, 
                    ident.1 = c("November 2021", "December 2021"), 
                    ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(m, paste0("gender/2023_11_30/", i, ".M.DEG.longterm.all.csv"))
}















sessionInfo()

rm(list=ls())

ls()




