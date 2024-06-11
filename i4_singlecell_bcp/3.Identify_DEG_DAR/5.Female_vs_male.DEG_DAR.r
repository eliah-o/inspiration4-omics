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



temp.sobj$gender <- ifelse(temp.sobj$ID == 'C001', 'M', 
                              ifelse(temp.sobj$ID == 'C004', 'M', 
                                     ifelse(temp.sobj$ID == 'C002', 'F', 'F')))



Idents(temp.sobj) <- "gender"


gender <- c("F", "M")



Idents(temp.sobj) <- 'timepoint'
temp.immediate.sobj <- subset(temp.sobj, idents = c("September Post-launch"))
temp.longterm.sobj <- subset(temp.sobj, idents = c("November 2021", "December 2021"))

pbmc.immediate.sobj <- temp.immediate.sobj
pbmc.longterm.sobj <- temp.longterm.sobj





Idents(temp.sobj) <- 'timepoint'

pbmc.tp5.sobj <- subset(temp.sobj, idents = c("November 2021"))
pbmc.tp6.sobj <- subset(temp.sobj, idents = c("December 2021"))





Idents(temp.sobj) <- "predicted.id"


cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "CD14_Mono", "CD16_Mono", "DC", "other")

all_list <- c('pbmc', cell_list)

pbmc.sobj <- temp.sobj



for (i in cell_list){
    
    F <- get("temp.immediate.sobj")
    M <- get("temp.longterm.sobj")
    
    Idents(F) <- 'predicted.id'
    Idents(M) <- 'predicted.id'
    
    f <- subset(F, idents = paste0(i))
    assign(paste0(i, '.immediate.sobj'), f)
    
    m <- subset(M, idents = paste0(i))
    assign(paste0(i, '.longterm.sobj'), m)
}



for (i in cell_list){
    
    a <- get("pbmc.tp5.sobj")
    b <- get("pbmc.tp6.sobj")
    
    Idents(a) <- 'predicted.id'
    Idents(b) <- 'predicted.id'
    
    A <- subset(a, idents = paste0(i))
    assign(paste0(i, '.tp5.sobj'), a)
    
    B <- subset(b, idents = paste0(i))
    assign(paste0(i, '.tp6.sobj'), b)
}





for (i in all_list){
    F <- get(paste0(i, '.immediate.sobj'))
    M <- get(paste0(i, '.longterm.sobj'))
    
    Idents(F) <- "gender"
    DefaultAssay(F) <- "SCT"
    
    Idents(M) <- "gender"
    DefaultAssay(M) <- "SCT"
    
    
    f = FindMarkers(F, 
                    ident.1 = c("F"), 
                    ident.2 = c("M"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(f, paste0("gender/2023_11_30/F_to_M/", i, ".immediate.F_to_M.DEG.csv"))
    
    
    m = FindMarkers(M, 
                    ident.1 = c("F"), 
                    ident.2 = c("M"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(m, paste0("gender/2023_11_30/F_to_M/", i, ".longterm.F_to_M.DEG.filter.csv"))
}











for (i in all_list){
    F <- get(paste0(i, '.tp5.sobj'))
    M <- get(paste0(i, '.tp6.sobj'))
    
    Idents(F) <- "gender"
    DefaultAssay(F) <- "SCT"
    
    Idents(M) <- "gender"
    DefaultAssay(M) <- "SCT"
    
    
    f = FindMarkers(F, 
                    ident.1 = c("F"), 
                    ident.2 = c("M"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(f, paste0("gender/2023_11_30/F_to_M/", i, ".tp5.F_to_M.DEG.filter.csv"))
    
    
    m = FindMarkers(M, 
                    ident.1 = c("F"), 
                    ident.2 = c("M"),
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    write.csv(m, paste0("gender/2023_11_30/F_to_M/", i, ".tp6.F_to_M.DEG.filter.csv"))
}










