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




table(temp.sobj$orig.ident, temp.sobj$timepoint)



temp.sobj$predicted.id <- ifelse(temp.sobj$celltype == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype == 'CD14 Mono', 'CD14_Mono',
                                                                 ifelse(temp.sobj$celltype == 'CD16 Mono', 'CD16_Mono',
                                                                         ifelse(temp.sobj$celltype == 'DC', 'DC','other'))))))))



cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "CD14_Mono", "CD16_Mono", "DC", "other")





T.sobj <- subset(temp.sobj, idents = c("CD4 T", "CD8 T", "other T") )



B.sobj <- subset(temp.sobj, idents = c("B") )



Mono.sobj <- subset(temp.sobj, idents = c("CD14 Mono", "CD16 Mono") )

NK.sobj <- subset(temp.sobj, idents = c("NK") )

DC.sobj <- subset(temp.sobj, idents = c("DC") )





cell_list = c("T", "B", "Mono", "NK", "DC")



for (i in cell_list){
    x <- get(paste0(i, ".sobj"))
    DefaultAssay(x) <- "SCT"
    Idents(x) <- "orig.ident"
    
    y <- as.data.frame(AverageExpression(x, verbose=TRUE)$SCT)
    
    y$gene <- row.names(y)

    write.csv(y, paste0("Resources_paper/update/", i, ".average.expression.SCT.orig.ident.csv"))
    
    Idents(x) <- "timepoint"    
    z <- as.data.frame(AverageExpression(x, verbose=TRUE)$SCT)
    
    z$gene <- row.names(z)
    
    
}

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
    
    
    
    
}







sessionInfo()

rm(list=ls())




