library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)

library(BiocParallel)
register(SerialParam())
options(Seurat.object.assay.version = 'v5')

library(data.table)



temp.sobj <- readRDS("path_to_seurat_object")

temp.sobj$predicted.id <- ifelse(temp.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(temp.sobj$celltype.l1 == 'DC', 'DC','other')))))))



cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "Mono", "DC", "other")

for (i in cell_list){
    Idents(temp.sobj) <- 'predicted.id'
    x <- subset(temp.sobj, idents = paste0(i))
    assign(paste0(i, '.sobj'), x)
}



all_list <- c('pbmc', cell_list)
pbmc.sobj <- temp.sobj





for (i in all_list){
    x <- get(paste0(i, '.sobj'))
    
    y <- AggregateExpression(DC.sobj, assays = "RNA", 
                                   return.seurat = TRUE, slot = 'counts',
                                 
                                   group.by = c("orig.ident", 'Day'))
    
    y$orig.ident<- sapply(strsplit(rownames(y[[]]), split='_', fixed=TRUE), `[`, 1)
    
    
    assign(paste0(i, ".pseudo"), y)
    
    
    
}





for (i in all_list){
    x <- get(paste0(i, '.pseudo'))
    Idents(x) <- "orig.ident"
    
    y1 = FindMarkers(x, ident.1 = c("Day2-Rep1", "Day2-Rep2", "Day2-Rep3"), 
                     ident.2 = c("Day1-Rep1", "Day1-Rep2", "Day1-Rep3"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    
}



for (i in all_list){
    x <- get(paste0(i, '.sobj'))
    
    y <- AggregateExpression(DC.sobj, assays = "ATAC", 
                                   return.seurat = TRUE, slot = 'counts',
                                 
                                   group.by = c("orig.ident", 'Day'))
    
    y$orig.ident<- sapply(strsplit(rownames(y[[]]), split='_', fixed=TRUE), `[`, 1)
    
    
    assign(paste0(i, ".pseudo"), y)
    
    
    
}



for (i in all_list){
    x <- get(paste0(i, '.pseudo'))
    Idents(x) <- "orig.ident"
    
    y1 = FindMarkers(x, ident.1 = c("Day2-Rep1", "Day2-Rep2", "Day2-Rep3"), 
                     ident.2 = c("Day1-Rep1", "Day1-Rep2", "Day1-Rep3"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    
    
    
    
}





sessionInfo()

rm(list=ls())




