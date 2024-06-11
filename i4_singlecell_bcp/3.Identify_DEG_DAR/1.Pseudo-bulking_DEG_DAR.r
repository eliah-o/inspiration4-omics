library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)

library(data.table)

temp.sobj <- readRDS("path_to_seurat_objct")

Idents(temp.sobj) <- "celltype"


temp.sobj$predicted.id <- ifelse(temp.sobj$celltype == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype == 'CD14 Mono', 'Mono',
                                                                 ifelse(temp.sobj$celltype == 'CD16 Mono', 'Mono',
                                                                         ifelse(temp.sobj$celltype == 'DC', 'DC','other'))))))))





temp.sobj$Timepoint <- ifelse(temp.sobj$timepoint == 'September Post-launch', 'R+1', 
                              ifelse(temp.sobj$timepoint == 'November 2021', 'R+45', 
                                     ifelse(temp.sobj$timepoint == 'December 2021', 'R+82', 'Pre-flight')))





cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "Mono", "DC", "other")

Idents(temp.sobj) <- "predicted.id"


for (i in cell_list){
    x <- subset(temp.sobj, idents = paste0(i))
    assign(paste0(i, '.sobj'), x)
}





all_list <- c('pbmc', cell_list)

pbmc.sobj <- temp.sobj

rm(list=c('temp.sobj'))



pseudo_B <- AggregateExpression(B.sobj, assays = "RNA", 
                                   return.seurat = TRUE,
                                 
                                   group.by = c("ID", "timepoint"))

pseudo_DC <- AggregateExpression(DC.sobj, assays = "RNA", 
                                   return.seurat = TRUE, slot = 'counts',
                                 
                                   group.by = c("ID", "timepoint"))

for (i in all_list){
    x <- get(paste0(i, '.sobj'))
    
    y <- AggregateExpression(x, assays = "RNA", 
                                   return.seurat = TRUE, slot = 'counts',
                                 
                                   group.by = c("ID", "timepoint"))
    
    y$timepoint <- sapply(strsplit(rownames(y[[]]), split='_', fixed=TRUE), `[`, 2)
    
    assign(paste0(i, ".pseudo"), y)
    
    
    
}



for (i in all_list){
    x <- get(paste0(i, ".pseudo"))
    
    Idents(x) <- 'timepoint'
    
    y1 = FindMarkers(x, ident.1 = c("September Post-launch"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y2 = FindMarkers(x, ident.1 = c("November 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y3 = FindMarkers(x, ident.1 = c("December 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y4 = FindMarkers(x, ident.1 = c('November 2021', "December 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y5 = FindMarkers(x, ident.1 = c('September Post-launch', 'November 2021', "December 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y6 = FindMarkers(x, ident.1 = c('September Post-launch'), 
                     ident.2 = c("August 2021"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y7 = FindMarkers(x, ident.1 = c('November 2021'), 
                     ident.2 = c("September Post-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y8 = FindMarkers(x, ident.1 = c('December 2021'), 
                     ident.2 = c("September Post-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y9 = FindMarkers(x, ident.1 = c('November 2021', 'December 2021'), 
                     ident.2 = c("September Post-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    write.csv(y1, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.FP1.csv'))
    write.csv(y2, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.LP1.csv'))
    write.csv(y3, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.LP2.csv'))
    write.csv(y4, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.LP3.csv'))
    write.csv(y5, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.FP2.csv'))
    write.csv(y6, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.FP4.csv'))
    write.csv(y7, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.RP1.csv'))
    write.csv(y8, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.RP2.csv'))
    write.csv(y9, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DEGs.', i, '.FP3.csv'))
    
    
    
}



for (i in all_list){
    x <- get(paste0(i, '.sobj'))
    
    y <- AggregateExpression(x, assays = "peaks", 
                                   return.seurat = TRUE, slot = 'counts',
                                 
                                   group.by = c("ID", "timepoint"))
    
    y$timepoint <- sapply(strsplit(rownames(y[[]]), split='_', fixed=TRUE), `[`, 2)
    
    assign(paste0(i, ".pseudo"), y)
    
    
    
}

for (i in all_list){
    x <- get(paste0(i, ".pseudo"))
    
    Idents(x) <- 'timepoint'
    
    y1 = FindMarkers(x, ident.1 = c("September Post-launch"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y2 = FindMarkers(x, ident.1 = c("November 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y3 = FindMarkers(x, ident.1 = c("December 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y4 = FindMarkers(x, ident.1 = c('November 2021', "December 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y5 = FindMarkers(x, ident.1 = c('September Post-launch', 'November 2021', "December 2021"), 
                     ident.2 = c("June 2021", "August 2021", "September Pre-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y6 = FindMarkers(x, ident.1 = c('September Post-launch'), 
                     ident.2 = c("August 2021"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y7 = FindMarkers(x, ident.1 = c('November 2021'), 
                     ident.2 = c("September Post-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y8 = FindMarkers(x, ident.1 = c('December 2021'), 
                     ident.2 = c("September Post-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    y9 = FindMarkers(x, ident.1 = c('November 2021', 'December 2021'), 
                     ident.2 = c("September Post-launch"),
                         test.use = "DESeq2", 
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.cells.feature = 1,
                     min.cells.group = 1)
    
    write.csv(y1, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.FP1.csv'))
    write.csv(y2, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.LP1.csv'))
    write.csv(y3, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.LP2.csv'))
    write.csv(y4, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.LP3.csv'))
    write.csv(y5, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.FP2.csv'))
    write.csv(y6, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.FP4.csv'))
    write.csv(y7, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.RP1.csv'))
    write.csv(y8, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.RP2.csv'))
    write.csv(y9, paste0('pval/update/all.list/update/DESeq2_pseudobulk/DARs.', i, '.FP3.csv'))
    
    
    
}



sessionInfo()

rm(list=ls())






