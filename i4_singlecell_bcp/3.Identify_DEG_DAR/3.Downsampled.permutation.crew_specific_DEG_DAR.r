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

temp.sobj$celltype.l1.label <- ifelse(temp.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(temp.sobj$celltype.l1 == 'DC', 'DC','other')))))))



cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", "Mono", "DC", "other")



Idents(temp.sobj) <- "celltype.l1.label"

for (i in cell_list){
    x <- subset(temp.sobj, idents = paste0(i))
    assign(paste0(i, '.sobj'), x)
}




pbmc.sobj <- temp.sobj

all_list <- c('pbmc', cell_list)



for (i in c('C001', 'C002', 'C003', 'C004')){
    x <- readRDS(paste0("../../../Motif_analysis/2022_10_05/", i, ".motif.chromvar.TP.downsample.1000.sobj.rds"))
    
    assign(paste0(i, ".sobj"), x)
}







C001.sobj$celltype.l1.label <- ifelse(C001.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(C001.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(C001.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(C001.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(C001.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(C001.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(C001.sobj$celltype.l1 == 'DC', 'DC','other')))))))



C002.sobj$celltype.l1.label <- ifelse(C002.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(C002.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(C002.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(C002.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(C002.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(C002.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(C002.sobj$celltype.l1 == 'DC', 'DC','other')))))))



C003.sobj$celltype.l1.label <- ifelse(C003.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(C003.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(C003.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(C003.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(C003.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(C003.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(C003.sobj$celltype.l1 == 'DC', 'DC','other')))))))



C004.sobj$celltype.l1.label <- ifelse(C004.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(C004.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(C004.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(C004.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(C004.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(C004.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(C004.sobj$celltype.l1 == 'DC', 'DC','other')))))))



for (j in c('C001', 'C002', 'C003', 'C004')){
    x1 <- get(paste0(j, '.sobj'))
    Idents(x1) <- 'celltype.l1.label'
    
    for (i in cell_list){
        y1 <- subset(x1, idents = paste0(i))
        assign(paste0(j, ".", i, '.sobj'), y1)
}}



C001.pbmc.sobj <- C001.sobj
C002.pbmc.sobj <- C002.sobj
C003.pbmc.sobj <- C003.sobj
C004.pbmc.sobj <- C004.sobj







for (i in c('B', 'CD8_T', 'DC', 'Mono', 'NK', 'other', 'other_T')){
    
    
    
    x <- readRDS(paste0("../../../Motif_analysis/2022_10_05/", i, ".down.max.sobj.rds"))
    
    
    assign(paste0(i, ".down.sobj"), x)
}



for (i in c('B', 'CD8_T', 'DC', 'Mono', 'NK', 'other', 'other_T')){
    
    y1 <- readRDS(paste0("../../../Motif_analysis/2022_10_05/C001.", i, ".down.max.sobj.rds"))
    y2 <- readRDS(paste0("../../../Motif_analysis/2022_10_05/C002.", i, ".down.max.sobj.rds"))
    y3 <- readRDS(paste0("../../../Motif_analysis/2022_10_05/C003.", i, ".down.max.sobj.rds"))
    y4 <- readRDS(paste0("../../../Motif_analysis/2022_10_05/C004.", i, ".down.max.sobj.rds"))
    
    assign(paste0("C001.", i, ".down.sobj"), y1)
    assign(paste0("C002.", i, ".down.sobj"), y2)
    assign(paste0("C003.", i, ".down.sobj"), y3)
    assign(paste0("C004.", i, ".down.sobj"), y4)
}



for (i in c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B', 'NK', 'Mono', 'DC', 'other')){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "TP"
    DefaultAssay(x) <- "SCT"
    
   
    y1 = FindMarkers(x, ident.1 = "R+1", 
                     ident.2 = c("Pre-flight"))
    y8 = FindMarkers(x, ident.1 = c("R+45&R+82"), 
                     ident.2 = c("Pre-flight"))

    assign(paste0("fcs_", i, '_1'), y1)
    assign(paste0("fcs_", i, '_8'), y8)
    
    write.csv(y1, file = paste0("pval/update/all.list/update/downsample/optimize/DEGs_",i,"_R+1_preflight.filtered.csv"))
    write.csv(y8, file = paste0("pval/update/all.list/update/downsample/optimize/DEGs_",i,"_R+45&R+82_preflight.filtered.csv"))
    
    
    
    
    
    
}





for (i in c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B', 'NK',  'Mono', 'DC', 'other')){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "TP"
    DefaultAssay(x) <- "peaks"
    
   
    y1 = FindMarkers(x, ident.1 = "R+1", 
                     ident.2 = c("Pre-flight"))
    y8 = FindMarkers(x, ident.1 = c("R+45&R+82"), 
                     ident.2 = c("Pre-flight"))

    assign(paste0("fcs_", i, '_1'), y1)
    assign(paste0("fcs_", i, '_8'), y8)
    
    write.csv(y1, file = paste0("pval/update/all.list/update/downsample/optimize/DARs_",i,"_R+1_preflight.filtered.csv"))
    write.csv(y8, file = paste0("pval/update/all.list/update/downsample/optimize/DARs_",i,"_R+45&R+82_preflight.filtered.csv"))
    
    
    
    
    
    
}





for (j in c('C001', 'C002', 'C003', 'C004')){
    
    

    for (i in c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B', 'NK', 'Mono', 'DC', 'other')){
        x <- get(paste0(j, ".", i, '.sobj'))
        Idents(x) <- "TP"
        DefaultAssay(x) <- "SCT"
    
   
        y1 = FindMarkers(x, ident.1 = "R+1", 
                         ident.2 = c("Pre-flight"))
        y8 = FindMarkers(x, ident.1 = c("R+45&R+82"), 
                         ident.2 = c("Pre-flight"))

        assign(paste0("fcs_", i, '_1'), y1)
        assign(paste0("fcs_", i, '_8'), y8)
    
        write.csv(y1, file = paste0("pval/update/all.list/update/downsample/optimize/sample.specific/DEGs_",j, ".", i,"_R+1_preflight.filtered.csv"))
        write.csv(y8, file = paste0("pval/update/all.list/update/downsample/optimize/sample.specific/DEGs_",j, ".", i,"_R+45&R+82_preflight.filtered.csv"))
    
    
    
    
    
    
}}





for (j in c('C001', 'C002', 'C003', 'C004')){
    
    

    for (i in c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B', 'NK', 'Mono', 'DC', 'other')){
        x <- get(paste0(j, ".", i, '.sobj'))
        Idents(x) <- "TP"
        DefaultAssay(x) <- "peaks"
    
   
        y1 = FindMarkers(x, ident.1 = "R+1", 
                         ident.2 = c("Pre-flight"))
        y8 = FindMarkers(x, ident.1 = c("R+45&R+82"), 
                         ident.2 = c("Pre-flight"))

        assign(paste0("fcs_", i, '_1'), y1)
        assign(paste0("fcs_", i, '_8'), y8)
    
        write.csv(y1, file = paste0("pval/update/all.list/update/downsample/optimize/sample.specific/DARs_",j, ".", i,"_R+1_preflight.filtered.csv"))
        write.csv(y8, file = paste0("pval/update/all.list/update/downsample/optimize/sample.specific/DARs_",j, ".", i,"_R+45&R+82_preflight.filtered.csv"))
    
    
    
    
    
    
}}








