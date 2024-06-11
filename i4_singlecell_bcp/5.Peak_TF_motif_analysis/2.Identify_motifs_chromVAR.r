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



temp.sobj$predicted.id <- ifelse(temp.sobj$celltype == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype == 'CD14 Mono', 'CD14_Mono',
                                                                 ifelse(temp.sobj$celltype == 'CD16 Mono', 'CD16_Mono',
                                                                         ifelse(temp.sobj$celltype == 'DC', 'DC','other'))))))))



cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", 
              "CD14_Mono", "CD16_Mono", "DC", 'other')

Idents(temp.sobj) <- 'predicted.id'

for (i in cell_list){
    x <- subset(temp.sobj, ident = i)
    assign(paste0(i,".sobj"), x)
}

pbmc.sobj <- temp.sobj

all_list <- c("pbmc", cell_list)









for (i in all_list){
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "chromvar"
    
    Idents(x) <- 'timepoint'
    y1 <- FindMarkers(
        object = x,
        ident.1 = c("June 2021", "August 2021", "September Pre-launch"),
        ident.2 = "September Post-launch")
    
    y2 <- FindMarkers(
        object = x,
        ident.1 = c("June 2021", "August 2021", "September Pre-launch"),
        ident.2 = "November 2021")
    
    y3 <- FindMarkers(
        object = x,
        ident.1 = c("June 2021", "August 2021", "September Pre-launch"),
        ident.2 = "December 2021")
    
    assign(paste0(i,".motifs_Sept_post"), y1)
    assign(paste0(i,".motifs_Nov"), y2)
    assign(paste0(i,".motifss_Dec"), y3)
    
    
    write.csv(y1, paste0("output/2023_08_08/", i, ".motifs.Sept_post.csv"))
    write.csv(y2, paste0("output/2023_08_08/", i, ".motifs.Nov.csv"))
    write.csv(y3, paste0("output/2023_08_08/", i, ".motifs.Dec.csv"))
    
    
    
}







sessionInfo()

rm(list=ls())


