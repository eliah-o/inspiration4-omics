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

for (i in c('B', 'CD8_T', 'DC', 'Mono', 'NK', 'other', 'other_T')){
    
    x <- readRDS(paste0("path_to_Seurat_object"))
    
    
    assign(paste0(i, '.down.sobj'), x)

}





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





B.sobj <- B.down.sobj
CD8_T.sobj <- CD8_T.down.sobj
DC.sobj <- DC.down.sobj
Mono.sobj <- Mono.down.sobj
NK.sobj <- NK.down.sobj
other.sobj <- other.down.sobj
other_T.sobj <- other_T.down.sobj




for (i in all_list){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "orig.ident"
    DefaultAssay(x) <- "SCT"
    
    
    x <- PrepSCTFindMarkers(object = x)
   
    y1 = FindMarkers(x, ident.1 = "Day1_Rep3", 
                     ident.2 = "Day1_Rep2")
    y2 = FindMarkers(x, ident.1 = "Day2_Rep1", 
                     ident.2 = "Day1_Rep2")
    y3 = FindMarkers(x, ident.1 = "Day2_Rep2", 
                     ident.2 = "Day1_Rep2")
    y4 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day1_Rep2")
    y5 = FindMarkers(x, ident.1 = "Day2_Rep1", 
                     ident.2 = "Day1_Rep3")
    y6 = FindMarkers(x, ident.1 = "Day2_Rep2", 
                     ident.2 = "Day1_Rep3")
    y7 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day1_Rep3")
    y8 = FindMarkers(x, ident.1 = "Day2_Rep2", 
                     ident.2 = "Day2_Rep1")
    y9 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day2_Rep1")
    y10 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day2_Rep2")


    
     
}





for (i in all_list){
    x <- get(paste0(i, '.sobj'))
    Idents(x) <- "orig.ident"
    DefaultAssay(x) <- "ATAC"
    
    
   
    y1 = FindMarkers(x, ident.1 = "Day1_Rep3", 
                     ident.2 = "Day1_Rep2")
    y2 = FindMarkers(x, ident.1 = "Day2_Rep1", 
                     ident.2 = "Day1_Rep2")
    y3 = FindMarkers(x, ident.1 = "Day2_Rep2", 
                     ident.2 = "Day1_Rep2")
    y4 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day1_Rep2")
    y5 = FindMarkers(x, ident.1 = "Day2_Rep1", 
                     ident.2 = "Day1_Rep3")
    y6 = FindMarkers(x, ident.1 = "Day2_Rep2", 
                     ident.2 = "Day1_Rep3")
    y7 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day1_Rep3")
    y8 = FindMarkers(x, ident.1 = "Day2_Rep2", 
                     ident.2 = "Day2_Rep1")
    y9 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day2_Rep1")
    y10 = FindMarkers(x, ident.1 = "Day2_Rep3", 
                     ident.2 = "Day2_Rep2")


    
    
    
    
    
    
    
}








