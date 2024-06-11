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









all_list <- c('pbmc', "CD4_T", "CD8_T", "other_T", "B", "NK", 
              "CD14_Mono", "CD16_Mono", "DC", 'other')

pbmc.sobj <- temp.sobj



for (i in all_list){
    
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "peaks"
    
    Idents(x) <- 'timepoint'
    
    y1 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/differential/DARs_", i, "_November_JAS.filter.csv"))
    y3 <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/differential/DARs_", i, "_December_JAS.filter.csv"))
    
    z1 <- y1[y1$p_val < 0.005, ]$X
    z2 <- y2[y2$p_val < 0.005, ]$X
    z3 <- y3[y3$p_val < 0.005, ]$X
    
    assign(paste0(i,".top.da.peaks.Sept_post"), z1)
    assign(paste0(i,".top.da.peaks.Nov"), z2)
    assign(paste0(i,".top.da.peaks.Dec"), z3)
    
    a1 <- FindMotifs(
        object = x,
        features = z1)
    
    a2 <- FindMotifs(
        object = x,
        features = z2)
    
    a3 <- FindMotifs(
        object = x,
        features = z3)
    
    write.csv(a1, paste0("2023_08_08_update/motif.csv/", i, ".enriched.motifs.Sept_post.csv"))
    write.csv(a2, paste0("2023_08_08_update/motif.csv/", i, ".enriched.motifs.Nov.csv"))
    write.csv(a3, paste0("2023_08_08_update/motif.csv/", i, ".enriched.motifs.Dec.csv"))
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), a1)
    assign(paste0(i, ".enriched.motifs.Nov"), a2)
    assign(paste0(i, ".enriched.motifs.Dec"), a3)
    
    
    
}









for (i in all_list){
    x <- get(paste0(i, ".sobj"))
    y <- get(paste0(i, ".enriched.motifs.Sept_post"))
    
    options(repr.plot.width=16, repr.plot.height=8)
    
    a <- MotifPlot(
      object = x,
      motifs = y$motif[1:12]
)+ 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size=15),
       strip.text = element_text(size =30))
    
    
    ggsave(paste0("2023_08_08_update/motif.image/", i, ".tiff"), width=16, height=8, dpi=300, compression='lzw')
    ggsave(paste0("2023_08_08_update/motif.image/", i, ".png"), width=16, height=8, dpi=300)
    ggsave(paste0("2023_08_08_update/motif.image/", i, ".svg"), width=16, height=8, dpi=300)
    ggsave(paste0("2023_08_08_update/motif.image/", i, ".pdf"), width=16, height=8, dpi=300)
    
    print(a)
}



B_cells <- c("B_Naive", "B_intermediate", "B_memory")

NK_cells <- c("NK_Proliferating", "NK_CD56bright")

temp.sobj$celltype.l3 <- ifelse(temp.sobj$celltype.l2 == "B naive", "B_Naive",
                               ifelse(temp.sobj$celltype.l2 == "B intermediate", "B_intermediate",
                                     ifelse(temp.sobj$celltype.l2 == "B memory", "B_memory",
                                           ifelse(temp.sobj$celltype.l2 == "NK Proliferating", "NK_Proliferating",
                                                 ifelse(temp.sobj$celltype.l2 == "NK_CD56bright", "NK_CD56bright",
                                                       ifelse(temp.sobj$celltype.l2 == "ILC", "ILC", "other"))))))



Idents(temp.sobj) <- 'celltype.l3'

for (i in B_cells){
    x <- subset(temp.sobj, ident = i)
    assign(paste0(i,".sobj"), x)
}





for (i in B_cells){
    a <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    
    x <- get(paste0(i, ".sobj"))
    
    a1 <- (a[a$p_val < 0.005,])$X

    DefaultAssay(x) <- 'peaks'

    z <- FindMotifs(object = x,
                                        features = a1)
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), z)



    
}

Idents(temp.sobj) <- 'celltype.l3'

for (i in NK_cells){
    x <- subset(temp.sobj, ident = i)
    assign(paste0(i,".sobj"), x)
}

for (i in "NK_CD56bright"){
    a <- read.csv(paste0("../../../Differential_expression/with_filter/celltype/pval/update/celltype.l2/differential/DARs_", i, "_September_Post_JAS.filter.csv"))
    
    x <- get(paste0(i, ".sobj"))
    
    a1 <- (a[a$p_val < 0.005,])$X

    DefaultAssay(x) <- 'peaks'

    z <- FindMotifs(object = x,
                                        features = a1)
    
    assign(paste0(i, ".enriched.motifs.Sept_post"), z)



    
}



for (i in c(B_cells, "NK_CD56bright")){
    x <- get(paste0(i, ".sobj"))
    y <- get(paste0(i, ".enriched.motifs.Sept_post"))
    
    DefaultAssay(x) <- "peaks"
    
    
    options(repr.plot.width=16, repr.plot.height=8)
    
    a <- MotifPlot(
      object = x,
      motifs = y$motif[1:12])+ 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size=15),
       strip.text = element_text(size =30))
    
    
    ggsave(paste0("2023_08_08_update/motif.image/B_NK/", i, ".tiff"), width=16, height=8, dpi=300, compression='lzw')
    ggsave(paste0("2023_08_08_update/motif.image/B_NK/", i, ".png"), width=16, height=8, dpi=300)
    ggsave(paste0("2023_08_08_update/motif.image/B_NK/", i, ".svg"), width=16, height=8, dpi=300)
    ggsave(paste0("2023_08_08_update/motif.image/B_NK/", i, ".pdf"), width=16, height=8, dpi=300)
    
    print(a)
}





getwd()





sessionInfo()

rm(list=ls())










