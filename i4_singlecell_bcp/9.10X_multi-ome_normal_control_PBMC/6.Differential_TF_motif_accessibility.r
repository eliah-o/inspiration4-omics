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

library(BSgenome.Hsapiens.UCSC.hg38)


library(patchwork)


library(patchwork)
library(TFBSTools)

library(JASPAR2020)
set.seed(1234)


library(BiocParallel)
register(SerialParam())

library(ggrepel)
library(stringr)
library(ggpubr)
library(patchwork)
theme_set(theme_classic(base_size = 20))



temp.sobj <- readRDS("path_to_seurat_object")

temp.sobj$predicted.id <- ifelse(temp.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(temp.sobj$celltype.l1 == 'DC', 'DC','other')))))))



temp.sobj$Day <- ifelse(temp.sobj$orig.ident == 'Day1_Rep1', 'Day1',
                       ifelse(temp.sobj$orig.ident == 'Day1_Rep2', 'Day1',
                             ifelse(temp.sobj$orig.ident == 'Day1_Rep3', 'Day1', 'Day2')))







cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", 
             "Mono", "DC", 'other')

Idents(temp.sobj) <- 'predicted.id'



for (i in cell_list){
    x <- subset(temp.sobj, ident = i)
    assign(paste0(i,".sobj"), x)
}

all_list <- c("pbmc", cell_list)

pbmc.sobj <- temp.sobj


all_list <- c("pbmc", cell_list)

rm(list=c('temp.sobj'))

all_list <- c("pbmc", cell_list)

id_to_name <- read.csv("path_to_motif.id.to.name.csv")

head(id_to_name)
id_to_name <- id_to_name[, 2:3]
colnames(id_to_name)[1] <- 'motif'
head(id_to_name)



for (i in all_list){
    x <- get(paste0(i,".sobj"))
    DefaultAssay(x) <- "chromvar"
    
    Idents(x) <- 'Day'
    y1 <- FindMarkers(
        object = x,
        ident.1 = "Day2",
        ident.2 = c("Day1"),
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    Idents(x) <- 'orig.ident'
    y2 <- FindMarkers(
        object = x,
        ident.1 = "Day1_Rep2",
        ident.2 = "Day1_Rep1",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y3 <- FindMarkers(
        object = x,
        ident.1 = "Day1_Rep3",
        ident.2 = "Day1_Rep1",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y4 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep1",
        ident.2 = "Day1_Rep1",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y5 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep2",
        ident.2 = "Day1_Rep1",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    
    
    y6 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep3",
        ident.2 = "Day1_Rep1",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    
    y7 <- FindMarkers(
        object = x,
        ident.1 = "Day1_Rep3",
        ident.2 = "Day1_Rep2",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y8 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep1",
        ident.2 = "Day1_Rep2",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y9 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep2",
        ident.2 = "Day1_Rep2",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y10 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep3",
        ident.2 = "Day1_Rep2",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    
    y11 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep1",
        ident.2 = "Day1_Rep3",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y12 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep2",
        ident.2 = "Day1_Rep3",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y13 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep3",
        ident.2 = "Day1_Rep3",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y14 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep2",
        ident.2 = "Day2_Rep1",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y15 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep3",
        ident.2 = "Day2_Rep1",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
    y16 <- FindMarkers(
        object = x,
        ident.1 = "Day2_Rep3",
        ident.2 = "Day2_Rep2",
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    logfc.threshold = 0)
    
 
    y1$motif <- rownames(y1)
    y2$motif <- rownames(y2)
    y3$motif <- rownames(y3)
    y4$motif <- rownames(y4)
    y5$motif <- rownames(y5)
    y6$motif <- rownames(y6)
    y7$motif <- rownames(y7)
    y8$motif <- rownames(y8)
    y9$motif <- rownames(y9)
    y10$motif <- rownames(y10)
    y11$motif <- rownames(y11)
    y12$motif <- rownames(y12)
    y13$motif <- rownames(y13)
    y14$motif <- rownames(y14)
    y15$motif <- rownames(y15)
    y16$motif <- rownames(y16)
    
    
    z1 <- y1 %>% left_join(id_to_name)
    z2 <- y2 %>% left_join(id_to_name)
    z3 <- y3 %>% left_join(id_to_name)
    z4 <- y4 %>% left_join(id_to_name)
    z5 <- y5 %>% left_join(id_to_name)
    z6 <- y6 %>% left_join(id_to_name)
    z7 <- y7 %>% left_join(id_to_name)
    z8 <- y8 %>% left_join(id_to_name)
    z9 <- y9 %>% left_join(id_to_name)
    z10 <- y10 %>% left_join(id_to_name)
    z11 <- y11 %>% left_join(id_to_name)
    z12 <- y12 %>% left_join(id_to_name)
    z13 <- y13 %>% left_join(id_to_name)
    z14 <- y14 %>% left_join(id_to_name)
    z15 <- y15 %>% left_join(id_to_name)
    z16 <- y16 %>% left_join(id_to_name)
    
    
    z1$fdr <- p.adjust(z1$p_val, method = 'fdr')
    z1$DIRECTION <- ifelse(z1$avg_diff > 0.25 & z1$p_val_adj < 0.05,'enriched',
                                            ifelse(z1$avg_diff < -0.25 & z1$p_val_adj < 0.05,'depleted',NA))
    z1$LABEL <- ifelse(z1$DIRECTION %in% c('enriched', 'depleted'), z1$name, NA)
    
    z2$fdr <- p.adjust(z2$p_val, method = 'fdr')
    z2$DIRECTION <- ifelse(z2$avg_diff > 0.25 & z2$p_val_adj < 0.05,'enriched',
                                            ifelse(z2$avg_diff < -0.25 & z2$p_val_adj < 0.05,'depleted',NA))
    z2$LABEL <- ifelse(z2$DIRECTION %in% c('enriched', 'depleted'), z2$name, NA)
    
    
    z3$fdr <- p.adjust(z3$p_val, method = 'fdr')
    z3$DIRECTION <- ifelse(z3$avg_diff > 0.25 & z3$p_val_adj < 0.05,'enriched',
                                            ifelse(z3$avg_diff < -0.25 & z3$p_val_adj < 0.05,'depleted',NA))
    z3$LABEL <- ifelse(z3$DIRECTION %in% c('enriched', 'depleted'), z3$name, NA)

    
    z4$fdr <- p.adjust(z4$p_val, method = 'fdr')
    z4$DIRECTION <- ifelse(z4$avg_diff > 0.25 & z4$p_val_adj < 0.05,'enriched',
                                            ifelse(z4$avg_diff < -0.25 & z4$p_val_adj < 0.05,'depleted',NA))
    z4$LABEL <- ifelse(z4$DIRECTION %in% c('enriched', 'depleted'), z4$name, NA)

    
    z5$fdr <- p.adjust(z5$p_val, method = 'fdr')
    z5$DIRECTION <- ifelse(z5$avg_diff > 0.25 & z5$p_val_adj < 0.05,'enriched',
                                            ifelse(z5$avg_diff < -0.25 & z5$p_val_adj < 0.05,'depleted',NA))
    z5$LABEL <- ifelse(z5$DIRECTION %in% c('enriched', 'depleted'), z5$name, NA)

    
    z6$fdr <- p.adjust(z6$p_val, method = 'fdr')
    z6$DIRECTION <- ifelse(z6$avg_diff > 0.25 & z6$p_val_adj < 0.05,'enriched',
                                            ifelse(z6$avg_diff < -0.25 & z6$p_val_adj < 0.05,'depleted',NA))
    z6$LABEL <- ifelse(z6$DIRECTION %in% c('enriched', 'depleted'), z6$name, NA)

    
    z7$fdr <- p.adjust(z7$p_val, method = 'fdr')
    z7$DIRECTION <- ifelse(z7$avg_diff > 0.25 & z7$p_val_adj < 0.05,'enriched',
                                            ifelse(z7$avg_diff < -0.25 & z7$p_val_adj < 0.05,'depleted',NA))
    z7$LABEL <- ifelse(z7$DIRECTION %in% c('enriched', 'depleted'), z7$name, NA)

    
    z8$fdr <- p.adjust(z8$p_val, method = 'fdr')
    z8$DIRECTION <- ifelse(z8$avg_diff > 0.25 & z8$p_val_adj < 0.05,'enriched',
                                            ifelse(z8$avg_diff < -0.25 & z8$p_val_adj < 0.05,'depleted',NA))
    z8$LABEL <- ifelse(z8$DIRECTION %in% c('enriched', 'depleted'), z8$name, NA)

    z9$fdr <- p.adjust(z9$p_val, method = 'fdr')
    z9$DIRECTION <- ifelse(z9$avg_diff > 0.25 & z9$p_val_adj < 0.05,'enriched',
                                            ifelse(z9$avg_diff < -0.25 & z9$p_val_adj < 0.05,'depleted',NA))
    z9$LABEL <- ifelse(z9$DIRECTION %in% c('enriched', 'depleted'), z9$name, NA)

    z10$fdr <- p.adjust(z10$p_val, method = 'fdr')
    z10$DIRECTION <- ifelse(z10$avg_diff > 0.25 & z10$p_val_adj < 0.05,'enriched',
                                            ifelse(z10$avg_diff < -0.25 & z10$p_val_adj < 0.05,'depleted',NA))
    z10$LABEL <- ifelse(z10$DIRECTION %in% c('enriched', 'depleted'), z10$name, NA)

    z11$fdr <- p.adjust(z11$p_val, method = 'fdr')
    z11$DIRECTION <- ifelse(z11$avg_diff > 0.25 & z11$p_val_adj < 0.05,'enriched',
                                            ifelse(z11$avg_diff < -0.25 & z11$p_val_adj < 0.05,'depleted',NA))
    z11$LABEL <- ifelse(z11$DIRECTION %in% c('enriched', 'depleted'), z11$name, NA)

    z12$fdr <- p.adjust(z12$p_val, method = 'fdr')
    z12$DIRECTION <- ifelse(z12$avg_diff > 0.25 & z12$p_val_adj < 0.05,'enriched',
                                            ifelse(z12$avg_diff < -0.25 & z12$p_val_adj < 0.05,'depleted',NA))
    z12$LABEL <- ifelse(z12$DIRECTION %in% c('enriched', 'depleted'), z12$name, NA)

    z13$fdr <- p.adjust(z13$p_val, method = 'fdr')
    z13$DIRECTION <- ifelse(z13$avg_diff > 0.25 & z13$p_val_adj < 0.05,'enriched',
                                            ifelse(z13$avg_diff < -0.25 & z13$p_val_adj < 0.05,'depleted',NA))
    z13$LABEL <- ifelse(z13$DIRECTION %in% c('enriched', 'depleted'), z13$name, NA)

    z14$fdr <- p.adjust(z14$p_val, method = 'fdr')
    z14$DIRECTION <- ifelse(z14$avg_diff > 0.25 & z14$p_val_adj < 0.05,'enriched',
                                            ifelse(z14$avg_diff < -0.25 & z14$p_val_adj < 0.05,'depleted',NA))
    z14$LABEL <- ifelse(z14$DIRECTION %in% c('enriched', 'depleted'), z14$name, NA)

    z15$fdr <- p.adjust(z15$p_val, method = 'fdr')
    z15$DIRECTION <- ifelse(z15$avg_diff > 0.25 & z15$p_val_adj < 0.05,'enriched',
                                            ifelse(z15$avg_diff < -0.25 & z15$p_val_adj < 0.05,'depleted',NA))
    z15$LABEL <- ifelse(z15$DIRECTION %in% c('enriched', 'depleted'), z15$name, NA)

    z16$fdr <- p.adjust(z16$p_val, method = 'fdr')
    z16$DIRECTION <- ifelse(z16$avg_diff > 0.25 & z16$p_val_adj < 0.05,'enriched',
                                            ifelse(z16$avg_diff < -0.25 & z16$p_val_adj < 0.05,'depleted',NA))
    z16$LABEL <- ifelse(z16$DIRECTION %in% c('enriched', 'depleted'), z16$name, NA)

    
    
    
    assign(paste0(i,".Day2.Day1.label"), z1)
    assign(paste0(i,".Day1.Rep2.Day1.Rep1.label"), z2)
    assign(paste0(i,".Day1.Rep3.Day1.Rep1.label"), z3)
    assign(paste0(i,".Day2.Rep1.Day1.Rep1.label"), z4)
    assign(paste0(i,".Day2.Rep2.Day1.Rep1.label"), z5)
    assign(paste0(i,".Day2.Rep3.Day1.Rep1.label"), z6)
    assign(paste0(i,".Day1.Rep3.Day1.Rep2.label"), z7)
    assign(paste0(i,".Day2.Rep1.Day1.Rep2.label"), z8)
    assign(paste0(i,".Day2.Rep2.Day1.Rep2.label"), z9)
    assign(paste0(i,".Day2.Rep3.Day1.Rep2.label"), z10)
    assign(paste0(i,".Day2.Rep1.Day1.Rep3.label"), z11)
    assign(paste0(i,".Day2.Rep2.Day1.Rep3.label"), z12)
    assign(paste0(i,".Day2.Rep3.Day1.Rep3.label"), z13)
    assign(paste0(i,".Day2.Rep2.Day2.Rep1.label"), z14)
    assign(paste0(i,".Day2.Rep3.Day2.Rep1.label"), z15)
    assign(paste0(i,".Day2.Rep3.Day2.Rep2.label"), z16)
    
    
    
   
    
    
}






