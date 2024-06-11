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


library(BiocParallel)
register(SerialParam())

library(ggrepel)
library(stringr)
library(ggpubr)
library(patchwork)
theme_set(theme_classic(base_size = 20))





cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", 
              "CD14_Mono", "CD16_Mono", "DC", 'other')

all_list <- c("pbmc", cell_list)

id_to_name <- read.csv("path_to_motif.id.to.name")

head(id_to_name)
id_to_name <- id_to_name[, 2:3]
colnames(id_to_name)[1] <- 'motif'
head(id_to_name)



for (i in all_list){
    x1 <- read.csv(paste0("path_to_TF.motifs.from.chromvar.R+45.vs.R+1"))
    x2 <- read.csv(paste0("path_to_TF.motifs.from.chromvar.R+82.vs.R+1"))
   
    
    
    assign(paste0(i,".motifs_R.45.vs.R.1.label"), x1)
    assign(paste0(i,".motifs_R.82.vs.R.1.label"), x2)
    
    
}



motifs_R.45.vs.R.1.shared.motifs <- Reduce(intersect, list(pbmc.motifs_R.45.vs.R.1.label$name, CD4_T.motifs_R.45.vs.R.1.label$name,
                        CD8_T.motifs_R.45.vs.R.1.label$name, other_T.motifs_R.45.vs.R.1.label$name,
                        B.motifs_R.45.vs.R.1.label$name, NK.motifs_R.45.vs.R.1.label$name,
                        CD14_Mono.motifs_R.45.vs.R.1.label$name, CD16_Mono.motifs_R.45.vs.R.1.label$name,
                        DC.motifs_R.45.vs.R.1.label$name, other.motifs_R.45.vs.R.1.label$name))


motifs_R.82.vs.R.1.shared.motifs <- Reduce(intersect, list(pbmc.motifs_R.82.vs.R.1.label$name, CD4_T.motifs_R.82.vs.R.1.label$name,
                        CD8_T.motifs_R.82.vs.R.1.label$name, other_T.motifs_R.82.vs.R.1.label$name,
                        B.motifs_R.82.vs.R.1.label$name, NK.motifs_R.82.vs.R.1.label$name,
                        CD14_Mono.motifs_R.82.vs.R.1.label$name, CD16_Mono.motifs_R.82.vs.R.1.label$name,
                        DC.motifs_R.82.vs.R.1.label$name, other.motifs_R.82.vs.R.1.label$name))


for (i in all_list){
    x1 <- get(paste0(i,".motifs_R.45.vs.R.1.label"))
    x2 <- get(paste0(i,".motifs_R.82.vs.R.1.label"))
    
    x1 <- x1 %>% dplyr::filter(name %in% motifs_R.45.vs.R.1.shared.motifs)
    x2<- x2 %>% dplyr::filter(name %in% motifs_R.82.vs.R.1.shared.motifs)
    
    z1 <- x1 %>% dplyr::filter(avg_diff > 0)
    z2 <- x1 %>% dplyr::filter(avg_diff < 0)
    
    z3 <- x2 %>% dplyr::filter(avg_diff > 0)
    z4 <- x2 %>% dplyr::filter(avg_diff < 0)
    
    y1 <- x1$name[1:5]
    y2 <- x2$name[1:5]
    
    b1 <- z1$name[1:5]
    b2 <- z2$name[1:5]
    
    b3 <- z3$name[1:5]
    b4 <- z4$name[1:5]
    
    a1 <- x1[order(x1$name),]
    a2 <- x2[order(x2$name),]
    
    
   
    assign(paste0(i,".motifs_R.45.vs.R.1.sort"), a1)
    assign(paste0(i,".motifs_R.82.vs.R.1.sort"), a2)
    
    assign(paste0(i,".R.45.vs.R.1.top5"), y1)
    assign(paste0(i,".R.82.vs.R.1.top5"), y2)
    
    
    assign(paste0(i,".R.45.vs.R.1.up.top5"), b1)
    assign(paste0(i,".R.45.vs.R.1.down.top5"), b2)

    assign(paste0(i,".R.82.vs.R.1.up.top5"), b3)
    assign(paste0(i,".R.82.vs.R.1.down.top5"), b4)
    
    print(paste0(i, ": R.45.vs.R.1 up ", b1))
    print(paste0(i, ": R.45.vs.R.1 down ", b2))
    print(paste0(i, ": R.82.vs.R.1 up ", b3))
    print(paste0(i, ": R.82.vs.R.1 down ", b4))
    
    
    
    
}



up.R.45.vs.R.1.motifs <- c(pbmc.R.45.vs.R.1.up.top5, CD4_T.R.45.vs.R.1.up.top5,
                        CD8_T.R.45.vs.R.1.up.top5, other_T.R.45.vs.R.1.up.top5,
                        B.R.45.vs.R.1.up.top5, NK.R.45.vs.R.1.up.top5,
                        CD14_Mono.R.45.vs.R.1.up.top5, CD16_Mono.R.45.vs.R.1.up.top5,
                        DC.R.45.vs.R.1.up.top5, other.R.45.vs.R.1.up.top5)


down.R.45.vs.R.1.motifs <- c(pbmc.R.45.vs.R.1.down.top5, CD4_T.R.45.vs.R.1.down.top5,
                        CD8_T.R.45.vs.R.1.down.top5, other_T.R.45.vs.R.1.down.top5,
                        B.R.45.vs.R.1.down.top5, NK.R.45.vs.R.1.down.top5,
                        CD14_Mono.R.45.vs.R.1.down.top5, CD16_Mono.R.45.vs.R.1.down.top5,
                        DC.R.45.vs.R.1.down.top5, other.R.45.vs.R.1.down.top5)





up.R.82.vs.R.1.motifs <- c(pbmc.R.82.vs.R.1.up.top5, CD4_T.R.82.vs.R.1.up.top5,
                        CD8_T.R.82.vs.R.1.up.top5, other_T.R.82.vs.R.1.up.top5,
                        B.R.82.vs.R.1.up.top5, NK.R.82.vs.R.1.up.top5,
                        CD14_Mono.R.82.vs.R.1.up.top5, CD16_Mono.R.82.vs.R.1.up.top5,
                        DC.R.82.vs.R.1.up.top5, other.R.82.vs.R.1.up.top5)


down.R.82.vs.R.1.motifs <- c(pbmc.R.82.vs.R.1.down.top5, CD4_T.R.82.vs.R.1.down.top5,
                        CD8_T.R.82.vs.R.1.down.top5, other_T.R.82.vs.R.1.down.top5,
                        B.R.82.vs.R.1.down.top5, NK.R.82.vs.R.1.down.top5,
                        CD14_Mono.R.82.vs.R.1.down.top5, CD16_Mono.R.82.vs.R.1.down.top5,
                        DC.R.82.vs.R.1.down.top5, other.R.82.vs.R.1.down.top5)








for (i in all_list){
    x1 <- get(paste0(i,".motifs_R.45.vs.R.1.sort"))
    x2 <- get(paste0(i,".motifs_R.82.vs.R.1.sort"))
    
    z1 <- x1 %>% dplyr::filter(name %in% unique(up.R.45.vs.R.1.motifs))
    z2 <- x1 %>% dplyr::filter(name %in% unique(down.R.45.vs.R.1.motifs))
    
    z3 <- x2 %>% dplyr::filter(name %in% unique(up.R.82.vs.R.1.motifs))
    z4 <- x2 %>% dplyr::filter(name %in% unique(down.R.82.vs.R.1.motifs))
    
    z1 <- z1[order(z1$name),]
    z2 <- z2[order(z2$name),]
    
    z3 <- z3[order(z3$name),]
    z4 <- z4[order(z4$name),]
    
    assign(paste0(i,".R.45.vs.R.1.up"), z1)
    assign(paste0(i,".R.45.vs.R.1.down"), z2)
    
    assign(paste0(i,".R.82.vs.R.1.up"), z3)
    assign(paste0(i,".R.82.vs.R.1.down"), z4)

    
    
    
    
    
}





up.R.45.vs.R.1 <- data.frame(PBMC = pbmc.R.45.vs.R.1.up$avg_diff, CD4_T = CD4_T.R.45.vs.R.1.up$avg_diff,
                        CD8_T = CD8_T.R.45.vs.R.1.up$avg_diff, other_T = other_T.R.45.vs.R.1.up$avg_diff,
                        B = B.R.45.vs.R.1.up$avg_diff, NK = NK.R.45.vs.R.1.up$avg_diff,
                        CD14_Mono = CD14_Mono.R.45.vs.R.1.up$avg_diff, CD16_Mono = CD16_Mono.R.45.vs.R.1.up$avg_diff,
                        DC = DC.R.45.vs.R.1.up$avg_diff, other = other.R.45.vs.R.1.up$avg_diff)


down.R.45.vs.R.1 <- data.frame(PBMC = pbmc.R.45.vs.R.1.down$avg_diff, CD4_T = CD4_T.R.45.vs.R.1.down$avg_diff,
                        CD8_T = CD8_T.R.45.vs.R.1.down$avg_diff, other_T = other_T.R.45.vs.R.1.down$avg_diff,
                        B = B.R.45.vs.R.1.down$avg_diff, NK = NK.R.45.vs.R.1.down$avg_diff,
                        CD14_Mono = CD14_Mono.R.45.vs.R.1.down$avg_diff, CD16_Mono = CD16_Mono.R.45.vs.R.1.down$avg_diff,
                        DC = DC.R.45.vs.R.1.down$avg_diff, other = other.R.45.vs.R.1.down$avg_diff)

up.R.82.vs.R.1 <- data.frame(PBMC = pbmc.R.82.vs.R.1.up$avg_diff, CD4_T = CD4_T.R.82.vs.R.1.up$avg_diff,
                        CD8_T = CD8_T.R.82.vs.R.1.up$avg_diff, other_T = other_T.R.82.vs.R.1.up$avg_diff,
                        B = B.R.82.vs.R.1.up$avg_diff, NK = NK.R.82.vs.R.1.up$avg_diff,
                        CD14_Mono = CD14_Mono.R.82.vs.R.1.up$avg_diff, CD16_Mono = CD16_Mono.R.82.vs.R.1.up$avg_diff,
                        DC = DC.R.82.vs.R.1.up$avg_diff, other = other.R.82.vs.R.1.up$avg_diff)


down.R.82.vs.R.1 <- data.frame(PBMC = pbmc.R.82.vs.R.1.down$avg_diff, CD4_T = CD4_T.R.82.vs.R.1.down$avg_diff,
                        CD8_T = CD8_T.R.82.vs.R.1.down$avg_diff, other_T = other_T.R.82.vs.R.1.down$avg_diff,
                        B = B.R.82.vs.R.1.down$avg_diff, NK = NK.R.82.vs.R.1.down$avg_diff,
                        CD14_Mono = CD14_Mono.R.82.vs.R.1.down$avg_diff, CD16_Mono = CD16_Mono.R.82.vs.R.1.down$avg_diff,
                        DC = DC.R.82.vs.R.1.down$avg_diff, other = other.R.82.vs.R.1.down$avg_diff)




row.names(up.R.45.vs.R.1) <- pbmc.R.45.vs.R.1.up$name
row.names(down.R.45.vs.R.1) <- pbmc.R.45.vs.R.1.down$name

row.names(up.R.82.vs.R.1) <- pbmc.R.82.vs.R.1.up$name
row.names(down.R.82.vs.R.1) <- pbmc.R.82.vs.R.1.down$name



up.down.R45.vs.R.1 <- rbind(up.R.45.vs.R.1, down.R.45.vs.R.1)
up.down.R82.vs.R.1 <- rbind(up.R.82.vs.R.1, down.R.82.vs.R.1)





library(ggplot2)

library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(pheatmap)

library(ggrepel)

library(tidyr)
library(ggsignif)
library(ggrepel)
library(tidyverse)

library(circlize)
library(ComplexHeatmap)




paletteLength <- 200
low = colorRampPalette(c("black","navy","white"))(100)
high = colorRampPalette(c("white","firebrick1","firebrick2","firebrick3","firebrick4"))(100)
combined = c(low,high)
myColor2 <- combined

myBreaks.R45.vs.R.1 <- c(seq(min(up.down.R45.vs.R.1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(up.down.R45.vs.R.1)/paletteLength, max(up.down.R45.vs.R.1), length.out=floor(paletteLength/2)))


myBreaks.R82.vs.R.1 <- c(seq(min(up.down.R82.vs.R.1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(up.down.R82.vs.R.1)/paletteLength, max(up.down.R82.vs.R.1), length.out=floor(paletteLength/2)))




col_fun <- colorRamp2(seq(-2,2,0.4),rev(brewer.pal(11, "RdBu")))



options(repr.plot.width=10, repr.plot.height=10)

Heatmap(up.down.R45.vs.R.1,
        column_labels = colnames(up.down.R45.vs.R.1),
        border = TRUE,
        use_raster = FALSE,
        row_split = ifelse(rownames(up.down.R45.vs.R.1)%in%up.R.45.vs.R.1.motifs,"Up-regulated","Down-regulated"),
        cluster_row_slices = FALSE,
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        heatmap_legend_param = list(direction = "vertical", title = "z-score", title_position = "topcenter"),
        ,
                                   col = col_fun
       )

options(repr.plot.width=10, repr.plot.height=10)

Heatmap(up.down.R82.vs.R.1,
        column_labels = colnames(up.down.R82.vs.R.1),
        border = TRUE,
        use_raster = FALSE,
        row_split = ifelse(rownames(up.down.R82.vs.R.1)%in%up.R.82.vs.R.1.motifs,"Up-regulated","Down-regulated"),
        cluster_row_slices = FALSE,
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        heatmap_legend_param = list(direction = "vertical", title = "z-score", title_position = "topcenter"),
                                   col = col_fun
       )







sessionInfo()

rm(list=ls())




