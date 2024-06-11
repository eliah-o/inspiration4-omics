library(ggplot2)
library(dplyr)
library(tidyr)


library(patchwork)


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
    x1 <- read.csv(paste0("path_to_motifs.Sept_post"))
    x2 <- read.csv(paste0("path_to_motifs.longterm"))
   
    
    
    assign(paste0(i,".motifs_FP1.label"), x1)
    assign(paste0(i,".motifs_LP3.label"), x2)
    
    
}

for (i in all_list){
    x1 <- read.csv(paste0("path_to_motifs.R+45.vs.R+1"))
    x2 <- read.csv(paste0("path_to_motifs.R+82.vs.R+1"))
   
    
    
    assign(paste0(i,".motifs_R.45.vs.R.1.label"), x1)
    assign(paste0(i,".motifs_R.82.vs.R.1.label"), x2)
    
    
}



motifs_FP1.shared.motifs <- Reduce(intersect, list(pbmc.motifs_FP1.label$name, CD4_T.motifs_FP1.label$name,
                        CD8_T.motifs_FP1.label$name, other_T.motifs_FP1.label$name,
                        B.motifs_FP1.label$name, NK.motifs_FP1.label$name,
                        CD14_Mono.motifs_FP1.label$name, CD16_Mono.motifs_FP1.label$name,
                        DC.motifs_FP1.label$name, other.motifs_FP1.label$name 
))


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
    x1 <- get(paste0(i,".motifs_FP1.label"))
    
    
    
    z1 <- x1 %>% dplyr::filter(p_val_adj < 0.05)

    assign(paste0(i,".FP1.filter"), z1)
    
    
    
    
}

motifs_FP1.filtered <- Reduce(intersect, list(pbmc.FP1.filter$name, CD4_T.FP1.filter$name,
                        CD8_T.FP1.filter$name, other_T.FP1.filter$name,
                        B.FP1.filter$name, NK.FP1.filter$name,
                        CD14_Mono.FP1.filter$name, CD16_Mono.FP1.filter$name,
                        DC.FP1.filter$name, other.FP1.filter$name
 
))




for (i in all_list){
    x1 <- get(paste0(i,".motifs_R.82.vs.R.1.label"))
    x2 <- get(paste0(i,".motifs_R.45.vs.R.1.label"))
    x3 <- get(paste0(i,".motifs_FP1.label"))
    
    y <- x3 %>% filter(p_val_adj < 0.05)
    
    
    
    z1 <- x1 %>% dplyr::filter(name %in% y$name)
    z2 <- x2 %>% dplyr::filter(name %in% y$name)
    z3 <- x3 %>% dplyr::filter(name %in% y$name)
    
    colnames(z1)[3] <-'avg_diff.RP2'
    colnames(z2)[3] <-'avg_diff.RP1'
    colnames(z3)[3] <-'avg_diff.FP1'
    
    z1 <- z1 %>% select(LABEL, avg_diff.RP2)
    z2 <- z2 %>% select(LABEL, avg_diff.RP1)
    z3 <- z3 %>% select(LABEL, avg_diff.FP1)

    assign(paste0(i,".RP2.filter"), z1)
    assign(paste0(i,".RP1.filter"), z2)
    assign(paste0(i,".FP1.filter"), z3)
    
    
    
    
}



for (i in all_list){
    x1 <- get(paste0(i,".RP2.filter"))
    x2 <- get(paste0(i,".RP1.filter"))
    x3 <- get(paste0(i,".FP1.filter"))
    
    z1 <- x3 %>% left_join(x1) 
    z2 <- x3 %>% left_join(x2)
    
    z1$label <- 'RP2'
    z2$label <- 'RP1'
    
    colnames(z1)[3] <- 'avg_diff'
    colnames(z2)[3] <- 'avg_diff'
    
    
    
    
    
    assign(paste0(i, ".FP1.RP2"), z1)
    assign(paste0(i, ".FP1.RP1"), z2)
    
    y <- rbind(z1, z2)
    
    assign(paste0(i, '.merge'), y)
    

    
    
    
}











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


options(repr.plot.width=7, repr.plot.height=6)

ggplot(pbmc.merge, aes(x = avg_diff.FP1,y = avg_diff))+
    geom_point(aes(color = label))+
    labs(color="Stage", x = "Delta z-score for motif accessibiility\nRecovery Profile (RP1 or RP2)",
        y = "Delta z-score for motif accessibiility\nFlight Profile (FP1)")+
    scale_color_manual(values = c("#74c476","#dd1c77"))+
    xlim(-3, 3)+
    ylim(-3, 3)+ 
    geom_hline(yintercept=0, size=1)+
    geom_vline(xintercept = 0, size=1)











for (i in all_list){
    
    x1 <- get(paste0(i,".merge"))
    
    options(repr.plot.width=7, repr.plot.height=6)

    p1 <- ggplot(x1, aes(x = avg_diff.FP1,y = avg_diff))+
    geom_point(aes(color = label))+
    labs(color="Stage", x = "Delta z-score for motif accessibiility\nRecovery Profile (RP1 or RP2)",
        y = "Delta z-score for motif accessibiility\nFlight Profile (FP1)")+
    scale_color_manual(values = c("#74c476","#dd1c77"))+
    xlim(-3, 3)+
    ylim(-3, 3)+ 
    geom_hline(yintercept=0, size=1)+
    geom_vline(xintercept = 0, size=1)
    
    
   
    print(p1)
    
    ggsave(paste0("output/2024_03_07/", i, ".motif.pdf"),width=7,height=7,limitsize = FALSE)
    ggsave(paste0("output/2024_03_07/", i, ".motif.png"),width=7,height=7,limitsize = FALSE)
    ggsave(paste0("output/2024_03_07/", i, ".motif.tiff"),width=7,height=7,
           dpi=300, compression='lzw', limitsize = FALSE)
    
}







getwd()

sessionInfo()

rm(list=ls())




