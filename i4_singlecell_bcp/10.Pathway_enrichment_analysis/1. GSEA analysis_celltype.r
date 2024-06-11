library(ggplot2)
library(dplyr)
library(tidyr)

library(fgsea)
library(ggplot2)

## Setup
### Bioconductor and CRAN libraries used
library(DESeq2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(pheatmap)


library(ggthemes)
library(ggtext)
library(stringr)
library(forcats)
library("RColorBrewer")
library(ggpubr)

library(msigdbr)

m_df.hu = msigdbr(species = "Homo sapiens")

m_df.hu.hallmark = m_df.hu %>% filter(gs_cat %in% c("H"))
m_df.hu.C2 = m_df.hu %>% filter(gs_cat %in% c("C2"))
m_df.hu.C5 = m_df.hu %>% filter(gs_cat %in% c("C5"))
m_df.hu.C7 = m_df.hu %>% filter(gs_cat %in% c("C7"))
m_df.hu.C8 = m_df.hu %>% filter(gs_cat %in% c("C8"))

m_list.hu = m_df.hu %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.hallmark = m_df.hu.hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C2 = m_df.hu.C2 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C5 = m_df.hu.C5 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C7 = m_df.hu.C7 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C8 = m_df.hu.C8 %>% split(x = .$gene_symbol, f = .$gs_name)

hallmark <- names(m_list.hu.hallmark)
C2 <- names(m_list.hu.C2)
C5 <- names(m_list.hu.C5)



cell_list <- c('pbmc', "CD4_T", "CD8_T", "other_T", "B", "NK", 
              "CD14_Mono", "CD16_Mono", "DC", 'other')



for (i in cell_list){
    
    y4 <- read.csv(paste0("path_to_DEG_R+1_preflight.csv"))
    y5 <- read.csv(paste0("path_to_DEG_R+45&R+82_preflight.csv"))
    
  
    colnames(y4)[1] <- "Gene"
    colnames(y5)[1] <- "Gene"
    
  
    y4$timepoint <- "Immediately Post-flight"
    y5$timepoint <- "Long-term Post-flight"
 
    
    
    assign(paste0(i,".immediately.FC"), y4)
    assign(paste0(i,".longterm.FC"), y5)
    
    
}



for (i in cell_list){
    
   
    
    x4 <- get(paste0(i,".immediately.FC"))
    x5 <- get(paste0(i,".longterm.FC"))
    
    
    x4 <- x4 %>% filter(!is.na(avg_log2FC))
    x4 <- x4 %>% filter(!is.infinite(avg_log2FC))
    
    x5 <- x5 %>% filter(!is.na(avg_log2FC))
    x5 <- x5 %>% filter(!is.infinite(avg_log2FC))
    
    
    y4 <- x4$avg_log2FC
    y5 <- x5$avg_log2FC
    
    
    names(y4) <- x4$Gene
    names(y5) <- x5$Gene
    
    
    assign(paste0(i,".gse.immediately.FC"), y4)
    assign(paste0(i,".gse.longterm.FC"), y5)
    
}

for (i in cell_list){
    
   
    
    x4 <- get(paste0(i,".immediately.FC"))
    x5 <- get(paste0(i,".longterm.FC"))
    
    
    x4 <- x4 %>% filter(!is.na(avg_log2FC))
    x4 <- x4 %>% filter(!is.infinite(avg_log2FC))
    
    x5 <- x5 %>% filter(!is.na(avg_log2FC))
    x5 <- x5 %>% filter(!is.infinite(avg_log2FC))
    
    
    y4 <- x4$p_val_adj
    y5 <- x5$p_val_adj
    
    
    names(y4) <- x4$Gene
    names(y5) <- x5$Gene
    
    
    assign(paste0(i,".gse.immediately.padj"), y4)
    assign(paste0(i,".gse.longterm.pajd"), y5)
    
}







for (i in cell_list){
    
   
    
    x4 <- get(paste0(i,".immediately.FC"))
    x5 <- get(paste0(i,".longterm.FC"))
    
    
    x4 <- x4 %>% filter(!is.na(avg_log2FC))
    x4 <- x4 %>% filter(!is.infinite(avg_log2FC))
    
    x5 <- x5 %>% filter(!is.na(avg_log2FC))
    x5 <- x5 %>% filter(!is.infinite(avg_log2FC))
    
    
    r4 <- sign(x4$avg_log2FC)*(-log10(x4$p_val_adj))
    names(r4) <- x4$Gene
    
    r5 <- sign(x5$avg_log2FC)*(-log10(x5$p_val_adj))
    names(r5) <- x5$Gene
    
    r4 <- sort(r4, decreasing = TRUE)
    r5 <- sort(r5, decreasing = TRUE)
    
    
    assign(paste0(i,".gse.immediately.rank"), r4)
    assign(paste0(i,".gse.longterm.rank"), r5)
    
}

for (i in cell_list){
    x1 <- get(paste0(i, ".gse.immediately.padj"))
    
    
    for (j in c(1:2)){
    
        
    
        options(repr.plot.width=8, repr.plot.height=8)

        plot(x1)
        print(length(x1))
        
        
    }
    

    
    
}





for (i in cell_list){
    
    
    x4 <- get(paste0(i,".gse.immediately.FC"))
    x5 <- get(paste0(i,".gse.longterm.FC"))
    
    
    
    a4 <- fgsea(pathways = m_list.hu.hallmark,
               stats = x4)
    a5 <- fgsea(pathways = m_list.hu.hallmark,
               stats = x5)
    
    
    
    
    b4 <- fgsea(pathways = m_list.hu.C2,
               stats = x4)
    b5 <- fgsea(pathways = m_list.hu.C2,
               stats = x5)
    
    
    
    c4 <- fgsea(pathways = m_list.hu.C5,
               stats = x4)
    c5 <- fgsea(pathways = m_list.hu.C5,
               stats = x5)
    
    
    
    
    assign(paste0(i, ".fgsea.hallmark.immediately"), a4)
    assign(paste0(i, ".fgsea.hallmark.longterm"), a5)
    
    
    
    assign(paste0(i, ".fgsea.C2.immediately"), b4)
    assign(paste0(i, ".fgsea.C2.longterm"), b5)
    
    
    
    assign(paste0(i, ".fgsea.C5.immediately"), c4)
    assign(paste0(i, ".fgsea.C5.longterm"), c5)
    
}





for (i in cell_list){
    
    
    a4 <- get(paste0(i,".fgsea.hallmark.immediately"))
    a5 <- get(paste0(i,".fgsea.hallmark.longterm"))
    
    b4 <- get(paste0(i,".fgsea.C2.immediately"))
    b5 <- get(paste0(i,".fgsea.C2.longterm"))
    
    c4 <- get(paste0(i,".fgsea.C5.immediately"))
    c5 <- get(paste0(i,".fgsea.C5.longterm"))
    
    
    a4$timepoint <- "Immediately Post-flight"
    a5$timepoint <- "Long-term Post-flight"
    
    
    b4$timepoint <- "Immediately Post-flight"
    b5$timepoint <- "Long-term Post-flight"
    
    
    c4$timepoint <- "Immediately Post-flight"
    c5$timepoint <- "Long-term Post-flight"
    
    a4$celltype <- i
    a5$celltype <- i
    
    
    b4$celltype <- i
    b5$celltype <- i
    
    
    c4$celltype <- i
    c5$celltype <- i
    
    
    a4$sample <- "PBMC"
    a5$sample <- "PBMC"
    
    b4$sample <- "PBMC"
    b5$sample <- "PBMC"
    
    
    c4$sample <- "PBMC"
    c5$sample <- "PBMC"

    
    a4$info <- "Hallmark"
    a5$info <- "Hallmark"
    
    b4$info <- "C2"
    b5$info <- "C2"
    
    
    c4$info <- "C5"
    c5$info <- "C5"

    
    H.2 <- rbind(a4, a5)
    C2.2 <- rbind(b4, b5)
    C5.2 <- rbind(c4, c5)
    
    H.2 <- H.2 %>% filter(padj <= 0.2)
    
    
    
    
    C2.2.up <- C2.2 %>% filter(NES > 0)
    C2.2.down <- C2.2 %>% filter(NES < 0)
    
    C2.2.up <- C2.2.up[order(padj,)][1:10]
    C2.2.down <- C2.2.down[order(padj,)][1:10]
    
    
    C5.2.up <- C5.2 %>% filter(NES > 0)
    C5.2.down <- C5.2 %>% filter(NES < 0)
    
    C5.2.up <- C5.2.up[order(padj,)][1:10]
    C5.2.down <- C5.2.down[order(padj,)][1:10]
    
    
     
    assign(paste0(i,".hallmark.2"), H.2)
    assign(paste0(i,".C2.2.up"), C2.2.up)
    assign(paste0(i,".C2.2.down"), C2.2.down)
    assign(paste0(i,".C5.2.up"), C5.2.up)
    assign(paste0(i,".C5.2.down"), C5.2.down)
    
    
}






hallmark.2 <- rbind(pbmc.hallmark.2, CD4_T.hallmark.2, CD8_T.hallmark.2, other_T.hallmark.2, B.hallmark.2,
                   NK.hallmark.2, CD14_Mono.hallmark.2, CD16_Mono.hallmark.2, DC.hallmark.2,
                   other.hallmark.2)


C2.2 <- rbind(pbmc.C2.2.up, CD4_T.C2.2.up, CD4_T.C2.2.down, CD8_T.C2.2.up, CD8_T.C2.2.down,
                other_T.C2.2.up, other_T.C2.2.down, B.C2.2.up, B.C2.2.down,
                   NK.C2.2.up, NK.C2.2.down, CD14_Mono.C2.2.up, CD14_Mono.C2.2.down,
                CD16_Mono.C2.2.up, CD16_Mono.C2.2.down, DC.C2.2.up, DC.C2.2.down,
                   other.C2.2.up, other.C2.2.down)



C5.2 <- rbind(pbmc.C5.2.up, CD4_T.C5.2.up, CD4_T.C5.2.down, CD8_T.C5.2.up, CD8_T.C5.2.down,
                other_T.C5.2.up, other_T.C5.2.down, B.C5.2.up, B.C5.2.down,
                   NK.C5.2.up, NK.C5.2.down, CD14_Mono.C5.2.up, CD14_Mono.C5.2.down,
                CD16_Mono.C5.2.up, CD16_Mono.C5.2.down, DC.C5.2.up, DC.C5.2.down,
                   other.C5.2.up, other.C5.2.down)






pbmc.sub.2 <- rbind(hallmark.2, C2.2, C5.2)


pbmc.sub.2$category_with_color2 <- pbmc.sub.2$pathway




pbmc.sub.2$info <- factor(pbmc.sub.2$info,
                            levels = c("Hallmark", "C2", "C5"))





pbmc.sub.2$celltype <- ifelse(pbmc.sub.2$celltype == "pbmc", "PBMC",
                              ifelse(pbmc.sub.2$celltype == "CD4_T", "CD4 T",
                                 ifelse(pbmc.sub.2$celltype == "CD8_T", "CD8 T",
                                       ifelse(pbmc.sub.2$celltype == "other_T", "other T",
                                             ifelse(pbmc.sub.2$celltype == "B", "B",
                                                   ifelse(pbmc.sub.2$celltype == "NK", "NK",
                                                         ifelse(pbmc.sub.2$celltype == "CD14_Mono", "CD14 Mono",
                                                               ifelse(pbmc.sub.2$celltype == "CD16_Mono", "CD16 Mono",
                                                                     ifelse(pbmc.sub.2$celltype == "DC", "DC", 'other')))))))))



pbmc.sub.2$celltype <- factor(pbmc.sub.2$celltype,
                             levels = c('PBMC', "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))




pbmc.2$celltype <- factor(pbmc.2$celltype,
                         levels = c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))





int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}

Category.labs <- c("Hallmark",
                   "C2",
                   "C5")

names(Category.labs) <- c("Hallmark",
                   "C2",
                   "C5")

sample.labs <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")

names(sample.labs) <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")




PBMC.2.hallmark <- pbmc.2 %>% filter(info == "Hallmark")



PBMC.2.hallmark <- pbmc.2 %>% filter(info == "Hallmark")
PBMC.2.C2 <- pbmc.2 %>% filter(info == "C2")
PBMC.2.C5 <- pbmc.2 %>% filter(info == "C5")






int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}



dGSEA.2.merge <- ggplot(pbmc.2 %>% filter(padj<0.05), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 18)
        )+ 
  scale_x_continuous(breaks = int_breaks) 


options(repr.plot.width=48, repr.plot.height=49)
dGSEA.2.merge

dGSEA.2.merge <- ggplot(pbmc.2 %>% filter(padj<0.05) %>% filter(info %in% c('Hallmark', 'C2')), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks) 


options(repr.plot.width=48, repr.plot.height=30)
dGSEA.2.merge




PBMC.2.hallmark <- pbmc.2 %>% filter(info == "Hallmark")
PBMC.2.C2 <- pbmc.2 %>% filter(info == "C2")
PBMC.2.C5 <- pbmc.2 %>% filter(info == "C5")






dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark

dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 28)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark





dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark %>% filter(padj < 0.05), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 28)
        )+ 
  scale_x_continuous(breaks = int_breaks)



options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark











dGSEA.2.merge.C2 <- ggplot(PBMC.2.C2, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)

options(repr.plot.width=35, repr.plot.height=20)
dGSEA.2.merge.C2

dGSEA.2.merge.C5 <- ggplot(PBMC.2.C5, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)

options(repr.plot.width=40, repr.plot.height=20)
dGSEA.2.merge.C5











for (i in cell_list){
    
    
    a4 <- get(paste0(i,".fgsea.hallmark.immediately"))
    a5 <- get(paste0(i,".fgsea.hallmark.longterm"))
    
    b4 <- get(paste0(i,".fgsea.C2.immediately"))
    b5 <- get(paste0(i,".fgsea.C2.longterm"))
    
    c4 <- get(paste0(i,".fgsea.C5.immediately"))
    c5 <- get(paste0(i,".fgsea.C5.longterm"))
    
    
    a4$timepoint <- "Immediately Post-flight"
    a5$timepoint <- "Long-term Post-flight"
    
    
    b4$timepoint <- "Immediately Post-flight"
    b5$timepoint <- "Long-term Post-flight"
    
    
    c4$timepoint <- "Immediately Post-flight"
    c5$timepoint <- "Long-term Post-flight"
    
    a4$celltype <- i
    a5$celltype <- i
    
    
    b4$celltype <- i
    b5$celltype <- i
    
    
    c4$celltype <- i
    c5$celltype <- i
    
    
    a4$sample <- "PBMC"
    a5$sample <- "PBMC"
    
    b4$sample <- "PBMC"
    b5$sample <- "PBMC"
    
    
    c4$sample <- "PBMC"
    c5$sample <- "PBMC"

    
    a4$info <- "Hallmark"
    a5$info <- "Hallmark"
    
    b4$info <- "C2"
    b5$info <- "C2"
    
    
    c4$info <- "C5"
    c5$info <- "C5"

    
    H.2 <- rbind(a4, a5)
    C2.2 <- rbind(b4, b5)
    C5.2 <- rbind(c4, c5)
    
    H.2 <- H.2 %>% filter(padj <= 0.2)
    
    
    
    
    C2.2.up <- C2.2 %>% filter(NES > 0)
    C2.2.down <- C2.2 %>% filter(NES < 0)
    
    C2.2.up <- C2.2.up[order(padj,)][1:7]
    C2.2.down <- C2.2.down[order(padj,)][1:7]
    
    
    C5.2.up <- C5.2 %>% filter(NES > 0)
    C5.2.down <- C5.2 %>% filter(NES < 0)
    
    C5.2.up <- C5.2.up[order(padj,)][1:7]
    C5.2.down <- C5.2.down[order(padj,)][1:7]
    
    
     
    assign(paste0(i,".hallmark.2"), H.2)
    assign(paste0(i,".C2.2.up"), C2.2.up)
    assign(paste0(i,".C2.2.down"), C2.2.down)
    assign(paste0(i,".C5.2.up"), C5.2.up)
    assign(paste0(i,".C5.2.down"), C5.2.down)
    
    
}






hallmark.2 <- rbind(pbmc.hallmark.2, CD4_T.hallmark.2, CD8_T.hallmark.2, other_T.hallmark.2, B.hallmark.2,
                   NK.hallmark.2, CD14_Mono.hallmark.2, CD16_Mono.hallmark.2, DC.hallmark.2,
                   other.hallmark.2)


C2.2 <- rbind(pbmc.C2.2.up, CD4_T.C2.2.up, CD4_T.C2.2.down, CD8_T.C2.2.up, CD8_T.C2.2.down,
                other_T.C2.2.up, other_T.C2.2.down, B.C2.2.up, B.C2.2.down,
                   NK.C2.2.up, NK.C2.2.down, CD14_Mono.C2.2.up, CD14_Mono.C2.2.down,
                CD16_Mono.C2.2.up, CD16_Mono.C2.2.down, DC.C2.2.up, DC.C2.2.down,
                   other.C2.2.up, other.C2.2.down)



C5.2 <- rbind(pbmc.C5.2.up, CD4_T.C5.2.up, CD4_T.C5.2.down, CD8_T.C5.2.up, CD8_T.C5.2.down,
                other_T.C5.2.up, other_T.C5.2.down, B.C5.2.up, B.C5.2.down,
                   NK.C5.2.up, NK.C5.2.down, CD14_Mono.C5.2.up, CD14_Mono.C5.2.down,
                CD16_Mono.C5.2.up, CD16_Mono.C5.2.down, DC.C5.2.up, DC.C5.2.down,
                   other.C5.2.up, other.C5.2.down)






pbmc.sub.2 <- rbind(hallmark.2, C2.2, C5.2)


pbmc.sub.2$category_with_color2 <- pbmc.sub.2$pathway




pbmc.sub.2$info <- factor(pbmc.sub.2$info,
                            levels = c("Hallmark", "C2", "C5"))





pbmc.sub.2$celltype <- ifelse(pbmc.sub.2$celltype == "pbmc", "PBMC",
                              ifelse(pbmc.sub.2$celltype == "CD4_T", "CD4 T",
                                 ifelse(pbmc.sub.2$celltype == "CD8_T", "CD8 T",
                                       ifelse(pbmc.sub.2$celltype == "other_T", "other T",
                                             ifelse(pbmc.sub.2$celltype == "B", "B",
                                                   ifelse(pbmc.sub.2$celltype == "NK", "NK",
                                                         ifelse(pbmc.sub.2$celltype == "CD14_Mono", "CD14 Mono",
                                                               ifelse(pbmc.sub.2$celltype == "CD16_Mono", "CD16 Mono",
                                                                     ifelse(pbmc.sub.2$celltype == "DC", "DC", 'other')))))))))



pbmc.sub.2$celltype <- factor(pbmc.sub.2$celltype,
                             levels = c('PBMC', "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))




pbmc.2 <- rbind(pbmc.sub.2)




pbmc.2$celltype <- factor(pbmc.2$celltype,
                         levels = c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))





int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}

Category.labs <- c("Hallmark",
                   "C2",
                   "C5")

names(Category.labs) <- c("Hallmark",
                   "C2",
                   "C5")

sample.labs <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")

names(sample.labs) <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")






PBMC.2.hallmark <- pbmc.2 %>% filter(info == "Hallmark")
PBMC.2.C2 <- pbmc.2 %>% filter(info == "C2")
PBMC.2.C5 <- pbmc.2 %>% filter(info == "C5")










int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}



dGSEA.2.merge <- ggplot(pbmc.2, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 18)
        )+ 
  scale_x_continuous(breaks = int_breaks) 


options(repr.plot.width=48, repr.plot.height=49)
dGSEA.2.merge






PBMC.2.hallmark <- pbmc.2 %>% filter(info == "Hallmark")
PBMC.2.C2 <- pbmc.2 %>% filter(info == "C2")
PBMC.2.C5 <- pbmc.2 %>% filter(info == "C5")






dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark

dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 28)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark

dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark %>% filter(padj < 0.05), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 28)
        )+ 
  scale_x_continuous(breaks = int_breaks)



options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark











dGSEA.2.merge.C2 <- ggplot(PBMC.2.C2, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=20)
dGSEA.2.merge.C2

dGSEA.2.merge.C5 <- ggplot(PBMC.2.C5, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=40, repr.plot.height=20)
dGSEA.2.merge.C5













for (i in cell_list){
    
    
    a4 <- get(paste0(i,".fgsea.hallmark.immediately"))
    a5 <- get(paste0(i,".fgsea.hallmark.longterm"))
    
    b4 <- get(paste0(i,".fgsea.C2.immediately"))
    b5 <- get(paste0(i,".fgsea.C2.longterm"))
    
    c4 <- get(paste0(i,".fgsea.C5.immediately"))
    c5 <- get(paste0(i,".fgsea.C5.longterm"))
    
    
    a4$timepoint <- "Immediately Post-flight"
    a5$timepoint <- "Long-term Post-flight"
    
    
    b4$timepoint <- "Immediately Post-flight"
    b5$timepoint <- "Long-term Post-flight"
    
    
    c4$timepoint <- "Immediately Post-flight"
    c5$timepoint <- "Long-term Post-flight"
    
    a4$celltype <- i
    a5$celltype <- i
    
    
    b4$celltype <- i
    b5$celltype <- i
    
    
    c4$celltype <- i
    c5$celltype <- i
    
    
    a4$sample <- "PBMC"
    a5$sample <- "PBMC"
    
    b4$sample <- "PBMC"
    b5$sample <- "PBMC"
    
    
    c4$sample <- "PBMC"
    c5$sample <- "PBMC"

    
    a4$info <- "Hallmark"
    a5$info <- "Hallmark"
    
    b4$info <- "C2"
    b5$info <- "C2"
    
    
    c4$info <- "C5"
    c5$info <- "C5"

    
    H.2 <- rbind(a4, a5)
    C2.2 <- rbind(b4, b5)
    C5.2 <- rbind(c4, c5)
    
    H.2 <- H.2 %>% filter(padj <= 0.2)
    
    
    
    
    C2.2.up <- C2.2 %>% filter(NES > 0)
    C2.2.down <- C2.2 %>% filter(NES < 0)
    
    C2.2.up <- C2.2.up[order(padj,)][1:5]
    C2.2.down <- C2.2.down[order(padj,)][1:5]
    
    
    C5.2.up <- C5.2 %>% filter(NES > 0)
    C5.2.down <- C5.2 %>% filter(NES < 0)
    
    C5.2.up <- C5.2.up[order(padj,)][1:5]
    C5.2.down <- C5.2.down[order(padj,)][1:5]
    
    
     
    assign(paste0(i,".hallmark.2"), H.2)
    assign(paste0(i,".C2.2.up"), C2.2.up)
    assign(paste0(i,".C2.2.down"), C2.2.down)
    assign(paste0(i,".C5.2.up"), C5.2.up)
    assign(paste0(i,".C5.2.down"), C5.2.down)
    
    
}






hallmark.2 <- rbind(pbmc.hallmark.2, CD4_T.hallmark.2, CD8_T.hallmark.2, other_T.hallmark.2, B.hallmark.2,
                   NK.hallmark.2, CD14_Mono.hallmark.2, CD16_Mono.hallmark.2, DC.hallmark.2,
                   other.hallmark.2)


C2.2 <- rbind(pbmc.C2.2.up, CD4_T.C2.2.up, CD4_T.C2.2.down, CD8_T.C2.2.up, CD8_T.C2.2.down,
                other_T.C2.2.up, other_T.C2.2.down, B.C2.2.up, B.C2.2.down,
                   NK.C2.2.up, NK.C2.2.down, CD14_Mono.C2.2.up, CD14_Mono.C2.2.down,
                CD16_Mono.C2.2.up, CD16_Mono.C2.2.down, DC.C2.2.up, DC.C2.2.down,
                   other.C2.2.up, other.C2.2.down)



C5.2 <- rbind(pbmc.C5.2.up, CD4_T.C5.2.up, CD4_T.C5.2.down, CD8_T.C5.2.up, CD8_T.C5.2.down,
                other_T.C5.2.up, other_T.C5.2.down, B.C5.2.up, B.C5.2.down,
                   NK.C5.2.up, NK.C5.2.down, CD14_Mono.C5.2.up, CD14_Mono.C5.2.down,
                CD16_Mono.C5.2.up, CD16_Mono.C5.2.down, DC.C5.2.up, DC.C5.2.down,
                   other.C5.2.up, other.C5.2.down)






pbmc.sub.2 <- rbind(hallmark.2, C2.2, C5.2)


pbmc.sub.2$category_with_color2 <- pbmc.sub.2$pathway




pbmc.sub.2$info <- factor(pbmc.sub.2$info,
                            levels = c("Hallmark", "C2", "C5"))





pbmc.sub.2$celltype <- ifelse(pbmc.sub.2$celltype == "pbmc", "PBMC",
                              ifelse(pbmc.sub.2$celltype == "CD4_T", "CD4 T",
                                 ifelse(pbmc.sub.2$celltype == "CD8_T", "CD8 T",
                                       ifelse(pbmc.sub.2$celltype == "other_T", "other T",
                                             ifelse(pbmc.sub.2$celltype == "B", "B",
                                                   ifelse(pbmc.sub.2$celltype == "NK", "NK",
                                                         ifelse(pbmc.sub.2$celltype == "CD14_Mono", "CD14 Mono",
                                                               ifelse(pbmc.sub.2$celltype == "CD16_Mono", "CD16 Mono",
                                                                     ifelse(pbmc.sub.2$celltype == "DC", "DC", 'other')))))))))



pbmc.sub.2$celltype <- factor(pbmc.sub.2$celltype,
                             levels = c('PBMC', "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))






pbmc.2 <- rbind(pbmc.sub.2)




pbmc.2$celltype <- factor(pbmc.2$celltype,
                         levels = c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))





int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}

Category.labs <- c("Hallmark",
                   "C2",
                   "C5")

names(Category.labs) <- c("Hallmark",
                   "C2",
                   "C5")

sample.labs <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")

names(sample.labs) <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")






PBMC.2.hallmark <- pbmc.2 %>% filter(info == "Hallmark")
PBMC.2.C2 <- pbmc.2 %>% filter(info == "C2")
PBMC.2.C5 <- pbmc.2 %>% filter(info == "C5")










int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}



dGSEA.2.merge <- ggplot(pbmc.2, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 18)
        )+ 
  scale_x_continuous(breaks = int_breaks) 


options(repr.plot.width=48, repr.plot.height=49)
dGSEA.2.merge






PBMC.2.hallmark <- pbmc.2 %>% filter(info == "Hallmark")
PBMC.2.C2 <- pbmc.2 %>% filter(info == "C2")
PBMC.2.C5 <- pbmc.2 %>% filter(info == "C5")






dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark

dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 28)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark

dGSEA.2.merge.hallmark <- ggplot(PBMC.2.hallmark %>% filter(padj < 0.05), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 28)
        )+ 
  scale_x_continuous(breaks = int_breaks)



options(repr.plot.width=35, repr.plot.height=15)
dGSEA.2.merge.hallmark











dGSEA.2.merge.C2 <- ggplot(PBMC.2.C2, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=35, repr.plot.height=20)
dGSEA.2.merge.C2

dGSEA.2.merge.C5 <- ggplot(PBMC.2.C5, aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,15,17,18,21))+
  ylab(NULL) +
  geom_vline(xintercept=0, col="black") +
  guides(shape = guide_legend(override.aes = list(size =10))) +
#  labs(title = "MitoPathways 3.0"
       #    ,subtitle = "FDR cutoff < 0.25"
       # ,caption = "FDR cutoff = 0.25"
 # ) +
  theme(plot.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=20)) +
  facet_grid(info ~ celltype, #space = "free", 
             scales = "free_y",
             labeller = labeller(info = Category.labs, sample = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 20)
        )+ 
  scale_x_continuous(breaks = int_breaks)


options(repr.plot.width=40, repr.plot.height=20)
dGSEA.2.merge.C5














