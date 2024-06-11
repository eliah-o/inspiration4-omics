library(ggplot2)
library(dplyr)
library(tidyr)

library(fgsea)
library(DESeq2)
library(tidyverse)

library(RColorBrewer)
library(pheatmap)

library(ggrepel)


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

cell_list <- c("CD4_T", "CD8_T", "other_T", "B", "NK", 
              "CD14_Mono", "CD16_Mono", "DC", 'other')

all_list <- c('pbmc', cell_list)

for (i in all_list){
    
    y1 <- read.csv(paste0("path_to.Female.DEG.all.csv"))
    y2 <- read.csv(paste0("path_to.Male.DEG.all.csv"))
    
    colnames(y1)[1] <- "Gene"
    colnames(y2)[1] <- "Gene"

    
    y1$timepoint <- "Immediately Post-flight Female"
    y2$timepoint <- "Immediately Post-flight Male"
    
   
    
    assign(paste0(i,".immediately.FC.F"), y1)
    assign(paste0(i,".immediately.FC.M"), y2)

}

for (i in all_list){
    
    x4 <- get(paste0(i,".immediately.FC.F"))
    x5 <- get(paste0(i,".immediately.FC.M"))
    
    
    x4 <- x4 %>% filter(!is.na(avg_log2FC))
    x4 <- x4 %>% filter(!is.infinite(avg_log2FC))

    x5 <- x5 %>% filter(!is.na(avg_log2FC))
    x5 <- x5 %>% filter(!is.infinite(avg_log2FC))

    
    y4 <- x4$avg_log2FC
    y5 <- x5$avg_log2FC
    
    names(y4) <- x4$Gene
    names(y5) <- x5$Gene
    
    assign(paste0(i,".gse.immediately.FC.F"), y4)
    assign(paste0(i,".gse.immediately.FC.M"), y5)
    
}

for (i in all_list){
    
    x4 <- get(paste0(i,".gse.immediately.FC.F"))
    x5 <- get(paste0(i,".gse.immediately.FC.M"))
    
    
    
    a4 <- fgsea(pathways = m_list.hu.hallmark,
               stats = x4)
    a5 <- fgsea(pathways = m_list.hu.hallmark,
               stats = x5)
    
    assign(paste0(i, ".fgsea.hallmark.immediately.F"), a4)
    assign(paste0(i, ".fgsea.hallmark.immediately.M"), a5)
    
}

for (i in all_list){
    
    a4 <- get(paste0(i,".fgsea.hallmark.immediately.F"))
    a5 <- get(paste0(i,".fgsea.hallmark.immediately.M"))
    
    a4$timepoint <- "Immediately Post-flight Female"
    a5$timepoint <- "Immediately Post-flight Male"
    
    
    a4$celltype <- i
    a5$celltype <- i

    a4$sample <- "PBMCs"
    a5$sample <- "PBMCs"

    a4$info <- "Hallmark"
    a5$info <- "Hallmark"
    
    a4$sex <- "Female"
    a5$sex <- "Male"
    
    
    H.2 <- rbind(a4, a5)
    
    
    assign(paste0(i,".hallmark.2"), H.2)
    
}


Hallmark.2 <- rbind(pbmc.hallmark.2, CD4_T.hallmark.2, CD8_T.hallmark.2, other_T.hallmark.2,
                   B.hallmark.2, NK.hallmark.2, CD14_Mono.hallmark.2, CD16_Mono.hallmark.2,
                   DC.hallmark.2, other.hallmark.2)



Hallmark.2$category_with_color2 <- Hallmark.2$pathway


Hallmark.2$celltype <- ifelse(Hallmark.2$celltype == "CD4_T", "CD4 T",
                                 ifelse(Hallmark.2$celltype == "CD8_T", "CD8 T",
                                       ifelse(Hallmark.2$celltype == "other_T", "other T",
                                             ifelse(Hallmark.2$celltype == "B", "B",
                                                   ifelse(Hallmark.2$celltype == "NK", "NK",
                                                         ifelse(Hallmark.2$celltype == "CD14_Mono", "CD14 Mono",
                                                               ifelse(Hallmark.2$celltype == "CD16_Mono", "CD16 Mono",
                                                                     ifelse(Hallmark.2$celltype == "DC", "DC", 
                                                                           ifelse(Hallmark.2$celltype == "pbmc", "PBMCs", "other")))))))))


Hallmark.2$celltype <- factor(Hallmark.2$celltype,
                             levels = c("PBMCs", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))








int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}

Category.labs <- c("Hallmark")

names(Category.labs) <- c("Hallmark")

sample.labs <- c("PBMCs", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")

names(sample.labs) <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")


dGSEA.1.merge <- ggplot(Hallmark.2 %>% filter(padj <=0.3), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
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
dGSEA.1.merge





dGSEA.1.merge <- ggplot(Hallmark.2 %>% filter(padj <=0.2), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25, 15,17,18,21))+
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
dGSEA.1.merge





dGSEA.1.merge <- ggplot(Hallmark.2 %>% filter(padj <=0.1), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25, 15,17,18,21))+
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
dGSEA.1.merge







dGSEA.1.merge <- ggplot(Hallmark.2 %>% filter(padj <=0.05), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25, 15,17,18,21))+
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
dGSEA.1.merge






