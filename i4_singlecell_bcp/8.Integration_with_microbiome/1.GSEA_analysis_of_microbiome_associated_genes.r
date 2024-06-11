library(ggplot2)
library(dplyr)
library(tidyr)

library(fgsea)

library(DESeq2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggrepel)
library(pheatmap)
library(AnnotationHub)
library(ensembldb)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)


library(ggthemes)
library(ggtext)
library(stringr)
library(forcats)
library("RColorBrewer")
library(ggpubr)
library(msigdbr)

m_df.hu = msigdbr(species = "Homo sapiens")

m_df.hu = msigdbr(species = "Homo sapiens")
m_df.hu.hallmark = m_df.hu %>% filter(gs_cat %in% c("H"))
m_df.hu.C2 = m_df.hu %>% filter(gs_cat %in% c("C2"))
m_df.hu.C5 = m_df.hu %>% filter(gs_cat %in% c("C5"))
m_df.hu.C6 = m_df.hu %>% filter(gs_cat %in% c("C6"))
m_df.hu.C7 = m_df.hu %>% filter(gs_cat %in% c("C7"))
m_df.hu.C8 = m_df.hu %>% filter(gs_cat %in% c("C8"))

m_list.hu = m_df.hu %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.hallmark = m_df.hu.hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C2 = m_df.hu.C2 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C5 = m_df.hu.C5 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C6 = m_df.hu.C6 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C7 = m_df.hu.C7 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C8 = m_df.hu.C8 %>% split(x = .$gene_symbol, f = .$gs_name)

target <- c("GOBP_LYMPHOCYTE_ACTIVATION","GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
            "GOBP_IMMUNE_RESPONSE","GOBP_T_CELL_ACTIVATION","GOBP_B_CELL_ACTIVATION",
            "GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_INTERSPECIES_INTERACTION_BETWEEN_ORGANISMS",
            "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE","GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY",
            "GOBP_CYTOKINE_PRODUCTION","GOBP_IMMUNE_EFFECTOR_PROCESS",
            "GOBP_DEFENSE_RESPONSE_TO_OTHER_ORGANISM","GOBP_INFLAMMATORY_RESPONSE",
            "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION","GOBP_VIRAL_PROCESS",
            "GOBP_RESPONSE_TO_VIRUS","GOBP_MACROPHAGE_ACTIVATION",
            "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
            "GOBP_NATURAL_KILLER_CELL_ACTIVATION","GOBP_RESPONSE_TO_BACTERIUM",
            "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN")





cell_list <- c('pbmc', 'CD4_T', 'CD8_T', 'other_T', 'B',
              'NK', 'CD14_Mono', 'CD16_Mono', 'DC', 'other')

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
    
    
    x4 <- get(paste0(i,".gse.immediately.FC"))
    x5 <- get(paste0(i,".gse.longterm.FC"))
    
    
    
    a4 <- fgsea(pathways = m_list.hu.C5[names(m_list.hu.C5) %in% target],
               stats = x4)
    a5 <- fgsea(pathways = m_list.hu.C5[names(m_list.hu.C5) %in% target],
               stats = x5)
    
    
    
    
    
    
    
    assign(paste0(i, ".fgsea.microbiome.immediately"), a4)
    assign(paste0(i, ".fgsea.microbiome.longterm"), a5)
    
    
    

    
}

for (i in cell_list){
    
    a4 <- get(paste0(i,".fgsea.microbiome.immediately"))
    a5 <- get(paste0(i,".fgsea.microbiome.longterm"))
    
    a4$timepoint <- "Immediately Post-flight"
    a5$timepoint <- "Long-term Post-flight"
    
    a4$celltype <- i
    a5$celltype <- i

    a4$sample <- "PBMC"
    a5$sample <- "PBMC"
    
    
    a4$info <- "i4 PBMC"
    a5$info <- "i4 PBMC"
    
    PBMC.sub.fine.1 <- rbind(a4, a5)
    
    assign(paste0(i,".PBMC.sub.fine.1"), PBMC.sub.fine.1)
    
    
}


PBMC.sub.2 <- rbind(pbmc.PBMC.sub.fine.1, CD4_T.PBMC.sub.fine.1, CD8_T.PBMC.sub.fine.1, other_T.PBMC.sub.fine.1, B.PBMC.sub.fine.1,
                   NK.PBMC.sub.fine.1, CD14_Mono.PBMC.sub.fine.1, CD16_Mono.PBMC.sub.fine.1,
                   DC.PBMC.sub.fine.1, other.PBMC.sub.fine.1)




PBMC_pathway <- rbind(PBMC.sub.2)

PBMC_pathway$category_with_color2 <- PBMC_pathway$pathway



PBMC_pathway$celltype <- ifelse(PBMC_pathway$celltype == "pbmc", "PBMC",
                                ifelse(PBMC_pathway$celltype == "CD4_T", "CD4 T",
                                 ifelse(PBMC_pathway$celltype == "CD8_T", "CD8 T",
                                       ifelse(PBMC_pathway$celltype == "other_T", "other T",
                                             ifelse(PBMC_pathway$celltype == "B", "B",
                                                   ifelse(PBMC_pathway$celltype == "NK", "NK",
                                                         ifelse(PBMC_pathway$celltype == "CD14_Mono", "CD14 Mono",
                                                               ifelse(PBMC_pathway$celltype == "CD16_Mono", "CD16 Mono",
                                                                     ifelse(PBMC_pathway$celltype == "DC", "DC", 
                                                                           ifelse(PBMC_pathway$celltype == 'other', 'other', 'PBMCs'))))))))))


PBMC_pathway$celltype <- factor(PBMC_pathway$celltype,
                             levels = c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other"))




PBMC_pathway$pathway <- factor(PBMC_pathway$pathway,
                              levels= target)




PBMC_pathway$category_with_color2 <- PBMC_pathway$pathway




int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}  



PBMC_pathway$info <- "Microbiome"

Category.labs <- c("Microbiome")

names(Category.labs) <- c("Microbiome")

sample.labs <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14\nMono", "CD16\nMono",
                                            "DC", "other")

names(sample.labs) <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")


PBMC_pathway$category_with_color2 <- factor(PBMC_pathway$category_with_color2,
                              levels= target)









dGSEA.1.merge <- ggplot(PBMC_pathway %>% filter(padj < 0.3), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,17,18,15,21))+
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
             labeller = labeller(info = Category.labs, celltype = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 30)
        )+ 
  scale_x_continuous(breaks = int_breaks) 

options(repr.plot.width=40, repr.plot.height=12)
dGSEA.1.merge







dGSEA.1.merge <- ggplot(PBMC_pathway %>% filter(padj < 0.2), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,17,18,15,21))+
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
             labeller = labeller(info = Category.labs, celltype = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 30)
        )+ 
  scale_x_continuous(breaks = int_breaks) 

options(repr.plot.width=40, repr.plot.height=12)
dGSEA.1.merge



dGSEA.1.merge <- ggplot(PBMC_pathway %>% filter(padj < 0.1), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,17,18,15,21))+
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
             labeller = labeller(info = Category.labs, celltype = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 30)
        )+ 
  scale_x_continuous(breaks = int_breaks) 

options(repr.plot.width=40, repr.plot.height=12)
dGSEA.1.merge



dGSEA.1.merge <- ggplot(PBMC_pathway %>% filter(padj < 0.05), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
  geom_point(aes(size = padj, color = NES, shape = timepoint, fill = NES)) +
  geom_segment( aes(x=0, xend=NES, y=category_with_color2, yend=category_with_color2, color= NES), size = 1.2) +
  theme_bw(base_size=20) +
  scale_color_gradient2(low="navy", mid="white", high="firebrick3", 
                        midpoint=0) +
  scale_fill_gradient2(low="navy", mid="white", high="firebrick3", 
                       midpoint=0) +
  scale_size(range = c(3,8), name="padj", trans = 'reverse') +
  scale_shape_manual(name = "Group", values=c(16,25,17,18,15,21))+
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
             labeller = labeller(info = Category.labs, celltype = sample.labs)) +
  theme(axis.text.x = element_text(size = 20),
        strip.text = element_text(face="bold", size=20),
        strip.background = element_rect(fill="lightblue", colour="black",size=1)
        ,axis.text.y = element_markdown(size = 30)
        )+ 
  scale_x_continuous(breaks = int_breaks) 

options(repr.plot.width=40, repr.plot.height=12)
dGSEA.1.merge






