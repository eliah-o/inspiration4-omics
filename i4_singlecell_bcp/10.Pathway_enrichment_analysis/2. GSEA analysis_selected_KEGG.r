library(ggplot2)
library(dplyr)
library(tidyr)

library(fgsea)
library(ggplot2)

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
    
    y4 <- y4 %>% filter(p_val_adj < 0.05)
    y5 <- y5 %>% filter(p_val_adj < 0.05)
    
    
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





pathway_fine <- readRDS("path_to_selected_KEGG_pathways")





for (i in cell_list){
    
    x4 <- get(paste0(i,".gse.immediately.FC"))
    x5 <- get(paste0(i,".gse.longterm.FC"))
    
    
    a4 <- fgsea(pathways = pathway_fine,
               stats = x4)
    a5 <- fgsea(pathways = pathway_fine,
               stats = x5)
    
    
   
    
    assign(paste0(i, ".fgsea.KEGG.immediately"), a4)
    assign(paste0(i, ".fgsea.KEGG.longterm"), a5)
    
   
    
}





for (i in cell_list){
    
    a4 <- get(paste0(i,".fgsea.KEGG.immediately"))
    a5 <- get(paste0(i,".fgsea.KEGG.longterm"))
    
  
    a4$timepoint <- "Immediately Post-flight"
    a5$timepoint <- "Long-term Post-flight"
  
    a4$celltype <- i
    a5$celltype <- i
    
    
    a4$info <- "C2"
    a5$info <- "C2"
    
    a4$sample <- "PBMC"
    a5$sample <- "PBMC"
    
 
    H.2 <- rbind(a4, a5)
    
    assign(paste0(i,".KEGG.2"), H.2)

    
}






KEGG.2 <- rbind(pbmc.KEGG.2, CD4_T.KEGG.2, CD8_T.KEGG.2, other_T.KEGG.2, B.KEGG.2,
                   NK.KEGG.2, CD14_Mono.KEGG.2, CD16_Mono.KEGG.2, DC.KEGG.2,
                   other.KEGG.2)







pbmc.sub.2 <- rbind(KEGG.2)


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





pbmc.2.KEGG.1 <- pbmc.2 %>% filter(pathway %in% c('KEGG_ALLOGRAFT_REJECTION', 'KEGG_ALZHEIMERS_DISEASE',"KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
                                                  'KEGG_ASTHMA', 'KEGG_AUTOIMMUNE_THYROID_DISEASE', 'KEGG_GRAFT_VERSUS_HOST_DISEASE',
                                                  'KEGG_LEISHMANIA_INFECTION', 'KEGG_PARKINSONS_DISEASE', "KEGG_PATHOGENIC_ESCHERICHIA_COLI_INFECTION",
                                                  'KEGG_PRION_DISEASES', 'KEGG_TYPE_I_DIABETES_MELLITUS', 'KEGG_VIRAL_MYOCARDITIS', names(pathway_fine)))





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


Category.labs <- c("Hallmark",
                   "C2",
                   "C5")

names(Category.labs) <- c("Hallmark",
                   "KEGG",
                   "C5")

sample.labs <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")

names(sample.labs) <- c("PBMC", "CD4 T", "CD8 T", "other T", "B", "NK", "CD14 Mono", "CD16 Mono",
                                            "DC", "other")






int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
}





dGSEA.1.merge <- ggplot(pbmc.2.KEGG.1 %>% filter(padj < 0.05), aes(x = NES, y = fct_reorder(category_with_color2, NES))) + 
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
        ,axis.text.y = element_markdown(size = 28)
        )+ 
  scale_x_continuous(breaks = int_breaks) 




options(repr.plot.width=40, repr.plot.height=10)
dGSEA.1.merge








