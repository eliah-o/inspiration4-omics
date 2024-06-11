library(fgsea)
library(ggplot2)

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

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)

library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(ggprism)
library(ggsignif)
library(ggrepel)
library(tidyverse)

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








temp.sobj <- readRDS("path_to_seurat_object")

temp.sobj$predicted.id <- ifelse(temp.sobj$celltype == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype == 'CD14 Mono', 'CD14_Mono',
                                                                 ifelse(temp.sobj$celltype == 'CD16 Mono', 'CD16_Mono',
                                                                         ifelse(temp.sobj$celltype == 'DC', 'DC','other'))))))))



temp.sobj$Timepoint <- ifelse(temp.sobj$timepoint == "June 2021", "Pre-flight",
                             ifelse(temp.sobj$timepoint == "August 2021", "Pre-flight",
                                   ifelse(temp.sobj$timepoint == "September Pre-launch", "Pre-flight",
                                         ifelse(temp.sobj$timepoint == "September Post-launch", "R+1",
                                               ifelse(temp.sobj$timepoint == "November 2021", "R+45", "R+82")))))



cell_list <- c('pbmc', "CD4_T", "CD8_T", "other_T", "B", "NK", 
              "CD14_Mono", "CD16_Mono", "DC", 'other')





for (i in cell_list){
    
    
    y4 <- read.csv(paste0("path_to_DEG_September_Post_JAS.filter.csv"))
    y5 <- read.csv(paste0("path_to_DEG_Nov+Dec_September_Post.filter.csv"))
    
    
    colnames(y4)[1] <- "Gene"
    colnames(y5)[1] <- "Gene"
    
    
    y4$timepoint <- "Immediately Post-flight"
    y5$timepoint <- "Long-term Post-flight"
    
    
    
    assign(paste0(i,".immediately.FC"), y4)
    assign(paste0(i,".longterm.FC"), y5)
    
    
}

T_cells <- c("CD4_Naive",  
             "CD4_TCM", "CD4_TEM", "Treg", "CD8_Naive",
           "CD8_TCM", "CD8_TEM", "dnT",
            "gdT", "MAIT")

for (i in T_cells){
    
    y1 <- read.csv(paste0("path_to_DEG_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("path_to_DEG_Nov+Dec_September_Post.filter.csv"))
    
    
    colnames(y1)[1] <- "Gene"
    colnames(y2)[1] <- "Gene"

    
    y1$timepoint <- "Immediately Post-flight"
    y2$timepoint <- "Long-term Post-flight"

    write.csv(y1, paste0(i,".immediately.FC.tp.csv"))
    write.csv(y2, paste0(i,".longterm.FC.tp.csv"))
    
    assign(paste0(i,".immediately.FC"), y1)
    assign(paste0(i,".longterm.FC"), y2)

}

B_cells <- c("B_Naive", "B_intermediate", "B_intermediate", "B_memory")

for (i in B_cells){
    
    y1 <- read.csv(paste0("path_to_DEG_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("path_to_DEG_Nov+Dec_September_Post.filter.csv"))    
    
    
    colnames(y1)[1] <- "Gene"
    colnames(y2)[1] <- "Gene"

    
    y1$timepoint <- "Immediately Post-flight"
    y2$timepoint <- "Long-term Post-flight"

    write.csv(y1, paste0(i,".immediately.FC.tp.csv"))
    write.csv(y2, paste0(i,".longterm.FC.tp.csv"))
    
    assign(paste0(i,".immediately.FC"), y1)
    assign(paste0(i,".longterm.FC"), y2)

}

for (i in 'NK_CD56bright'){
    
    y1 <- read.csv(paste0("path_to_DEG_September_Post_JAS.filter.csv"))
    y2 <- read.csv(paste0("path_to_DEG_Nov+Dec_September_Post.filter.csv"))    
    
    
    colnames(y1)[1] <- "Gene"
    colnames(y2)[1] <- "Gene"

    
    y1$timepoint <- "Immediately Post-flight"
    y2$timepoint <- "Long-term Post-flight"

    write.csv(y1, paste0(i,".immediately.FC.tp.csv"))
    write.csv(y2, paste0(i,".longterm.FC.tp.csv"))
    
    assign(paste0(i,".immediately.FC"), y1)
    assign(paste0(i,".longterm.FC"), y2)

}







FOXP3 <- m_list.hu$ZHENG_BOUND_BY_FOXP3

Tcf21 <- m_list.hu$CUI_TCF21_TARGETS_2_DN







BCR <- m_list.hu$SIG_BCR_SIGNALING_PATHWAY



temp2.sobj <- temp.sobj



DefaultAssay(temp2.sobj) <- 'SCT'

temp2.sobj <- AddModuleScore(
    object = temp2.sobj,
    features = list(FOXP3),
    name = "FOXP3 targets")

DefaultAssay(temp2.sobj) <- 'SCT'

temp2.sobj <- AddModuleScore(
    object = temp2.sobj,
    features = list(Tcf21),
    name = "Tcf21 targets")



Idents(temp.sobj) <- 'predicted.id'
B.sobj <- subset(temp.sobj, ident="B")


Idents(temp.sobj) <- 'celltype.l2'
B_intermediate.sobj <- subset(temp.sobj, ident="B intermediate")

Idents(temp.sobj) <- 'celltype.l2'
B_memory.sobj <- subset(temp.sobj, ident="B memory")

Idents(temp.sobj) <- 'celltype.l2'
B_Naive.sobj <- subset(temp.sobj, ident="B naive")


Idents(temp.sobj) <- 'celltype.l2'
plasmablast.sobj <- subset(temp.sobj, ident="Plasmablast")


Idents(temp.sobj) <- 'predicted.id'
NK.sobj <- subset(temp.sobj, ident="NK")

Idents(temp.sobj) <- 'celltype.l2'
NK_CD56bright.sobj <- subset(temp.sobj, ident="NK_CD56bright")


Idents(temp.sobj) <- 'predicted.id'
CD4_T.sobj <- subset(temp.sobj, ident="CD4_T")

Idents(temp.sobj) <- 'predicted.id'
CD8_T.sobj <- subset(temp.sobj, ident="CD8_T")

Idents(temp.sobj) <- 'predicted.id'
other_T.sobj <- subset(temp.sobj, ident="other_T")

Idents(temp.sobj) <- 'celltype.l2'
CD4_Naive.sobj <- subset(temp.sobj, ident="CD4 Naive")

Idents(temp.sobj) <- 'celltype.l2'
CD4_TCM.sobj <- subset(temp.sobj, ident="CD4 TCM")

Idents(temp.sobj) <- 'celltype.l2'
CD4_TEM.sobj <- subset(temp.sobj, ident="CD4 TEM")

Idents(temp.sobj) <- 'celltype.l2'
Treg.sobj <- subset(temp.sobj, ident="Treg")

Idents(temp.sobj) <- 'celltype.l2'
CD8_Naive.sobj <- subset(temp.sobj, ident="CD8 Naive")

Idents(temp.sobj) <- 'celltype.l2'
CD8_TCM.sobj <- subset(temp.sobj, ident="CD8 TCM")

Idents(temp.sobj) <- 'celltype.l2'
CD8_TEM.sobj <- subset(temp.sobj, ident="CD8 TEM")

Idents(temp.sobj) <- 'celltype.l2'
dnT.sobj <- subset(temp.sobj, ident="dnT")

Idents(temp.sobj) <- 'celltype.l2'
gdT.sobj <- subset(temp.sobj, ident="gdT")


Idents(temp.sobj) <- 'celltype.l2'
MAIT.sobj <- subset(temp.sobj, ident="MAIT")





Idents(temp2.sobj) <- 'predicted.id'
B2.sobj <- subset(temp2.sobj, ident="B")


Idents(temp2.sobj) <- 'celltype.l2'
B_intermediate2.sobj <- subset(temp2.sobj, ident="B intermediate")

Idents(temp2.sobj) <- 'celltype.l2'
B_memory2.sobj <- subset(temp2.sobj, ident="B memory")

Idents(temp2.sobj) <- 'celltype.l2'
B_Naive2.sobj <- subset(temp2.sobj, ident="B naive")


Idents(temp2.sobj) <- 'celltype.l2'
plasmablast2.sobj <- subset(temp2.sobj, ident="Plasmablast")


Idents(temp2.sobj) <- 'predicted.id'
NK2.sobj <- subset(temp2.sobj, ident="NK")

Idents(temp2.sobj) <- 'celltype.l2'
NK_CD56bright2.sobj <- subset(temp2.sobj, ident="NK_CD56bright")


Idents(temp2.sobj) <- 'predicted.id'
CD4_T2.sobj <- subset(temp2.sobj, ident="CD4_T")

Idents(temp2.sobj) <- 'predicted.id'
CD8_T2.sobj <- subset(temp2.sobj, ident="CD8_T")

Idents(temp2.sobj) <- 'predicted.id'
other_T2.sobj <- subset(temp2.sobj, ident="other_T")

Idents(temp2.sobj) <- 'celltype.l2'
CD4_Naive2.sobj <- subset(temp2.sobj, ident="CD4 Naive")

Idents(temp2.sobj) <- 'celltype.l2'
CD4_TCM2.sobj <- subset(temp2.sobj, ident="CD4 TCM")

Idents(temp2.sobj) <- 'celltype.l2'
CD4_TEM2.sobj <- subset(temp2.sobj, ident="CD4 TEM")

Idents(temp2.sobj) <- 'celltype.l2'
Treg2.sobj <- subset(temp2.sobj, ident="Treg")

Idents(temp2.sobj) <- 'celltype.l2'
CD8_Naive2.sobj <- subset(temp2.sobj, ident="CD8 Naive")

Idents(temp2.sobj) <- 'celltype.l2'
CD8_TCM2.sobj <- subset(temp2.sobj, ident="CD8 TCM")

Idents(temp2.sobj) <- 'celltype.l2'
CD8_TEM2.sobj <- subset(temp2.sobj, ident="CD8 TEM")

Idents(temp2.sobj) <- 'celltype.l2'
dnT2.sobj <- subset(temp2.sobj, ident="dnT")

Idents(temp2.sobj) <- 'celltype.l2'
gdT2.sobj <- subset(temp2.sobj, ident="gdT")


Idents(temp2.sobj) <- 'celltype.l2'
MAIT2.sobj <- subset(temp2.sobj, ident="MAIT")



T <- c("CD4_T", "CD8_T", "other_T", T_cells)







for (i in T){
    x <- get(paste0(i,"2.sobj"))
    
    
    
    colnames(x@meta.data)[33] <- 'FOXP3 targets'
    colnames(x@meta.data)[34] <- 'Tcf21 targets'
    
    DefaultAssay(x) <- "SCT"
    Idents(x) <- "Timepoint"
    
    options(repr.plot.width=8, repr.plot.height=5)
    a1<-DotPlot(x, 
      features=c("FOXP3 targets", "Tcf21 targets"))+ 
  theme(plot.title = element_text(size = 30), 
  axis.text = element_text(size=18),
  axis.title = element_text(size=18),
  axis.text.x = element_text(size=18, angle=-45),
  axis.text.y = element_text(size=18),
  legend.text = element_text(size=18)
)+
  coord_flip()+
    labs(x='', y='Timepoint')+
scale_color_gradient2(low="navy", mid="white", high="firebrick3", midpoint=0)
    
    

}

for (i in T){
    x <- get(paste0(i,"2.sobj"))
    
    
    
    colnames(x@meta.data)[33] <- 'FOXP3 targets'
    colnames(x@meta.data)[34] <- 'Tcf21 targets'
    
    DefaultAssay(x) <- "SCT"
    Idents(x) <- "Timepoint"
    
    options(repr.plot.width=8, repr.plot.height=4.5)
    a1<-DotPlot(x, 
      features=c("FOXP3 targets"))+ 
  theme(plot.title = element_text(size = 30), 
  axis.text = element_text(size=18),
  axis.title = element_text(size=18),
  axis.text.x = element_text(size=18, angle=-45),
  axis.text.y = element_text(size=18),
  legend.text = element_text(size=18)
)+
  coord_flip()+
    labs(x='', y='Timepoint')+
scale_color_gradient2(low="navy", mid="white", high="firebrick3", midpoint=0)
    
    

}

for (i in T){
    x <- get(paste0(i,"2.sobj"))
    
    
    
    colnames(x@meta.data)[33] <- 'FOXP3 targets'
    colnames(x@meta.data)[34] <- 'Tcf21 targets'
    
    DefaultAssay(x) <- "SCT"
    Idents(x) <- "Timepoint"
    
    options(repr.plot.width=8, repr.plot.height=4.5)
    a1<-DotPlot(x, 
      features=c("Tcf21 targets"))+ 
  theme(plot.title = element_text(size = 30), 
  axis.text = element_text(size=18),
  axis.title = element_text(size=18),
  axis.text.x = element_text(size=18, angle=-45),
  axis.text.y = element_text(size=18),
  legend.text = element_text(size=18)
)+
  coord_flip()+
    labs(x='', y='Timepoint')+
scale_color_gradient2(low="navy", mid="white", high="firebrick3", midpoint=0)
    
    

}







for (i in c('B', B_cells)){
    x <- get(paste0(i,"2.sobj"))
    
    
    
    colnames(x@meta.data)[33] <- 'FOXP3 targets'
    colnames(x@meta.data)[34] <- 'Tcf21 targets'
    
    DefaultAssay(x) <- "SCT"
    Idents(x) <- "Timepoint"
    
    options(repr.plot.width=8, repr.plot.height=5)
    a1<-DotPlot(x, 
      features=c("FOXP3 targets", "Tcf21 targets"))+ 
  theme(plot.title = element_text(size = 30), 
  axis.text = element_text(size=18),
  axis.title = element_text(size=18),
  axis.text.x = element_text(size=18, angle=-45),
  axis.text.y = element_text(size=18),
  legend.text = element_text(size=18)
)+
  coord_flip()+
    labs(x='', y='Timepoint')+
scale_color_gradient2(low="navy", mid="white", high="firebrick3", midpoint=0)
    
    

}

for (i in c('B', B_cells)){
    x <- get(paste0(i,"2.sobj"))
    
    
    
    colnames(x@meta.data)[33] <- 'FOXP3 targets'
    colnames(x@meta.data)[34] <- 'Tcf21 targets'
    
    DefaultAssay(x) <- "SCT"
    Idents(x) <- "Timepoint"
    
    options(repr.plot.width=8, repr.plot.height=4.5)
    a1<-DotPlot(x, 
      features=c("FOXP3 targets"))+ 
  theme(plot.title = element_text(size = 30), 
  axis.text = element_text(size=18),
  axis.title = element_text(size=18),
  axis.text.x = element_text(size=18, angle=-45),
  axis.text.y = element_text(size=18),
  legend.text = element_text(size=18)
)+
  coord_flip()+
    labs(x='', y='Timepoint')+
scale_color_gradient2(low="navy", mid="white", high="firebrick3", midpoint=0)
    
    

}

for (i in c('B', B_cells)){
    x <- get(paste0(i,"2.sobj"))
    
    
    
    colnames(x@meta.data)[33] <- 'FOXP3 targets'
    colnames(x@meta.data)[34] <- 'Tcf21 targets'
    
    DefaultAssay(x) <- "SCT"
    Idents(x) <- "Timepoint"
    
    options(repr.plot.width=8, repr.plot.height=4.5)
    a1<-DotPlot(x, 
      features=c("Tcf21 targets"))+ 
  theme(plot.title = element_text(size = 30), 
  axis.text = element_text(size=18),
  axis.title = element_text(size=18),
  axis.text.x = element_text(size=18, angle=-45),
  axis.text.y = element_text(size=18),
  legend.text = element_text(size=18)
)+
  coord_flip()+
    labs(x='', y='Timepoint')+
scale_color_gradient2(low="navy", mid="white", high="firebrick3", midpoint=0)
    
    

}










