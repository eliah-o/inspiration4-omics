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




library("GO.db")

library(GSVA)
library(GSEABase)

tbl=toTable(GOTERM)

T.list.GO <- c("GO:2000566","GO:0043372","GO:0043378","GO:0043382","GO:0042104","GO:0045585","GO:0045588","GO:2000516","GO:2000568","GO:1903905","GO:2000563","GO:1900281","GO:2001187","GO:2001193","GO:0033091","GO:0045591","GO:2000451","GO:2000454","GO:1905404")
T.list.name <- c("POSITIVE_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION","POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_MEMORY_T_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_ACTIVATED_T_CELL_PROLIFERATION","POSITIVE_REGULATION_OF_CYTOTOXIC_T_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_GAMMA_DELTA_T_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_CD4_POSITIVE_ ALPHA_BETA_T_CELL_ACTIVATION","POSITIVE_REGULATION_OF_MEMORY_T_CELL_ACTIVATION","POSITIVE_REGULATION_OF_ESTABLISHMENT_OF_T_CELL_POLARITY","POSITIVE_REGULATION_OF_CD4_POSITIVE_ ALPHA_BETA_T_CELL_PROLIFERATION","POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_COSTIMULATION","POSITIVE_REGULATION_OF_CD8_POSITIVE_ ALPHA_BETA_T_CELL_ACTIVATION","POSITIVE_REGULATION_OF_GAMMA_DELTA_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE","POSITIVE_REGULATION_OF_IMMATURE_T_CELL_PROLIFERATION","POSITIVE_REGULATION_OF_REGULATORY_T_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_CD8_POSITIVE,_ALPHA_BETA_T_CELL_EXTRAVASATION","POSITIVE_REGULATION_OF_CD8_POSITIVE_ ALPHA_BETA_CYTOTOXIC_T_CELL_EXTRAVASATION","POSITIVE_REGULATION_OF_ACTIVATED_CD8_POSITIVE_ALPHA_BETA_T_CELL_APOPTOTIC_PROCESS")

B.list.GO <- c("GO:0050871","GO:2000538","GO:0002663","GO:0002904","GO:0030890","GO:0045579","GO:0050861")
B.list.name <- c("POSITIVE_REGULATION_OF_B_CELL_ACTIVATION","POSITIVE_REGULATION_OF_B_CELL_CHEMOTAXIS","POSITIVE_REGULATION_OF_B_CELL_TOLERANCE_INDUCTION","POSITIVE_REGULATION_OF_B_CELL_APOPTOTIC_PROCESS","POSITIVE REGULATION OF B CELL PROLIFERATION","POSITIVE_REGULATION_OF_B_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_B_CELL_RECEPTOR_SIGNALING_PATHWAY")

NK.list.GO <- c("GO:2000503","GO:0043323","GO:0032816","GO:0032819","GO:0032825","GO:0045954","GO:0002717","GO:0002729","GO:0002857","GO:0002860")
NK.list.name <- c("POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_CHEMOTAXIS","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_DEGRANULATION","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_ACTIVATION","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_PROLIFERATION","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION","POSITIV_ REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL","POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY_DIRECTED_AGAINST_TUMOR_CELL_TARGET")

DC.list.GO <- c("GO:0030887","GO:2001200","GO:2000670","GO:2000549","GO:2000529","GO:2000510","GO:0002732","GO:0002735","GO:0002606")
DC.list.name <- c("POSITIVE_REGULATION_OF_MYELOID_DENDRITIC_CELL_ACTIVATION","POSITIVE_REGULATION_OF_DENDRITIC_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_DENDRITIC_CELL_APOPTOTIC_PROCESS","POSITIVE_REGULATION_OF_DENDRITIC_CELL_DENDRITE_ASSEMBLY","POSITIVE_REGULATION_OF_MYELOID_DENDRITIC_CELL_CHEMOTAXIS","POSITIVE_REGULATION_OF_DENDRITIC_CELL_CHEMOTAXIS","POSITIVE_REGULATION_OF_DENDRITIC_CELL_CYTOKINE_PRODUCTION","POSITIVE_REGULATION_OF_MYELOID_DENDRITIC_CELL_CYTOKINE_PRODUCTION","POSITIVE_REGULATION_OF_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION")

Mono.list.GO <- c("GO:0042117","GO:0045657","GO:2000439","GO:1900625","GO:0090026","GO:0071639")
Mono.list.name <- c("MONOCYTE_ACTIVATION","POSITIVE_REGULATION_OF_MONOCYTE_DIFFERENTIATION","POSITIVE_REGULATION_OF_MONOCYTE_EXTRAVASATION","POSITIVE_REGULATION_OF_MONOCYTE_AGGREGATION","POSITIVE_REGULATION_OF_MONOCYTE_CHEMOTAXIS","POSITIVE_REGULATION_OF_MONOCYTE_CHEMOTACTIC_PROTEIN-1_PRODUCTION")







T.db <- readRDS("path_to_T_cell_object")

B.db <- readRDS("path_to_B_cell_object")

NK.db <- readRDS("path_to_NK_cell_object")

DC.db <- readRDS("path_to_DC_cell_object")

Mono.db <- readRDS("path_to_Mono_cell_object")







cell_list <- c('pbmc', "CD4_T", "CD8_T", "other_T", "B", "NK", "CD14_Mono", "CD16_Mono", "DC", "other")





for (i in cell_list){
    x <- read.csv(paste0("path_to_gene_expression_file"))
    
    assign(paste0(i), x)
}







for (i in cell_list){
    x <- get(i)
    
    y <- x %>% select(gene, June.2021, August.2021, September.Pre.launch, 
                      September.Post.launch, November.2021, December.2021)
    y$Pre <- (y$June.2021+y$August.2021+y$September.Pre.launch)/3
    
    
    z1 <- y %>% select(June.2021, August.2021, September.Pre.launch, 
                      September.Post.launch, November.2021, December.2021)
    rownames(z1) <- y$gene
    
    z2 <- y %>% select(Pre, 
                      September.Post.launch, November.2021, December.2021)
    rownames(z2) <- y$gene
    
    
    
    colnames(z1) <- c("L-92", "L-44", "L-3", "R+1", "R+45", "R+82")
    
    colnames(z2) <- c("Pre-flight", "R+1", "R+45", "R+82")
    
    
    assign(paste0(i,".gsea.input.all"), z1)
    assign(paste0(i,".gsea.input.simple"), z2)
    
    }



CD4_T.1 <- gsva(as.matrix(CD4_T.gsea.input.all),gset.idx.list=T.db,method="ssgsea")
CD4_T.2 <- gsva(as.matrix(CD4_T.gsea.input.simple),gset.idx.list=T.db,method="ssgsea")

CD8_T.1 <- gsva(as.matrix(CD8_T.gsea.input.all),gset.idx.list=T.db,method="ssgsea")
CD8_T.2 <- gsva(as.matrix(CD8_T.gsea.input.simple),gset.idx.list=T.db,method="ssgsea")

other_T.1 <- gsva(as.matrix(other_T.gsea.input.all),gset.idx.list=T.db,method="ssgsea")
other_T.2 <- gsva(as.matrix(other_T.gsea.input.simple),gset.idx.list=T.db,method="ssgsea")

B.1 <- gsva(as.matrix(B.gsea.input.all),gset.idx.list=B.db,method="ssgsea")
B.2 <- gsva(as.matrix(B.gsea.input.simple),gset.idx.list=B.db,method="ssgsea")

NK.1 <- gsva(as.matrix(NK.gsea.input.all),gset.idx.list=NK.db,method="ssgsea")
NK.2 <- gsva(as.matrix(NK.gsea.input.simple),gset.idx.list=NK.db,method="ssgsea")

CD14_Mono.1 <- gsva(as.matrix(CD14_Mono.gsea.input.all),gset.idx.list=Mono.db,method="ssgsea")
CD14_Mono.2 <- gsva(as.matrix(CD14_Mono.gsea.input.simple),gset.idx.list=Mono.db,method="ssgsea")


CD16_Mono.1 <- gsva(as.matrix(CD16_Mono.gsea.input.all),gset.idx.list=Mono.db,method="ssgsea")
CD16_Mono.2 <- gsva(as.matrix(CD16_Mono.gsea.input.simple),gset.idx.list=Mono.db,method="ssgsea")

DC.1 <- gsva(as.matrix(DC.gsea.input.all),gset.idx.list=DC.db,method="ssgsea")
DC.2 <- gsva(as.matrix(DC.gsea.input.simple),gset.idx.list=DC.db,method="ssgsea")














paletteLength <- 200
low = colorRampPalette(c("black","navy","white"))(100)
high = colorRampPalette(c("white","firebrick1","firebrick2","firebrick3","firebrick4"))(100)
combined = c(low,high)
myColor2 <- combined





cell_list2 <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD4_Proliferating", "CD4_CTL", "CD8_Naive", "CD8_TCM",
               "CD8_TEM", "CD8_Proliferating", "dnT", "MAIT", "gdT", "Treg", "NK_CD56bright", "NK_Proliferating",
               "B_Naive", "B_intermediate", "B_memory", "cDC1", "cDC2", "pDC",
               "ASDC", "Plasmablast")



for (i in cell_list2){
    x <- read.csv(paste0("path_to_gene_exprssion_file"))
    
    assign(paste0(i), x)
}

cell_list3 <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD4_Proliferating", "CD8_Naive", "CD8_TCM",
               "CD8_TEM", "dnT", "MAIT", "gdT", "Treg", "NK_CD56bright", "NK_Proliferating",
               "B_Naive", "B_intermediate", "B_memory", "cDC2", "pDC",
               "ASDC", "Plasmablast")



for (i in cell_list3){
    x <- get(i)
    
    y <- x %>% select(gene, June.2021, August.2021,  
                      September.Post.launch, November.2021, December.2021)
    y$Pre <- (y$June.2021+y$August.2021)/2
    
    
    z1 <- y %>% select(June.2021, August.2021, 
                      September.Post.launch, November.2021, December.2021)
    rownames(z1) <- y$gene
    
    z2 <- y %>% select(Pre, 
                      September.Post.launch, November.2021, December.2021)
    rownames(z2) <- y$gene
    
    
    
    colnames(z1) <- c("L-92", "L-44", "R+1", "R+45", "R+82")
    
    colnames(z2) <- c("Pre-flight", "R+1", "R+45", "R+82")
    
    
    assign(paste0(i,".gsea.input.all"), z1)
    assign(paste0(i,".gsea.input.simple"), z2)
    
    }

cell_list3 <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD4_Proliferating", "CD8_Naive", "CD8_TCM",
               "CD8_TEM", "dnT", "MAIT", "gdT", "Treg", "NK_CD56bright", "NK_Proliferating",
               "B_Naive", "B_intermediate", "B_memory", "cDC2", "pDC",
               "ASDC", "Plasmablast")



T.cells <- c("CD4_Naive", "CD4_TCM", "CD4_TEM", "CD4_Proliferating", "CD8_Naive", "CD8_TCM","CD8_TEM", "dnT", "MAIT", "gdT", "Treg")

B.cells <- c("B_Naive", "B_intermediate", "B_memory", "Plasmablast")
NK.cells <- c("NK_CD56bright", "NK_Proliferating")
DC.cells <- c("cDC2", "pDC","ASDC")







CD4_T_GO <- c("POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION", "POSITIVE_REGULATION_OF_ACTIVATED_T_CELL_PROLIFERATION",
             "POSITIVE_REGULATION_OF_CD4_POSITIVE_ ALPHA_BETA_T_CELL_ACTIVATION", "POSITIVE_REGULATION_OF_ESTABLISHMENT_OF_T_CELL_POLARITY",
             "POSITIVE_REGULATION_OF_CD4_POSITIVE_ ALPHA_BETA_T_CELL_PROLIFERATION", "POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_COSTIMULATION")

CD8_T_GO <- c("POSITIVE_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION", "POSITIVE_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION",
             "POSITIVE_REGULATION_OF_ACTIVATED_T_CELL_PROLIFERATION", "POSITIVE_REGULATION_OF_ESTABLISHMENT_OF_T_CELL_POLARITY",
             "POSITIVE_REGULATION_OF_CD8_POSITIVE_ ALPHA_BETA_T_CELL_ACTIVATION", "POSITIVE_REGULATION_OF_CD8_POSITIVE,_ALPHA_BETA_T_CELL_EXTRAVASATION",
             "POSITIVE_REGULATION_OF_CD8_POSITIVE_ ALPHA_BETA_CYTOTOXIC_T_CELL_EXTRAVASATION", "POSITIVE_REGULATION_OF_ACTIVATED_CD8_POSITIVE_ALPHA_BETA_T_CELL_APOPTOTIC_PROCESS")

T_memory_GO <- c("POSITIVE_REGULATION_OF_MEMORY_T_CELL_DIFFERENTIATION", "POSITIVE_REGULATION_OF_MEMORY_T_CELL_ACTIVATION")

gdT_GO <- c("POSITIVE_REGULATION_OF_GAMMA_DELTA_T_CELL_DIFFERENTIATION", "POSITIVE_REGULATION_OF_GAMMA_DELTA_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")

CD4_Naive_GO <- c("POSITIVE_REGULATION_OF_IMMATURE_T_CELL_PROLIFERATION")

DC_GO <- c("POSITIVE_REGULATION_OF_MYELOID_DENDRITIC_CELL_ACTIVATION","POSITIVE_REGULATION_OF_DENDRITIC_CELL_DIFFERENTIATION","POSITIVE_REGULATION_OF_DENDRITIC_CELL_APOPTOTIC_PROCESS","POSITIVE_REGULATION_OF_DENDRITIC_CELL_DENDRITE_ASSEMBLY","POSITIVE_REGULATION_OF_MYELOID_DENDRITIC_CELL_CHEMOTAXIS","POSITIVE_REGULATION_OF_DENDRITIC_CELL_CHEMOTAXIS","POSITIVE_REGULATION_OF_DENDRITIC_CELL_CYTOKINE_PRODUCTION","POSITIVE_REGULATION_OF_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION")

cDC_GO <- c("POSITIVE_REGULATION_OF_MYELOID_DENDRITIC_CELL_CYTOKINE_PRODUCTION")



CD4_T.2.filter <- as.matrix(as.data.frame(CD4_T.2) %>% filter(row.names(as.data.frame(CD4_T.2)) %in% CD4_T_GO))

CD8_T.2.filter <- as.matrix(as.data.frame(CD8_T.2) %>% filter(row.names(as.data.frame(CD8_T.2)) %in% CD8_T_GO))

CD4_TCM.2.filter <- as.matrix(as.data.frame(CD4_TCM.2) %>% filter(row.names(as.data.frame(CD4_TCM.2)) %in% T_memory_GO))
CD4_TEM.2.filter <- as.matrix(as.data.frame(CD4_TEM.2) %>% filter(row.names(as.data.frame(CD4_TEM.2)) %in% T_memory_GO))

CD8_TCM.2.filter <- as.matrix(as.data.frame(CD8_TCM.2) %>% filter(row.names(as.data.frame(CD8_TCM.2)) %in% T_memory_GO))
CD8_TEM.2.filter <- as.matrix(as.data.frame(CD8_TEM.2) %>% filter(row.names(as.data.frame(CD8_TEM.2)) %in% T_memory_GO))

gdT.2.filter <- as.matrix(as.data.frame(gdT.2) %>% filter(row.names(as.data.frame(gdT.2)) %in% gdT_GO))

CD4_Naive.2.filter <- as.matrix(as.data.frame(CD4_Naive.2) %>% filter(row.names(as.data.frame(CD4_Naive.2)) %in% CD4_Naive_GO))



T.merge.2 <- rbind(CD4_T.2.filter, CD8_T.2.filter, CD4_TCM.2.filter, CD4_TEM.2.filter, CD8_TCM.2.filter,
                  CD8_TEM.2.filter, gdT.2.filter)

T <- c('CD4_T', 'CD8_T', 'CD4_TCM', 'CD4_TEM',
      'CD8_TCM', 'CD8_TEM', 'gdT')

for (i in T){
    x <- get(paste0(i, ".2.filter"))
    x <- as.data.frame(x)
    x$celltype <- i
    
    assign(paste0(i, ".2.filter"), x)
}



T.merge.2 <- rbind(CD4_T.2.filter, CD8_T.2.filter, CD4_TCM.2.filter, CD4_TEM.2.filter, CD8_TCM.2.filter,
                  CD8_TEM.2.filter, gdT.2.filter)







DC.2.filter <- as.matrix(as.data.frame(DC.2) %>% filter(row.names(as.data.frame(DC.2)) %in% DC_GO))

cDC2.2.filter <- as.matrix(as.data.frame(cDC2.2) %>% filter(row.names(as.data.frame(cDC2.2)) %in% cDC_GO))

DC.merge.2 <- rbind(DC.2.filter, cDC2.2.filter)

anno_time.dc <- data.frame(Timepoint = c("Pre-flight", "R+1",  
                                         "R+45", "R+82"))

row.names(anno_time.dc) <- colnames(DC.merge.2)



DC <- c('DC', 'cDC2')

for (i in DC){
    x <- get(paste0(i, ".2.filter"))
    x <- as.data.frame(x)
    x$celltype <- i
    
    assign(paste0(i, ".2.filter"), x)
}

DC.merge.2 <- rbind(DC.2.filter, cDC2.2.filter)



B.2.filter <- B.2

B <- c("B", 'B_intermediate', 'B_memory', 'B_Naive', 'Plasmablast')

for (i in B){
    x <- get(paste0(i, ".2"))
    x <- as.data.frame(x)
    x$celltype <- i
    
    assign(paste0(i, ".2"), x)
}

B.merge.2 <- rbind(B.2, B_intermediate.2, B_memory.2, B_Naive.2, Plasmablast.2)

mono <- c("CD14_Mono", 'CD16_Mono')

for (i in mono){
    x <- get(paste0(i, ".2"))
    x <- as.data.frame(x)
    x$celltype <- i
    
    assign(paste0(i, ".2.2"), x)
}



Mono.merge.2 <- rbind(CD14_Mono.2.2, CD16_Mono.2.2)





NK <- c("NK", 'NK_Proliferating', 'NK_CD56bright')

for (i in NK){
    x <- get(paste0(i, ".2"))
    x <- as.data.frame(x)
    x$celltype <- i
    
    assign(paste0(i, ".2"), x)
}

NK.merge.2 <- rbind(NK.2, NK_Proliferating.2, NK_CD56bright.2)



merge.2 <- rbind(T.merge.2, B.merge.2, DC.merge.2, Mono.merge.2, NK.merge.2)








