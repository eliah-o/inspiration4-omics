# All codes are written by Matias Fuentealba Valenzuela (mfuentealba@buckinstitute.org).

library(tidyverse)
library(GeneOverlap)
library(ggpubr)

calc_fisher <- function(genelist1,genelist2,universe){
  genelist1 <- intersect(genelist1,universe)
  genelist2 <- intersect(genelist2,universe)
  fisher <- testGeneOverlap(newGeneOverlap(genelist1,genelist2, genome.size = length(universe)))
  return(tibble(overlap = length(intersect(genelist1,genelist2)),or = fisher@odds.ratio,p = fisher@pval))
}

#Define genelab up- and down-regulated genes
detected_i4 <- read_rds("./data/genes_detected_i4.rds")
genelab <- read_tsv("./data/genelab_meta_analysis.tsv")
ortologues <- read_tsv("./data/human2mouse_symbols.txt") %>% set_names("human","mouse") %>% na.omit
background <- intersect(detected_i4,unique(ortologues$human[ortologues$mouse%in%genelab$gene]))
genelab_dn <- intersect(unique(ortologues$human[ortologues$mouse%in%c(genelab$gene[genelab$log2FoldChange<0&genelab$padj<0.05])]),background)
genelab_up <- intersect(unique(ortologues$human[ortologues$mouse%in%c(genelab$gene[genelab$log2FoldChange>0&genelab$padj<0.05])]),background)

#Fisher's test - level 1 cell annotation
i4_l1 <- do.call("rbind",lapply(list.files("./data/DEGs_inspiration4/celltype.l1/", full.names = TRUE), function(x) read_csv(x) %>% mutate(cell = gsub(".FP1.csv|DEGs_","",basename(x)))))
i4_l1 <- rbind(i4_l1 %>% filter(p_val_adj<0.05&avg_log2FC>0) %>% group_by(cell) %>% summarise(i4 = list(intersect(X,background)), dir = "up"),
               i4_l1 %>% filter(p_val_adj<0.05&avg_log2FC<0) %>% group_by(cell) %>% summarise(i4 = list(intersect(X,background)), dir = "dn"))
i4_l1 <- i4_l1 %>% mutate(genelab = if_else(dir == "dn", list(genelab_dn), list(genelab_up))) 
i4_l1$background <- list(background)
i4_l1 <- i4_l1 %>% mutate(fisher = pmap(list(i4,genelab,background),calc_fisher)) %>% unnest(fisher)
i4_l1$fdr <- p.adjust(i4_l1$p)
i4_l1$i4 <- lapply(i4_l1$i4, length) %>% unlist
i4_l1$genelab <- lapply(i4_l1$genelab, length) %>% unlist
i4_l1$background <- lapply(i4_l1$background, length) %>% unlist
i4_l1$dir <- ifelse(i4_l1$dir=="up","Up-regulated","Down-regulated")
i4_l1$cell <- gsub("_"," ",i4_l1$cell)
i4_l1$cell <- gsub("other","Other",i4_l1$cell)
i4_l1$cell <- gsub("pbmc","PBMC",i4_l1$cell)
i4_l1$cell <- factor(i4_l1$cell, levels = sort(unique(i4_l1$cell)) %>% rev)
mean(i4_l1$fdr<0.05)

p1 <- ggplot(data=i4_l1, aes(x=cell, y=-log10(p), fill=dir)) +
  geom_bar(aes(alpha=ifelse(fdr < 0.05, 1, 0.4)), stat="identity", position=position_dodge())+
  scale_alpha_continuous(limits = c(0, 1), range = c(0, 1))+
  coord_flip()+
  geom_vline(xintercept = seq(0.5, length(unique(i4_l1$cell)), by=1), linetype="dashed", color="black", lwd = 0.2) + # Add lines
  scale_fill_manual(values = c("#4393c3","#d6604d"))+
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")+
  labs(x = "Cell type", y = "-log10(p)", fill = "")

#Fisher's test - level 2 cell annotation
i4_l2 <- do.call("rbind",lapply(list.files("./data/DEGs_inspiration4/celltype.l2/", full.names = TRUE), function(x) read_csv(x) %>% mutate(cell = gsub("_September_Post_JAS.csv|DEGs_","",basename(x)))))
i4_l2 <- rbind(i4_l2 %>% filter(p_val_adj<0.05&avg_log2FC>0) %>% group_by(cell) %>% summarise(i4 = list(intersect(...1,background)), dir = "up"),
               i4_l2 %>% filter(p_val_adj<0.05&avg_log2FC<0) %>% group_by(cell) %>% summarise(i4 = list(intersect(...1,background)), dir = "dn"))
i4_l2 <- i4_l2 %>% mutate(genelab = if_else(dir == "dn", list(genelab_dn), list(genelab_up))) 
i4_l2$background <- list(background)
i4_l2 <- i4_l2 %>% mutate(fisher = pmap(list(i4,genelab,background),calc_fisher)) %>% unnest(fisher)
i4_l2$fdr <- p.adjust(i4_l2$p)
i4_l2$i4 <- lapply(i4_l2$i4, length) %>% unlist
i4_l2$genelab <- lapply(i4_l2$genelab, length) %>% unlist
i4_l2$background <- lapply(i4_l2$background, length) %>% unlist
i4_l2$dir <- ifelse(i4_l2$dir=="up","Up-regulated","Down-regulated")
i4_l2$cell <- gsub("_"," ",i4_l2$cell)
i4_l2$cell <- factor(i4_l2$cell, levels = sort(unique(i4_l2$cell)) %>% rev)
mean(i4_l2$fdr<0.05)

p2 <- ggplot(data=i4_l2, aes(x=cell, y=-log10(p), fill=dir)) +
  geom_bar(aes(alpha=ifelse(fdr < 0.05, 1, 0.4)), stat="identity", position=position_dodge())+
  scale_alpha_continuous(limits = c(0, 1), range = c(0, 1))+
  coord_flip()+
  geom_vline(xintercept = seq(0.5, length(unique(i4_l2$cell)), by=1), linetype="dashed", color="black", lwd = 0.2) + # Add lines
  scale_fill_manual(values = c("#4393c3","#d6604d"))+
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")+
  labs(x = "Cell type", y = "-log10(p)", fill = "")

#Print figure
pdf(file = "./output/02_fisher.pdf",width=9.00,height=4.00)
ggarrange(p1,p2, nrow = 1, widths = c(1,1))
dev.off()
