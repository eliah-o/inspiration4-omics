library(tidyverse)
library(fgsea)
library(ggpubr)

#Define genelab rank ordered list
ortologues <- read_tsv("./data/human2mouse_symbols.txt") %>% set_names("human_gene","gene") %>% na.omit
genelab <- read_tsv("./data/genelab_meta_analysis.tsv") %>% 
  dplyr::select(gene, pvalue, log2FoldChange) %>% 
  left_join(ortologues) %>% 
  na.omit() %>% 
  mutate(score = -log10(pvalue)*sign(log2FoldChange)) %>% 
  group_by(human_gene) %>%
  summarise(score = mean(score)) %>%
  set_names("gene","score") %>% 
  arrange(score)
genelab <- setNames(genelab$score,genelab$gene)

#GSEA - level 1 cell annotation
i4_l1 <- do.call("rbind",lapply(list.files("./data/DEGs_inspiration4/celltype.l1/", full.names = TRUE), function(x) read_csv(x) %>% mutate(cell = gsub(".FP1.csv|DEGs_","",basename(x)))))
i4_l1 <- rbind(i4_l1 %>% filter(p_val_adj<0.05&avg_log2FC>0) %>% group_by(cell) %>% summarise(genes = list(X), dir = "up"),
               i4_l1 %>% filter(p_val_adj<0.05&avg_log2FC<0) %>% group_by(cell) %>% summarise(genes = list(X), dir = "dn"))
i4_l1$label <- interaction(i4_l1$dir,i4_l1$cell) %>% as.character()
i4_l1 <- setNames(i4_l1$genes, i4_l1$label)

gsea_l1 <- fgsea(i4_l1, genelab, maxSize=10000) %>% 
  dplyr::select(pathway, pval, padj, ES, NES) %>% 
  separate(col = "pathway", into = c("dir","cell"), sep = "\\.")
gsea_l1$dir <- ifelse(gsea_l1$dir=="up","Up-regulated","Down-regulated")
gsea_l1$cell <- gsub("_"," ",gsea_l1$cell)
gsea_l1$cell <- gsub("other","Other",gsea_l1$cell)
gsea_l1$cell <- gsub("pbmc","PBMC",gsea_l1$cell)
gsea_l1$cell <- factor(gsea_l1$cell, levels = sort(unique(gsea_l1$cell)) %>% rev)

p1 <- ggplot(data=gsea_l1, aes(x=cell, y=NES, fill=dir)) +
  geom_bar(aes(alpha=ifelse(padj < 0.05, 1, 0.4)), stat="identity", position=position_dodge())+
  scale_alpha_continuous(limits = c(0, 1), range = c(0, 1))+
  coord_flip()+
  ylim(-3,3)+
  geom_vline(xintercept = seq(0.5, length(unique(gsea_l1$cell)), by=1), linetype="dashed", color="black", lwd = 0.2) + # Add lines
  scale_fill_manual(values = c("#4393c3","#d6604d"))+
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")+
  labs(x = "Cell type", y = "Normalized enrichment score", fill = "")

#GSEA - level 2 cell annotation
i4_l2 <- do.call("rbind",lapply(list.files("./data/DEGs_inspiration4/celltype.l2/", full.names = TRUE), function(x) read_csv(x) %>% mutate(cell = gsub("_September_Post_JAS.csv|DEGs_","",basename(x)))))
i4_l2 <- rbind(i4_l2 %>% filter(p_val_adj<0.05&avg_log2FC>0) %>% group_by(cell) %>% summarise(genes = list(...1), dir = "up"),
               i4_l2 %>% filter(p_val_adj<0.05&avg_log2FC<0) %>% group_by(cell) %>% summarise(genes = list(...1), dir = "dn"))
i4_l2$label <- interaction(i4_l2$dir,i4_l2$cell) %>% as.character()
i4_l2 <- setNames(i4_l2$genes, i4_l2$label)

gsea_l2 <- fgsea(i4_l2, genelab, maxSize=10000) %>% 
  dplyr::select(pathway, pval, padj, ES, NES) %>% 
  separate(col = "pathway", into = c("dir","cell"), sep = "\\.")
gsea_l2$dir <- ifelse(gsea_l2$dir=="up","Up-regulated","Down-regulated")
gsea_l2$cell <- gsub("_"," ",gsea_l2$cell)
gsea_l2$cell <- factor(gsea_l2$cell, levels = sort(unique(gsea_l2$cell)) %>% rev)

p2 <- ggplot(data=gsea_l2, aes(x=cell, y=NES, fill=dir)) +
  geom_bar(aes(alpha=ifelse(padj < 0.05, 1, 0.4)), stat="identity", position=position_dodge())+
  scale_alpha_continuous(limits = c(0, 1), range = c(0, 1))+
  coord_flip()+
  ylim(-3,3)+
  geom_vline(xintercept = seq(0.5, length(unique(gsea_l2$cell)), by=1), linetype="dashed", color="black", lwd = 0.2) + # Add lines
  scale_fill_manual(values = c("#4393c3","#d6604d"))+
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")+
  labs(x = "Cell type", y = "Normalized enrichment score", fill = "")

#Print figure
pdf(file = "./output/03_gsea.pdf",width=9.00,height=4.00)
ggarrange(p1,p2, nrow = 1, widths = c(1,1))
dev.off()
