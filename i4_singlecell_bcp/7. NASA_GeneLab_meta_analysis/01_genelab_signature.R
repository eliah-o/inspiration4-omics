library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(tidyverse)

#load genelab signature fold changes
genelab_log2fc <- read_tsv("./data/genelab_signature_log2fc.tsv") %>% column_to_rownames(var = "gene") %>% as.matrix

#color palette heatmap
col_fun <- colorRamp2(seq(-1,1,0.2),rev(brewer.pal(11, "RdBu")))

#define annotations heatmap
ha = HeatmapAnnotation(Tissue = sapply(colnames(genelab_log2fc), function(x) strsplit(x, split = "\\.")[[1]][4]),
                       Sex = str_to_title(sapply(colnames(genelab_log2fc), function(x) strsplit(x, split = "\\.")[[1]][3])),
                       Age = sapply(colnames(genelab_log2fc), function(x) strsplit(x, split = "\\.")[[1]][2]) %>% as.numeric,
                       Duration = sapply(colnames(genelab_log2fc), function(x) strsplit(x, split = "\\.")[[1]][6]) %>% as.numeric,
                       col = list(Age = colorRamp2(c(0, 40), c("white", "#009988")),
                                  Duration = colorRamp2(c(0, 40), c("white", "#ee7733")),
                                  Tissue = setNames(brewer.pal(10, "Set3"),unique(sapply(colnames(genelab_log2fc), function(x) strsplit(x, split = "\\.")[[1]][4]))),
                                  Sex = setNames(brewer.pal(3, "Set1")[1:2],c("Male","Female"))))

#load genelab signature
genelab <- read_tsv("./data/genelab_meta_analysis.tsv")
genelab_up <- genelab %>% filter(log2FoldChange>0) %>% pull(gene)
genelab_dn <- genelab %>% filter(log2FoldChange<0) %>% pull(gene)
  
#plot heatmap
pdf(file = paste0("./output/01_genelab_signature.pdf"), width=8,height=5.87)
Heatmap(genelab_log2fc,
        column_labels = sapply(colnames(genelab_log2fc), function(x) strsplit(x, split = "\\.")[[1]][1]),
        bottom_annotation = ha,
        border = TRUE,
        use_raster = FALSE,
        row_split = ifelse(rownames(genelab_log2fc)%in%genelab_up,"Up-regulated","Down-regulated"),
        cluster_row_slices = FALSE,
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        heatmap_legend_param = list(direction = "vertical", title = "Log2\nfold\nchange", title_position = "topcenter"),
        col = col_fun)
dev.off()
