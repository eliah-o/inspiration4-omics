library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(BiocParallel)
register(SerialParam())
options(Seurat.object.assay.version = 'v5')

library(data.table)

library(BSgenome.Hsapiens.UCSC.hg38)


library(patchwork)
library(TFBSTools)

library(JASPAR2020)


set.seed(1234)

temp.sobj <- readRDS("path_to_Seurat_object")

temp.sobj$predicted.id <- ifelse(temp.sobj$celltype.l1 == 'CD4 T', 'CD4_T', 
                              ifelse(temp.sobj$celltype.l1 == 'CD8 T', 'CD8_T', 
                                     ifelse(temp.sobj$celltype.l1 == 'other T', 'other_T', 
                                            ifelse(temp.sobj$celltype.l1 == 'B', 'B', 
                                                   ifelse(temp.sobj$celltype.l1 == 'NK', 'NK',
                                                          ifelse(temp.sobj$celltype.l1 == 'Mono', 'Mono',
                                                                         ifelse(temp.sobj$celltype.l1 == 'DC', 'DC','other')))))))



DefaultAssay(temp.sobj) <- "ATAC"

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
temp.sobj <- AddMotifs(
  object = temp.sobj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)



temp.sobj <- RunChromVAR(
  object = temp.sobj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)



sessionInfo()

rm(list=ls())


