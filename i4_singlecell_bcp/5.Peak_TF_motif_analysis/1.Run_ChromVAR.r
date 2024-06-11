library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)


library(patchwork)
library(TFBSTools)

library(JASPAR2020)


set.seed(1234)


library(BiocParallel)


register(SerialParam())





temp.sobj <- readRDS("path_to_Seurat_object")

DefaultAssay(temp.sobj) <- "peaks"

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




