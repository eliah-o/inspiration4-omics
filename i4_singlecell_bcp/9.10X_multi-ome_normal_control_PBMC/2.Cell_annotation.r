library(data.table)
library(Signac)
library(Seurat)
library(BiocParallel)
register(SerialParam())
options(Seurat.object.assay.version = 'v5')

ref.sobj.path <- "path_to_reference"

ref.sobj <- readRDS(ref.sobj.path)

temp.sobj <- readRDS("path_to_seurat_object")



assay.name <- "SCT"
label.names <- "celltype.l1"

print("Setting default assays")
Idents(ref.sobj) <- "celltype.l1"

DefaultAssay(temp.sobj) <- assay.name

DefaultAssay(ref.sobj) <- assay.name

transfer_anchors <- FindTransferAnchors(
  reference = ref.sobj,
  query = temp.sobj,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)



predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = ref.sobj$celltype.l1,
  weight.reduction = temp.sobj[['pca']],
  dims = 1:50
)



temp.sobj <- AddMetaData(
  object = temp.sobj,
  metadata = predictions
)

temp.sobj$celltype.l1 <- temp.sobj$predicted.id



Idents(temp.sobj) <- "predicted.id"




assay.name <- "SCT"
label.names <- "celltype.l2"

print("Setting default assays")
Idents(ref.sobj) <- "celltype.l2"

DefaultAssay(temp.sobj) <- assay.name

DefaultAssay(ref.sobj) <- assay.name

transfer_anchors <- FindTransferAnchors(
  reference = ref.sobj,
  query = temp.sobj,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)



predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = ref.sobj$celltype.l2,
  weight.reduction = temp.sobj[['pca']],
  dims = 1:50
)



temp.sobj <- AddMetaData(
  object = temp.sobj,
  metadata = predictions
) 







sessionInfo()

rm(list=ls())


