library(data.table)


library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(doMC)
library(ggplot2)



ref.sobj.path <- "path_to_reference"

print("Loading reference object")
ref.sobj <- readRDS(ref.sobj.path)

temp.sobj <- readRDS("path_to_Seurat_object")

assay.name <- "SCT"
label.names <- "celltype.l1"

print("Setting default assays")
Idents(ref.sobj) <- "celltype.l1"

DefaultAssay(temp.sobj) <- assay.name

DefaultAssay(ref.sobj) <- assay.name

temp.sobj$timepoint <- ifelse(temp.sobj$orig.ident == 'C001_1', 'June 2021', 
                              ifelse(temp.sobj$orig.ident == 'C001_2', 'August 2021', 
                                     ifelse(temp.sobj$orig.ident == 'C001_3', 'September Pre-launch', 
                                            ifelse(temp.sobj$orig.ident == 'C001_4', 'September Post-launch', 
                                                   ifelse(temp.sobj$orig.ident == 'C001_5', 'November 2021',
                                                          ifelse(temp.sobj$orig.ident == 'C001_6', 'December 2021',
                                                                ifelse(temp.sobj$orig.ident == 'C002_1', 'June 2021',
                                                                      ifelse(temp.sobj$orig.ident == 'C002_2', 'August 2021',
                                                                            ifelse(temp.sobj$orig.ident == 'C002_3', 'September Pre-launch',
                                                                                  ifelse(temp.sobj$orig.ident == 'C02_4', 'September Post-launch',
                                                                                        ifelse(temp.sobj$orig.ident == 'C002_5', 'November 2021',
                                                                                              ifelse(temp.sobj$orig.ident == 'C002_6', 'Decmeber 2021',
                                                                                                    ifelse(temp.sobj$orig.ident == 'C003_1', 'June 2021',
                                                                                                          ifelse(temp.sobj$orig.ident == 'C003_2', 'August 2021',
                                                                                                                ifelse(temp.sobj$orig.ident == 'C003_3', 'September Pre-launch',
                                                                                                                      ifelse(temp.sobj$orig.ident == 'C003_4', 'September Post-launch',
                                                                                                                            ifelse(temp.sobj$orig.ident == 'C003_5', 'November 2021',
                                                                                                                                  ifelse(temp.sobj$orig.ident == 'C003_6', 'December 2021',
                                                                                                                                        ifelse(temp.sobj$orig.ident == 'C004_1', 'June 2021',
                                                                                                                                              ifelse(temp.sobj$orig.ident == 'C004_2', 'August 2021',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_3', 'September Pre-launch',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_4', 'September Post-launch',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_5', 'November 2021', 'December 2021')))))))))))))))))))))))



temp.sobj$ID <- ifelse(temp.sobj$orig.ident == 'C001_1', 'C001', 
                              ifelse(temp.sobj$orig.ident == 'C001_2', 'C001', 
                                     ifelse(temp.sobj$orig.ident == 'C001_3', 'C001', 
                                            ifelse(temp.sobj$orig.ident == 'C001_4', 'C001', 
                                                   ifelse(temp.sobj$orig.ident == 'C001_5', 'C001',
                                                          ifelse(temp.sobj$orig.ident == 'C001_6', 'C001',
                                                                ifelse(temp.sobj$orig.ident == 'C002_1', 'C002',
                                                                      ifelse(temp.sobj$orig.ident == 'C002_2', 'C002',
                                                                            ifelse(temp.sobj$orig.ident == 'C002_3', 'C002',
                                                                                  ifelse(temp.sobj$orig.ident == 'C002_4', 'C002',
                                                                                        ifelse(temp.sobj$orig.ident == 'C002_5', 'C002',
                                                                                              ifelse(temp.sobj$orig.ident == 'C002_6', 'C002',
                                                                                                    ifelse(temp.sobj$orig.ident == 'C003_1', 'C003',
                                                                                                          ifelse(temp.sobj$orig.ident == 'C003_2', 'C003',
                                                                                                                ifelse(temp.sobj$orig.ident == 'C003_3', 'C003',
                                                                                                                      ifelse(temp.sobj$orig.ident == 'C003_4', 'C003',
                                                                                                                            ifelse(temp.sobj$orig.ident == 'C003_5', 'C003',
                                                                                                                                  ifelse(temp.sobj$orig.ident == 'C003_6', 'C003',
                                                                                                                                        ifelse(temp.sobj$orig.ident == 'C004_1', 'C004',
                                                                                                                                              ifelse(temp.sobj$orig.ident == 'C004_2', 'C004',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_3', 'C004',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_4', 'C004',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_5', 'C004', 'C004')))))))))))))))))))))))



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

Idents(temp.sobj) <- "predicted.id"


temp.sobj$timepoint <- ifelse(temp.sobj$orig.ident == 'C001_1', 'June 2021', 
                              ifelse(temp.sobj$orig.ident == 'C001_2', 'August 2021', 
                                     ifelse(temp.sobj$orig.ident == 'C001_3', 'September Pre-launch', 
                                            ifelse(temp.sobj$orig.ident == 'C001_4', 'September Post-launch', 
                                                   ifelse(temp.sobj$orig.ident == 'C001_5', 'November 2021',
                                                          ifelse(temp.sobj$orig.ident == 'C001_6', 'December 2021',
                                                                ifelse(temp.sobj$orig.ident == 'C002_1', 'June 2021',
                                                                      ifelse(temp.sobj$orig.ident == 'C002_2', 'August 2021',
                                                                            ifelse(temp.sobj$orig.ident == 'C002_3', 'September Pre-launch',
                                                                                  ifelse(temp.sobj$orig.ident == 'C002_4', 'September Post-launch',
                                                                                        ifelse(temp.sobj$orig.ident == 'C002_5', 'November 2021',
                                                                                              ifelse(temp.sobj$orig.ident == 'C002_6', 'December 2021',
                                                                                                    ifelse(temp.sobj$orig.ident == 'C003_1', 'June 2021',
                                                                                                          ifelse(temp.sobj$orig.ident == 'C003_2', 'August 2021',
                                                                                                                ifelse(temp.sobj$orig.ident == 'C003_3', 'September Pre-launch',
                                                                                                                      ifelse(temp.sobj$orig.ident == 'C003_4', 'September Post-launch',
                                                                                                                            ifelse(temp.sobj$orig.ident == 'C003_5', 'November 2021',
                                                                                                                                  ifelse(temp.sobj$orig.ident == 'C003_6', 'December 2021',
                                                                                                                                        ifelse(temp.sobj$orig.ident == 'C004_1', 'June 2021',
                                                                                                                                              ifelse(temp.sobj$orig.ident == 'C004_2', 'August 2021',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_3', 'September Pre-launch',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_4', 'September Post-launch',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_5', 'November 2021', 'December 2021')))))))))))))))))))))))



saveRDS(temp.sobj, file = "celltype.l1.filtered.multiome.sobj.rds")

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

Idents(temp.sobj) <- "predicted.id"


temp.sobj$timepoint <- ifelse(temp.sobj$orig.ident == 'C001_1', 'June 2021', 
                              ifelse(temp.sobj$orig.ident == 'C001_2', 'August 2021', 
                                     ifelse(temp.sobj$orig.ident == 'C001_3', 'September Pre-launch', 
                                            ifelse(temp.sobj$orig.ident == 'C001_4', 'September Post-launch', 
                                                   ifelse(temp.sobj$orig.ident == 'C001_5', 'November 2021',
                                                          ifelse(temp.sobj$orig.ident == 'C001_6', 'December 2021',
                                                                ifelse(temp.sobj$orig.ident == 'C002_1', 'June 2021',
                                                                      ifelse(temp.sobj$orig.ident == 'C002_2', 'August 2021',
                                                                            ifelse(temp.sobj$orig.ident == 'C002_3', 'September Pre-launch',
                                                                                  ifelse(temp.sobj$orig.ident == 'C002_4', 'September Post-launch',
                                                                                        ifelse(temp.sobj$orig.ident == 'C002_5', 'November 2021',
                                                                                              ifelse(temp.sobj$orig.ident == 'C002_6', 'December 2021',
                                                                                                    ifelse(temp.sobj$orig.ident == 'C003_1', 'June 2021',
                                                                                                          ifelse(temp.sobj$orig.ident == 'C003_2', 'August 2021',
                                                                                                                ifelse(temp.sobj$orig.ident == 'C003_3', 'September Pre-launch',
                                                                                                                      ifelse(temp.sobj$orig.ident == 'C003_4', 'September Post-launch',
                                                                                                                            ifelse(temp.sobj$orig.ident == 'C003_5', 'November 2021',
                                                                                                                                  ifelse(temp.sobj$orig.ident == 'C003_6', 'December 2021',
                                                                                                                                        ifelse(temp.sobj$orig.ident == 'C004_1', 'June 2021',
                                                                                                                                              ifelse(temp.sobj$orig.ident == 'C004_2', 'August 2021',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_3', 'September Pre-launch',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_4', 'September Post-launch',
                                                                                                                                                    ifelse(temp.sobj$orig.ident == 'C004_5', 'November 2021', 'December 2021')))))))))))))))))))))))



saveRDS(temp.sobj, file = "celltype.l2.filtered.multiome.sobj.rds")



temp.l1.sobj <- readRDS("celltype.l1.filtered.multiome.sobj.rds")

temp.l2.sobj <- readRDS("celltype.l2.filtered.multiome.sobj.rds")

temp.sobj <- temp.l1.sobj

colnames(temp.sobj@meta.data)[19] <- 'celltype.l1'

temp.sobj$celltype.l2 <- temp.l2.sobj$predicted.id

temp.sobj$celltype <- ifelse(temp.sobj$celltype.l1 == 'B', 'B',
                            ifelse(temp.sobj$celltype.l1 == 'CD4 T', 'CD4 T',
                                  ifelse(temp.sobj$celltype.l1 == 'CD8 T', 'CD8 T',
                                        ifelse(temp.sobj$celltype.l1 == 'DC', 'DC',
                                                ifelse(temp.sobj$celltype.l1 == 'NK', 'NK',
                                                        ifelse(temp.sobj$celltype.l1 == 'other T', 'other T',
                                                                ifelse(temp.sobj$celltype.l2 == 'CD14 Mono', 'CD14 Mono',
                                                                    ifelse(temp.sobj$celltype.l2 == 'CD16 Mono', 'CD16 Mono', 'other'))))))))



saveRDS(temp.sobj, file = 'celltype.filtered.multiome.sobj.rds')





sessionInfo()

rm(list=ls())










