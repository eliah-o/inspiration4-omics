library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
register(SerialParam())
options(Seurat.object.assay.version = 'v5')

library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
ref.blacklist <- blacklist_hg38_unified
ref.genome <- BSgenome.Hsapiens.UCSC.hg38



if (!dir.exists("work")){dir.create("work")}
work <- "work"



samples <- c("Day1_Rep1", 'Day1_Rep2', 'Day1_Rep3',
            'Day2_Rep1', 'Day2_Rep2', 'Day2_Rep3')


for (i in samples){
    x <- readRDS(paste0('/athena/masonlab/scratch/projects/spacex/data/single_cell/raw/10X_multiome_control_PBMC/Donor1_', i, '/Seurat/Donor1_', i, '.rds'))
    
    assign(paste0(i), x)
}







print("Merging multiome objects:")
msobj.list <- paste0(samples)

print(msobj.list)



for (i in msobj.list){
  temp.sobj <- get(i)
  temp.sobj$orig.ident <- i
  Idents(temp.sobj) <- temp.sobj$orig.ident
  assign(i, temp.sobj)
}




if (length(samples) > 2){
  combined.msobj <- merge(x = get(paste0(samples[1])), y = mget(paste0(samples[2:length(samples)])), add.cell.ids = samples)
}else{
  combined.msobj <- merge(x = get(paste0(samples[1])), y = get(paste0(samples[2])), add.cell.ids = samples)
}



print("Saving Combined multiome object to")
print(file.path(getwd(), "work", "combined.multiome.merge.sobj.rds"))

saveRDS(combined.msobj, file = "work/combined.multiome.merge.sobj.rds")
print("Combined multiome signac object saved!")




DefaultAssay(combined.msobj) <- "RNA"
combined.msobj <- SCTransform(combined.msobj)
combined.msobj <- RunPCA(combined.msobj, npcs = 50, verbose = T)



combined.msobj

print("Processing combined ATAC data")
DefaultAssay(combined.msobj) <- "ATAC"
combined.msobj <- FindTopFeatures(combined.msobj, min.cutoff = "q50")
combined.msobj <- RunTFIDF(combined.msobj)
combined.msobj <- RunSVD(combined.msobj)



# build a joint neighbor graph using both assays
combined.msobj <- FindMultiModalNeighbors(
  object = combined.msobj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)





# build a joint UMAP visualization
combined.msobj <- RunUMAP(
  object = combined.msobj,
  nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
  verbose = TRUE
)





combined.msobj <- FindClusters(combined.msobj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)



combined.msobj

print("Saving Combined multiome object to")
print(file.path(getwd(), "work", "combined.multiome.before.sobj.rds"))

saveRDS(combined.msobj, file = "work/combined.multiome.before.sobj.rds")
print("Combined multiome signac object saved!")

library(ggplot2)

options(repr.plot.width=18, repr.plot.height=14)
DimPlot(combined.msobj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', 
        label.size = 16, reduction = "umap", raster = FALSE)+
  theme(axis.text = element_text(size = 48),
    axis.title = element_text(size=48),
    legend.text=element_text(size=48),
        title = element_text(size=0)
)

options(repr.plot.width=18, repr.plot.height=14)
DimPlot(combined.msobj, label = TRUE, repel = TRUE, group.by = 'orig.ident', 
        label.size = 16, reduction = "umap", raster = FALSE)+
  theme(axis.text = element_text(size = 48),
    axis.title = element_text(size=48),
    legend.text=element_text(size=48),
        title = element_text(size=0)
)



sessionInfo()

rm(list=ls())


