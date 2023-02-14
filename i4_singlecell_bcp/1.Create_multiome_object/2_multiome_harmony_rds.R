#! /bin/R
# by Christopher Chin
# 5/2021
# Use harmony on RNA and combined ATAC peaks for normalization
# Use RDS file as input

# Get runsheet path from arguments
args <- commandArgs(trailingOnly = T)
# Object path
msobj.path <- args[1]
# name of object
msobj.name <- args[2]
# Name of metadata to correct by
batch.name <- args[3]
# number of cores to use
if(is.na(args[4])){
  cores <- 12
}else{cores <- as.numeric(args[4])}

# Load libraries
library(future)
plan("multiprocess", workers = cores)
print(paste("Using", cores, "cores"))

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
# Running into memory issues
library(doMC)
registerDoMC(cores)
library(harmony)

# Create work directory
print(paste("Creating work directory at", file.path(getwd(),"work")))
if(!dir.exists("work")){dir.create("work")}

# Load data

temp.msobj.pca.path <- file.path("work","temp.msobj.pca.Rdata")
temp.msobj.lsi.path <- file.path("work","temp.msobj.lsi.Rdata")

if (file.exists(temp.msobj.lsi.path)){
  print("Found data from previous run")
  print(paste("Loading harmony corrected pca and lsi object from", temp.msobj.lsi.path))
  load(temp.msobj.lsi.path)
}else if (file.exists(temp.msobj.pca.path)){
  print("Found data from previous run")
  print(paste("Loading harmony corrected pca object from", temp.msobj.pca.path))
  load(temp.msobj.pca.path)
  temp.msobj <- RunHarmony(
    object = temp.msobj,
    group.by.vars = batch.name,
    reduction = 'lsi',
    assay.use = 'peaks',
    project.dim = FALSE
  )
  temp.msobj[["harmony.lsi"]] <- temp.msobj[["harmony"]]
  print(paste("Saving harmony corrected lsi to", temp.msobj.lsi.path))
  save(temp.msobj, file = temp.msobj.lsi.path)
}else {
  print("Loading data")
  temp.msobj <- readRDS(msobj.path)
  temp.msobj <- RunHarmony(
    object = temp.msobj,
    group.by.vars = batch.name,
    reduction = 'pca',
    assay.use = 'RNA',
    project.dim = FALSE
  )
  temp.msobj[["harmony_pca"]] <- temp.msobj[["harmony"]]
  print(paste("Saving harmony corrected pca to", temp.msobj.pca.path))
  save(temp.msobj, file = temp.msobj.pca.path)
}

if (!("harmony_pca" %in% names(temp.msobj@reductions))){
  temp.msobj <- RunHarmony(
    object = temp.msobj,
    group.by.vars = batch.name,
    reduction = 'pca',
    assay.use = 'RNA',
    project.dim = FALSE
  )
  temp.msobj[["harmony_pca"]] <- temp.msobj[["harmony"]]
  save(temp.msobj, file = temp.msobj.pca.path)
}

if (!("harmony.lsi" %in% names(temp.msobj@reductions))){
  temp.msobj <- RunHarmony(
    object = temp.msobj,
    group.by.vars = batch.name,
    reduction = 'lsi',
    assay.use = 'peaks',
    project.dim = FALSE
  )
  temp.msobj[["harmony.lsi"]] <- temp.msobj[["harmony"]]
  save(temp.msobj, file = temp.msobj.lsi.path)
}

print ("Building joint neighbors based on both RNA and ATAC")
temp.msobj <- FindMultiModalNeighbors(
  object = temp.msobj,
  reduction.list = list("harmony_pca", "harmony.lsi"), 
  dims.list = list(1:50, 2:50),
  modality.weight.name = c("RNA.weight","peak.weight"),
  verbose = TRUE
)

print("Running UMAP on RNA and ATAC data")
temp.msobj <- RunUMAP(
  object = temp.msobj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

print(paste("Saving final object to",file.path(getwd(), paste0(msobj.name, ".harmony.integrated.Rdata"))))

saveRDS(temp.msobj, file = file.path(getwd(), paste0(msobj.name, ".harmony.integrated.rds")))
