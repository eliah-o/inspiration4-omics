library(data.table)

runsheet.path <- "to_runsheet.path"



species <- "human"
cores <- as.numeric("1")

print("Loading runsheet")
runsheet <- fread(runsheet.path, header = T)
if (all(colnames(runsheet) == c("sample","path","timepoint","ID")) == F){
  print("Runsheet colnames do not match")
  stop()
}else{print("Runsheet accepted!")}

samples <- runsheet$sample
timepoints <- runsheet$timepoint
IDs <- runsheet$ID

library(future)
plan("multiprocess", workers = cores)
print(paste("Using", cores, "cores"))

library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(doMC)
registerDoMC(cores)

set.seed(42)

options(future.globals.maxSize= 53687091200000000)

if (!dir.exists("work")){dir.create("work")}
work <- "work"

if (species %in% c("Human", "human")){
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  ref.blacklist <- blacklist_hg38_unified
  ref.genome <- BSgenome.Hsapiens.UCSC.hg38
}else if (species %in% c("Mouse","mouse")){
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "mm10"
  ref.blacklist <- blacklist_mm10
  ref.genome <- BSgenome.Mmusculus.UCSC.mm10
}else {
  print("species not found")
  stop()
}



print("Creating seurat objects with RNA data")

for (i in samples){
  temp.obj.path <- file.path(work, paste0(i,".sobj.Rdata"))
  temp.outpath <- runsheet$path[runsheet$sample == i]
  if(!file.exists(temp.obj.path)){
    counts <- Read10X_h5(file.path(temp.outpath, "filtered_feature_bc_matrix.h5"))
    print(paste("Creating",i,"seurat object"))
    temp.sobj <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA", project = i
    )
    assign(paste0(i,".sobj"), temp.sobj)
    save(list = paste0(i,".sobj"), file = temp.obj.path)
  }else{
    print(paste("Loading seurat object from", temp.obj.path))
    load(temp.obj.path)
    temp.sobj <- get(paste0(i,".sobj"))
  }
  temp.msobj.path <- file.path(work, paste0(i,".indv.peaks.msobj.Rdata"))
  if(!file.exists(temp.msobj.path)){
    print("Creating ATAC assay and adding to object")
    temp.sobj[["ATAC"]] <- CreateChromatinAssay(
      counts = counts$Peaks,
      sep = c(":", "-"),
      fragments = file.path(temp.outpath, "atac_fragments.tsv.gz"),
      annotation = annotation
    )
    assign(paste0(i,".msobj"), temp.sobj)
    save(list = paste0(i,".msobj"), file = temp.msobj.path)
  }else{
    print(paste("Loading multiome seurat object (individual peaks) from", temp.msobj.path))
    load(temp.msobj.path)
  }
}



peak.file.path <- file.path(work, "combined.peaks.Rdata")
if(!file.exists(peak.file.path)){
  for (i in samples){
    temp.sobj <- get(paste0(i,".msobj"))
    print("Calling peaks with MACS2")
    DefaultAssay(temp.sobj) <- "ATAC"
    temp.peaks <- CallPeaks(temp.sobj, macs2.path = "/home/jak4013/.conda/envs/cellranger/bin/macs2")
    temp.peaks <- keepStandardChromosomes(temp.peaks, pruning.mode = "coarse")
    temp.peaks <- subsetByOverlaps(x = temp.peaks, ranges = ref.blacklist, invert = TRUE)
    assign(paste0(i,".peaks"), temp.peaks)
  }
  print("Combining peaks")
  all.peaks <- paste0(samples, ".peaks")
  combined.peaks <- get(all.peaks[1])
  for (i in all.peaks[2:length(all.peaks)]){
    combined.peaks <- c(combined.peaks,get(i))
  }
  combined.peaks <- reduce(combined.peaks)
  save(combined.peaks, file = peak.file.path)
}else{
  print(paste("Loading peaks from", peak.file.path))
  load(peak.file.path)
}





print("Flitering peaks")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]



for (i in samples){
  temp.frag.path <- file.path(work, paste0(i,".frag.Rdata"))
  if(!file.exists(temp.frag.path)){
    temp.sobj <- get(paste0(i,".sobj"))
    temp.outpath <- runsheet$path[runsheet$sample == i]
    print(paste("Creating",i, "fragment object"))
    temp.frag <- CreateFragmentObject(path = file.path(temp.outpath, "atac_fragments.tsv.gz"), cells = colnames(temp.sobj))
    assign(paste0(i,".frag"), temp.frag)
    save(list = paste0(i,".frag"), file = file.path(work, paste0(i,".frag.Rdata")))
  }else{
    print(paste("Loading fragment object from", temp.frag.path))
    load(file.path(work, paste0(i,".frag.Rdata")))
  }
}

rm(list = ls(pattern = "temp"))

for (i in samples){
  temp.msobj <- get(paste0(i,".msobj"))
  temp.sobj <- get(paste0(i,".sobj"))
  temp.msobj.path <- file.path(work, paste0(i,".combined.peaks.msobj.Rdata"))
  if(!file.exists(temp.msobj.path)){
    print("Creating ATAC assay and adding to object")
    DefaultAssay(temp.msobj) <- "ATAC"
    temp.frag <- get(paste0(i,".frag"))
    temp.counts <- FeatureMatrix(
      fragments = Fragments(temp.msobj),
      features = combined.peaks,
      cells = colnames(temp.msobj)
    )
    temp.sobj[["peaks"]] <- CreateChromatinAssay(
      counts = temp.counts,
      fragments = temp.frag,
      annotation = annotation,
      genome = "hg38"
    )
    print("Adding Nucleosome signal and TSS enrichment to object")
    DefaultAssay(temp.sobj) <- "peaks"
    temp.sobj <- NucleosomeSignal(temp.sobj)
    temp.sobj <- TSSEnrichment(temp.sobj)
    
    temp.sobj <- RenameCells(temp.sobj, new.names = paste(i, colnames(temp.sobj), sep = "_"))
    
    assign(paste0(i,".msobj"), temp.sobj)
    save(list = paste0(i,".msobj"), file = temp.msobj.path)
  }else{
    print(paste("Loading multiome seurat object (combined peaks) from", temp.msobj.path))
    load(temp.msobj.path)
  }
}


print("Merging multiome objects:")
msobj.list <- paste0(samples, ".msobj")
print(msobj.list)

for (i in msobj.list){
  temp.sobj <- get(i)
  Idents(temp.sobj) <- temp.sobj$orig.ident
  assign(i, temp.sobj)
}

if (length(samples) > 2){
  combined.msobj <- merge(x = get(paste0(samples[1],".msobj")), y = mget(paste0(samples[2:length(samples)], ".msobj")), add.cell.ids = samples)
}else{
  combined.msobj <- merge(x = get(paste0(samples[1],".msobj")), y = get(paste0(samples[2], ".msobj")), add.cell.ids = samples)
}


print("Saving Combined multiome object to")
print(file.path(getwd(), "work", "combined.multiome.merge.sobj.rds"))

saveRDS(combined.msobj, file = "work/combined.multiome.merge.sobj.rds")
print("Combined multiome signac object saved!")



combined.msobj <- readRDS("work/combined.multiome.merge.sobj.rds")



print("QC filtering")

DefaultAssay(combined.msobj) <- "RNA"
combined.msobj[["percent.mt"]] <- PercentageFeatureSet(combined.msobj, pattern = "^MT-")
combined.subset.msobj <- subset(
    x = combined.msobj, 
    subset = nFeature_RNA > 200 &
    nFeature_RNA < 4500 &
    percent.mt < 20 &
    nCount_peaks < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

print("Saving QC filtered combined multiome object to")
print(file.path(getwd(), "work", "combined.multiome.qc.sobj.rds"))

saveRDS(combined.subset.msobj, file = "work/combined.multiome.qc.sobj.rds")
print("Combined multiome signac object saved!")


DefaultAssay(combined.subset.msobj) <- "RNA"
combined.subset.msobj <- SCTransform(combined.subset.msobj)
combined.subset.msobj <- RunPCA(combined.subset.msobj, npcs = 50, verbose = T)


print("Processing combined ATAC data")
DefaultAssay(combined.subset.msobj) <- "peaks"
combined.subset.msobj <- FindTopFeatures(combined.subset.msobj, min.cutoff = "q50")
combined.subset.msobj <- RunTFIDF(combined.subset.msobj)
combined.subset.msobj <- RunSVD(combined.subset.msobj)

print ("Building joint neighbors based on both RNA and ATAC")
combined.subset.msobj <- FindMultiModalNeighbors(
  object = combined.subset.msobj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
print("Running UMAP on RNA and ATAC data")
combined.subset.msobj <- RunUMAP(
  object = combined.subset.msobj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)



print("Saving Combined multiome object to")
print(file.path(getwd(), "work", "combined.subset.multiome.before.sobj.rds"))

saveRDS(combined.subset.msobj, file = "work/combined.subset.multiome.before.sobj.rds")
print("Combined multiome signac object saved!")



print("Linking peaks to genes")
DefaultAssay(combined.subset.msobj) <- "peaks"




combined.subset.msobj <- RegionStats(combined.subset.msobj, genome = ref.genome)


combined.subset.msobj <- LinkPeaks(
  object = combined.subset.msobj,
  peak.assay = "peaks",
  expression.assay = "SCT"
)

print("Saving Combined multiome object to")
print(file.path(getwd(), "work", "combined.subset.multiome.sobj.rds"))

saveRDS(combined.subset.msobj, file = "work/combined.subset.multiome.sobj.rds")
print("Combined multiome signac object saved!")



sessionInfo()

rm(list=ls())
