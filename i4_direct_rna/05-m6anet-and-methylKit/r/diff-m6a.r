#!/usr/bin/env -S Rscript --vanilla
library("zeallot")
library("methylKit")

c(ann_fn, subjects, timepoints, output_fn) %<-% commandArgs(trailingOnly=T)
annotation = read.table(ann_fn, header=T)

if (subjects != "ALL") {
    annotation = annotation[
        (annotation$subject %in% unlist(strsplit(subjects, ","))),
    ]
}
if (timepoints != "ALL") {
    annotation = annotation[
        (annotation$timepoint %in% unlist(strsplit(timepoints, ","))),
    ]
}
if (min(annotation$pre_post_rec) == 0) {
    annotation$treatment = sapply(
        annotation$pre_post_rec, function (x) {min(x, 1)}
    )
} else {
    annotation$treatment = annotation$pre_post_rec - 1
}

input = methRead(
    as.list(annotation$input), sample.id=as.list(annotation$sample),
    treatment=as.vector(annotation$treatment), mincov=10,
    assembly="gencode.v41", context="m6A", pipeline=list(
        fraction=T, chr.col=1, start.col=2, end.col=2,
        coverage.col=3, freqC.col=4, strand.col=5
    )
)

mdiff = getMethylDiff(calculateDiffMeth(unite(input, destrand=F)), difference=0, qvalue=1)
write.csv(mdiff, output_fn, row.names=F, quote=F)
