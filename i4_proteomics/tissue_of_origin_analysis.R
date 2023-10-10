#Tissue of origin deconvolution based on fgsea
rm(list=ls())
setwd("/./i4_proteomics_analysis/")
library(fgsea)
library(enrichplot)

#Load datasets compiled from the Human Protein Atlas
tissue_elevated <- gmtPathways("~/Downloads/tissue_elevated.gmt")
tissue_specific <- gmtPathways("~/Downloads/tissue_specific.gmt")
exosomes_immPostVpre <- read_csv("./results/exosomes_limma_DE_immPostVpre_arrayWeighted.csv")

#Rank by t statistic
exosomes_immPostVpre <- exosomes_immPostVpre[order(exosomes_immPostVpre$t),]
ranks <- as.vector(exosomes_immPostVpre$t)
names(ranks) <- exosomes_immPostVpre$Gene

#fgsea analysis
fgseaRes_elevated <- fgsea(tissue_elevated, ranks, minSize=5, maxSize=500)
write_csv(fgseaRes_elevated, "./tissue_enriched_exosomes.csv")

fgseaRes_elevated <- fgsea(tissue_specific, ranks, minSize=5, maxSize=500)
write_csv(fgseaRes_elevated, "./tissue_specific_exosomes.csv")
