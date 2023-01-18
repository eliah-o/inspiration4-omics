# Inspiration4 metagenomics and metatranscriptomics analysis

You can find here the scripts used to run regressions and generate figures for the analysis used in the i4 metagenomics an metatranscriptomics manuscript, which is currently submitted and being reviewed. This README will be updated with the title and abstract upon eventual publication.

### i4_swab_metagenomics.Rmd

This R markdown contains the code used to generate all the figures presented in the manuscript. It also categorizes microbial features into the groups used in the paper (e.g., transiently increased).

### run_regressions.R

This was the script used to run the mixed model associations for all taxonomic data (e.g., bacteria, viruses).

### run_regressions_gene.R

This was the script used to run the mixed model gene catalog associations, which were batched to allow for parallel computing on a slurm cluster (e.g., take 5000 genes at a time, run the regressions, and later combine them and adjust p-values).

### run_immune_lasso_associations.R

This is the script used to compute associations between significant microbial features and human immune cell subtype gene expression.
