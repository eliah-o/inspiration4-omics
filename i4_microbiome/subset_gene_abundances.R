# get gene subset

library(tidyverse)

parsed = read.delim('/athena/masonlab/scratch/projects/spacex/analysis/i4_bodyswabs/immune_analysis/regression_data_timetrends_species_decontam.csv',sep=',')

siggenes = parsed %>% filter(dset == 'GENE-CATALOG') %>% select(yvar) %>% unlist %>% unname %>% unique

metag = readRDS('gene_metagenomics_decontam.rds')
metag2 = metag[siggenes,]
saveRDS(metag2,'gene_catalog90_filtered_relab_metagenomics_subset.rds')

metat = readRDS('gene_metatranscriptomics_decontam.rds')
metat2 = metat[siggenes,]
saveRDS(metat2,'gene_catalog90_filtered_relab_metatranscriptomics_subset.rds',)



























Amazing -- so these are all abundance tables (or metadata) for different microbial taxa or genes. 

Gene catalog info:

THIS FOLDER:

/athena/masonlab/scratch/projects/spacex/analysis/i4_bodyswabs/data_packet/gene_catalog/

THESE FILES:

gene_catalog90_filtered_relab_metagenomics.tsv -- Non-redundant gene relative abundances (90% identity clustering, normalized for sequencing depth and read length), metagenomic sample alignments.
gene_catalog90_filtered_relab_metatranscriptomics.tsv -- Non-redundant gene relative abundances (90% identity clustering, normalized for sequencing depth and read length), metatranscriptomic sample alignments.
all_prot_seqs_filteredpseudo_genecat30_rep_seq.fasta -- Consensus gene sequences, 30% identity clustering.
all_prot_seqs_filteredpseudo_genecat50_rep_seq.fasta -- Consensus gene sequences, 50% identity clustering.
all_prot_seqs_filteredpseudo_genecat70_rep_seq.fasta -- Consensus gene sequences, 70% identity clustering.
all_prot_seqs_filteredpseudo_genecat90_rep_seq.fasta -- Consensus gene sequences, 90% identity clustering. 
i4_calledORFs_full_annotation_data.tsv -- Annotation data for all called ORFs.

THIS FOLDER:

THESE FILES:







