library(tidyverse)

# remove non significant findings to make the whole genecat processing faster

a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_oral_merged_gene_metag.tsv",sep=' ') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% filter(grepl('Timepoint',term)) %>% mutate(dset = 'GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_oral_merged_gene_metag_FILTERED.rds')

a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_oral_merged_gene_metat.tsv",sep=' ') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% filter(grepl('Timepoint',term)) %>% mutate(dset = 'GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_oral_merged_gene_metat_FILTERED.rds')

a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_nasal_merged_gene_metag.tsv",sep=' ') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal') %>% filter(grepl('Timepoint',term))   %>% mutate(dset = 'GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_nasal_merged_gene_metag_FILTERED.rds')

a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_nasal_merged_gene_metat.tsv",sep=' ') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% filter(grepl('Timepoint',term)) %>% mutate(dset = 'GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_nasal_merged_gene_metat_FILTERED.rds')


a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skin_merged_gene_metag.tsv",sep=' ') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% filter(grepl('Timepoint',term)) %>% mutate(dset = '
	GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skin_merged_gene_metag_FILTERED.rds')


a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skin_merged_gene_metat.tsv",sep=' ') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% filter(grepl('Timepoint',term)) %>% mutate(dset = 'GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skin_merged_gene_metat_FILTERED.rds')


a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skinseparates_merged_gene_metag.tsv",sep=' ') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% filter(grepl('Timepoint',term)) %>% mutate(dset = 'GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skinseparates_merged_gene_metag_FILTERED.rds')

a = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skinseparates_merged_gene_metat.tsv",sep=' ') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% filter(grepl('Timepoint',term)) %>% mutate(dset = 'GENE-CATALOG') %>% filter(BH_adjusted<0.05 | BY_adjusted<0.05 )
saveRDS(a,'/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skinseparates_merged_gene_metat_FILTERED.rds')

