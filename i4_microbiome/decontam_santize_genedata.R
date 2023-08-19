#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(lme4)
library(broom)
library(broom.mixed)
library(reshape2)
library(lmerTest)
library(decontam)
library(phyloseq)

# load metadata
meta = read.csv('../../i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
meta$SeqID = gsub('ELMB_','',meta$SeqID)
meta$SeqID = gsub('SW_','',meta$SeqID)
meta$Timepoint_Recode = factor(meta$Timepoint)
levels(meta$Timepoint_Recode) = c(NA,'PRE-LAUNCH','POST-LAUNCH','PRE-LAUNCH','MID-FLIGHT','MID-FLIGHT','POST-LAUNCH','POST-LAUNCH','PRE-LAUNCH')

meta = meta %>% distinct %>% mutate(Timepoint_Recode2 = if_else(as.character(Timepoint_Recode) == 'MID-FLIGHT',Timepoint,as.character(Timepoint_Recode)))
meta$Timepoint_Recode2 = factor(meta$Timepoint_Recode2,levels = c('PRE-LAUNCH','Flight 1','Flight 2','POST-LAUNCH'))
meta$Timepoint = factor(meta$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec',NA))
meta$Timepoint_Numeric = as.numeric(meta$Timepoint)

convert_to_phyloseq <- function(otu_mat,tax_mat,df.meta){
  phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
  phylo_TAX<- tax_table(tax_mat)
  phylo_samples<- sample_data(df.meta)
  return(phyloseq(phylo_OTU, phylo_TAX, phylo_samples))
}

sanitize_sample_names <- function(data){
  temp = data %>% t %>% as.data.frame %>% rownames_to_column('temp') %>% mutate(namelengths = nchar(temp))
  temp = temp %>% mutate(temp = if_else(namelengths>=3 & str_sub(temp,nchar(temp),nchar(temp))=='D',str_sub(temp,1,nchar(temp)-1),temp))
  return(temp %>% column_to_rownames('temp') %>% select(-namelengths) %>% t %>% data.frame(check.names=F))
}

# FIND NEGATIVE CONTROLS
oac = meta %>% filter(location == 'Open Air Control' | Body.Location == 'Control Swab (0)' | Body.Location == 'Swab Water')
oac = oac %>% select(SeqID) %>% unlist %>% unname 

print('Loading WGS data')
data = read.delim('gene_catalog90_filtered_relab_metagenomics.tsv',header=T,sep='\t',check.names=F)
colnames(data)[1]='X'
data = data%>% as.data.frame(check.names=F)%>% column_to_rownames('X') 
data = sanitize_sample_names(data)
# get all negative controls into their own df
datanegs = data %>% select(any_of(oac)) %>% rownames_to_column('name')
datanegs[is.na(datanegs)] = 0
# separate out metagenomic and metatranscriptomic data and write to file
data_metag = data %>% select(-any_of(oac)) %>% select(-all_of(grep('CEM',colnames(.))))
saveRDS(data_metag,paste('gene_metagenomics_nodecontam.rds',sep=''))
# decontaminate
data_metag_w_negs = full_join(data_metag %>% rownames_to_column('name'),datanegs,by='name') %>% column_to_rownames('name')
data_metag_w_negs[is.na(data_metag_w_negs)] = 0
otus_metag = data.frame(rownames(data_metag_w_negs),row.names=rownames(data_metag_w_negs))
colnames(otus_metag) = 'Species'
otus_metag=otus_metag[rownames(data_metag_w_negs),,drop=F]
mdatframe = bind_rows(as.data.frame(colnames(datanegs),check.names=F) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),as.data.frame(colnames(data_metag),check.names=F) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
ps = convert_to_phyloseq(as.matrix(data_metag_w_negs),as.matrix(otus_metag),(mdatframe))
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = .1)
saveRDS(contamdf,paste('data_metag_contam_gene.rds',sep=''))
filtered_df <- contamdf %>% filter(contaminant == FALSE)
data_metag_decontam <- data_metag[rownames(data_metag) %in% rownames(filtered_df),]
saveRDS(data_metag_decontam,paste('gene_metagenomics_decontam.rds',sep=''))

data_metag_decontam_list = split(data_metag_decontam, rep(1:ceiling(nrow(data_metag_decontam)/15000), each=15000, length.out=nrow(data_metag_decontam)))

count = 0
for(val in data_metag_decontam_list){
  count = count + 1
  print(count)
  write.table(file = paste('gene_metag_',count,'.tsv',sep=''),x=val,quote=F,sep='\t')
  saveRDS(file = paste('gene_metag_',count,'.rds',sep=''),object = val)
}

data = read.delim('gene_catalog90_filtered_relab_metatranscriptomics.tsv',header=T,sep='\t',check.names=F)
colnames(data)[1]='X'
data = data%>% as.data.frame(check.names=F)%>% column_to_rownames('X') 
data = sanitize_sample_names(data)
# get all negative controls into their own df
datanegs = data %>% select(any_of(oac)) %>% rownames_to_column('name')
datanegs[is.na(datanegs)] = 0
# separate out metatenomic and metatranscriptomic data and write to file
data_metat = data %>% select(-any_of(oac)) %>% select(all_of(grep('CEM',colnames(.))))
saveRDS(data_metat,paste('gene_metatranscriptomics_nodecontam.rds',sep=''))
# decontaminate
data_metat_w_negs = full_join(data_metat %>% rownames_to_column('name'),datanegs,by='name') %>% column_to_rownames('name')
data_metat_w_negs[is.na(data_metat_w_negs)] = 0
otus_metat = as.data.frame(rownames(data_metat_w_negs),row.names=rownames(data_metat_w_negs),check.names=F)
colnames(otus_metat) = 'Species'
otus_metat=otus_metat[rownames(data_metat_w_negs),,drop=F]
mdatframe = bind_rows(as.data.frame(colnames(datanegs),check.names=F) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),as.data.frame(colnames(data_metat),check.names=F) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
ps = convert_to_phyloseq(as.matrix(data_metat_w_negs),as.matrix(otus_metat),(mdatframe))
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf <- isContaminant(ps, method="prevalence", neg="is.neg",threshold = .1)
saveRDS(contamdf,paste('data_metat_contam_gene.rds',sep=''))
filtered_df <- contamdf %>% filter(contaminant == FALSE)
data_metat_decontam <- data_metat[rownames(data_metat) %in% rownames(filtered_df),]
saveRDS(data_metat_decontam,paste('gene_metatranscriptomics_decontam.rds',sep=''))

data_metat_decontam_list = split(data_metat_decontam, rep(1:ceiling(nrow(data_metat_decontam)/15000), each=15000, length.out=nrow(data_metat_decontam)))

count = 0
for(val in data_metat_decontam_list){
  count = count + 1
  print(count)
  write.table(file = paste('gene_metat_',count,'.tsv',sep=''),x=val,quote=F,sep='\t')
  saveRDS(file = paste('gene_metat_',count,'.rds',sep=''),object=val)
}
  
  