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

org=args[[1]] 
taxlevel=args[[2]]
filepath=args[[3]]
algorithm=args[[4]]
dataframedescr=args[[5]]
metasubnegativepath=args[[6]]
outname = paste(org,algorithm,dataframedescr,taxlevel,sep='_')


#org='TEST' #args[[1]] #
#taxlevel='s__' #a rgs[[2]]
#filepath='xtree_merged/GTDB_.005_.0025_s_metagenomics_ra.tsv,xtree_merged/GTDB_.005_.0025_s_metatranscriptomics_ra.tsv'#args[[3]]
#algorithm='TEST'#args[[4]]
#dataframedescr='TEST'#args[[5]]
#metasubnegativepath='metasub_negatives/abundance_tables/GTDB_.005_.0025_s_metagenomics_ra.tsv'#args[[6]]
#outname = paste(org,algorithm,dataframedescr,taxlevel,sep='_')

# load metadata
meta = read.csv('~/Dropbox (Mason Lab)/i4/i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
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

# load abundance tables
if(algorithm=='phanta'){
  print('Loading phanta data')
  phanta = read.delim(filepath,check.names=F)  %>% filter(grepl('superkingdom_Viruses',Taxon_Lineage_with_Names)) %>% select(-Taxon_Lineage_with_IDs) %>% filter(grepl(taxlevel,Taxon_Lineage_with_Names)) %>% mutate(Taxon_Lineage_with_Names = strsplit(Taxon_Lineage_with_Names,taxlevel) %>% map_chr(2) %>% strsplit('\\|') %>% map_chr(1))
  phanta = phanta %>% melt %>% group_by(Taxon_Lineage_with_Names,variable) %>% summarise(value = sum(value)) %>% dcast(Taxon_Lineage_with_Names ~ variable,value.var = 'value')
  phanta = phanta%>% as.data.frame%>% column_to_rownames('Taxon_Lineage_with_Names') 
  phanta = sanitize_sample_names(phanta)
  # get all negative controls into their own df
  phantanegs = read.delim(metasubnegativepath,check.names=F)  %>% filter(grepl('superkingdom_Viruses',Taxon_Lineage_with_Names)) %>% filter(grepl(taxlevel,Taxon_Lineage_with_Names)) %>% select(-Taxon_Lineage_with_IDs)%>% mutate(Taxon_Lineage_with_Names = strsplit(Taxon_Lineage_with_Names,taxlevel) %>% map_chr(2) %>% strsplit('\\|') %>% map_chr(1))
  phantanegs = phantanegs%>% as.data.frame#%>% column_to_rownames('name') 
  phantanegs = phantanegs %>% melt %>% group_by(Taxon_Lineage_with_Names,variable) %>% summarise(value = sum(value)) %>% dcast(Taxon_Lineage_with_Names ~ variable,value.var = 'value')
  phantanegs2 = phanta %>% select(any_of(oac)) %>% rownames_to_column('Taxon_Lineage_with_Names')
  phantanegs = full_join(phantanegs,phantanegs2,by='Taxon_Lineage_with_Names') %>% column_to_rownames('Taxon_Lineage_with_Names')
  phantanegs[is.na(phantanegs)] = 0
  # separate out metagenomic and metatranscriptomic data and write to file
  phanta_metag = phanta %>% select(-any_of(oac)) %>% select(-all_of(grep('CEM',colnames(.))))
  phanta_metat = phanta %>% select(-any_of(oac)) %>% select(all_of(grep('CEM',colnames(.))))
  saveRDS(phanta_metag,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/phanta/',outname,'_metagenomics_nodecontam.rds',sep=''))
  saveRDS(phanta_metat,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/phanta/',outname,'_metatranscriptomics_nodecontam.rds',sep=''))
  # decontaminate
  phanta_metag_w_negs = full_join(phanta_metag %>% rownames_to_column('name'),phantanegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  phanta_metag_w_negs[is.na(phanta_metag_w_negs)] = 0
  otus_metag = data.frame(rownames(phanta_metag_w_negs),row.names=rownames(phanta_metag_w_negs))
  colnames(otus_metag) = 'Species'
  otus_metag=otus_metag[rownames(phanta_metag_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(phantanegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(phanta_metag)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(phanta_metag_w_negs),as.matrix(otus_metag),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/phanta/phanta_metag_contam_ ',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  phanta_metag_decontam <- phanta_metag[rownames(phanta_metag) %in% rownames(filtered_df),]
  phanta_metat_w_negs = full_join(phanta_metat %>% rownames_to_column('name'),phantanegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  phanta_metat_w_negs[is.na(phanta_metat_w_negs)] = 0
  otus_metat = data.frame(rownames(phanta_metat_w_negs),row.names=rownames(phanta_metat_w_negs))
  colnames(otus_metat) = 'Species'
  otus_metat=otus_metat[rownames(phanta_metat_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(phantanegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(phanta_metat)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(phanta_metat_w_negs),as.matrix(otus_metat),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/phanta/phanta_metat_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  phanta_metat_decontam <- phanta_metat[rownames(phanta_metat) %in% rownames(filtered_df),]
  # write new dfs to file
  saveRDS(phanta_metag_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/phanta/',outname,'_metagenomics_decontam.rds',sep=''))
  saveRDS(phanta_metat_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/phanta/',outname,'_metatranscriptomics_decontam.rds',sep=''))
}


# load abundance tables
if(algorithm=='kraken2'){
  print('Loading kraken2 data')
  kraken = read.delim(filepath,check.names=F)
  kraken = kraken %>% select(name,all_of(grep('bracken_frac',colnames(.)))) 
  colnames(kraken) = purrr::map(colnames(kraken),function(x) x %>% strsplit('\\.') %>% map_chr(1))
  kraken = kraken%>% as.data.frame%>% column_to_rownames('name') 
  kraken = sanitize_sample_names(kraken)
  # get all negative controls into their own df
  krakennegs = read.delim(metasubnegativepath,check.names=F)
  krakennegs = krakennegs %>% select(name,all_of(grep('bracken_frac',colnames(.)))) 
  colnames(krakennegs) = purrr::map(colnames(krakennegs),function(x) x %>% strsplit('\\.') %>% map_chr(1))
  krakennegs = krakennegs%>% as.data.frame#%>% column_to_rownames('name') 
  krakennegs2 = kraken %>% select(any_of(oac)) %>% rownames_to_column('name')
  krakennegs = full_join(krakennegs,krakennegs2,by='name') %>% column_to_rownames('name')
  krakennegs[is.na(krakennegs)] = 0
  # separate out metagenomic and metatranscriptomic data and write to file
  kraken_metag = kraken %>% select(-any_of(oac)) %>% select(-all_of(grep('CEM',colnames(.))))
  kraken_metat = kraken %>% select(-any_of(oac)) %>% select(all_of(grep('CEM',colnames(.))))
  saveRDS(kraken_metag,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/kraken/',outname,'_metagenomics_nodecontam.rds',sep=''))
  saveRDS(kraken_metat,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/kraken/',outname,'_metatranscriptomics_nodecontam.rds',sep=''))
  # decontaminate
  kraken_metag_w_negs = full_join(kraken_metag %>% rownames_to_column('name'),krakennegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  kraken_metag_w_negs[is.na(kraken_metag_w_negs)] = 0
  otus_metag = data.frame(rownames(kraken_metag_w_negs),row.names=rownames(kraken_metag_w_negs))
  colnames(otus_metag) = 'Species'
  otus_metag=otus_metag[rownames(kraken_metag_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(krakennegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(kraken_metag)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(kraken_metag_w_negs),as.matrix(otus_metag),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/kraken/kraken_metag_contam_ ',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  kraken_metag_decontam <- kraken_metag[rownames(kraken_metag) %in% rownames(filtered_df),]
  kraken_metat_w_negs = full_join(kraken_metat %>% rownames_to_column('name'),krakennegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  kraken_metat_w_negs[is.na(kraken_metat_w_negs)] = 0
  otus_metat = data.frame(rownames(kraken_metat_w_negs),row.names=rownames(kraken_metat_w_negs))
  colnames(otus_metat) = 'Species'
  otus_metat=otus_metat[rownames(kraken_metat_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(krakennegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(kraken_metat)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(kraken_metat_w_negs),as.matrix(otus_metat),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/kraken/kraken_metat_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  kraken_metat_decontam <- kraken_metat[rownames(kraken_metat) %in% rownames(filtered_df),]
  # write new dfs to file
  saveRDS(kraken_metag_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/kraken/',outname,'_metagenomics_decontam.rds',sep=''))
  saveRDS(kraken_metat_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/kraken/',outname,'_metatranscriptomics_decontam.rds',sep=''))
}

if(org=='bacteria' & algorithm=='xtree'){
  print('Loading xtree bacterial data')
  ### LOAD IN EACH FILE
  wgs_bacterial = read.csv(strsplit(filepath,',') %>% map_chr(1),sep='\t',check.names=F) %>% rownames_to_column('temp') %>% filter(temp!='Unknown') %>% column_to_rownames('temp')
  mtx_bacterial = read.csv(strsplit(filepath,',') %>% map_chr(2),sep='\t',check.names=F) %>% rownames_to_column('temp') %>% filter(temp!='Unknown') %>% column_to_rownames('temp')
  wgs_bacterial = sanitize_sample_names(wgs_bacterial)
  mtx_bacterial = sanitize_sample_names(mtx_bacterial)
  # get all negative controls into their own df
  xtreenegs = read.delim(metasubnegativepath,check.names=F)
  colnames(xtreenegs) = purrr::map(colnames(xtreenegs),function(x) x %>% strsplit('\\.') %>% map_chr(1))
  xtreenegs = xtreenegs%>% as.data.frame %>% rownames_to_column('name') 
  xtreenegs2 = wgs_bacterial %>% select(any_of(oac)) %>% rownames_to_column('name')
  xtreenegs3 = mtx_bacterial %>% select(any_of(oac)) %>% rownames_to_column('name')
  xtreenegs = full_join(xtreenegs,xtreenegs2,by='name') 
  xtreenegs = full_join(xtreenegs,xtreenegs3,by='name') %>% column_to_rownames('name')
  xtreenegs[is.na(xtreenegs)] = 0
  # separate out metagenomic and metatranscriptomic data and write to file
  wgs_bacterial = wgs_bacterial %>% select(-any_of(oac))
  mtx_bacterial = mtx_bacterial %>% select(-any_of(oac)) 
  saveRDS(wgs_bacterial,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/gtdb/',outname,'_metagenomics_nodecontam.rds',sep=''))
  saveRDS(mtx_bacterial,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/gtdb/',outname,'_metatranscriptomics_nodecontam.rds',sep=''))
  # decontaminate
  wgs_bacterial_w_negs = full_join(wgs_bacterial %>% rownames_to_column('name'),xtreenegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  wgs_bacterial_w_negs[is.na(wgs_bacterial_w_negs)] = 0
  otus_metag = data.frame(rownames(wgs_bacterial_w_negs),row.names=rownames(wgs_bacterial_w_negs))
  colnames(otus_metag) = 'Species'
  otus_metag=otus_metag[rownames(wgs_bacterial_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(xtreenegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(wgs_bacterial)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(wgs_bacterial_w_negs),as.matrix(otus_metag),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/gtdb/gtdb_metag_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  gtdb_metag_decontam <- wgs_bacterial[rownames(wgs_bacterial) %in% rownames(filtered_df),]
  mtx_bacterial_w_negs = full_join(mtx_bacterial %>% rownames_to_column('name'),xtreenegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  mtx_bacterial_w_negs[is.na(mtx_bacterial_w_negs)] = 0
  otus_metat = data.frame(rownames(mtx_bacterial_w_negs),row.names=rownames(mtx_bacterial_w_negs))
  colnames(otus_metat) = 'Species'
  otus_metat=otus_metat[rownames(mtx_bacterial_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(xtreenegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(mtx_bacterial)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(mtx_bacterial_w_negs),as.matrix(otus_metat),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/gtdb/gtdb_metat_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  gtdb_metat_decontam <- mtx_bacterial[rownames(mtx_bacterial) %in% rownames(filtered_df),]
  # write new dfs to file
  saveRDS(gtdb_metag_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/gtdb/',outname,'_metagenomics_decontam.rds',sep=''))
  saveRDS(gtdb_metat_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/gtdb/',outname,'_metatranscriptomics_decontam.rds',sep=''))
}

if(org=='virus' & algorithm=='xtree'){
  ### LOAD IN EACH FILE
  wgs_viral <- read.delim(strsplit(filepath,',') %>% map_chr(1), sep='\t', check.names=F, na.strings="NOT_NA")
  wgs_viral <- wgs_viral[rownames(wgs_viral) != "NA", ]
  wgs_viral <- wgs_viral[rownames(wgs_viral) != "Unknown", ]
  mtx_viral <- read.delim(strsplit(filepath,',') %>% map_chr(2), sep='\t', check.names=F, na.strings="NOT_NA")
  mtx_viral <- mtx_viral[rownames(mtx_viral) != "NA", ]
  mtx_viral <- mtx_viral[rownames(mtx_viral) != "Unknown", ]
  wgs_viral = sanitize_sample_names(wgs_viral)
  mtx_viral = sanitize_sample_names(mtx_viral)
  # get all negative controls into their own df
  xtreenegs <- read.delim(metasubnegativepath,check.names=F, na.strings="NOT_NA")
  xtreenegs <- xtreenegs[rownames(xtreenegs) != "NA", ]
  xtreenegs <- xtreenegs[rownames(xtreenegs) != "Unknown", ]
  colnames(xtreenegs) = purrr::map(colnames(xtreenegs),function(x) x %>% strsplit('\\.') %>% map_chr(1))
  xtreenegs = xtreenegs%>% as.data.frame %>% rownames_to_column('name') 
  xtreenegs2 = wgs_viral %>% select(any_of(oac)) %>% rownames_to_column('name')
  xtreenegs3 = mtx_viral %>% select(any_of(oac)) %>% rownames_to_column('name')
  xtreenegs = full_join(xtreenegs,xtreenegs2,by='name') 
  xtreenegs = full_join(xtreenegs,xtreenegs3,by='name') %>% column_to_rownames('name')
  xtreenegs[is.na(xtreenegs)] = 0
  # separate out metagenomic and metatranscriptomic data and write to file
  wgs_viral = wgs_viral %>% select(-any_of(oac))
  mtx_viral = mtx_viral %>% select(-any_of(oac)) 
  saveRDS(wgs_viral,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/genbank/',outname,'_metagenomics_nodecontam.rds',sep=''))
  saveRDS(mtx_viral,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/genbank/',outname,'_metatranscriptomics_nodecontam.rds',sep=''))
  # decontaminate
  wgs_viral_w_negs = full_join(wgs_viral %>% rownames_to_column('name'),xtreenegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  wgs_viral_w_negs[is.na(wgs_viral_w_negs)] = 0
  otus_metag = data.frame(rownames(wgs_viral_w_negs),row.names=rownames(wgs_viral_w_negs))
  colnames(otus_metag) = 'Species'
  otus_metag=otus_metag[rownames(wgs_viral_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(xtreenegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(wgs_viral)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(wgs_viral_w_negs),as.matrix(otus_metag),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg",threshold=.1)
  saveRDS(contamdf,paste('revisions/genbank/genbank_metag_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  genbank_metag_decontam <- wgs_viral[rownames(wgs_viral) %in% rownames(filtered_df),]
  mtx_viral_w_negs = full_join(mtx_viral %>% rownames_to_column('name'),xtreenegs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  mtx_viral_w_negs[is.na(mtx_viral_w_negs)] = 0
  otus_metat = data.frame(rownames(mtx_viral_w_negs),row.names=rownames(mtx_viral_w_negs))
  colnames(otus_metat) = 'Species'
  otus_metat=otus_metat[rownames(mtx_viral_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(xtreenegs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(mtx_viral)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(mtx_viral_w_negs),as.matrix(otus_metat),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg",threshold=.1)
  saveRDS(contamdf,paste('revisions/genbank/genbank_metat_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  genbank_metat_decontam <- mtx_viral[rownames(mtx_viral) %in% rownames(filtered_df),]
  # write new dfs to file
  saveRDS(genbank_metag_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/genbank/',outname,'_metagenomics_decontam.rds',sep=''))
  saveRDS(genbank_metat_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/genbank/',outname,'_metatranscriptomics_decontam.rds',sep=''))
}

if(algorithm=='metaphlan4'){
  print('Loading metaphlan4 data')
  if(taxlevel == 's__'){
    nextone = 't__'
  }
  if(taxlevel == 'g__'){
    nextone = 's__'
  }
  if(taxlevel == 'f__'){
    nextone = 'g__'
  }
  if(taxlevel == 'o__'){
    nextone = 'f__'
  }
  if(taxlevel == 'c__'){
    nextone = 'o__'
  }
  if(taxlevel == 'p__'){
    nextone = 'c__'
  }
  if(taxlevel == 't__'){
    nextone = 'nullnullnullnullnull'  
  }
  metaphlan4 = read.delim(filepath,check.names=F,sep='\t') %>% dplyr::rename(name=clade_name) %>% filter(grepl(taxlevel,name)) %>% filter(!grepl(nextone,name))
  colnames(metaphlan4) = purrr::map(colnames(metaphlan4),function(x) x %>% gsub('_metaphlan','',.))
  metaphlan4 = metaphlan4 %>% as.data.frame %>% column_to_rownames('name') 
  metaphlan4 = sanitize_sample_names(metaphlan4)
  # get all negative controls into their own df
  metaphlan4negs = read.delim(metasubnegativepath,check.names=F)%>% dplyr::rename(name=clade_name)%>% filter(grepl(taxlevel,name)) %>% filter(!grepl(nextone,name))
  metaphlan4negs = metaphlan4negs %>% select(name,all_of(grep('bracken_frac',colnames(.)))) 
  colnames(metaphlan4) = purrr::map(colnames(metaphlan4),function(x) x %>% gsub('_metaphlan','',.))
  metaphlan4negs = metaphlan4negs%>% as.data.frame#%>% column_to_rownames('name') 
  metaphlan4negs2 = metaphlan4 %>% select(any_of(oac)) %>% rownames_to_column('name')
  metaphlan4negs = full_join(metaphlan4negs,metaphlan4negs2,by='name') %>% column_to_rownames('name')
  metaphlan4negs[is.na(metaphlan4negs)] = 0
  # separate out metagenomic and metatranscriptomic data and write to file
  metaphlan4_metag = metaphlan4 %>% select(-any_of(oac)) %>% select(-all_of(grep('CEM',colnames(.))))
  metaphlan4_metat = metaphlan4 %>% select(-any_of(oac)) %>% select(all_of(grep('CEM',colnames(.))))
  saveRDS(metaphlan4_metag,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/metaphlan4/',outname,'metagenomics_nodecontam.rds',sep=''))
  saveRDS(metaphlan4_metat,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/metaphlan4/',outname,'metatranscriptomics_nodecontam.rds',sep=''))
  # decontaminate
  metaphlan4_metag_w_negs = full_join(metaphlan4_metag %>% rownames_to_column('name'),metaphlan4negs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  metaphlan4_metag_w_negs[is.na(metaphlan4_metag_w_negs)] = 0
  otus_metag = data.frame(rownames(metaphlan4_metag_w_negs),row.names=rownames(metaphlan4_metag_w_negs))
  colnames(otus_metag) = 'Species'
  otus_metag=otus_metag[rownames(metaphlan4_metag_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(metaphlan4negs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(metaphlan4_metag)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(metaphlan4_metag_w_negs),as.matrix(otus_metag),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/metaphlan4/metaphlan4_metag_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  metaphlan4_metag_decontam <- metaphlan4_metag[rownames(metaphlan4_metag) %in% rownames(filtered_df),]  
  metaphlan4_metat_w_negs = full_join(metaphlan4_metat %>% rownames_to_column('name'),metaphlan4negs%>% rownames_to_column('name'),by='name') %>% column_to_rownames('name')
  metaphlan4_metat_w_negs[is.na(metaphlan4_metat_w_negs)] = 0
  otus_metat = data.frame(rownames(metaphlan4_metat_w_negs),row.names=rownames(metaphlan4_metat_w_negs))
  colnames(otus_metat) = 'Species'
  otus_metat=otus_metat[rownames(metaphlan4_metat_w_negs),,drop=F]
  mdatframe = bind_rows(data.frame(colnames(metaphlan4negs)) %>% mutate(Sample_or_Control = 'Control') %>% column_to_rownames(colnames(.)[1]),data.frame(colnames(metaphlan4_metat)) %>% mutate(Sample_or_Control = 'True')%>% column_to_rownames(colnames(.)[1]))
  ps = convert_to_phyloseq(as.matrix(metaphlan4_metat_w_negs),as.matrix(otus_metat),(mdatframe))
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
  contamdf <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=.1)
  saveRDS(contamdf,paste('revisions/metaphlan4/metaphlan4_metat_contam',outname,'.rds',sep=''))
  filtered_df <- contamdf %>% filter(contaminant == FALSE)
  metaphlan4_metat_decontam <- metaphlan4_metat[rownames(metaphlan4_metat) %in% rownames(filtered_df),]
  # write new dfs to file
  saveRDS(metaphlan4_metag_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/metaphlan4/',outname,'metagenomics_decontam.rds',sep=''))
  saveRDS(metaphlan4_metat_decontam,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/metaphlan4/',outname,'metatranscriptomics_decontam.rds',sep=''))
}

if(org=='bacteria' & algorithm=='MAG'){
  print('Loading bacterial MAG data')
  ### LOAD IN EACH FILE
  wgs_bacterial = read.csv(strsplit(filepath,',') %>% map_chr(1),sep='\t',check.names=F) %>% rownames_to_column('temp') %>% filter(temp!='Unknown') %>% column_to_rownames('temp')
  mtx_bacterial = read.csv(strsplit(filepath,',') %>% map_chr(2),sep='\t',check.names=F) %>% rownames_to_column('temp') %>% filter(temp!='Unknown') %>% column_to_rownames('temp')
  wgs_bacterial = sanitize_sample_names(wgs_bacterial)
  mtx_bacterial = sanitize_sample_names(mtx_bacterial)
  wgs_bacterial = wgs_bacterial %>% select(-any_of(oac))
  mtx_bacterial = mtx_bacterial %>% select(-any_of(oac)) 
  saveRDS(wgs_bacterial,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/bacterial_MAG/',outname,'_metagenomics_nodecontam.rds',sep=''))
  saveRDS(mtx_bacterial,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/bacterial_MAG/',outname,'_metatranscriptomics_nodecontam.rds',sep=''))
}

if(org=='virus' & algorithm=='MAG'){
  print('Loading viral MAG data')
  ### LOAD IN EACH FILE
  wgs_bacterial = read.csv(strsplit(filepath,',') %>% map_chr(1),sep='\t',check.names=F) %>% rownames_to_column('temp') %>% filter(temp!='Unknown') %>% column_to_rownames('temp')
  mtx_bacterial = read.csv(strsplit(filepath,',') %>% map_chr(2),sep='\t',check.names=F) %>% rownames_to_column('temp') %>% filter(temp!='Unknown') %>% column_to_rownames('temp')
  wgs_bacterial = sanitize_sample_names(wgs_bacterial)
  mtx_bacterial = sanitize_sample_names(mtx_bacterial)
  wgs_bacterial = wgs_bacterial %>% select(-any_of(oac))
  mtx_bacterial = mtx_bacterial %>% select(-any_of(oac)) 
  saveRDS(wgs_bacterial,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/viral_MAG/',outname,'_metagenomics_nodecontam.rds',sep=''))
  saveRDS(mtx_bacterial,paste('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/viral_MAG/',outname,'_metatranscriptomics_nodecontam.rds',sep=''))
}


