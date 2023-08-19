library(tidyverse)
library(ggplot2)
library(umap)
library(lme4)
library(broom.mixed)
library(gridExtra)
library(broom)
library(ape)
library(cowplot)
library(reshape2)
library(taxonomizr)
library(circlize)
library(ggtree)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggbeeswarm)
library(lmerTest)
library(ComplexUpset)

args = commandArgs(trailingOnly=TRUE)

setwd('~/Dropbox (Mason Lab)/i4/i4_data_packet/')

sanitize_sample_names <- function(data){
  temp = data%>% rownames_to_column('temp') %>% mutate(namelengths = nchar(temp))
  temp = temp %>% mutate(temp = if_else(namelengths>=3 & str_sub(temp,nchar(temp),nchar(temp))=='D',str_sub(temp,1,nchar(temp)-1),temp))
  return(temp %>% column_to_rownames('temp') %>% select(-namelengths))
}

# load metadata
meta = read.csv('../i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
meta$SeqID = gsub('ELMB_','',meta$SeqID)
meta$SeqID = gsub('SW_','',meta$SeqID)
meta$Timepoint_Recode = factor(meta$Timepoint)
levels(meta$Timepoint_Recode) = c(NA,'PRE-LAUNCH','POST-LAUNCH','PRE-LAUNCH','MID-FLIGHT','MID-FLIGHT','POST-LAUNCH','POST-LAUNCH','PRE-LAUNCH')

meta = meta %>% distinct %>% mutate(Timepoint_Recode2 = if_else(as.character(Timepoint_Recode) == 'MID-FLIGHT',Timepoint,as.character(Timepoint_Recode)))
meta$Timepoint_Recode2 = factor(meta$Timepoint_Recode2,levels = c('PRE-LAUNCH','Flight 1','Flight 2','POST-LAUNCH'))
meta$Timepoint = factor(meta$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec',NA))
meta$Timepoint_Numeric = as.numeric(meta$Timepoint)
meta$Timepoint_Recode = factor(meta$Timepoint_Recode,levels = c('MID-FLIGHT','PRE-LAUNCH','POST-LAUNCH'))

#meta = meta %>% filter(Body.Location != "Swab Water",location!='Capsule',Body.Location != 'Open Air Control', Body.Location != "Deltoid - Pre-Biospy", !is.na(Body.Location))
meta = meta %>% mutate(isoral = if_else(Body.Location == 'Oral',1,0))
meta = meta %>% mutate(isnasal = if_else(Body.Location == 'Nasal',1,0))
meta = meta %>% mutate(isskin = if_else(Body.Location != 'Oral' & Body.Location != 'Nasal',1,0))
meta = meta %>% mutate(Armpit = if_else(Body.Location == 'Armpit',1,0))
meta = meta %>% mutate(web = if_else(Body.Location == 'Toe Web Space',1,0))
meta = meta %>% mutate(nape = if_else(Body.Location == 'Nape of Neck',1,0))
meta = meta %>% mutate(postauric = if_else(Body.Location == 'Post-Auricular',1,0))
meta = meta %>% mutate(fore = if_else(Body.Location == 'Forearm',1,0))
meta = meta %>% mutate(bb = if_else(Body.Location == 'Belly Button',1,0))
meta = meta %>% mutate(gc = if_else(Body.Location == 'Gluteal Crease',1,0))
meta = meta %>% mutate(nasal = if_else(Body.Location == 'Nasal',1,0))
meta = meta %>% mutate(Tzone = if_else(Body.Location == 'T-Zone',1,0))
meta = meta %>% mutate(Oral = if_else(Body.Location == 'Oral',1,0))

# kraken metag
kraken_metag = read.delim('kraken2_bracken_abundances/kraken2_bracken_metagenomics.tsv',check.names=F)
kraken_tax1 = kraken_metag %>% select(name,taxonomy_id)
kraken_metag_num = kraken_metag %>% select(name,all_of(grep('bracken_num',colnames(.))))%>% melt
counts = kraken_metag_num %>% ungroup %>% select(name,value) %>% group_by(name)  %>% mutate(count = if_else(value>0,1,0))  %>% summarise(s = sum(count))
tokeep = counts%>% select(name) %>% unlist %>% unname
kraken_metag = kraken_metag %>% select(name,all_of(grep('bracken_frac',colnames(.)))) 
colnames(kraken_metag) = map(colnames(kraken_metag),function(x) gsub('.qc.bracken_frac','',x))
kraken_metag = kraken_metag%>%  data.frame(check.names=F)%>% column_to_rownames('name') %>% t %>% data.frame(check.names=F)
kraken_metag = sanitize_sample_names(kraken_metag)

# kraken metat
kraken_metat = read.delim('kraken2_bracken_abundances/kraken2_bracken_metatranscriptomics.tsv',check.names=F)
kraken_tax2 = kraken_metat %>% select(name,taxonomy_id)
kraken_metat_num = kraken_metat %>% select(name,all_of(grep('bracken_num',colnames(.))))%>% melt
counts = kraken_metat_num %>% ungroup %>% select(name,value) %>% group_by(name) %>% mutate(count = if_else(value>0,1,0))  %>% summarise(s = sum(count))
tokeep = counts  %>% select(name) %>% unlist %>% unname
kraken_metat = kraken_metat %>% select(name,all_of(grep('bracken_frac',colnames(.))))%>% data.frame(check.names=F)
colnames(kraken_metat) = map(colnames(kraken_metat),function(x) gsub('.qc.bracken_frac','',x))
kraken_metat = kraken_metat%>% as.data.frame%>% column_to_rownames('name') %>% t%>% data.frame(check.names=F)

# bac GTDB metag/metat
gtdb_metag = read.table("GTDB_r207_xtree_abundances/GTDB_.005_.0025_metagenomics_species_ra.tsv",check.names=F,sep='\t')
gtdb_metat = read.table('GTDB_r207_xtree_abundances/GTDB_.005_.0025_metatranscriptomics_species_ra.tsv',check.names=F,sep='\t') %>% t %>% data.frame(check.names=F)

gtdb_metag = sanitize_sample_names(gtdb_metag  %>% t %>% as.data.frame(check.names=F))

# bac assembled metag/metat
bass_metag = read.table("ASSEMBLED-BACTERIAL-MAGS_.01_.005_metagenomics_species_ra.tsv",check.names=F,sep='\t')
bass_metat = read.table('ASSEMBLED-BACTERIAL-MAGS_.01_.005_metatranscriptomics_species_ra.tsv',check.names=F,sep='\t')

bass_metag = sanitize_sample_names(bass_metag %>% t %>% as.data.frame(check.names=F))

# viral metag/metat
viral_metag = read.csv("~/Dropbox (Mason Lab)/i4/i4_data_packet/GenBank_viral_xtree_abundances/genbank-viral_.02_.0.015_metagenomics_species_ra.tsv",check.names=F,sep='\t')
viral_metat = read.csv('~/Dropbox (Mason Lab)/i4/i4_data_packet/GenBank_viral_xtree_abundances/genbank-viral_.02_.0.015_metatranscriptomics_species_ra.tsv',check.names=F,sep='\t')

viral_metag = sanitize_sample_names(viral_metag  %>% t %>% as.data.frame(check.names=F))

# viral assembled metag/metat
viral_ass_metag = read.table("ASSEMBLED-VIRAL-GENOMES_.01_.0.005_metagenomics_species_ra.tsv",check.names=F,sep='\t')
viral_ass_metat = read.table('ASSEMBLED-VIRAL-GENOMES_.01_.0.005_metatranscriptomics_species_ra.tsv',check.names=F,sep='\t')

viral_ass_metag = sanitize_sample_names(viral_ass_metag %>% t %>% as.data.frame(check.names=F))

# Gene metag/metat

#### NEED TO UPDATE TO LOAD ONLY SIGNIFICANT GENES
gene_metag = read.table('~/Dropbox (Mason Lab)/i4/i4_data_packet/gene_catalog/gene_catalog90_filtered_relab_metagenomics_subset.tsv',sep='\t',check.names=F,header = T,row.names = 1) 
gene_metag = sanitize_sample_names(gene_metag %>% t %>% as.data.frame(check.names=F))

gene_metat = read.table('~/Dropbox (Mason Lab)/i4/i4_data_packet/gene_catalog/gene_catalog90_filtered_relab_metatranscriptomics_subset.tsv',sep='\t',check.names=F,header = T,row.names = 1)
gene_metat = gene_metat %>% t %>% as.data.frame(check.names=F)

# get the regression data
fortrendline = read.csv('~/Dropbox (Mason Lab)/i4/random_tables_supp/SUPPTABLE_regression_output.csv') 
associated_feat_mg = fortrendline %>% filter(grepl('Increased in/after flight',timetrend)) %>% filter(dset == 'GTDB' | dset == 'GENE-CATALOG' | dset == 'GENBANK-VIRUSES',seqtype == 'METAGENOMICS') %>% mutate(regressionclass = if_else(regclass2 == 'Oral' | regclass2 == 'Nasal',regclass2,'Skin')) %>% select(regressionclass,feature)
associated_feat_mt = fortrendline %>% filter(grepl('Increased in/after flight',timetrend)) %>% filter(dset == 'GTDB' | dset == 'GENE-CATALOG' | dset == 'GENBANK-VIRUSES',seqtype == 'METATRANSCRIPTOMICS')%>% mutate(regressionclass = if_else(regclass2 == 'Oral' | regclass2 == 'Nasal',regclass2,'Skin')) %>% select(regressionclass,feature)

bodylocsskin = fortrendline %>% filter(regclass2!='Oral',regclass2!='Skin',regclass2!='Overall',regclass2!='Nasal') %>% select(regclass2) %>% unlist %>% unname %>% unique 

# load expression data for a single cell type
immunedata = read.csv('~/Dropbox (Mason Lab)/i4/singlecell/Gene_expression/average.B.gene.expression.orig.ident.csv') %>% select(-X) %>% melt %>% mutate(Timepoint_Numeric = strsplit(as.character(variable),'_') %>% map_chr(2) %>% as.numeric,Crew.ID = strsplit(as.character(variable),'_') %>% map_chr(1)) 

immunedata$Timepoint_Numeric[immunedata$Timepoint_Numeric==6]=8
immunedata$Timepoint_Numeric[immunedata$Timepoint_Numeric==5]=7
immunedata$Timepoint_Numeric[immunedata$Timepoint_Numeric==4]=6

rawdata1 =  gtdb_metag %>% select(any_of(unique(associated_feat_mg$feature))) %>% rownames_to_column('SeqID') %>% melt
rawdata2 =  gtdb_metat %>% select(any_of(unique(associated_feat_mt$feature))) %>% rownames_to_column('SeqID') %>% melt
rawdata3 =  viral_metag %>% select(any_of(unique(associated_feat_mg$feature))) %>% rownames_to_column('SeqID') %>% melt
rawdata4 =  viral_metat %>% select(any_of(unique(associated_feat_mt$feature))) %>% rownames_to_column('SeqID') %>% melt
rawdata5 =  gene_metag %>% select(any_of(unique(associated_feat_mg$feature))) %>% rownames_to_column('SeqID') %>% melt
rawdata6 =  gene_metat %>% select(any_of(unique(associated_feat_mt$feature))) %>% rownames_to_column('SeqID') %>% melt

rawdata_mg = bind_rows(rawdata1,rawdata3,rawdata5)
rawdata_mt = bind_rows(rawdata2,rawdata4,rawdata6)

rawdata_meta_mg = inner_join(rawdata_mg,meta,by='SeqID')
rawdata_meta_mt = inner_join(rawdata_mt,meta,by='SeqID')

m = x
data_p = immunedata %>% rename(PATHWAY = gene)
data_m = rawdata_meta_mg
regressionoutput = associated_feat_mg



regress_microbe_immune <- function(m,data_p,data_m,regressionoutput,bodylocsskin){
  output=list()
  m = as.character(m)
  print(m)
  data_p_sub = data_p %>% select(PATHWAY,value,Crew.ID,Timepoint_Numeric) %>% dcast(Crew.ID + Timepoint_Numeric ~ PATHWAY,value.var ='value')
  data_p_sub[is.na(data_p_sub)]=0
  regoutsub = regressionoutput %>% filter(feature == m) %>% select(regressionclass) %>% unlist %>% unname %>% unique
  for(r in regoutsub){
    if(r == 'Skin'){
      data_m_sub = data_m %>% filter(variable == m,location %in% bodylocsskin) %>% group_by(Timepoint_Numeric,Crew.ID)  %>% summarise(value = mean(value))
    }
    if(r != 'Skin'){
      data_m_sub = data_m %>% filter(variable == m,location == r) %>% group_by(Timepoint_Numeric,Crew.ID) %>% group_by(Timepoint_Numeric,Crew.ID)  %>% summarise(value = mean(value))
    }
    for_regression = inner_join(data_m_sub,data_p_sub,by=c('Timepoint_Numeric','Crew.ID'))
    X <- for_regression %>% ungroup %>% select(any_of(unique(data_p$PATHWAY)))
    X = X[, which(colSums(X) != 0)]
    y <- log(for_regression$value + min(for_regression$value[for_regression$value>0]))
    preprocessParams<-preProcess(X, method = c("center", "scale"))
    X <- predict(preprocessParams, X)
    lasso <- 'none'
    try({lasso<-coef(cv.glmnet(as.matrix(X),y)) %>% as.matrix %>% as.data.frame %>% dplyr::rename(!!m := s1) %>% rownames_to_column('features')},silent=TRUE)
    output[[r]] = lasso
  }
  return(output)
}

date()
output = map(unique(rawdata_mg$variable), function(x) regress_microbe_immune(x,immunedata %>% rename(PATHWAY = gene),rawdata_meta_mg,associated_feat_mg,bodylocsskin))
output2 = map(output, function(x) if(typeof(x)!='character'){return(x)})
output2 = output2[which(!sapply(output2, is.null))]
# to anyone reading this, yes, I know this is ugly code, I'm tired and in a hurry.
count = 0
output3 =list()
for(o in output2){
  count = count + 1
  temp = list()
  for(i in names(o)){
    foo = o[[i]] %>% mutate(regressionclass = i,microbialfeaturename = colnames(o[[i]])[2])
    colnames(foo)[2] = 'coefficient'
    temp[[i]] = foo
  }
  temp2= bind_rows(temp)
  output3[[count]] = temp2
}
date()

print('Finished metagenomic modeling...')
write.csv(bind_rows(output3),paste('LASSO_OUT_METAGENOMIC','_','average.B.gene.expression.orig.ident.csv',sep=''))


date()
output = map(unique(rawdata_mt$variable), function(x) regress_microbe_immune(x,immunedata %>% rename(PATHWAY = gene),rawdata_meta_mt,associated_feat_mt,bodylocsskin))
output2 = map(output, function(x) if(typeof(x)!='character'){return(x)})
output2 = output2[which(!sapply(output2, is.null))]
# to anyone reading this, yes, I know this is ugly code, I'm tired and in a hurry.
count = 0
output3 =list()
for(o in output2){
  count = count + 1
  temp = list()
  for(i in names(o)){
    foo = o[[i]] %>% mutate(regressionclass = i,microbialfeaturename = colnames(o[[i]])[2])
    colnames(foo)[2] = 'coefficient'
    temp[[i]] = foo
  }
  temp2= bind_rows(temp)
  output3[[count]] = temp2
}
date()

output3 = bind_rows(output3)

print('Finished metatranscriptomic modeling...')
write.csv(bind_rows(output3),paste('LASSO_OUT_METATRANSCRIPTOMIC','_','average.B.gene.expression.orig.ident.csv',sep=''))










