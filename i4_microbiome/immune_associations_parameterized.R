#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(glmnet)
library(furrr)
library(caret)
library(progressr)

args = commandArgs(trailingOnly=TRUE)

# load metadata
meta = read.csv('i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
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

# bac GTDB metag/metat
xtreebacmetag = readRDS('revisions/gtdb/bacteria_xtree_05-0025_s___metagenomics_decontam.rds')
xtreebacmetat = readRDS('revisions/gtdb/bacteria_xtree_05-0025_s___metatranscriptomics_decontam.rds')

# bac metaphlan metag/metat
metaphbacmetag = readRDS('revisions/metaphlan4/bacteria_metaphlan4_default_s__metagenomics_decontam.rds')
metaphbacmetat = readRDS('revisions/metaphlan4/bacteria_metaphlan4_default_s__metatranscriptomics_decontam.rds')

# genbank viral metag/metat
xtreegenmetag = readRDS('revisions/genbank/virus_xtree_01-005_g___metagenomics_decontam.rds') 
xtreegenmetat = readRDS('revisions/genbank/virus_xtree_01-005_g___metatranscriptomics_decontam.rds')

# Gene metag/metat (load significant genes only)
gene_metag = readRDS('gene_catalog90_filtered_relab_metagenomics_subset.rds') 
gene_metat = readRDS('gene_catalog90_filtered_relab_metatranscriptomics_subset.rds')

# get the regression data
fortrendline = read.csv('regression_data_timetrends_species_decontam.csv') 
associated_feat_mg = fortrendline %>% filter(grepl('Increased in/after flight',timetrend)) %>% filter(dset == 'GTDB' | dset == 'GENE-CATALOG' | dset == 'METAPHLAN4' ,seqtype == 'METAGENOMICS') %>% mutate(regressionclass = if_else(regclass2 == 'Oral' | regclass2 == 'Nasal',regclass2,'Skin')) %>% select(regressionclass,yvar)
associated_feat_mt = fortrendline %>% filter(grepl('Increased in/after flight',timetrend)) %>% filter(dset == 'GTDB' | dset == 'GENE-CATALOG' | dset == 'METAPHLAN4' ,seqtype == 'METATRANSCRIPTOMICS')%>% mutate(regressionclass = if_else(regclass2 == 'Oral' | regclass2 == 'Nasal',regclass2,'Skin')) %>% select(regressionclass,yvar)

# get the regression data genbank virus
fortrendline2 = read.csv('regression_data_timetrends_genus_decontam.csv') 
associated_feat_mg2 = fortrendline2 %>% filter(grepl('Decreased in/after flight',timetrend)) %>% filter(dset == 'GENBANK-VIRUSES' ,seqtype == 'METAGENOMICS') %>% mutate(regressionclass = if_else(regclass2 == 'Oral' | regclass2 == 'Nasal',regclass2,'Skin')) %>% select(regressionclass,yvar)
associated_feat_mt2 = fortrendline2 %>% filter(grepl('Decreased in/after flight',timetrend)) %>% filter(dset == 'GENBANK-VIRUSES' ,seqtype == 'METATRANSCRIPTOMICS')%>% mutate(regressionclass = if_else(regclass2 == 'Oral' | regclass2 == 'Nasal',regclass2,'Skin')) %>% select(regressionclass,yvar)

associated_feat_mg = bind_rows(associated_feat_mg,associated_feat_mg2)
associated_feat_mt = bind_rows(associated_feat_mt,associated_feat_mt2)

bodylocsskin = fortrendline %>% filter(regclass2!='Oral',regclass2!='Skin',regclass2!='Overall',regclass2!='Nasal') %>% select(regclass2) %>% unlist %>% unname %>% unique 

# load expression data for a single cell type
#immunedata = read.csv(args[[1]]) %>% select(-X) %>% melt %>% mutate(Timepoint_Numeric = strsplit(as.character(variable),'_') %>% map_chr(2) %>% as.numeric,Crew.ID = strsplit(as.character(variable),'_') %>% map_chr(1)) 
immunedata = read.csv(args[[1]]) %>% select(-X) %>% melt %>% mutate(Timepoint_Numeric = strsplit(as.character(variable),'_') %>% map_chr(2) %>% as.numeric,Crew.ID = strsplit(as.character(variable),'_') %>% map_chr(1)) 
immunedata$Timepoint_Numeric[immunedata$Timepoint_Numeric==6]=8
immunedata$Timepoint_Numeric[immunedata$Timepoint_Numeric==5]=7
immunedata$Timepoint_Numeric[immunedata$Timepoint_Numeric==4]=6

rawdata1 =  xtreebacmetag %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mg$yvar))) %>% rownames_to_column('SeqID') %>% melt
rawdata2 =  xtreebacmetat %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mt$yvar))) %>% rownames_to_column('SeqID') %>% melt
rawdata3 =  xtreegenmetag %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mg$yvar))) %>% rownames_to_column('SeqID') %>% melt
rawdata4 =  xtreegenmetat %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mg$yvar))) %>% select(any_of(unique(associated_feat_mt$feature))) %>% rownames_to_column('SeqID') %>% melt
rawdata5 =  gene_metag %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mg$yvar))) %>% rownames_to_column('SeqID') %>% melt
rawdata6 =  gene_metat %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mt$yvar))) %>% rownames_to_column('SeqID') %>% melt
rawdata7 =  metaphbacmetag %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mg$yvar))) %>% rownames_to_column('SeqID') %>% melt
rawdata8 =  metaphbacmetat %>% t %>% as.data.frame(check.names=F)%>% select(any_of(unique(associated_feat_mt$yvar))) %>% rownames_to_column('SeqID') %>% melt

rawdata_mg = bind_rows(rawdata1,rawdata3,rawdata5,rawdata7)
rawdata_mt = bind_rows(rawdata2,rawdata4,rawdata6,rawdata8)

rawdata_meta_mg = inner_join(rawdata_mg,meta,by='SeqID')
rawdata_meta_mt = inner_join(rawdata_mt,meta,by='SeqID')

regress_microbe_immune <- function(m,data_p,data_m,regressionoutput,bodylocsskin){
  output=list()
  m = as.character(m)
  print(m)
  data_p_sub = data_p %>% select(PATHWAY,value,Crew.ID,Timepoint_Numeric) %>% dcast(Crew.ID + Timepoint_Numeric ~ PATHWAY,value.var ='value')
  data_p_sub[is.na(data_p_sub)]=0
  regoutsub = regressionoutput %>% filter(yvar == m) %>% select(regressionclass) %>% unlist %>% unname %>% unique
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
    try({lasso<-coef(cv.glmnet(as.matrix(X),y,grouped=F)) %>% as.matrix %>% as.data.frame %>% dplyr::rename(!!m := s1) %>% rownames_to_column('features')},silent=TRUE)
    output[[r]] = lasso
  }
  return(output)
}

plan(multisession, workers = 24)

date()
output2 = future_map(unique(rawdata_mg$variable), function(x) regress_microbe_immune(x,immunedata %>% rename(PATHWAY = gene),rawdata_meta_mg,associated_feat_mg,bodylocsskin), .options = furrr_options(seed = TRUE))
count = 0
output3 =list()
for(o in output2){
  count = count + 1
  temp = list()
  for(i in names(o)){
    if(typeof(o[[i]])=='character'){
      next
    }
    foo = o[[i]] %>% mutate(regressionclass = i,microbialfeaturename = colnames(o[[i]])[2])
    colnames(foo)[2] = 'coefficient'
    temp[[i]] = foo
  }
  temp2= bind_rows(temp)
  output3[[count]] = temp2
}
date()

print('Finished metagenomic modeling...')
write.csv(bind_rows(output3)%>% filter(coefficient!=0),paste('LASSO_OUT_METAGENOMIC','_',args[[1]],sep=''))


date()
output2 = future_map(unique(rawdata_mt$variable), function(x) regress_microbe_immune(x,immunedata %>% rename(PATHWAY = gene),rawdata_meta_mt,associated_feat_mt,bodylocsskin), .options = furrr_options(seed = TRUE))
count = 0
output3 =list()
for(o in output2){
  count = count + 1
  temp = list()
  for(i in names(o)){
    if(typeof(o[[i]])=='character'){
      next
    }
    foo = o[[i]] %>% mutate(regressionclass = i,microbialfeaturename = colnames(o[[i]])[2])
    colnames(foo)[2] = 'coefficient'
    temp[[i]] = foo
  }
  temp2= bind_rows(temp)
  output3[[count]] = temp2
}
date()


print('Finished metatranscriptomic modeling...')
write.csv(bind_rows(output3) %>% filter(coefficient!=0),paste('LASSO_OUT_METATRANSCRIPTOMIC','_',args[[1]],sep=''))










