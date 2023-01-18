#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(lme4)
library(broom)
library(broom.mixed)
library(reshape2)
library(lmerTest)
library(purrr)

dtype=strsplit(args[[1]],',') %>% map_chr(2)
filepath=strsplit(args[[1]],',') %>% map_chr(1)
outname = paste(strsplit(args[[1]],'/') %>% map_chr(4) %>% strsplit(',') %>% map_chr(1),sep='_')
print(outname)
# load metadata
meta = read.csv('i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
meta$SeqID = gsub('SW_','',meta$SeqID)
meta$Timepoint_Recode = factor(meta$Timepoint)
levels(meta$Timepoint_Recode) = c(NA,'PRE-LAUNCH','POST-LAUNCH','PRE-LAUNCH','MID-FLIGHT','MID-FLIGHT','POST-LAUNCH','POST-LAUNCH','PRE-LAUNCH')

meta = meta %>% distinct %>% mutate(Timepoint_Recode2 = if_else(as.character(Timepoint_Recode) == 'MID-FLIGHT',Timepoint,as.character(Timepoint_Recode)))
meta$Timepoint_Recode2 = factor(meta$Timepoint_Recode2,levels = c('PRE-LAUNCH','Flight 1','Flight 2','POST-LAUNCH'))
meta$Timepoint = factor(meta$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec',NA))
meta$Timepoint_Numeric = as.numeric(meta$Timepoint)

# load abundance tables

sanitize_sample_names <- function(data){
  temp = data %>% t %>% as.data.frame %>% rownames_to_column('temp') %>% mutate(namelengths = nchar(temp))
  temp = temp %>% mutate(temp = if_else(namelengths>=3 & str_sub(temp,nchar(temp),nchar(temp))=='D',str_sub(temp,1,nchar(temp)-1),temp))
  return(temp %>% column_to_rownames('temp') %>% select(-namelengths) %>% t %>% data.frame(check.names=F))
}

oac = meta %>% filter(location == 'Open Air Control' | Body.Location == 'Control Swab (0)' | Body.Location == 'Swab Water')
oac = oac %>% select(SeqID) %>% unlist %>% unname 

remove_potential_contamination <- function(data,oac){
  fpg = data %>% select(any_of(oac)) %>% rownames_to_column('microbe') %>% melt 
  todrop = fpg %>% filter(value>quantile(fpg$value,na.rm=T,.75)) %>% select(microbe) %>% unlist %>% unname %>% unique
  data = data[setdiff(rownames(data),todrop),]
  return(data)
}

if(dtype == 'dna'){
  print('Loading WGS data')
  wgs = read.table(filepath,sep='\t',check.names=F,header = T,row.names = 1) 
  abdata = remove_potential_contamination(sanitize_sample_names(wgs),oac)
}

if(dtype == 'rna'){
  print('Loading MTX data')
  mtx = read.table(filepath,sep='\t',check.names=F,header = T,row.names = 1) 
  abdata = remove_potential_contamination(mtx,oac)
}

# run the regressions

metasub = meta %>% filter(Body.Location != "Swab Water",location!='Capsule',Body.Location != 'Open Air Control', Body.Location != "Deltoid - Pre-Biospy", !is.na(Body.Location))

metasub$Timepoint_Recode = factor(metasub$Timepoint_Recode,levels = c('PRE-LAUNCH','MID-FLIGHT','POST-LAUNCH'))
metasub = metasub %>% mutate(Timepoint = if_else(Timepoint == 'Flight 1' | Timepoint == 'Flight 2','Mid-Flight',as.character(Timepoint)))
metasub$Timepoint = factor(metasub$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Mid-Flight','Sept Post-Return','November','21-Dec',NA))
metasub$Timepoint_Recode = factor(metasub$Timepoint_Recode,levels=c('MID-FLIGHT','PRE-LAUNCH','POST-LAUNCH'))
metasub = metasub %>% mutate(isoral = if_else(Body.Location == 'Oral',1,0))
metasub = metasub %>% mutate(isnasal = if_else(Body.Location == 'Nasal',1,0))
metasub = metasub %>% mutate(isskin = if_else(Body.Location != 'Oral' & Body.Location != 'Nasal',1,0))
metasub = metasub %>% mutate(Armpit = if_else(Body.Location == 'Armpit',1,0))
metasub = metasub %>% mutate(web = if_else(Body.Location == 'Toe Web Space',1,0))
metasub = metasub %>% mutate(nape = if_else(Body.Location == 'Nape of Neck',1,0))
metasub = metasub %>% mutate(postauric = if_else(Body.Location == 'Post-Auricular',1,0))
metasub = metasub %>% mutate(fore = if_else(Body.Location == 'Forearm',1,0))
metasub = metasub %>% mutate(bb = if_else(Body.Location == 'Belly Button',1,0))
metasub = metasub %>% mutate(gc = if_else(Body.Location == 'Gluteal Crease',1,0))
metasub = metasub %>% mutate(nasal = if_else(Body.Location == 'Nasal',1,0))
metasub = metasub %>% mutate(Tzone = if_else(Body.Location == 'T-Zone',1,0))
metasub = metasub %>% mutate(Oral = if_else(Body.Location == 'Oral',1,0))

microbesofinterest = rownames(abdata)
minval = min(abdata[abdata>0])

abdata_t = abdata %>% t %>% data.frame(check.names=F) %>% rownames_to_column('SeqID')

abdata_meta = inner_join(abdata_t,metasub,by='SeqID')

regression_output_skin = list()
regression_output_overall = list()
regression_output_skinseparates = list()
regression_output_nasal = list()
regression_output_oral = list()
for(m in microbesofinterest){
  regression_output_overall[[m]] = try(lmer(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode  + (1|Crew.ID)) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
  regression_output_skin[[m]] = try(lmer(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*isskin  + (1|Crew.ID)) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
  regression_output_nasal[[m]] = try(lmer(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*isnasal + (1|Crew.ID)) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
  regression_output_oral[[m]] = try(lmer(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*isoral + (1|Crew.ID)) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
  regression_output_skinseparates[[m]] = try(lmer(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*Armpit +Timepoint_Recode*web + Timepoint_Recode*nape + Timepoint_Recode*postauric + Timepoint_Recode*fore + Timepoint_Recode*bb + Timepoint_Recode*gc + Timepoint_Recode*Tzone + (1|Crew.ID)) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
}

regression_output_overall = bind_rows(regression_output_overall)
regression_output_overall = regression_output_overall %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_overall,paste('regression_output_overall_',outname,'.tsv',sep=''),quote=F,sep='\t')

regression_output_skin = bind_rows(regression_output_skin)
regression_output_skin = regression_output_skin %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_skin,paste('regression_output_skin_',outname,'.tsv',sep=''),quote=F,sep='\t')

regression_output_nasal = bind_rows(regression_output_nasal)
regression_output_nasal = regression_output_nasal %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_nasal,paste('regression_output_nasal_',outname,'.tsv',sep=''),quote=F,sep='\t')

regression_output_oral = bind_rows(regression_output_oral)
regression_output_oral = regression_output_oral %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_oral,paste('regression_output_oral_',outname,'.tsv',sep=''),quote=F,sep='\t')

regression_output_skinseparates = bind_rows(regression_output_skinseparates)
regression_output_skinseparates = regression_output_skinseparates  %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_skinseparates,paste('regression_output_skin_site_by_site_',outname,'.tsv',sep=''),quote=F,sep='\t')



