#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(lme4)
library(broom)
library(broom.mixed)
library(reshape2)
library(lmerTest)

dtype=args[[1]]
org=args[[2]]
taxlevel=args[[3]]
filepath=args[[4]]
algorithm=args[[5]]
cutoffs=args[[6]]
dataframedescr=args[[7]]
outname = paste(org,dtype,taxlevel,algorithm,cutoffs,dataframedescr,sep='_')

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

abdata = readRDS(filepath)

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

regression_output_overall = list()
regression_output_skin = list()
regression_output_skinseparates = list()
regression_output_nasal = list()
regression_output_oral = list()
for(m in microbesofinterest){
  # regression_output_overall[[m]] = try(lmer(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode  + (1|Crew.ID)) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
  regression_output_skin[[m]] = try(lm(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*isskin) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
  regression_output_nasal[[m]] = try(lm(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*isnasal) %>% tidy %>% mutate(yvar = m) %>%filter(term!='(Intercept)'),silent=T)
  regression_output_oral[[m]] = try(lm(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*isoral) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
  regression_output_skinseparates[[m]] = try(lm(data = abdata_meta,log(abdata_meta[,m] + minval) ~ Timepoint_Recode*Armpit +Timepoint_Recode*web + Timepoint_Recode*nape + Timepoint_Recode*postauric + Timepoint_Recode*fore + Timepoint_Recode*bb + Timepoint_Recode*gc + Timepoint_Recode*Tzone) %>% tidy %>% mutate(yvar = m) %>% filter(term!='(Intercept)'),silent=T)
}

regression_output_skin = bind_rows(regression_output_skin)
regression_output_skin = regression_output_skin %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_skin,paste('regression_output_skin_NON-LMER_',outname,'.tsv',sep=''),quote=F,sep='\t')

regression_output_nasal = bind_rows(regression_output_nasal)
regression_output_nasal = regression_output_nasal %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_nasal,paste('regression_output_nasal_NON-LMER_',outname,'.tsv',sep=''),quote=F,sep='\t')

regression_output_oral = bind_rows(regression_output_oral)
regression_output_oral = regression_output_oral %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_oral,paste('regression_output_oral_NON-LMER_',outname,'.tsv',sep=''),quote=F,sep='\t')

regression_output_skinseparates = bind_rows(regression_output_skinseparates)
regression_output_skinseparates = regression_output_skinseparates  %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni'))
write.table(regression_output_skinseparates,paste('regression_output_skin_site_by_site_NON-LMER_',outname,'.tsv',sep=''),quote=F,sep='\t')



