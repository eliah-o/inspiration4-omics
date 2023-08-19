# correlation between all lmer and non lmer results

library(tidyverse)
library(reshape2)

setwd('~/Dropbox (Mason Lab)/i4/revisions/')

dat = read.delim("outlist_nonlmer.csv",header=F) %>% filter(grepl('_decontam',V1)) %>% unlist %>% unname

overlapdata = list()
corrout = list()

for(filename in dat){
  print(filename)
  site = strsplit(filename,'/') %>% map_chr(4) %>% strsplit('_NON-LMER') %>% map_chr(1) %>% strsplit('regression_output_') %>% map_chr(2)
  seqtype = strsplit(filename,'/') %>% map_chr(4) %>% strsplit('_NON-LMER') %>% map_chr(2) %>% strsplit('_') %>%map_chr(3)
  orgtype = strsplit(filename,'/') %>% map_chr(4) %>% strsplit('_NON-LMER') %>% map_chr(2) %>% strsplit('_') %>%map_chr(2)
  rank = strsplit(filename,'/') %>% map_chr(4) %>% strsplit('_NON-LMER') %>% map_chr(2) %>% strsplit('_') %>%map_chr(4)
  algo = strsplit(filename,'/') %>% map_chr(4) %>% strsplit('_NON-LMER') %>% map_chr(2) %>% strsplit('_') %>%map_chr(5)
  thresholds = strsplit(filename,'/') %>% map_chr(4) %>% strsplit('_NON-LMER') %>% map_chr(2) %>% strsplit('_') %>%map_chr(7) %>% strsplit('\\.') %>% map_chr(1)
  nonlmer = read.delim(filename,sep='\t') 
  lmer = read.delim(gsub('_NON-LMER','',filename),sep='\t') 
  nonlmer = nonlmer  %>% mutate(regressiontype = 'GLM') #%>% rename(term = GLM_term,yvar = GLM_yvar)
  lmer = lmer  %>% mutate(regressiontype = 'LMER') #%>% rename(term = LMER_term,yvar = LMER_yvar)
  lmer = lmer %>% select(all_of(colnames(nonlmer)))
  nonlmer = nonlmer %>% select(all_of(colnames(lmer)))
  regdata_widesig = full_join(lmer %>% filter(BH_adjusted<0.05)%>% select(yvar,term,regressiontype) %>% rename(LMER = regressiontype) ,nonlmer%>% filter(BH_adjusted<0.05)%>% select(yvar,term,regressiontype) %>% rename(GLM = regressiontype) ,by=c('yvar','term'))
  overlapregdata = regdata_widesig %>% mutate(val = paste(yvar,term)) %>% select(-yvar,-term) %>% column_to_rownames("val")
  overlapregdata$LMER[is.na(overlapregdata$LMER)]=0
  overlapregdata$GLM[is.na(overlapregdata$GLM)]=0
  overlapregdata$LMER[overlapregdata$LMER=='LMER']=1
  overlapregdata$GLM[overlapregdata$GLM=='GLM']=1
  overlapregdata$LMER = as.numeric(overlapregdata$LMER)
  overlapregdata$GLM = as.numeric(overlapregdata$GLM)
  overlapdata[[filename]] = overlapregdata %>% mutate(site = site, rank = rank, orgtype = orgtype, algo = algo, thresholds = thresholds, seqtype = seqtype) %>% mutate(category = 'SHARED') %>% mutate(category = if_else(GLM == 0 & LMER == 1,'LMER-SPECIFIC',category))%>% mutate(category = if_else(GLM == 1 & LMER == 0,'GLM-SPECIFIC',category))%>% mutate(category = if_else(GLM == 1 & LMER == 1,'SHARED',category))
  overlapregdata2 = overlapregdata %>% filter(GLM==1 & LMER == 1)
  if(nrow(overlapregdata2)<10){
    next
  }
  corrout[[filename]] = cor.test(overlapregdata2$GLM,overlapregdata2$LMER) %>% tidy %>% mutate(site = site, rank = rank, orgtype = orgtype, algo = algo, thresholds = thresholds, seqtype = seqtype)
}
 
overlapdata = bind_rows(overlapdata)
corrout = bind_rows(corrout)

overlapdata_wide = overlapdata %>% group_by(site,rank,orgtype,algo,thresholds,seqtype) %>% count(category) %>% dcast(site+rank+orgtype+algo+thresholds+seqtype ~ category,value.var = 'n')
overlapdata_wide[is.na(overlapdata_wide)] = 0

write.csv(overlapdata_wide,'~/Dropbox (Mason Lab)/i4/revisions/regression_comparison/glm_vs_lmer_comparison.csv')

overlapdata_wide_sub = overlapdata_wide %>% filter(rank == 's' & algo == 'metaphlan4' | rank == 'g' & algo == 'xtree' & orgtype == 'virus' & thresholds == '01-005' | rank == 's' & algo == 'xtree' & orgtype == 'bacteria' |  rank == 's' & algo == 'phanta' | rank == 'species') %>% filter(site != 'skin')

ggplot(overlapdata_wide_sub %>% filter(seqtype == 'metagenomics')%>% melt(id.vars = c('site','rank','orgtype','algo','thresholds','seqtype')),aes(y = value, x = variable)) + geom_bar(stat = 'identity') + facet_grid(algo + thresholds ~ site,scales='free_y') + theme_minimal() + theme(axis.text.x = element_text(angle = 60,hjust = 1)) + ylab('Feature count') + xlab('') + ggtitle('Metagenomics')
ggsave('~/Dropbox (Mason Lab)/i4/revisions/regression_comparison/metagenomics_regression_comp.pdf',width=4.5,height=9)

ggplot(overlapdata_wide_sub %>% filter(seqtype == 'metatranscriptomics')%>% melt(id.vars = c('site','rank','orgtype','algo','thresholds','seqtype')),aes(y = value, x = variable)) + geom_bar(stat = 'identity') + facet_grid(algo + thresholds ~ site,scales='free_y') + theme_minimal() + theme(axis.text.x = element_text(angle = 60,hjust = 1)) + ylab('Feature count') + xlab('') + ggtitle('Metatranscriptomics')
ggsave('~/Dropbox (Mason Lab)/i4/revisions/regression_comparison/metatranscriptomics_regression_comp.pdf',width=4.5,height=9)

write.csv(corrout,'~/Dropbox (Mason Lab)/i4/revisions/regression_comparison/correlation_glm_lmer.csv')

overlapdata_wide_sub2 = overlapdata_wide_sub %>% mutate(total = `GLM-SPECIFIC` + `LMER-SPECIFIC` + `SHARED`,`GLM-SPECIFIC` = `GLM-SPECIFIC`/total,`SHARED` = `SHARED`/total,`LMER-SPECIFIC` = `LMER-SPECIFIC`/total) %>% select(seqtype,`GLM-SPECIFIC` ,`LMER-SPECIFIC`, `SHARED`) %>% melt

ggplot(overlapdata_wide_sub2,aes(x = variable,y = value))+theme_minimal() +facet_grid(.~seqtype)+ ggbeeswarm::geom_quasirandom() + theme(axis.text.x = element_text(angle=60,hjust =1))












