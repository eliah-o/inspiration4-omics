
library(tidyverse)
### LOAD IN REGRESSION OUTPUT AND MERGE
regdata = list()

### BACTERIAL ASSEMBLY
regdata[[1]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral_bacteria_metagenomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
regdata[[2]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral_bacteria_metatranscriptomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
regdata[[3]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal_bacteria_metagenomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'ASSEMBLED-BACTERIA')
regdata[[4]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal_bacteria_metatranscriptomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
regdata[[5]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_bacteria_metagenomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
regdata[[6]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_bacteria_metatranscriptomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin')  %>% mutate(dset = 'ASSEMBLED-BACTERIA')
regdata[[7]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site_bacteria_metagenomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components')  %>% mutate(dset = 'ASSEMBLED-BACTERIA')
regdata[[8]] = read.table("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site_bacteria_metatranscriptomics_p_MAG_nodecontam_05-0025.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'ASSEMBLED-BACTERIA')

### VIRAL ASSEMBLED

regdata[[9]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral_virus_metagenomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
regdata[[10]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral_virus_metatranscriptomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
regdata[[11]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal_virus_metagenomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'ASSEMBLED-VIRUSES')
regdata[[12]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal_virus_metatranscriptomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
regdata[[13]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_virus_metagenomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
regdata[[14]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_virus_metatranscriptomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
regdata[[15]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site_virus_metagenomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
regdata[[16]] = read.csv("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site_virus_metatranscriptomics_p_MAG_nodecontam_01-005.tsv",sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'ASSEMBLED-VIRUSES')

regression_output = bind_rows(regdata)

regression_output$term = gsub('Timepoint_Recode','',regression_output$term)
regression_output$term = gsub(':','',regression_output$term)
regression_output$term = gsub(':Body.Location','',regression_output$term)
regression_output$term = gsub('PRE-LAUNCH','PRE-LAUNCH --- ',regression_output$term)
regression_output$term = gsub('POST-LAUNCH','POST-LAUNCH --- ',regression_output$term)
regression_output$term = gsub('isskin','Skin',regression_output$term)
regression_output$term = gsub('isoral','Oral',regression_output$term)
regression_output$term = gsub('nape','Nape of Neck',regression_output$term)
regression_output$term = gsub('postauric','Post-Auricular',regression_output$term)
regression_output$term = gsub('bb','Belly Button',regression_output$term)
regression_output$term = gsub('gc','Gluteal Crease',regression_output$term)
regression_output$term = gsub('Tzone','T-Zone',regression_output$term)
regression_output$term = gsub('web','Toe Web Space',regression_output$term)
regression_output$term = gsub('fore','Forearm',regression_output$term)
regression_output$term = gsub(' ---  --- ',' --- ',regression_output$term)
regression_output$term = gsub('isnasal','Nasal',regression_output$term)
regression_output = regression_output %>% mutate(term = if_else(term == 'PRE-LAUNCH --- ','PRE-LAUNCH',term,))
regression_output = regression_output %>% mutate(term = if_else(term == 'POST-LAUNCH --- ','POST-LAUNCH',term,))
regression_output = regression_output %>% mutate(launch = if_else(grepl('PRE',term),'PRE-LAUNCH','POST-LAUNCH')) %>% mutate(direction = if_else(estimate>0,'positive','negative')) %>% filter(term != 'POST-LAUNCH --- OVERALL')%>% filter(term != 'PRE-LAUNCH --- OVERALL')
#temp0 =regression_output %>% filter(grepl('---',term)) %>% mutate(regclass2 = strsplit(as.character(term),' --- ') %>% map_chr(2))
#temp1 = regression_output %>% filter(!grepl('---',term),grepl('LAUNCH',term)) %>% mutate(regclass2 = 'Main effects (pre/post launch)')
#temp2 = regression_output %>% filter(!grepl('---',term),!grepl('LAUNCH',term)) %>% mutate(regclass2 = paste0('Main effects (',term,')'))
#regression_output = bind_rows(temp0,temp1,temp2) %>% group_by(regclass2) %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni')) %>% ungroup

#### LOAD IN BACTERIAL ANNOTATIONS
bacanno = read.delim('~/Dropbox (Mason Lab)/i4/i4_data_packet/gtdbtk_classification_bac120.tsv',sep='\t') %>% select(contig_identifier,classification) %>% rename(BACTERIAL_MAG_CLASSIFICATION = classification)

regression_output = left_join(regression_output,bacanno,by = c('yvar'='contig_identifier'))

saveRDS(regression_output,paste0('~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/parsed/regression_data_timetrends_MAGS.rds'))

regression_output = regression_output  %>% filter(grepl('---',term)) 
regression_output$term = factor(regression_output$term,levels = c("PRE-LAUNCH --- Armpit","PRE-LAUNCH --- Belly Button","PRE-LAUNCH --- Forearm", "PRE-LAUNCH --- Gluteal Crease","PRE-LAUNCH --- Nape of Neck","PRE-LAUNCH --- Nasal","PRE-LAUNCH --- Post-Auricular","PRE-LAUNCH --- T-Zone","PRE-LAUNCH --- Toe Web Space","PRE-LAUNCH --- Skin","PRE-LAUNCH --- Oral","POST-LAUNCH --- Armpit","POST-LAUNCH --- Belly Button", "POST-LAUNCH --- Forearm","POST-LAUNCH --- Gluteal Crease" ,"POST-LAUNCH --- Nape of Neck", "POST-LAUNCH --- Nasal", "POST-LAUNCH --- Post-Auricular","POST-LAUNCH --- T-Zone", "POST-LAUNCH --- Toe Web Space","POST-LAUNCH --- Skin","POST-LAUNCH --- Oral"))
regression_output = regression_output %>% mutate(launch = if_else(grepl('PRE',term),'PRE-LAUNCH','POST-LAUNCH')) %>% mutate(direction = if_else(estimate>0,'positive','negative')) %>% filter(term != 'POST-LAUNCH --- OVERALL')%>% filter(term != 'PRE-LAUNCH --- OVERALL')

a =regression_output %>% filter(launch == 'PRE-LAUNCH') %>% mutate(regclass2 = strsplit(as.character(term),' --- ') %>% map_chr(2))  %>% dplyr::rename(PRE=launch,adj_pre = BH_adjusted,est_pre = estimate) %>% select(yvar,dset,seqtype,regclass2,PRE,adj_pre,est_pre,BACTERIAL_MAG_CLASSIFICATION)
#a2 =regression_output %>% filter(launch == 'PRE-LAUNCH')%>% filter(term == 'PRE-LAUNCH --- OVERALL') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(PRE=launch,adj_pre = BH_adjusted,est_pre = estimate) %>% select(yvar,dset,seqtype,regclass2,PRE,adj_pre,est_pre)
#a3 =regression_output %>% filter(regressionclass == 'Overall') %>% filter(launch == 'PRE-LAUNCH')%>% filter(term == 'PRE-LAUNCH') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(PRE=launch,adj_pre = BH_adjusted,est_pre = estimate) %>% select(yvar,dset,seqtype,regclass2,PRE,adj_pre,est_pre)

#a = bind_rows(a,a3)

b =regression_output %>% filter(launch == 'POST-LAUNCH') %>% filter(grepl('---',term)) %>% mutate(regclass2 = strsplit(as.character(term),' --- ') %>% map_chr(2))%>% dplyr::rename(POST=launch,adj_post = BH_adjusted,est_post = estimate) %>% select(yvar,dset,seqtype,regclass2,POST,adj_post,est_post,BACTERIAL_MAG_CLASSIFICATION)
#b2 =regression_output %>% filter(launch == 'POST-LAUNCH') %>% filter(term == 'POST-LAUNCH --- OVERALL') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(POST=launch,adj_post = BH_adjusted,est_post = estimate) %>% select(yvar,dset,seqtype,regclass2,POST,adj_post,est_post)
#b3 =regression_output %>% filter(regressionclass == 'Overall')%>% filter(launch == 'POST-LAUNCH') %>% filter(term == 'POST-LAUNCH') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(POST=launch,adj_post = BH_adjusted,est_post = estimate) %>% select(yvar,dset,seqtype,regclass2,POST,adj_post,est_post)
#b = bind_rows(b,b3)

c = inner_join(a,b)

c = c %>% filter(adj_pre<0.05 |  adj_post<0.05)

### MODERATE CONSERVATISM
#c = c %>% mutate(timetrend = 'NONE')%>% mutate(timetrend = if_else(est_pre<0 & est_post>0 & adj_pre <0.05 |est_pre<0 & est_post>0 & adj_post <0.05,'Increased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post<0 & adj_pre <0.05 | est_pre>0 & est_post<0 & adj_post <0.05,'Decreased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 | est_pre<0 & est_post<0 & adj_post <0.05,'Increased in flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 | est_pre>0 & est_post>0  & adj_post <0.05,'Decreased in flight',timetrend)) 

### HIGH CONSERVATISM
#c = c %>% mutate(timetrend = 'NONE')%>% mutate(timetrend = if_else(est_pre<0 & est_post>0 & adj_pre <0.05 |est_pre<0 & est_post>0 & adj_post <0.05,'Increased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post<0 & adj_pre <0.05 | est_pre>0 & est_post<0 & adj_post <0.05,'Decreased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 & adj_post <0.05,'Increased in flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 & adj_post <0.05,'Decreased in flight',timetrend)) 

### BOTH
c = c %>% mutate(timetrendrank = 'NONE')%>% mutate(timetrendrank = if_else(est_pre<0 & est_post>0 & adj_pre <0.05 |est_pre<0 & est_post>0 & adj_post <0.05,'Persistent increase',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre>0 & est_post<0 & adj_pre <0.05 | est_pre>0 & est_post<0 & adj_post <0.05,'Persistent decrease',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 | est_pre<0 & est_post<0 & adj_post <0.05,'Transient increase -- LP',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 | est_pre>0 & est_post>0  & adj_post <0.05,'Transient decrease -- LP',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 & adj_post <0.05,'Transient increase',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 & adj_post <0.05,'Transient decrease',timetrendrank)) 
c = c %>% mutate(timetrend = gsub(' -- LP','',timetrendrank))

fortrendline = c %>% filter(regclass2 != 'OVERALL') %>% mutate(mergeid = paste(yvar,dset,seqtype,regclass2))

write.csv(fortrendline,paste0('~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/parsed/regression_data_timetrends_MAGS.csv'))

