# generate simple summary plots of taxa per sample

library(tidyverse)
library(ggplot2)
library(reshape2)
library(cowplot)
library(tidytext)

####
args = commandArgs(trailingOnly=TRUE)

filepath = '~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/gtdb/bacteria_xtree_005-0025_g___metagenomics_decontam.rds'

if(grepl('metaphlan',filepath)){
  taxsplitchar  = "\\|"
  if(grepl('_p__',filepath)){
    num=2
  }
  if(grepl('_c__',filepath)){
    num=3
  }
  if(grepl('_o__',filepath)){
    num=4
  }
  if(grepl('_f__',filepath)){
    num=5
  }
  if(grepl('_g__',filepath)){
    num=6
  }
  if(grepl('_s__',filepath)){
    num=7
  }
  if(grepl('_t__',filepath)){
    num=8
  }
}

if(grepl('xtree',filepath)){
  taxsplitchar  = ";"
  if(grepl('_p__',filepath)){
    num=2
  }
  if(grepl('_c__',filepath)){
    num=3
  }
  if(grepl('_o__',filepath)){
    num=4
  }
  if(grepl('_f__',filepath)){
    num=5
  }
  if(grepl('_g__',filepath)){
    num=6
  }
  if(grepl('_s__',filepath)){
    num=7
  }
}

if(grepl('MAG',filepath)){
  taxsplitchar  = 'none'
}

if(grepl('kraken',filepath)){
  taxsplitchar = 'none'
}

if(grepl('genbank',filepath)){
  taxsplitchar  = ";"
  if(grepl('_p__',filepath)){
    num=1
  }
  if(grepl('_c__',filepath)){
    num=2
  }
  if(grepl('_o__',filepath)){
    num=3
  }
  if(grepl('_f__',filepath)){
    num=4
  }
  if(grepl('_g__',filepath)){
    num=5
  }
  if(grepl('_s__',filepath)){
    num=6
  }
}

outname = gsub('.rds','',filepath)

# load metadata
meta = read.csv('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
meta$SeqID = gsub('SW_','',meta$SeqID)
meta$Timepoint_Recode = factor(meta$Timepoint)
levels(meta$Timepoint_Recode) = c(NA,'PRE-LAUNCH','POST-LAUNCH','PRE-LAUNCH','MID-FLIGHT','MID-FLIGHT','POST-LAUNCH','POST-LAUNCH','PRE-LAUNCH')

meta = meta %>% distinct %>% mutate(Timepoint_Recode2 = if_else(as.character(Timepoint_Recode) == 'MID-FLIGHT',Timepoint,as.character(Timepoint_Recode)))
meta$Timepoint_Recode2 = factor(meta$Timepoint_Recode2,levels = c('PRE-LAUNCH','Flight 1','Flight 2','POST-LAUNCH'))
meta$Timepoint = factor(meta$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec',NA))
meta$Timepoint_Numeric = as.numeric(meta$Timepoint)

# load abundance tables
abdata = readRDS(filepath)

# merge
metasub = meta %>% filter(Body.Location != "Swab Water",Body.Location != 'Open Air Control', Body.Location != "Deltoid - Pre-Biospy", !is.na(Body.Location))
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
abdata_t = abdata_t %>% melt

abdata_meta = left_join(abdata_t,metasub,by='SeqID')
if(taxsplitchar != 'none'){
  abdata_meta$variable = strsplit(as.character(abdata_meta$variable),taxsplitchar) %>% map_chr(num)
}
abdata_meta$location = factor(abdata_meta$location,levels = c('Oral','Nasal','Armpit','Belly Button','Forearm','Gluteal Crease','Nape of Neck','Post-Auricular','T-Zone','Toe Web Space','Stool','Capsule'))

### BY BODY SITE
#barplotdat = abdata_meta  %>% filter(!is.na(location),!is.na(Timepoint_Recode)) %>% dplyr::group_by(location,variable,Timepoint_Recode) %>% summarise(m = mean(value),sd = sd(value))%>% ungroup %>% dplyr::group_by(location,Timepoint_Recode) %>% arrange(desc(m)) %>% slice_max(m,n=10) 

#all = ggplot(barplotdat %>% ungroup,aes(x = reorder_within(variable,m,location), y = m)) + theme_cowplot() + geom_bar(stat='identity') + facet_wrap(Timepoint_Recode ~ location,scales='free',ncol= 3 ,nrow=10) + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") + scale_x_reordered() + ylab('')

barplotdat_pre = abdata_meta %>% filter(!is.na(location),!is.na(Timepoint_Recode),Timepoint_Recode == 'PRE-LAUNCH') %>% group_by(location,variable) %>% summarise(m = mean(value),sd = sd(value)) %>% arrange(desc(m)) %>% slice_max(m,n=10) %>% ungroup %>% filter(m>0)
barplotdat_mid = abdata_meta %>% filter(!is.na(location),!is.na(Timepoint_Recode),Timepoint_Recode == 'MID-FLIGHT') %>% group_by(location,variable) %>% summarise(m = mean(value),sd = sd(value)) %>% arrange(desc(m)) %>% slice_max(m,n=10) %>% ungroup %>% filter(m>0)
barplotdat_post = abdata_meta %>% filter(!is.na(location),!is.na(Timepoint_Recode),Timepoint_Recode == 'POST-LAUNCH') %>% group_by(location,variable) %>% summarise(m = mean(value),sd = sd(value)) %>% arrange(desc(m)) %>% slice_max(m,n=10) %>% ungroup %>% filter(m>0)

pre = ggplot(barplotdat_pre %>% ungroup,aes(x = reorder_within(variable,m,location), y = m)) + theme_cowplot() + geom_bar(stat='identity') + facet_wrap(. ~ location,scales='free_x',nrow=2,ncol=6) + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('')+ scale_x_reordered() + ylab('') #+ geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") + scale_x_reordered() + ylab('')
#ggsave(paste(outname,'_pre.pdf',sep=''),width=14,height=8)
mid = ggplot(barplotdat_mid %>% ungroup,aes(x = reorder_within(variable,m,location), y = m)) + theme_cowplot() + geom_bar(stat='identity') + facet_wrap(. ~ location,scales='free_x',nrow=2,ncol=6) + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('')+ scale_x_reordered() + ylab('') #+ geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") + scale_x_reordered() + ylab('')
#ggsave(paste(outname,'_mid.pdf',sep=''),width=14,height=8)
post =  ggplot(barplotdat_post %>% ungroup,aes(x = reorder_within(variable,m,location), y = m)) + theme_cowplot() + geom_bar(stat='identity') + facet_wrap(. ~ location,scales='free_x',nrow=2,ncol=6) + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('') + scale_x_reordered() + ylab('')#+ geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") 
#ggsave(paste(outname,'_post.pdf',sep=''),width=14,height=8)


#### OVERALL

barplotdat = abdata_meta %>% filter(!is.na(location),!is.na(Timepoint_Recode)) %>% group_by(location,variable) %>% summarise(m = mean(value),sd = sd(value)) %>% arrange(desc(m)) %>% slice_max(m,n=10) %>% ungroup %>% filter(m>0)

overall = ggplot(barplotdat %>% ungroup,aes(x = reorder_within(variable,m,location), y = m)) + theme_cowplot() + geom_bar(stat='identity') + facet_wrap(. ~ location,scales='free_x',nrow=2,ncol=6) + theme(axis.text.x = element_text(angle=45,hjust=1,size=6)) + xlab('')+ scale_x_reordered() + ylab('') #+ geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") + scale_x_reordered() + ylab('')
ggsave(paste(outname,'_overall.pdf',sep=''),width=8.5,height=5.5)

### BODY SITE GROUPED, SPECIFIC SITE AGNOSTIC

abdata_meta = abdata_meta %>% mutate(skin_oral_nasal = if_else(isnasal == 1,'NASAL',if_else(isoral == 1,'ORAL',if_else(isskin == 1,'SKIN',if_else(location == 'Capsule','CAPSULE',if_else(location == 'Stool','Stool',NA))))))

abdata_meta2 = abdata_meta %>% group_by(skin_oral_nasal,variable) %>% summarise(m = mean(value),sd = sd(value)) %>% arrange(desc(m)) %>% slice_max(m,n=5) %>% ungroup %>% filter(m>0) %>% filter(!is.na(skin_oral_nasal))

overall = ggplot(abdata_meta2 %>% ungroup,aes(x = reorder_within(variable,m,skin_oral_nasal), y = m)) + theme_cowplot() + geom_bar(stat='identity') + facet_wrap(. ~ skin_oral_nasal,scales='free_x',nrow=2,ncol=6) + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('')+ scale_x_reordered() + ylab('') #+ geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") + scale_x_reordered() + ylab('')
ggsave(paste(outname,'_overall_acrosssites.pdf',sep=''),width=6,height=3)


#CAPSULE -- changed filenames manually as needed for plots

barplotdat_pre = abdata_meta %>% filter(location == 'Capsule',!is.na(location),!is.na(Timepoint_Recode),Timepoint_Recode == 'PRE-LAUNCH') %>% group_by(location,variable) %>% summarise(m = mean(value),sd = sd(value)) %>% arrange(desc(m)) %>% slice_max(m,n=25) %>% ungroup %>% filter(m>0)
barplotdat_mid = abdata_meta %>% filter(location == 'Capsule',!is.na(location),!is.na(Timepoint_Recode),Timepoint_Recode == 'MID-FLIGHT') %>% group_by(location,variable) %>% summarise(m = mean(value),sd = sd(value)) %>% arrange(desc(m)) %>% slice_max(m,n=25) %>% ungroup %>% filter(m>0)

pre = ggplot(barplotdat_pre %>% ungroup,aes(x = reorder_within(variable,m,location), y = m)) + theme_cowplot() + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('')+ scale_x_reordered() + ylab('') +ggtitle('Capsule -- Pre-flight') + ylab('Relative abundance')#+ geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") + scale_x_reordered() + ylab('')
ggsave('~/Dropbox (Mason Lab)/i4/revisions/new_revisions_plots/capsule_pre_metag.pdf',width=8,height=8)
#ggsave(paste(outname,'_pre.pdf',sep=''),width=14,height=8)
mid = ggplot(barplotdat_mid %>% ungroup,aes(x = reorder_within(variable,m,location), y = m)) + theme_cowplot() + geom_bar(stat='identity') + theme(axis.text.x = element_text(angle=45,hjust=1)) + xlab('')+ scale_x_reordered() + ylab('') +ggtitle('Capsule -- Mid-flight') + ylab('Relative abundance')#+ geom_errorbar(aes(ymin=m-sd, ymax=m+sd), width=0.2, size=1, color="black") + scale_x_reordered() + ylab('')
ggsave('~/Dropbox (Mason Lab)/i4/revisions/new_revisions_plots/capsule_mid_metag.pdf',width=8,height=8)

