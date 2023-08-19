#### REVISIONS ANALYSIS

library(tidyverse)


# load metadata
meta = read.csv('~/Dropbox (Mason Lab)/i4/i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
meta$SeqID = gsub('SW_','',meta$SeqID)
meta$Timepoint_Recode = factor(meta$Timepoint)
levels(meta$Timepoint_Recode) = c(NA,'PRE-LAUNCH','POST-LAUNCH','PRE-LAUNCH','MID-FLIGHT','MID-FLIGHT','POST-LAUNCH','POST-LAUNCH','PRE-LAUNCH')

meta = meta %>% distinct %>% mutate(Timepoint_Recode2 = if_else(as.character(Timepoint_Recode) == 'MID-FLIGHT',Timepoint,as.character(Timepoint_Recode)))
meta$Timepoint_Recode2 = factor(meta$Timepoint_Recode2,levels = c('PRE-LAUNCH','Flight 1','Flight 2','POST-LAUNCH'))
meta$Timepoint = factor(meta$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec',NA))
meta$Timepoint_Numeric = as.numeric(meta$Timepoint)

metasub = meta %>% filter(Body.Location != "Swab Water",location!='Capsule',Body.Location != 'Open Air Control', Body.Location != "Deltoid - Pre-Biospy", !is.na(Body.Location))

metasub$Timepoint_Recode = factor(metasub$Timepoint_Recode,levels = c('PRE-LAUNCH','MID-FLIGHT','POST-LAUNCH'))
metasub = metasub %>% mutate(Timepoint = if_else(Timepoint == 'Flight 1' | Timepoint == 'Flight 2','Mid-Flight',as.character(Timepoint)))
metasub$Timepoint = factor(metasub$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Mid-Flight','Sept Post-Return','November','21-Dec',NA))
metasub$Timepoint_Recode = factor(metasub$Timepoint_Recode,levels=c('MID-FLIGHT','PRE-LAUNCH','POST-LAUNCH'))
metasub = metasub %>% mutate(isoral = if_else(Body.Location == 'Oral',1,0))
metasub = metasub %>% mutate(isnasal = if_else(Body.Location == 'Nasal',1,0))
metasub = metasub %>% mutate(isstool = if_else(Body.Location == 'Stool',1,0))
metasub = metasub %>% mutate(isskin = if_else(Body.Location != 'Oral' & Body.Location != 'Nasal' & Body.Location != 'Stool',1,0))
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


sanitize_sample_names1 <- function(data){
  temp = data %>% mutate(namelengths = nchar(SeqID))
  temp = temp %>% mutate(SeqID = if_else(namelengths>=3 & str_sub(SeqID,nchar(SeqID),nchar(SeqID))=='D',str_sub(SeqID,1,nchar(SeqID)-1),SeqID))
  return(temp %>% select(-namelengths))
}

sanitize_sample_names <- function(data){
  temp = data %>% t %>% as.data.frame %>% rownames_to_column('temp') %>% mutate(namelengths = nchar(temp))
  temp = temp %>% mutate(temp = if_else(namelengths>=3 & str_sub(temp,nchar(temp),nchar(temp))=='D',str_sub(temp,1,nchar(temp)-1),temp))
  return(temp %>% column_to_rownames('temp') %>% select(-namelengths) %>% t %>% data.frame(check.names=F))
}

setwd('~/Dropbox (Mason Lab)/i4/revisions/human_read_alignment/')

t2t = read.delim('t2t_alignment_data',sep=',',header=F) %>% distinct %>% mutate(PERCENTAGE_T2T = (.01*as.numeric(gsub('%','',V3))),totalreads  =round(V2/PERCENTAGE_T2T)) %>% select(-V3) %>% mutate(PERCENTAGE_T2T = 100*(1 - PERCENTAGE_T2T))
colnames(t2t) = c('sampleid','COUNTS_T2T','PERCENTAGE_T2T','totalreads')
t2t$COUNTS_T2T = t2t$totalreads - t2t$COUNTS_T2T

hg38 = read.delim('hg38_alignment_data',sep=',',header=F) %>% distinct %>% mutate(PERCENTAGE_HG38 = (.01*as.numeric(gsub('%','',V3))),totalreads  =round(V2/PERCENTAGE_HG38)) %>% select(-V3) %>% mutate(PERCENTAGE_HG38 = 100*(1 - PERCENTAGE_HG38))
colnames(hg38) = c('sampleid','COUNTS_HG38','PERCENTAGE_HG38','totalreads')
hg38$COUNTS_HG38 = hg38$totalreads - hg38$COUNTS_HG38
hg38 = hg38 %>% select(-totalreads)

merged = full_join(t2t,hg38,by='sampleid') %>% filter(!is.na(PERCENTAGE_HG38) & !is.na(PERCENTAGE_T2T))%>% melt(id.vars = c('sampleid','totalreads')) %>% mutate(reference = strsplit(as.character(variable),'_') %>% map_chr(2)) %>% mutate(type = strsplit(as.character(variable),'_') %>% map_chr(1)) 

merged = merged %>% mutate(SeqID = strsplit(sampleid,'_') %>% map_chr(1)) %>% sanitize_sample_names1()

merged_meta = inner_join(merged,metasub %>% select(SeqID,location))

ggplot(merged_meta %>% filter(reference == 'HG38'),aes(x = reorder(sampleid,value),fill=location,y = value)) + theme(axis.text.x = element_blank()) + geom_bar(stat='identity') + facet_grid(type  ~ .,scales = 'free_y') + scale_fill_brewer(palette = 'Set3') + xlab('') + ylab('')
ggsave('~/Dropbox (Mason Lab)/i4/revisions/human_read_alignment/human_read_alignment.pdf',width=6,height=3.5)

ggplot(merged %>% filter(type == 'COUNTS'),aes(x = reference, y= value)) + theme_cowplot()+ ggbeeswarm::geom_quasirandom() + ggpubr::stat_compare_means()  + ggtitle('# READS ALIGNED, T2T vs HG38') + geom_line(aes(group = sampleid),color='gray',alpha=.6) + ylab('# Aligned Reads') + xlab('Reference Genome')
ggsave('~/Dropbox (Mason Lab)/i4/revisions/human_read_alignment/t2t_v_hg38.png',width=6,height=6)

### LOAD IN AND GET PERCENT ALIGNED FOR OTHER TAX CALLS

# load in ACTUAL READCOUNTS
readcounts = read.delim('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/read_counts_i4_data.tsv',sep='\t',header=F)%>% mutate(V1 = gsub('Sample_','',V1)) %>% mutate(V1 = gsub('ELMB_','',V1)) 
readcounts = readcounts %>% mutate(V1 = if_else(!grepl('SPX',V1),strsplit(V1,'_') %>% map_chr(1),V1))%>% column_to_rownames('V1') %>% t

colnames(readcounts) = gsub('_metaphlan','',colnames(readcounts))
readcounts = sanitize_sample_names(readcounts) %>% melt
colnames(readcounts) = c('variable','total_read_counts')
readcounts$total_read_counts = readcounts$total_read_counts*2

#### METAGENOMICS
# mag alignments
a = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/ASSEMBLED-BACTERIAL-MAGS_.05_.0025_metagenomics_species_counts.rds') %>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp')%>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- BACTERIAL_MAGS')
b = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/ASSEMBLED-BACTERIAL-MAGS_.05_.0025_metatranscriptomics_species_counts.rds')%>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp') %>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- BACTERIAL_MAGS')

bacmagalign = bind_rows(a,b)

# viral bin alignments
a = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/ASSEMBLED-VIRAL-GENOMES_.1_.0.05_metagenomics_species_counts.rds') %>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp')%>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- VIRAL_MAGS')
b = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/ASSEMBLED-VIRAL-GENOMES_.1_.0.05_metatranscriptomics_species_counts.rds') %>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp')%>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- VIRAL_MAGS')

viralmagalign = bind_rows(a,b)

# genbank alignments
a = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/genbank-viral_.01_.0.005_metagenomics_species_counts.rds')%>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp') %>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- GENBANK_VIRUSES')
b = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/genbank-viral_.01_.0.005_metatranscriptomics_species_counts.rds') %>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp') %>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- GENBANK_VIRUSES')

genbankalign = bind_rows(a,b)

# gtdb alignments
a = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/GTDB_.005_.0025_metagenomics_species_counts.rds') %>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp') %>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- GTDB')
b = readRDS('~/Dropbox (Mason Lab)/i4/revisions/xtree_read_alignment_data/GTDB_.005_.0025_metatranscriptomics_species_counts.rds') %>% rownames_to_column('temp') %>% filter(temp!='Unknown')%>% column_to_rownames('temp') %>% sanitize_sample_names() %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'XTREE -- GTDB')

gtdbalign = bind_rows(a,b)

# metaphlan4 alignments
metaphlan4 = read.delim('~/Dropbox (Mason Lab)/i4/revisions/metaphlan_readcounts.tsv',check.names=F,header=T)
colnames(metaphlan4) = c('V1','V2')
metaphlan4 = metaphlan4 %>% column_to_rownames('V1') %>% t %>% sanitize_sample_names %>% t %>% data.frame %>% dplyr::rename(valsum = V2) %>% rownames_to_column('variable')%>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS'))%>% mutate(method = 'METAPHLAN4')


# phanta

phanta = read.delim('~/Dropbox (Mason Lab)/i4/revisions/phanta_total_reads.tsv',sep='\t',header=T)%>% select(Samp_Name,Assigned_Step_Three) %>% dplyr::rename(SeqID = Samp_Name) %>%  sanitize_sample_names1 %>% rename(variable = SeqID)%>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS'))%>% mutate(method = 'PHANTA') %>% rename(valsum = Assigned_Step_Three)

# kraken2 masked conf
kraken2maskedconf = read.delim('~/Dropbox (Mason Lab)/i4/i4_data_packet/kraken2_bracken_abundances/i4_wgs_mtx_bracken_conf_mask_species.tsv',check.names=F)  %>% column_to_rownames('name') %>% select(all_of(grep('bracken_num',colnames(.)))) 
colnames(kraken2maskedconf) = gsub('.S.maskedCONF.bracken_num','',colnames(kraken2maskedconf))
kraken2maskedconf = sanitize_sample_names(kraken2maskedconf) %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS'))%>% mutate(method = 'KRAKEN2 -- MASK/CONF-0.2 -- REFSEQ')

# kraken2 nomasked conf
kraken2nomaskedconf = read.delim('~/Dropbox (Mason Lab)/i4/i4_data_packet/kraken2_bracken_abundances/i4_wgs_mtx_bracken_conf_nomask_species.tsv',check.names=F)  %>% column_to_rownames('name') %>% select(all_of(grep('bracken_num',colnames(.)))) 
colnames(kraken2nomaskedconf) = gsub('.S.CONF.bracken_num','',colnames(kraken2nomaskedconf))
kraken2nomaskedconf = sanitize_sample_names(kraken2nomaskedconf) %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS'))%>% mutate(method = 'KRAKEN2 -- CONF-0.2 -- REFSEQ')

# kraken2 masked noconf
kraken2maskednoconf = read.delim('~/Dropbox (Mason Lab)/i4/i4_data_packet/kraken2_bracken_abundances/i4_wgs_mtx_bracken_noconf_mask_species.tsv',check.names=F)  %>% column_to_rownames('name') %>% select(all_of(grep('bracken_num',colnames(.)))) 
colnames(kraken2maskednoconf) = gsub('.S.masked.bracken_num','',colnames(kraken2maskednoconf))
kraken2maskednoconf = sanitize_sample_names(kraken2maskednoconf)%>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS'))%>% mutate(method = 'KRAKEN2 -- MASK -- REFSEQ')

# kraken2 nomasked noconf
kraken2nomaskednoconf = read.delim('~/Dropbox (Mason Lab)/i4/i4_data_packet/kraken2_bracken_abundances/i4_wgs_mtx_bracken_noconf_nomask_species.tsv',check.names=F)  %>% column_to_rownames('name') %>% select(all_of(grep('bracken_num',colnames(.)))) 
colnames(kraken2nomaskednoconf) = gsub('.S.qc.bracken_num','',colnames(kraken2nomaskednoconf))
kraken2nomaskednoconf = kraken2nomaskednoconf %>% melt %>% filter(value>0) %>% dplyr::group_by(variable) %>% dplyr::summarise(valsum = sum(value)) %>% mutate(type = if_else(grepl('CEM',variable),'METATRANSCRIPTOMICS','METAGENOMICS')) %>% mutate(method = 'KRAKEN2 -- REFSEQ')

alldata = bind_rows(bacmagalign,gtdbalign,genbankalign,viralmagalign,kraken2nomaskednoconf,kraken2maskednoconf,kraken2nomaskedconf,kraken2maskedconf,metaphlan4,phanta)

alldata = full_join(alldata,readcounts)
alldata$readfrac = alldata$valsum/alldata$total_read_counts
alldata$valsum = log10(alldata$valsum)

### ADD IN METADATA 

# load metadata
meta = read.csv('~/Dropbox (Mason Lab)/i4/i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
meta$SeqID = gsub('SW_','',meta$SeqID)
meta$Timepoint_Recode = factor(meta$Timepoint)
levels(meta$Timepoint_Recode) = c(NA,'PRE-LAUNCH','POST-LAUNCH','PRE-LAUNCH','MID-FLIGHT','MID-FLIGHT','POST-LAUNCH','POST-LAUNCH','PRE-LAUNCH')

meta = meta %>% distinct %>% mutate(Timepoint_Recode2 = if_else(as.character(Timepoint_Recode) == 'MID-FLIGHT',Timepoint,as.character(Timepoint_Recode)))
meta$Timepoint_Recode2 = factor(meta$Timepoint_Recode2,levels = c('PRE-LAUNCH','Flight 1','Flight 2','POST-LAUNCH'))
meta$Timepoint = factor(meta$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec',NA))
meta$Timepoint_Numeric = as.numeric(meta$Timepoint)

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

alldata_meta = inner_join(alldata,metasub,by=c('variable'='SeqID'))%>% dplyr::rename(SeqID = variable, aligned_read_counts = valsum, aligned_read_fraction = readfrac) %>% select(SeqID,type,method,location,aligned_read_counts,aligned_read_fraction) %>% melt(id.vars = c('SeqID','type','method','location'))

p1 = ggplot(alldata_meta %>% filter(variable == 'aligned_read_counts'),aes(x=SeqID,y=value,fill=location)) + theme_minimal() + geom_bar(stat='identity') + theme_cowplot() + theme(axis.text.x = element_blank())  + facet_wrap(method ~ type + variable ,scales= 'free_x',ncol = 4) + scale_fill_brewer(palette = 'Set3')
ggsave(plot=p1,'~/Dropbox (Mason Lab)/i4/revisions/human_read_alignment/read_alignment_by_algorithm_counts.pdf',width=14,height=10)

p2 = ggplot(alldata_meta %>% filter(variable != 'aligned_read_counts'),aes(x=SeqID,y=value,fill=location))+ theme_minimal() + geom_bar(stat='identity') + theme_cowplot() + theme(axis.text.x = element_blank()) +ylim(0,1)+ facet_wrap(method ~ type + variable ,scales= 'free_x',ncol = 4) + scale_fill_brewer(palette = 'Set3')
ggsave(plot = p2,'~/Dropbox (Mason Lab)/i4/revisions/human_read_alignment/read_alignment_by_algorithm_percent.pdf',width=14,height=10)

p3 = ggplot(alldata_meta,aes(y=value,x=method)) + geom_boxplot() + theme_minimal()+ theme(axis.text.x = element_text(angle=45,hjust=1))+ facet_grid(variable ~ type ,scales= 'free')
ggsave(plot = p3,'~/Dropbox (Mason Lab)/i4/revisions/human_read_alignment//overall_read_alignment.pdf',width=10,height=8)







