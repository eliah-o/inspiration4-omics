# plot out the strain sharing data

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(reshape2)


metaphlan4 = read.delim('~/Dropbox (Mason Lab)/i4/i4_data_packet/metaphlan4/metaphlan_i4.tsv',check.names=F) %>% filter(grepl('t__',clade_name))
names = metaphlan4 %>% select(clade_name) %>% mutate(strain = strsplit(clade_name,'\\|') %>% map_chr(8)) %>% mutate(species = strsplit(clade_name,'\\|s__') %>% map_chr(2))%>% select(clade_name,species,strain)

setwd('~/Dropbox (Mason Lab)/i4/revisions/strain_sharing/')

strainshare = read.delim('~/Dropbox (Mason Lab)/i4/revisions/strain_sharing/i4_wgs_shared_strains.tsv')

negatives = strainshare %>% filter(bodyloc1 == 'Control Swab (0)'|bodyloc2 == 'Control Swab (0)')

strainshare = strainshare %>% filter(!(strain %in% negatives$strain)) %>% filter(bodyloc1 != 'Swab Water', bodyloc2 != 'Swab Water',bodyloc1 != 'Control Swab (0)',bodyloc2 != 'Control Swab (0)',bodyloc1 != "Deltoid - Pre-Biospy",bodyloc2 != "Deltoid - Pre-Biospy")
strainshare = strainshare %>% mutate(bodyloc1 = if_else(subjectID_1 == 'Capsule',gsub("^", "CAPSULE: ", bodyloc1),gsub("^", "HUMAN: ",bodyloc1)))
strainshare = strainshare %>% mutate(bodyloc2 = if_else(subjectID_2 == 'Capsule',gsub("^", "CAPSULE: ", bodyloc2),gsub("^", "HUMAN: ",bodyloc2)))
                                                                          
strainshare$timepoint1_3 = factor(strainshare$timepoint1_3,levels = c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec'))
strainshare$timepoint2_3 = factor(strainshare$timepoint2_3,levels = c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec'))

strainshare$timepoint1[strainshare$timepoint1 == 'POST-LAUNCH'] = 'POST-FLIGHT'
strainshare$timepoint2[strainshare$timepoint2 == 'POST-LAUNCH'] = 'POST-FLIGHT'

strainshare$timepoint1 = factor(strainshare$timepoint1,levels = c('PRE-LAUNCH','MID-FLIGHT','POST-FLIGHT'))
strainshare$timepoint2 = factor(strainshare$timepoint2,levels = c('PRE-LAUNCH','MID-FLIGHT','POST-FLIGHT'))

# number shared between individuals across time
strainshare_sub_between = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual!='same_individual',timepoint1_3 == timepoint2_3) %>% select(subjectID_1,subjectID_2,timepoint1_3,strain) %>% distinct 
strainshare_sub_between = strainshare_sub_between%>%mutate(subjectID_1 = pmin(subjectID_1, subjectID_2), subjectID_2 = pmax(subjectID_1, subjectID_2)) %>%distinct(subjectID_1, subjectID_2, strain, .keep_all = TRUE)
strainshare_sub_between = strainshare_sub_between %>% select(timepoint1_3,strain) %>% group_by(timepoint1_3) %>% count %>% mutate(type = 'Between astronauts')

strainshare_sub_within = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual=='same_individual',timepoint1_3 == timepoint2_3) %>% select(subjectID_1,subjectID_2,timepoint1_3,strain) %>% distinct 
strainshare_sub_within = strainshare_sub_within%>%mutate(subjectID_1 = pmin(subjectID_1, subjectID_2), subjectID_2 = pmax(subjectID_1, subjectID_2)) %>%distinct(subjectID_1, subjectID_2, strain, .keep_all = TRUE)
strainshare_sub_within = strainshare_sub_within %>% select(timepoint1_3,strain) %>% group_by(timepoint1_3) %>% count %>% mutate(type = 'Within a single astronaut')


strainshare_sub_capsule = strainshare %>% filter(subjectID_1 == 'Capsule' | subjectID_2 == 'Capsule',same_individual!='same_individual',timepoint1_3 == timepoint2_3) %>% select(subjectID_1,subjectID_2,timepoint1_3,strain) %>% distinct
strainshare_sub_capsule = strainshare_sub_capsule%>%mutate(subjectID_1 = pmin(subjectID_1, subjectID_2), subjectID_2 = pmax(subjectID_1, subjectID_2)) %>%distinct(subjectID_1, subjectID_2, strain, .keep_all = TRUE)
strainshare_sub_capsule = strainshare_sub_capsule %>% select(timepoint1_3,strain) %>% group_by(timepoint1_3) %>% count %>% mutate(type = 'Capsule and astronauts')

strainsharecount = bind_rows(strainshare_sub_between,strainshare_sub_within,strainshare_sub_capsule)

ggplot(strainsharecount,aes(x = timepoint1_3, y = n,color = type, group = type,fill=type)) + geom_point(size=6,shape='diamond',alpha=.8) + geom_line(aes(x = timepoint1_3, y = n)) + theme_cowplot() + ggtitle('Strain sharing events') + ylab('Count') + xlab('') + theme(axis.text.x=element_text(angle=45,hjust = 1),legend.title = element_blank()) 
ggsave('strains_shared_over_time.pdf',width=6,height=4)


# body site stransit within individuals
strainshare_sub_between = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual!='same_individual',timepoint1 == timepoint2) %>% select(bodyloc1,bodyloc2,timepoint1,strain) %>% distinct 
strainshare_sub_between = strainshare_sub_between%>%mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, strain, .keep_all = TRUE)%>%mutate(Pair = paste(sep='---',pmin(bodyloc1, bodyloc2), pmax(bodyloc1, bodyloc2))) %>% arrange(Pair) %>% select(-Pair)
strainshare_sub_between = strainshare_sub_between %>% dcast(timepoint1 + bodyloc1 ~ bodyloc2,value.var ='strain') %>% arrange(timepoint1)

strainshare_sub_between$bodyloc1 = strsplit(strainshare_sub_between$bodyloc1,': ') %>% map_chr(2)
colnames(strainshare_sub_between) = gsub('HUMAN: ','',colnames(strainshare_sub_between))
colnames(strainshare_sub_between) = gsub('CAPSULE: ','',colnames(strainshare_sub_between))


pdf('between_ind_heatmap.pdf',width=8,height=6)
Heatmap(strainshare_sub_between %>% select(-bodyloc1,-timepoint1),cluster_columns=T,cluster_rows=T,row_split = strainshare_sub_between$timepoint1,row_labels = strainshare_sub_between$bodyloc1,cluster_row_slices = F,name = 'Between astronauts')
dev.off()
# body site transit within individuals
strainshare_sub_within = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual=='same_individual',timepoint1 == timepoint2) %>% select(bodyloc1,bodyloc2,timepoint1,strain) %>% distinct 
strainshare_sub_within = strainshare_sub_within%>%mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, strain, .keep_all = TRUE)%>%mutate(Pair = paste(sep='---',pmin(bodyloc1, bodyloc2), pmax(bodyloc1, bodyloc2))) %>% arrange(Pair) %>% select(-Pair)
strainshare_sub_within = strainshare_sub_within %>% dcast(timepoint1 + bodyloc1 ~ bodyloc2,value.var ='strain') %>% arrange(timepoint1)# %>% melt
#strainshare_sub_within[colnames(strainshare_sub_within) == strainshare_sub_within$bodyloc1] = 0

strainshare_sub_within$bodyloc1 = strsplit(strainshare_sub_within$bodyloc1,': ') %>% map_chr(2)
colnames(strainshare_sub_within) = gsub('HUMAN: ','',colnames(strainshare_sub_within))
colnames(strainshare_sub_within) = gsub('CAPSULE: ','',colnames(strainshare_sub_within))

pdf('within_ind_heatmap.pdf',width=8,height=6)
Heatmap(strainshare_sub_within %>% select(-bodyloc1,-timepoint1),cluster_columns=T,cluster_rows=T,row_split = strainshare_sub_within$timepoint1,row_labels = strainshare_sub_within$bodyloc1,cluster_row_slices = F,name = 'Within a single astronaut',col = colorRamp2(c(0,15,30), c("black","red", "cyan")))
dev.off()

# body site transit within capsule
strainshare_sub_capsule = strainshare %>% filter(subjectID_1 == 'Capsule' | subjectID_2 == 'Capsule',!(subjectID_1 == 'Capsule' &subjectID_2 == 'Capsule'),timepoint1 == timepoint2) %>% select(bodyloc1,bodyloc2,timepoint1,strain) %>% distinct 
strainshare_sub_capsule = strainshare_sub_capsule %>% mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, timepoint1,strain, .keep_all = TRUE) %>%mutate(Pair = paste(sep='---',pmin(bodyloc1, bodyloc2), pmax(bodyloc1, bodyloc2))) %>% arrange(Pair) %>% select(Pair,timepoint1,strain) 
 
strainshare_sub_capsule = strainshare_sub_capsule %>% group_by(Pair,timepoint1) %>% count() %>% mutate(bodyloc1 = strsplit(Pair,'---') %>% map_chr(1),bodyloc2 = strsplit(Pair,'---') %>% map_chr(2))%>% ungroup %>% select(-Pair) %>% filter(bodyloc1!=bodyloc2)# %>% filter(grepl('\\(',bodyloc1))%>% filter(!grepl('\\(',bodyloc2))

strainshare_sub_capsule = strainshare_sub_capsule %>% dcast(timepoint1 + bodyloc1 ~ bodyloc2,value.var ='n') %>% arrange(timepoint1)
strainshare_sub_capsule[is.na(strainshare_sub_capsule)] = 0
strainshare_sub_capsule$Oral = 0
strainshare_sub_capsule$`Gluteal Crease` = 0
strainshare_sub_capsule$`Belly Button` = 0
strainshare_sub_capsule$`Armpit` = 0
strainshare_sub_capsule$Oral = 0

strainshare_sub_capsule$bodyloc1 = strsplit(strainshare_sub_capsule$bodyloc1,': ') %>% map_chr(2) 
strainshare_sub_capsule$bodyloc1 = strsplit(strainshare_sub_capsule$bodyloc1 ,' \\(') %>% map_chr(1)
colnames(strainshare_sub_capsule) = gsub('HUMAN: ','',colnames(strainshare_sub_capsule))
colnames(strainshare_sub_capsule) = gsub('CAPSULE: ','',colnames(strainshare_sub_capsule))

pdf('capsule_sharing_heatmap.pdf',width=7,height=3)
Heatmap(strainshare_sub_capsule %>% select(-bodyloc1,-timepoint1),cluster_columns=T,cluster_rows=T,row_split = strainshare_sub_capsule$timepoint1,row_labels = strainshare_sub_capsule$bodyloc1,'Between capsule and astronauts',col = colorRamp2(c(0,2,4), c("black","red", "cyan")))
dev.off()

### GET PROMISCUOUS STRAINS
strainshare_sub_between = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual!='same_individual',timepoint1 == timepoint2) %>% select(bodyloc1,bodyloc2,timepoint1,strain) %>% distinct 
strainshare_sub_between = strainshare_sub_between%>%mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, strain, .keep_all = TRUE) %>% select(strain,timepoint1) %>% mutate(type = 'Between\nastronauts')

strainshare_sub_within = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual=='same_individual',timepoint1 == timepoint2) %>% select(bodyloc1,bodyloc2,timepoint1,strain) %>% distinct 
strainshare_sub_within = strainshare_sub_within%>%mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, strain, .keep_all = TRUE) %>% select(strain,timepoint1) %>% mutate(type = 'Within a\nsingle astronaut')

strainshare_sub_capsule = strainshare %>% filter(subjectID_1 == 'Capsule' | subjectID_2 == 'Capsule',!(subjectID_1 == 'Capsule' &subjectID_2 == 'Capsule'),timepoint1 == timepoint2) %>% select(bodyloc1,bodyloc2,timepoint1,strain) %>% distinct 
strainshare_sub_capsule = strainshare_sub_capsule %>% mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, timepoint1,strain, .keep_all = TRUE) %>%mutate(Pair = paste(sep='---',pmin(bodyloc1, bodyloc2), pmax(bodyloc1, bodyloc2))) %>% arrange(Pair) %>% select(strain,timepoint1) %>% mutate(type = 'Between capsule\nand astronauts')

promdata = bind_rows(strainshare_sub_between,strainshare_sub_within,strainshare_sub_capsule)
promdata = promdata %>% group_by(type,timepoint1,strain) %>% count() 

promdata$type = factor(promdata$type, levels = c('Within a\nsingle astronaut','Between\nastronauts','Between capsule\nand astronauts'))

promdata = inner_join(promdata,names) #%>% dplyr::rename('Strain sharing events\nper taxon')

ggplot(promdata %>% filter(n>1),aes(x = timepoint1,y=reorder(species,n),size=n)) + geom_point() + theme_cowplot() + tidytext::scale_y_reordered() + ylab('') + xlab('') + facet_grid(. ~ type)+ theme(axis.text.x = element_text(angle = 45,hjust =1))

#### PROM BUGS BY TIMEPOINT

strainshare_sub_between = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual!='same_individual',timepoint1_3 == timepoint2_3) %>% select(bodyloc1,bodyloc2,timepoint1_3,strain) %>% distinct 
strainshare_sub_between = strainshare_sub_between%>%mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, strain, .keep_all = TRUE) %>% select(strain,timepoint1_3) %>% mutate(type = 'Between\nastronauts')

strainshare_sub_within = strainshare %>% filter(subjectID_1 != 'Capsule',subjectID_2 != 'Capsule',same_individual=='same_individual',timepoint1_3 == timepoint2_3) %>% select(bodyloc1,bodyloc2,timepoint1_3,strain) %>% distinct 
strainshare_sub_within = strainshare_sub_within%>%mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, strain, .keep_all = TRUE) %>% select(strain,timepoint1_3) %>% mutate(type = 'Within a\nsingle astronaut')

strainshare_sub_capsule = strainshare %>% filter(subjectID_1 == 'Capsule' | subjectID_2 == 'Capsule',!(subjectID_1 == 'Capsule' &subjectID_2 == 'Capsule'),timepoint1_3 == timepoint2_3) %>% select(bodyloc1,bodyloc2,timepoint1_3,strain) %>% distinct 
strainshare_sub_capsule = strainshare_sub_capsule %>% mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2, timepoint1_3,strain, .keep_all = TRUE) %>%mutate(Pair = paste(sep='---',pmin(bodyloc1, bodyloc2), pmax(bodyloc1, bodyloc2))) %>% arrange(Pair) %>% select(strain,timepoint1_3) %>% mutate(type = 'Between capsule\nand astronauts')

promdata = bind_rows(strainshare_sub_between,strainshare_sub_within,strainshare_sub_capsule)
promdata = promdata %>% group_by(type,timepoint1_3,strain) %>% count() 

promdata$type = factor(promdata$type, levels = c('Within a\nsingle astronaut','Between\nastronauts','Between capsule\nand astronauts'))

promdata = inner_join(promdata,names) #%>% dplyr::rename(n = 'Strain sharing events\nper taxon')

promdata = promdata %>% filter(type != 'Between capsule\nand astronauts' | type == 'Between capsule\nand astronauts' & timepoint1_3 == '21-Jun' |type == 'Between capsule\nand astronauts' & timepoint1_3 == 'Flight 2')


ggplot(promdata %>% filter(n>1),aes(x = timepoint1_3,y=reorder(species,n),size=n)) + geom_point() + theme_cowplot() + tidytext::scale_y_reordered() + ylab('') + xlab('') + facet_grid(. ~ type)+ theme(axis.text.x = element_text(angle = 45,hjust =1))
ggsave('promiscuous_strains.pdf',width=10,height=6)


#### BIG CHORD DIAGRAM

strainshare_chord = strainshare %>%filter(bodyloc1 != bodyloc2)%>% filter(timepoint1 == 'MID-FLIGHT',timepoint2 == 'MID-FLIGHT') %>% select(subjectID_1,subjectID_2,bodyloc1,bodyloc2,strain)  %>% mutate(bodyloc1 = gsub('HUMAN: ','',bodyloc1))%>% mutate(bodyloc2 = gsub('HUMAN: ','',bodyloc2)) %>% mutate(bodyloc1 = if_else(grepl('CAPSULE: ',bodyloc1),'CAPSULE',bodyloc1))%>% mutate(bodyloc2 = if_else(grepl('CAPSULE: ',bodyloc2),'CAPSULE',bodyloc2))
strainshare_chord = strainshare_chord %>% mutate(bodyloc1 = pmin(bodyloc1, bodyloc2), bodyloc2 = pmax(bodyloc1, bodyloc2)) %>%distinct(bodyloc1, bodyloc2,strain, .keep_all = TRUE) %>%mutate(Pair = paste(sep='---',pmin(bodyloc1, bodyloc2), pmax(bodyloc1, bodyloc2))) %>% arrange(Pair)  %>% select(-Pair)

strainshare_chord = strainshare_chord %>%  group_by(subjectID_1,subjectID_2,bodyloc1,bodyloc2) %>% count  %>% mutate(a = paste0(subjectID_1,' --- ',bodyloc1),b = paste0(subjectID_2,' --- ',bodyloc2)) %>%ungroup


grouping = setNames(c(strainshare_chord$subjectID_1,strainshare_chord$subjectID_2),c(strainshare_chord$a,strainshare_chord$b))


pdf('~/Dropbox (Mason Lab)/i4/revisions/strain_sharing/circos_plot.odf',width=5,height=5)
circos.clear()
chordDiagramFromDataFrame(df = strainshare_chord %>% select(a,b,n),group = grouping,directional = 0,annotationTrack = "grid")
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0,.5))
}, bg.border = NA)
circos.clear()
dev.off()

