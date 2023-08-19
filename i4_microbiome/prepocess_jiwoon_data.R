# preprocessing for jiwoon's data

library(tidyverse)
library(reshape2)

i4swabdata = read.delim('i4_swab_metadata.csv',sep=',',header=T) %>% filter(grepl('Deltoid',Body.Location))%>% mutate(Timepoint = if_else(Timepoint == '21-Aug','Pre','Post')) %>% rename(SlideName = Crew.ID) %>% select(SlideName,Timepoint,SeqID)

jiwoondata = read.delim('SkinSpatial_Counts.csv',sep=',',header=T)
datacols = colnames(jiwoondata[,36:ncol(jiwoondata)])

jiwoon_swabdata_mapped = left_join(jiwoondata,i4swabdata) 

bodylocs = unique(unlist(unname(jiwoon_swabdata_mapped$ROIType)))

for(b in bodylocs){
	sub = jiwoon_swabdata_mapped  %>% filter(ROIType == b)%>% select(SeqID,SlideName,all_of(datacols)) %>% melt(id.vars = c('SeqID','SlideName')) %>% group_by(SeqID,SlideName,variable) %>% summarise(value = mean(value))%>% ungroup %>% dcast(SeqID ~ variable,value.var = 'value')
	sub = sub %>% column_to_rownames('SeqID')
	sub = sub %>% select(all_of(datacols))
	saveRDS(sub,paste0('skin_swab_spatial_t_',b,'.rds'))
}

### outputs should be:
# jiwoon raw data with the correct sample ids -- one matrix for every layer