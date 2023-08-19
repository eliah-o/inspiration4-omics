library(tidyverse)
library(reshape2)

metaphlan = readRDS('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/metaphlan4/bacteria_metaphlan4_default_t__metagenomics_decontam.rds')
metaphlan = metaphlan %>% rownames_to_column('clade_name')%>% melt %>% filter(value>0) %>% group_by(clade_name) %>% count %>% arrange(desc(n))

a = metaphlan %>% filter(n>=3) %>% select(clade_name) %>% unlist %>% unname %>% strsplit('\\|') %>% map_chr(8) %>% data.frame %>% mutate(type  = 'decontaminated')
  
#metaphlan = readRDS('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/metaphlan4/bacteria_metaphlan4_default_t__metagenomics_nodecontam.rds')
#metaphlan = metaphlan %>% rownames_to_column('clade_name')%>% melt %>% filter(value>0) %>% group_by(clade_name) %>% count %>% arrange(desc(n))

#b= metaphlan %>% filter(n>10) %>% select(clade_name) %>% filter(!(clade_name %in% a$.)) %>% unlist %>% unname %>% strsplit('\\|') %>% map_chr(8) %>% data.frame %>% mutate(type  = 'not decontaminated')

#c = bind_rows(a,b)
write.table(a,'~/Dropbox (Mason Lab)/i4/i4_data_packet/bugs_for_strainphlan.tsv',quote=F,sep='\t')