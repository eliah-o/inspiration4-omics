# get removed taxa

setwd('~/Dropbox (Mason Lab)/i4/i4_data_packet/')

clist = read.delim('contamlist.tsv',sep='\t',header=T) %>% unlist %>% unname

output = list()
for(file in clist){
  f = readRDS(file) %>% filter(contaminant == TRUE) 
  file = gsub(' ','',file)
  file = gsub('contam_bac','contambac',file)
  algo = strsplit(file,'/')%>% map_chr(5) %>% strsplit('_') %>% map_chr(1) 
  seqtype = strsplit(file,'/')%>% map_chr(5) %>% strsplit('_') %>% map_chr(2) 
  org = strsplit(file,'/')%>% map_chr(5) %>% strsplit('_') %>% map_chr(3) %>% gsub('contam','',.) 
  params = strsplit(file,'/')%>% map_chr(5) %>% strsplit('_') %>% map_chr(5) 
  rank = strsplit(file,'/')%>% map_chr(5) %>% strsplit('_') %>% map_chr(6) 
  f =f %>% rownames_to_column('tax') %>% select(tax) %>% mutate(algo = algo,seqtype = seqtype,org=org,params=params,rank=rank)
  output[[file]] = f
}

output = bind_rows(output)
output = output %>% mutate(org = if_else(org == '' & algo == 'phanta','virus',org))

write.csv(output,'~/Dropbox (Mason Lab)/i4/revisions/Figures and Tables (I4WGSMTX)/revisions/SUPPTABLE8_taxa_removed_decontamination.csv')



