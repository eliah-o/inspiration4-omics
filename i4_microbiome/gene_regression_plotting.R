# parse gene catalog analysis output

library(ggplot2)
library(tidyverse)

setwd('~/Dropbox (Mason Lab)/i4/regression_output/gene_regressions/')

hyps = readRDS('hypothetical_protein_outputs.rds') 
hyps = hyps %>% filter(grepl('Timepoint',term))%>% mutate(direction = if_else(estimate>0,'Decreased in launch','Increased in launch'))
nonhyps = readRDS('nonhypothetical_protein_outputs.rds') %>% filter(grepl('Timepoint',term))
nonhyps = nonhyps  %>% filter(grepl('Timepoint',term))%>% mutate(direction = if_else(estimate>0,'Decreased in launch','Increased in launch'))

regressiondata = bind_rows(hyps,nonhyps)
regressiondata$term = gsub('Timepoint_Recode','',regressiondata$term)
regressiondata$term = gsub(':','',regressiondata$term)
regressiondata$term = gsub(':Body.Location','',regressiondata$term)
regressiondata$term = gsub('PRE-LAUNCH','PRE-LAUNCH --- ',regressiondata$term)
regressiondata$term = gsub('POST-LAUNCH','POST-LAUNCH --- ',regressiondata$term)
regressiondata$term = gsub('isskin','Skin',regressiondata$term)
regressiondata$term = gsub('isoral','Oral',regressiondata$term)
regressiondata$term = gsub('nape','Nape of Neck',regressiondata$term)
regressiondata$term = gsub('postauric','Post-Auricular',regressiondata$term)
regressiondata$term = gsub('bb','Belly Button',regressiondata$term)
regressiondata$term = gsub('gc','Gluteal Crease',regressiondata$term)
regressiondata$term = gsub('Tzone','T-Zone',regressiondata$term)
regressiondata$term = gsub('web','Toe Web Space',regressiondata$term)
regressiondata$term = gsub('fore','Forearm',regressiondata$term)
regressiondata$term = gsub('isnasal','Nasal',regressiondata$term)

regressiondata = regressiondata %>% mutate(term = if_else(term == 'PRE-LAUNCH --- ','PRE-LAUNCH',term,))
regressiondata = regressiondata %>% mutate(term = if_else(term == 'POST-LAUNCH --- ','POST-LAUNCH',term,))

regressiondata$term = factor(regressiondata$term,levels = c("PRE-LAUNCH","POST-LAUNCH","PRE-LAUNCH --- Armpit","PRE-LAUNCH --- Belly Button","PRE-LAUNCH --- Forearm", "PRE-LAUNCH --- Gluteal Crease","PRE-LAUNCH --- Nape of Neck","PRE-LAUNCH --- Nasal","PRE-LAUNCH --- Post-Auricular","PRE-LAUNCH --- T-Zone","PRE-LAUNCH --- Toe Web Space","PRE-LAUNCH --- Oral","POST-LAUNCH --- Armpit","POST-LAUNCH --- Belly Button", "POST-LAUNCH --- Forearm","POST-LAUNCH --- Gluteal Crease" ,"POST-LAUNCH --- Nape of Neck", "POST-LAUNCH --- Nasal", "POST-LAUNCH --- Post-Auricular","POST-LAUNCH --- T-Zone", "POST-LAUNCH --- Toe Web Space","POST-LAUNCH --- Oral"))

# plot total up/down -
sigcountsdna = regressiondata %>% mutate(direction = if_else(estimate>0,'Decreased in launch','Increased in launch')) %>% filter(tag == 'Metagenomics',adjusted<0.05) %>% select(term,direction) %>% group_by(term) %>% count(direction) %>% mutate(type = 'Metagenomics')

sigcountsrna = regressiondata  %>% mutate(direction = if_else(estimate>0,'Decreased in launch','Increased in launch')) %>% filter(tag == 'Metatranscriptomics',adjusted<0.05) %>% select(term,direction) %>% group_by(term) %>% count(direction) %>% mutate(type = 'Metatranscriptomics')

sigcounts = bind_rows(sigcountsdna,sigcountsrna)
sigcounts = sigcounts %>% mutate(term = if_else(!grepl('---',as.character(term)),paste(as.character(term),'--- Overall'),as.character(term))) %>% mutate(`Location` = strsplit(term,' --- ') %>% map_chr(2),`Time` = strsplit(term,' --- ') %>% map_chr(1)) %>% arrange(desc(n))

sigcounts$Location = fct_reorder(as.factor(sigcounts$Location),sigcounts$n)
sigcounts$Time = factor(sigcounts$Time,levels = c("PRE-LAUNCH","POST-LAUNCH"))

ggplot(data = sigcounts %>% filter(Location != 'Overall'),aes(x=n,y=`Location`,fill=direction)) +geom_bar(stat='identity',position='dodge',width=.95) + theme_bw() + facet_grid(rows = vars(type),cols = vars(Time)) + scale_fill_viridis_d()
ggsave(paste('sig_launch_by_body_site_','genes','.pdf',sep=''))

#### HOW TO HANDLE WHICH TERMS TO INCLUDE
# build concordance barplots
concordanceplot = regressiondata %>% filter(tag == 'Metatranscriptomics' & grepl('PRE-LAUNCH',term)| tag == 'Metagenomics',tag == 'Metatranscriptomics' & adjusted<0.05 | tag == 'Metagenomics' ) %>% select(term,estimate,yvar,regressionclass,tag,adjusted,Product) %>% group_by(term,tag) %>% arrange(desc(estimate))

concordanceplot2 = inner_join(concordanceplot %>% filter(tag == 'Metatranscriptomics'),concordanceplot%>% filter(tag == 'Metagenomics') %>% select(-Product)%>% mutate(estimate = if_else(adjusted>0.05,0,estimate)),by=c('term','yvar'))
concordanceplot2 = concordanceplot2 %>% mutate(mean = (estimate.x + estimate.y)/2,concordance = if_else(estimate.x<0 & estimate.y<0 | estimate.x>0 & estimate.y>0,'Concordant','Discordant'))

p1 = ggplot(concordanceplot2,aes(x=estimate.x,y=estimate.y,color=concordance)) + geom_point(alpha=.5) + theme_bw() + xlab('Estimate, pre-flight, MTX') + ylab('Estimate, pre-flight, WGS') + scale_color_manual(values = c('black','blue')) + theme(legend.position = 'bottom') 
sub = concordanceplot2%>% ungroup %>% filter(concordance == 'Discordant',Product != 'hypothetical protein')%>% slice_max(abs(estimate.x),n=20) %>% arrange(desc(estimate.x))
sub$Product = fct_reorder(sub$Product,sub$estimate.x)
p2 = ggplot(sub,aes(x=estimate.x,y=Product)) + geom_bar(stat='identity') + theme_bw() + xlab('Estimate, pre-flight, MTX') + ylab('')+ scale_y_discrete(position = "right")
p1|p2
ggsave(paste('gene_concordance.pdf',sep=''),width=15,height=6)

# get top negative positive hits by overall, oral, nasal
regressiondata_nonhyps = regressiondata %>% filter(effect == 'fixed',Product != 'hypothetical protein',adjusted<0.05)%>% mutate(direction = if_else(estimate>0,'Increased','Decreased'))
regressiondata_nonhyps = regressiondata_nonhyps %>% mutate(`Time` = strsplit(as.character(term),' --- ') %>% map_chr(1)) %>% filter(Time == 'PRE-LAUNCH')
regressiondata_nonhyps2 = regressiondata_nonhyps%>% filter(tag == 'Metagenomics') %>% group_by(Product,tag,regressionclass,direction) %>% count %>% slice_max(n,n=10) %>% filter(n>15)
regressiondata_nonhyps3 = regressiondata_nonhyps%>% filter(tag == 'Metatranscriptomics') %>% group_by(Product,tag,regressionclass,direction) %>% count %>% slice_max(n,n=5) %>% filter(n>1)
#regressiondata_nonhyps3 = bind_rows(regressiondata_nonhyps2,regressiondata_nonhyps3)
regressiondata_nonhyps2$Product = fct_reorder(regressiondata_nonhyps2$Product,regressiondata_nonhyps2$n)
regressiondata_nonhyps3$Product = fct_reorder(regressiondata_nonhyps3$Product,regressiondata_nonhyps3$n)
p1 = ggplot(data = regressiondata_nonhyps2 %>% filter(direction == 'Decreased'),aes(y=Product, x = n)) + geom_bar(stat='identity') + facet_grid(rows=vars(regressionclass),cols = vars(tag),scale='free') + xlab('') + ylab('')
p2 = ggplot(data = regressiondata_nonhyps3 %>% filter(direction == 'Decreased'),aes(y=Product, x = n)) + geom_bar(stat='identity') + facet_grid(rows=vars(regressionclass),cols = vars(tag),scale='free') + xlab('') + ylab('')
p1|p2
ggsave('top_increased_changed_functions.pdf',width=8,height=12)











