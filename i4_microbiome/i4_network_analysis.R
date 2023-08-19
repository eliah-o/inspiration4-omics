library(tidyverse)
library(ggraph)
library(tidygraph)
library(reshape2)
library(corrplot)
library(ComplexHeatmap)

setwd('~/Dropbox (Mason Lab)/i4/revisions/network_analysis/')

corrmapplot <- function(data,tokeep,t){
  #dates = c('PRE-LAUNCH','MID-FLIGHT','POST-LAUNCH')
  heatlist = list()
  heatlistlist = list()
  dates = c('Mid-Flight','21-Jun','21-Aug','Sept Pre-Launch','Sept Post-Return','November','21-Dec')
 # for(l in c('Oral','Skin','Nasal')){
    #regout_sub = regout %>% filter(grepl(tolower(l),tolower(regressionclass)))
    data_sub = data %>% rownames_to_column('yvar') %>% filter(yvar %in% unique(tokeep)) %>% column_to_rownames('yvar')
    #rows_to_keep <- apply(data, 1, function(row) sum(row != 0) > 10)
    #data_sub <- data_sub[rows_to_keep, ]
    count=0
  #  if(l == 'Skin'){
   #   l = metasub %>% filter(isskin == 1) %>% select(Body.Location) %>% unique %>% unlist %>% unname
  #  }
    for(d in dates){
      data_sub2 = data_sub %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(any_of(metasub %>% filter(Timepoint == d) %>% select(SeqID) %>% unlist %>% unname)) %>% t %>% cor()
      data_sub2[is.na(data_sub2)] = 0
      if(d == 'Mid-Flight'){
        tokeep2 = rowSums((data_sub2)) %>% data.frame %>% filter(.!=1) %>% rownames
      }
      data_sub2=data_sub2[tokeep2,tokeep2]
      count = count+1
      data_sub_annot = data_sub2 %>%data.frame %>% mutate(phylum = strsplit(rownames(data_sub2),'p__') %>% map_chr(2) %>% strsplit(';') %>% map_chr(1)) %>% arrange(phylum) %>% select(phylum)
      tokeep = table(data_sub_annot$phylum) %>% data.frame %>% arrange(desc(Freq)) %>% head(7) %>% select(Var1) %>% unlist %>% unname
      data_sub_annot = data_sub_annot %>% mutate(phylum = if_else(phylum %in% tokeep,phylum,'Other'))
      color_scale <- setNames(RColorBrewer::brewer.pal('Set1',n=8),unique(data_sub_annot$phylum))
      color_scale = color_scale[!is.na(names(color_scale))] 
      anno = ComplexHeatmap::rowAnnotation(phylum = data_sub_annot$phylum,col = list(phylum = color_scale))
      data_sub2 = data_sub2[,rownames(data_sub2)]
      #data_sub2[abs(data_sub2)<.7] = 0
      if(count == 1){
        dissimilarity_matrix <- as.dist(1 - abs(data_sub2))
        column_order <- order.dendrogram(as.dendrogram(hclust(dissimilarity_matrix)))
        data_sub2 = data_sub2[column_order, column_order]
        rowstokeep = colnames(data_sub2)
        heatlist[[d]] = Heatmap(data_sub2,cluster_rows= T, cluster_columns=F,show_row_names=F, show_column_names = F,name=d,show_heatmap_legend = F,show_row_dend = F)
      }
      if(count>1){
        data_sub2 = data_sub2[rowstokeep, rowstokeep]
        if(d == '21-Jun'){
          heatlist[[d]] = Heatmap(data_sub2,cluster_rows= T, cluster_columns=F,show_row_names=F, show_column_names = F,name=d,show_heatmap_legend = T,left_annotation = anno,heatmap_legend_param = list(title = "Pearson correlation", at = c(-1, 0, 1))) 
        }
        if(d!='21-Jun'){
          heatlist[[d]] = Heatmap(data_sub2,cluster_rows= T, cluster_columns=F,show_row_names=F, show_column_names = F,name=d,show_heatmap_legend = F) 
        }
      }
    }
    pdf(paste('~/Dropbox (Mason Lab)/i4/revisions/network_analysis/network_timeplot_',t,'.pdf',sep=''),width = 10,height=4)
    draw( heatlist[[2]] + heatlist[[3]] + heatlist [[4]] +heatlist[[1]] + heatlist[[5]] + heatlist[[6]] + heatlist[[7]], main_heatmap = "Mid-Flight")
    #draw(heatlist[[1]] + heatlist[[2]] + heatlist[[3]], main_heatmap = 'MID-FLIGHT')
    dev.off()
   # heatlistlist[[l]] = heatlist
 # }
  return(heatlistlist)
}

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


#### BACTERIAL
# load data
data = readRDS('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/gtdb/bacteria_xtree_005-0025_f___metagenomics_decontam.rds')
metag_filt = corrmapplot(data,rownames(data),'wgs')

# load data
data = readRDS('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/gtdb/bacteria_xtree_005-0025_f___metatranscriptomics_decontam.rds')
metat_filt = corrmapplot(data,rownames(data),'mtx')

#### VIRAL
# load data
data = readRDS('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/genbank/virus_xtree_01-005_f___metagenomics_decontam.rds')
metag_filt = corrmapplot(data,rownames(data),'wgs-viral')

# load data
data = readRDS('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/genbank/virus_xtree_01-005_f___metatranscriptomics_decontam.rds')
metat_filt = corrmapplot(data,rownames(data),'mtx-viral')









