library(dplyr)
library(maplet)
library(ggplot2)


load("./preprocessed_i4_metabolomics.RData")

#### Figure 4 ####
#Metabolites of interest
mets<-c("LysoPC(15:0)", "LysoPC(14:0)", "LysoPC(17:0)","LysoPC(16:1)", "Sphingosine 1-phosphate", "Protoporphyrin", "4-Aminobutanoate (GABA)")

# Combine metabolite data with time data
violin_dat <- data.frame(t(assay(i4_D)[match(mets,rowData(i4_D)$Name ),]), colData(i4_D)$time)
names(violin_dat)<-c(mets, "time")

# Label pre- and post-  flight
violin_dat$group <- unlist(lapply(violin_dat$time, function(i){
  if (i<4){
    return("L-92; L-44; L-3")
  }
  else{
    return("R+1; R+45; R+82")
  }
}))

# Label immediate post-flight
addit <- violin_dat[violin_dat$time==4,]
addit$group <- rep("R+1", nrow(addit))
violin_dat <- rbind(violin_dat, addit)

# For each metabolite, make a violin plot
plot_list <- list()
tp_colors = c("#63297A", "#E41A1C", "#057652")
for (met in mets){
  plot_list[[met]]<-ggplot(violin_dat, aes(x=group, y=.data[[met]], fill=group)) + 
    geom_violin() + 
    scale_fill_manual(values=tp_colors) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1, fill="black") + 
    theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
    ylab("Abundance")+
    ggtitle(paste0("Metabolite: ", met))
}


#### Extended Data Figure 6 #### 
# Build bar plot df
bar_plot_df <- data.frame(rowData(i4_D)) %>% 
  dplyr::group_by(SUPER_PATHWAY, SUB_PATHWAY) %>% 
  dplyr::add_count() %>% # Get count of unique sub pathways
  dplyr::ungroup() %>% 
  dplyr::group_by(SUPER_PATHWAY) %>% 
  dplyr::add_count() %>% # Get count of unique super pathways
  dplyr::mutate(SUPER_PATHWAY = ifelse(nn >5, SUPER_PATHWAY, "Other")) %>% # Make small pathways other
  dplyr::mutate(SUB_PATHWAY = ifelse(n >5, SUB_PATHWAY, "Other")) %>% # Make small pathways other
  dplyr::rename(`Metabolite Class`="SUPER_PATHWAY") %>% 
  dplyr::rename(`Metabolite Subclass` = "SUB_PATHWAY")

### Plot ggplot bar plot with main category super pathway stacked with subpathway
order<-bar_plot_df %>% 
  dplyr::select(`Metabolite Class`,nn) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(`Metabolite Class`) %>% 
  dplyr::summarize(sum(nn)) %>% 
  arrange(`sum(nn)`) %>% 
  .$`Metabolite Class`

ggplot(bar_plot_df, aes(y = `Metabolite Class`, fill=`Metabolite Subclass`))+
  geom_bar()+
  scale_y_discrete(limits = order)+
  theme(legend.position=c(0.65,0.3))+
  scale_fill_discrete(name=NULL)+
  ggtitle("Metabolite Class and Subclass Annotations")

