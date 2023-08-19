# loading packages
library(readr)
library(dplyr)
library(ggplot2)
library(cutpointr)
library(tidyverse)

# load metadata
meta = read.csv('i4_swab_metadata.csv') %>% mutate(location = if_else(Crew.ID == 'Capsule','Capsule',Body.Location))
meta$SeqID = gsub('ELMB_','',meta$SeqID)
meta$SeqID = gsub('SW_','',meta$SeqID)
meta$Timepoint_Recode = factor(meta$Timepoint)
levels(meta$Timepoint_Recode) = c(NA,'PRE-LAUNCH','POST-LAUNCH','PRE-LAUNCH','MID-FLIGHT','MID-FLIGHT','POST-LAUNCH','POST-LAUNCH','PRE-LAUNCH')

meta = meta %>% distinct %>% mutate(Timepoint_Recode2 = if_else(as.character(Timepoint_Recode) == 'MID-FLIGHT',Timepoint,as.character(Timepoint_Recode)))
meta$Timepoint_Recode2 = factor(meta$Timepoint_Recode2,levels = c('PRE-LAUNCH','Flight 1','Flight 2','POST-LAUNCH'))
meta$Timepoint = factor(meta$Timepoint,levels=c('21-Jun','21-Aug','Sept Pre-Launch','Flight 1','Flight 2','Sept Post-Return','November','21-Dec',NA))
meta$Timepoint_Numeric = as.numeric(meta$Timepoint)
meta$Timepoint_Recode = factor(meta$Timepoint_Recode,levels = c('MID-FLIGHT','PRE-LAUNCH','POST-LAUNCH'))

sanitize_sample_names <- function(data){
  temp = data %>% mutate(namelengths = nchar(X1))
  temp = temp %>% mutate(X1 = if_else(namelengths>=3 & str_sub(X1,nchar(X1),nchar(X1))=='D',str_sub(X1,1,nchar(X1)-1),X1))
  temp = temp %>% mutate(namelengths = nchar(X2))
  temp = temp %>% mutate(X2 = if_else(namelengths>=3 & str_sub(X2,nchar(X2),nchar(X2))=='D',str_sub(X2,1,nchar(X2)-1),X2))
  return(temp %>% select(-namelengths))
}



filelist = read.delim('treedists',header=F) %>% unlist %>% unname

cutoffs = list()
sharedata = list()
for(fileloc in filelist){
	#print(fileloc)
	# reading the tsv table of pairwise distances
	nGD <- read_tsv(file = fileloc, col_names = F, show_col_types = F)
	nGD <- sanitize_sample_names(nGD)

	# Adding the metadata to the table
	nGD <- left_join(nGD %>% select(sampleID_1 = X1, everything()),
	                 meta %>% select(sampleID_1 = SeqID,
	                               subjectID_1 = Crew.ID,
	                               bodyloc1 = Body.Location,
	                               days_from_first_collection_1 = Timepoint_Numeric,timepoint1 = Timepoint_Recode,timepoint1_2 = Timepoint_Recode2,loc1_2 = location,timepoint1_3 = Timepoint))
	nGD <- left_join(nGD %>% select(sampleID_2 = X2, everything()),
	                 meta %>% select(sampleID_2 = SeqID,
	                               subjectID_2 = Crew.ID,
	                               bodyloc2 = Body.Location,
	                               days_from_first_collection_2 = Timepoint_Numeric,timepoint2 = Timepoint_Recode,timepoint2_2 = Timepoint_Recode2,loc2_2 = location,timepoint2_3 = Timepoint))

	# Computing time difference between sample (important for longitudinal samples)
	nGD$time_day_diff <- abs(nGD$days_from_first_collection_1 - nGD$days_from_first_collection_2)

	# Annotating pairs of samples. Are they related? Are they from the same individual?
	nGD$same_individual <- ifelse(nGD$subjectID_1 == nGD$subjectID_2, "same_individual", "different_individual")
	#nGD$related <- ifelse(nGD$bodyloc1 == nGD$bodyloc2, "same_swabloc", "different_swabloc")

	# Keeping only the training data
	nGD_training <- rbind(nGD %>% filter(same_individual == "same_individual") %>% group_by(subjectID_1) %>% arrange(subjectID_1, time_day_diff) %>% slice_head(n = 1) %>% ungroup(),nGD %>% filter(same_individual != "same_individual") %>% group_by(subjectID_1, subjectID_2) %>% slice_head(n = 1) %>% ungroup())

	sameindcount = table(nGD_training$same_individual) %>% data.frame %>% filter(Var1 == 'same_individual')%>% select(Freq) %>% unlist %>% unname

	if(length(sameindcount)==0){
		next
	}

	if(sameindcount>=50){
		print('here')
		res_youden <- cutpointr(data = nGD_training, x = X3, class = same_individual, method = maximize_metric, metric = youden) %>% select(optimal_cutpoint) %>% unlist %>% unname
		quantile_5pc <- nGD_training %>% select(X3) %>% quantile(0.05) %>% unlist %>% unname
		if(res_youden<quantile_5pc){
			out = res_youden
		}
		if(res_youden>=quantile_5pc){
			out = quantile_5pc
		}
	}

	if(sameindcount<50){
		out <- nGD_training %>% pull(X3) %>% quantile(0.03) %>% unlist %>% unname
	}
	cutoffs[[fileloc]] = c(strsplit(fileloc,'/') %>% map_chr(2),out)

	sharedata[[fileloc]] = nGD %>% filter(X3 < out) %>% mutate(strain = strsplit(fileloc,'/') %>% map_chr(2))
}

cutoffsout <- as.data.frame(do.call(rbind, cutoffs))
rownames(cutoffsout) = NULL

sharedata = bind_rows(sharedata)

write.table(cutoffsout,'strainphlan_sharing_cutoffs.tsv',sep='\t',quote=F)


write.table(sharedata,'i4_wgs_shared_strains.tsv',quote=F,sep='\t')




