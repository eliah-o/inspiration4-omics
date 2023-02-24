### Load i4 metabolomics data ####
library(dplyr)
library(maplet)
library(readxl)


# Preprocess data
#' Preprocessing steps for metabolomics data
#'
#' Filters out missingness, quotient normalizes, logs, scales
#' filters outliers, and knn imputes
#'
#' @param D SummarizedExperiment
#'
#' @return D with preprocessed data
#'
#' @noRd
preprocess_data <- function(D){
  
  processed_dat <- D %>% 
    
    # Filter out metabolites with over 20% missingness
    mt_pre_filter_missingness(feat_max=0.2) %>%
    
    #Filter out samples with over 10% missingness
    mt_pre_filter_missingness(samp_max=0.1) %>%
    
    #Quotient Normalize
    mt_pre_norm_quot() %>%
    
    #Log transform
    mt_pre_trans_log() %>%
    
    #Scale transform
    mt_pre_trans_scale() %>%
    
    #Outlier correction
    mt_pre_outlier_to_na() %>%
    
    # Filter out metabolites with over 20% missingness
    mt_pre_filter_missingness(feat_max=0.2) %>%
    
    #Impute NAs
    mt_pre_impute_knn()
  
  
  processed_dat
}


##### Summarized Experiment creation ####
time <- c(rep(1,4),rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4)) # Create Timepoints
subject <- rep(1:4,6) # Create subject ids

# Read C18 Hydrophobic LC data
RPPOS <- read_excel("./1024-2022 Space X plasma RPPOS-NEG extracted data.xlsx") #Read raw metabolomics data
rppos_metinfo  <- RPPOS[,c(1,3:7, 38:42),] # Pull out relevant sample information
colnames(rppos_metinfo) <-c("Name", 
                            "Confidence",
                            "Formula",
                            "Mass",
                            "RT",
                            "MODE",
                            "Compound_ID",
                            "SUPER_PATHWAY",
                            "SUB_PATHWAY",
                            "KEGG",
                            "HMDB")
rppos_data <-  RPPOS[,c(8:31)] # Pull out metabolomics information

# Read Aqueous Neutral Phase hydrophilic LC data
ANPPOS <- read_excel("./1024-2022 Space X plasma ANPPOS-NEG extracted data.xlsx")

anppos_metinfo <- ANPPOS[,c(1:3,3:7,9, 34, 35)]
colnames(anppos_metinfo)<-c("SUPER_PATHWAY", 
                            "SUB_PATHWAY",
                            "Name","Compound_ID",
                            "Confidence",
                            "Formula",
                            "Mass",
                            "RT", 
                            "MODE",
                            "KEGG",
                            "HMDB")
# Match sample info from two methods
anppos_metinfo <- anppos_metinfo[,match(names(rppos_metinfo) , names(anppos_metinfo))] 
anppos_data <-  ANPPOS[,c(10:33),] # Pull out metabolomics information
names(anppos_data) <- names(rppos_data)

# Sample information
sampinfo <- data.frame(subj_time = colnames(rppos_data), time, subject)

# Summarized Experiment of all data
all_met_SE <- SummarizedExperiment(assays = rbind(rppos_data, anppos_data), 
                                   rowData = rbind(rppos_metinfo, anppos_metinfo),
                                   colData = sampinfo)

# Run preprocessing steps
i4_D <- preprocess_data(all_met_SE)

save(i4_D, file = "./preprocessed_i4_metabolomics.RData")

