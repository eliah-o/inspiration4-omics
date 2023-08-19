# Required Libraries
#install.packages(c("readr", "writexl"))
library(readr)
library(writexl)
library(tidyverse)

setwd('~/Dropbox (Mason Lab)/i4/i4_data_packet/revisions/')

# 1. Load list of filenames
filename_list <- read_lines("mag")  

# Load dataframes from rds files
df_list <- lapply(filename_list, function(file){
  readRDS(file) %>% rownames_to_column('taxon')
})

# Associate names (for tabs) with loaded dataframes
names(df_list) <- tools::file_path_sans_ext(basename(filename_list))
names(df_list) = map(names(df_list), function(x) gsub('bacteria_xtree_005-0025_','',x))
names(df_list) = map(names(df_list), function(x) gsub('bacteria_metaphlan4_default_','',x))
names(df_list) = map(names(df_list), function(x) gsub('bac-vir-fung_kraken2_','',x))
names(df_list) = map(names(df_list), function(x) gsub('noconf','',x))
names(df_list) = map(names(df_list), function(x) gsub('conf','-0.2',x))
names(df_list) = map(names(df_list), function(x) gsub('virus_xtree_01-005_','',x))
names(df_list) = map(names(df_list), function(x) gsub('virus_phanta_default_','',x))
names(df_list) = map(names(df_list), function(x) gsub('_metatranscriptomics_decontam','META-T',x))
names(df_list) = map(names(df_list), function(x) gsub('_metagenomics_decontam','META-G',x))
names(df_list) = map(names(df_list), function(x) gsub('_metatranscriptomics_nodecontam','META-T',x))
names(df_list) = map(names(df_list), function(x) gsub('_metagenomics_nodecontam','META-G',x))

# 2. Write out to single excel spreadsheet
write_xlsx(df_list, "~/Dropbox (Mason Lab)/i4/revisions/Figures and Tables (I4WGSMTX)/revisions/SUPPTABLE_MAG_abundances.xlsx")

