############### LOAD LIBRARIES
library(DEP)
library(dplyr)
library(tidyverse)
set.seed(seed=0)

############### DATA PRE-PROCESSING
# Read EXOSOMES TABLE, AND METDADA
# Set working directory 
#setwd("..")
Protein_exosomes <- read.delim("./Exosomes_i4/EVPs_proteomics_preprocessed_data.txt")
metadata <- read_csv("./Exosomes_i4/EVPs_sample_metadata.csv")

# Extract gene name "GN=.."
temp <- as.data.frame(Protein_exosomes) %>% tidyr::separate(col=Description, into=c("V1","V2"),sep="GN=" )
temp <- as.data.frame(temp) %>% tidyr::separate(col=V2, into=c("gene","other"),sep=" " )
temp$Gene.names <- ifelse(is.na(temp$gene), temp$Accession, temp$gene)
# Adding gene names
df <- cbind(Protein_exosomes,Gene.names=temp$Gene.names)
df[df == 0] <- NA

# Filter based on coefficient of variation (CV) of all timepoints
Average_Sept_Post <- rowMeans(as.matrix(df[,grep("_4", colnames(df))]), na.rm=TRUE)
SD_Sept_Post <- rowSds(as.matrix(df[,grep("_4", colnames(df))]), na.rm=TRUE)
CV_Sept_Post <- SD_Sept_Post/Average_Sept_Post
Average_Sept_Pre <- rowMeans(as.matrix(df[,grep("_3", colnames(df))]), na.rm=TRUE)
SD_Sept_Pre <- rowSds(as.matrix(df[,grep("_3", colnames(df))]), na.rm=TRUE)
CV_Sept_Pre <- SD_Sept_Pre/Average_Sept_Pre
Average_Nov_Post <- rowMeans(as.matrix(df[,grep("_5", colnames(df))]), na.rm=TRUE)
SD_Nov_Post <- rowSds(as.matrix(df[,grep("_5", colnames(df))]), na.rm=TRUE)
CV_Nov_Post <- SD_Nov_Post/Average_Nov_Post
Average_Aug_Pre <- rowMeans(as.matrix(df[,grep("_2", colnames(df))]), na.rm=TRUE)
SD_Aug_Pre <- rowSds(as.matrix(df[,grep("_2", colnames(df))]), na.rm=TRUE)
CV_Aug_Pre <- SD_Aug_Pre/Average_Aug_Pre
Average_Dec_Post <- rowMeans(as.matrix(df[,grep("_6", colnames(df))]), na.rm=TRUE)
SD_Dec_Post <- rowSds(as.matrix(df[,grep("_6", colnames(df))]), na.rm=TRUE)
CV_Dec_Post <- SD_Dec_Post/Average_Dec_Post
Average_June_Pre <- rowMeans(as.matrix(df[,grep("_1", colnames(df))]), na.rm=TRUE)
SD_June_Pre <- rowSds(as.matrix(df[,grep("_1", colnames(df))]), na.rm=TRUE)
CV_June_Pre <- SD_June_Pre/Average_June_Pre
df[is.na(df)] <- 0
df <- df[(CV_Sept_Post <0.5) | (CV_Sept_Pre <0.5) | (CV_Nov_Post <0.5) | (CV_Aug_Pre <0.5) | (CV_Dec_Post <0.5) | (CV_June_Pre <0.5),]
df <- df[!is.na(df$Accession),]
write.csv(df, file="./Exosomes_i4/exosomes_CV_filt_data.csv")

# Read CV filtered file
df <- read_csv("./Exosomes_i4/exosomes_CV_filt_data.csv", 
               col_types = cols(...1 = col_skip()))
metadata <- read_csv("./Exosomes_i4/Area_sample_code.csv")
# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(df, "Gene.names", "Accession", delim = ";")
# Format metadata
columns_assay <- grep("C00", colnames(data_unique))
metadata <- metadata[1:24,]
# Create Post_V_Pre SE object
experimental_design <- data.frame(label = paste0(metadata$Sample,"_",metadata$PreVpost), condition =metadata$PreVpost, replicate = metadata$Sample)
colnames(data_unique) <- c(colnames(data_unique)[1:2],experimental_design$label, colnames(data_unique)[27:length(colnames(data_unique))])
data_se <- make_se(data_unique, columns_assay, experimental_design)
save(data_se, file = "./RData_objects/exosomes_postVpre_CV_filt.RData")

# Create R+1_v_Pre SE object 
experimental_design <- data.frame(label = paste0(metadata$Sample,"_",metadata$ImmediateVLongTerm), condition =metadata$ImmediateVLongTerm, replicate = metadata$Sample)
colnames(data_unique) <- c(colnames(data_unique)[1:2],experimental_design$label, colnames(data_unique)[27:length(colnames(data_unique))])
data_se <- make_se(data_unique, columns_assay, experimental_design)
save(data_se, file = "./RData_objects/exosomes_immPostVpre_CV_filt.RData")

# Filter based on number of NAs
load("./final_code/RData_objects/exosomes_postVpre_CV_filt.RData")
# Filter for proteins that are identified in all replicates of at least one condition
data_filt0 <- filter_missval(data_se, thr = 0)
save(data_filt0, file = "./RData_objects/exosomes_postVpre_filt_CV_filt.RData")
load("./final_code/RData_objects/exosomes_immPostVpre_CV_filt.RData")
# Filter for proteins that are identified in all replicates of at least one condition
data_filt0 <- filter_missval(data_se, thr = 0)
save(data_filt0, file = "./RData_objects/exosomes_immPostVpre_filt_CV_filt.RData")

# Normalization and Imputation
load("./final_code/RData_objects/exosomes_postVpre_filt_CV_filt.RData")
data_filt <- data_filt0
# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Impute with MinProb
MinProb_imputation <- impute(data_norm, fun = "MinProb", q = 0.01)
save(MinProb_imputation, file = "./RData_objects/exosomes_postVpre_norm_imputed.RData")
# Normalization and Imputation
load("./final_code/RData_objects/exosomes_immPostVpre_filt_CV_filt.RData")
data_filt <- data_filt0
# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Impute with MinProb
MinProb_imputation <- impute(data_norm, fun = "MinProb", q = 0.01)
save(MinProb_imputation, file = "./RData_objects/exosomes_immPostVpre_norm_imputed.RData")


######### DE ANALYSIS
# AllPost_V_AllPre
load("./final_code/RData_objects/exosomes_postVpre_norm_imputed.RData")
metadata <- read_csv("./Exosomes_i4/Area_sample_code.csv")
metadata <- metadata[1:24,]
metadata$PreVpost <- stringr::str_replace(metadata$PreVpost,"-","_")
experimental_design <- data.frame(label = paste0(metadata$Sample,"_",metadata$PreVpost), condition =metadata$PreVpost, replicate = metadata$Astronaut, timepoint = metadata$Timepoint)
design_matrix <- model.matrix(~ 0+ condition + replicate, data = experimental_design)
df_wide <- get_df_wide(MinProb_imputation)
data_matrix <- df_wide[2:25]
row.names(data_matrix) <- df_wide$name
row.names(design_matrix) <- colnames(data_matrix)
data_matrix <- as.matrix(data_matrix)
arrayw<-arrayWeights(data_matrix,design=design_matrix)
fit <- lmFit(data_matrix, design_matrix,weights=arrayw)
contr <- makeContrasts(PostVPre=conditionPost_flight - conditionPre_flight, levels = design_matrix)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit)
results_table <- topTable(fit, coef="PostVPre", n=1000)
View(results_table)
results_table$Gene <- rownames(results_table)
write.csv(results_table,"./final_code/results/exosomes_postVpre_limma_arrayWeighted.csv")


# R+1_V_AllPre
load("./final_code/RData_objects/exosomes_immPostVpre_norm_imputed.RData")
metadata <- read_csv("./Exosomes_i4/Area_sample_code.csv")
metadata <- metadata[1:24,]
metadata$ImmediateVLongTerm <- stringr::str_replace(metadata$ImmediateVLongTerm," ","_")
experimental_design <- data.frame(label = paste0(metadata$Sample,"_",metadata$PreVpost), condition =metadata$ImmediateVLongTerm, replicate = metadata$Astronaut, timepoint = metadata$Timepoint)
design_matrix <- model.matrix(~ 0+ condition + replicate, data = experimental_design)
df_wide <- get_df_wide(MinProb_imputation)
write_csv(df_wide, "./final_code/exosomes_normalized_imputed_data.csv")
data_matrix <- df_wide[2:25]
row.names(data_matrix) <- df_wide$name
row.names(design_matrix) <- colnames(data_matrix)
arrayw<-arrayWeights(data_matrix,design=design_matrix)
fit <- lmFit(data_matrix, design_matrix,weights=arrayw)
contr <- makeContrasts(ImmPostVSPre=conditionImmediate_Postflight - conditionPreflight,
                       levels = design_matrix)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit)
results <- decideTests(fit)
results_table <- topTable(fit, coef="ImmPostVSPre", n=1000)
results_table$Gene <- rownames(results_table)
write_csv(results_table,"./final_code/results/exosomes_limma_DE_immPostVpre_arrayWeighted.csv")

