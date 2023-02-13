############### LOAD LIBRARIES
library(DEP)
library(dplyr)
library(tidyverse)
set.seed(seed=0)

############### DATA PRE-PROCESSING
# Read PLASMA TABLE, AND METDADA
# Set working directory 
#setwd("..")
Protein_seer <- read.delim("./Seer_complete_data_Library_based/Library-based_PAS/Protein_Group_Panel.tsv")
metadata <- read_csv("./Seer_complete_data_Library_based/metadata_all_samples_collapsed.csv")

# Format data and collapse technical replicates (2 per sample)
df <- data.frame(Intensity = 10^(Protein_seer$Intensity..Log10.), Protein.IDs = Protein_seer$Protein.Group, name = Protein_seer$Sample.Name, Gene.names = Protein_seer$Gene.Names)
df <- df %>% pivot_wider(values_from = "Intensity", values_fill=NA)
data_unique <- make_unique(df, "Gene.names", "Protein.IDs", delim = ";")
columns_assay <- grep("C00", colnames(data_unique))
biological_reps <- unique(separate(data.frame(V1=colnames(data_unique)[columns_assay]),into=c("Sample","replicate"),col=V1, sep="-")$Sample)
data_collapsed <- data_unique[,c("name","ID")]
for(i in 1:length(biological_reps))
{
  rep1 <- paste0(biological_reps[i],"-1")
  rep2 <- paste0(biological_reps[i],"-2")
  cols <- data_unique[,c(rep1,rep2)]
  mean_cols <- rowMeans(cols,na.rm = TRUE)
  data_collapsed <- cbind(data_collapsed, mean_cols)
  colnames(data_collapsed)[i+2] <- biological_reps[i]
}
# Replace NaN with NA
data_collapsed <- data_collapsed %>% mutate_all(~ifelse(is.nan(.), NA, .))

# Filter based on coefficient of variation (CV) of all timepoints
Average_Sept_Post <- rowMeans(as.matrix(data_collapsed[,grep("Sept_Post", colnames(data_collapsed))]), na.rm=TRUE)
SD_Sept_Post <- rowSds(as.matrix(data_collapsed[,grep("Sept_Post", colnames(data_collapsed))]), na.rm=TRUE)
CV_Sept_Post <- SD_Sept_Post/Average_Sept_Post
Average_Sept_Pre <- rowMeans(as.matrix(data_collapsed[,grep("Sept_Pre", colnames(data_collapsed))]), na.rm=TRUE)
SD_Sept_Pre <- rowSds(as.matrix(data_collapsed[,grep("Sept_Pre", colnames(data_collapsed))]), na.rm=TRUE)
CV_Sept_Pre <- SD_Sept_Pre/Average_Sept_Pre
Average_Nov_Post <- rowMeans(as.matrix(data_collapsed[,grep("Nov_Post", colnames(data_collapsed))]), na.rm=TRUE)
SD_Nov_Post <- rowSds(as.matrix(data_collapsed[,grep("Nov_Post", colnames(data_collapsed))]), na.rm=TRUE)
CV_Nov_Post <- SD_Nov_Post/Average_Nov_Post
Average_Aug_Pre <- rowMeans(as.matrix(data_collapsed[,grep("Aug_Pre", colnames(data_collapsed))]), na.rm=TRUE)
SD_Aug_Pre <- rowSds(as.matrix(data_collapsed[,grep("Aug_Pre", colnames(data_collapsed))]), na.rm=TRUE)
CV_Aug_Pre <- SD_Aug_Pre/Average_Aug_Pre
Average_Dec_Post <- rowMeans(as.matrix(data_collapsed[,grep("Dec_Post", colnames(data_collapsed))]), na.rm=TRUE)
SD_Dec_Post <- rowSds(as.matrix(data_collapsed[,grep("Dec_Post", colnames(data_collapsed))]), na.rm=TRUE)
CV_Dec_Post <- SD_Dec_Post/Average_Dec_Post
Average_June_Pre <- rowMeans(as.matrix(data_collapsed[,grep("June_Pre", colnames(data_collapsed))]), na.rm=TRUE)
SD_June_Pre <- rowSds(as.matrix(data_collapsed[,grep("June_Pre", colnames(data_collapsed))]), na.rm=TRUE)
CV_June_Pre <- SD_June_Pre/Average_June_Pre

data_collapsed <- data_collapsed[(CV_Sept_Post <0.5) | (CV_Sept_Pre <0.5) | (CV_Nov_Post <0.5) | (CV_Aug_Pre <0.5) | (CV_Dec_Post <0.5) | (CV_June_Pre <0.5),]
data_collapsed <- data_collapsed[!is.na(data_collapsed$ID),]
write.csv(data_collapsed, file="./Seer_complete_data_Library_based/CV_filt_collapsed_data.csv")

#Read CV filtered file
data_collapsed <- read_csv("./Seer_complete_data_Library_based/CV_filt_collapsed_data.csv")

# Create Post_V_Pre SE object
experimental_design <- data.frame(label = metadata$Sample_name, condition = metadata$PreVpost, replicate = metadata$Sample_name, timepoint= metadata$Timepoint)
experimental_design <- experimental_design[order(experimental_design$label),]
data_ordered_collapsed <- data_collapsed[,order(colnames(data_collapsed))]
columns_assay <- grep("C00", colnames(data_ordered_collapsed))
data_se_collapsed <- make_se(data_ordered_collapsed, columns_assay, experimental_design)
save(data_se_collapsed, file = "./RData_objects/plasma_postVpre_collapsed_CV_filt.RData")

# Create R+1_v_Pre SE object 
experimental_design <- data.frame(label = metadata$Sample_name, condition = metadata$ImmediateVLongTerm, replicate = metadata$Sample_name)
experimental_design <- experimental_design[order(experimental_design$label),]
data_ordered_collapsed <- data_collapsed[,order(colnames(data_collapsed))]
columns_assay <- grep("C00", colnames(data_ordered_collapsed))
data_se_collapsed <- make_se(data_ordered_collapsed, columns_assay, experimental_design)
save(data_se_collapsed, file = "./RData_objects/plasma_immPostVpre_collapsed_CV_filt.RData")

# Filter based on number of NAs
load("./final_code/RData_objects/plasma_postVpre_collapsed_CV_filt.RData")
# Filter for proteins that are identified in all replicates of at least one condition
data_filt0 <- filter_missval(data_se, thr = 0)
save(data_filt0, file = "./RData_objects/plasma_postVpre_filt_collapsed_CV_filt.RData")
load("./final_code/RData_objects/plasma_immPostVpre_collapsed_CV_filt.RData")
# Filter for proteins that are identified in all replicates of at least one condition
data_filt0 <- filter_missval(data_se, thr = 0)
save(data_filt0, file = "./RData_objects/plasma_immPostVpre_filt_collapsed_CV_filt.RData")

# Normalization and Imputation
load("./final_code/RData_objects/plasma_postVpre_filt_collapsed_CV_filt.RData")
data_filt <- data_filt0
# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Impute with MinProb
MinProb_imputation <- impute(data_norm, fun = "MinProb", q = 0.01)
save(MinProb_imputation, file = "./RData_objects/plasma_postVpre_norm_imputed.RData")
# Normalization and Imputation
load("./final_code/RData_objects/plasma_immPostVpre_filt_collapsed_CV_filt.RData")
data_filt <- data_filt0
# Normalize the data
data_norm <- normalize_vsn(data_filt)
# Impute with MinProb
MinProb_imputation <- impute(data_norm, fun = "MinProb", q = 0.01)
save(MinProb_imputation, file = "./RData_objects/plasma_immPostVpre_norm_imputed.RData")

######### DE ANALYSIS
# AllPost_V_AllPre
load("./final_code/RData_objects/plasma_postVpre_norm_imputed.RData")
metadata <- read_csv("./Seer_complete_data_Library_based/metadata_all_samples_collapsed.csv")
experimental_design <- data.frame(label = metadata$Sample_name, condition = as.factor(metadata$PreVpost), 
                                  replicate = metadata$Sample_name, timepoint= metadata$Timepoint,
                                  astronaut = metadata$Astronaut)
experimental_design <- experimental_design[order(experimental_design$label),]
design_matrix <- model.matrix(~ 0 + condition + astronaut, data = experimental_design)
df_wide <- get_df_wide(MinProb_imputation)
data_matrix <- df_wide[2:22]
row.names(data_matrix) <- df_wide$name
row.names(design_matrix) <- colnames(data_matrix)
arrayw<-arrayWeights(data_matrix,design=design_matrix)
fit <- lmFit(data_matrix, design_matrix,weights=arrayw)
contr <- makeContrasts(PostVPre= conditionPostflight - conditionPreflight, levels = design_matrix)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit)
results_table <- topTable(fit, coef="PostVPre", n=2000)
results_table$Gene <- rownames(results_table)
setwd("./final_code/results/")
write_csv(results_table, "./plasma_postVpre_limma_DE_results_CV_filt_arrayWeighted.csv")

# R+1_V_AllPre
load("./final_code/RData_objects/plasma_immPostVpre_norm_imputed.RData")
metadata <- read_csv("./Seer_complete_data_Library_based/metadata_all_samples_collapsed.csv")
experimental_design <- data.frame(label = metadata$Sample_name, condition = as.factor(metadata$ImmediateVLongTerm), 
                                  replicate = metadata$Sample_name, timepoint= metadata$Timepoint,
                                  astronaut = metadata$Astronaut)
experimental_design <- experimental_design[order(experimental_design$label),]
design_matrix <- model.matrix(~ 0 + condition+astronaut, data = experimental_design)
df_wide <- get_df_wide(MinProb_imputation)
write_csv(df_wide, "./final_code/plasma_normalized_imputed_data.csv")
data_matrix <- df_wide[2:22]
row.names(data_matrix) <- df_wide$name
row.names(design_matrix) <- colnames(data_matrix)
arrayw<-arrayWeights(data_matrix,design=design_matrix)
fit <- lmFit(data_matrix, design_matrix,weights=arrayw)
contr <- makeContrasts(ImmPostVSPre=conditionImmediate_Postflight - conditionPreflight,
                       levels = design_matrix)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit)
results <- decideTests(fit)
results_table <- topTable(fit, coef="ImmPostVSPre", n=2000)
results_table$Gene <- rownames(results_table)
write_csv(results_table,"./final_code/results/plasma_immPostVpre_limma_DE_analysis_arrayWeighted.csv")
