library(dplyr)
library(maplet)
library(ggplot2)
library(matrixStats)
library(limma)
library(gridExtra)

#### Functions ####

# Metabolomics boxplotting
#' i4 boxplots
#'
#' plots box plots of all significant hits in metabolomics data
#'
#' @param sig_table dataframe with metabolite ID
#' @param test Comparison name for plot pdf filename 
#'
#' @return Saves PDF file with boxplots of all significant hits
#'
#' @noRd
i4_boxplots <- function(sig_table, test){
  
  color_palette <- c("#68BBDA","#377EB8","#3E4D8A", "#E41A1C","#D9742C","#ECB03B")
  
  plot_list <- list()
  for (i in 1:nrow(sig_table)){
    met_ind <- which(rowData(i4_D)$Name==sig_table$ID[i])
    plot_df <- data.frame(met=assay(i4_D)[met_ind,], 
                          time = as.factor(colData(i4_D)$time))
    
    plot_list[[i]]<-ggplot(plot_df, aes(x=time, y=met, color=time))+
      geom_boxplot() + 
      scale_color_manual(values=color_palette) +
      ggtitle(label = rowData(i4_D)$Name[met_ind],
              subtitle =  paste0( " (Sub-pathway: ",
                                  rowData(i4_D)$SUB_PATHWAY[met_ind],
                                  ")" ))
  }
  
  pdf(height = 9, width=6, file=sprintf("%s_boxes.pdf", test))
  for(i in seq(1, length(plot_list), by=3)){
    if((i+5)<length(plot_list)){
      grid.arrange(grobs=plot_list[i:(i+2)], ncol=1)
    } else{
      grid.arrange(grobs=plot_list[i:length(plot_list)], ncol=1)
    }
  }
  dev.off()
}


#### Limma DE Analysis ####

load("./preprocessed_i4_metabolomics.RData")

# Metabolomics data
data_matrix <- assay(i4_D)
row.names(data_matrix) <- rowData(i4_D)$Name

# Add pre and post flight information
colData(i4_D)$PreVpost <- lapply(colData(i4_D)$time, 
                                 function(i) if (i<4) "Preflight" else "Postflight") %>% 
  unlist()

# Add immediate vs. long term post information
colData(i4_D)$ImmediateVLongTerm <- lapply(colData(i4_D)$time, 
                                           function(i) if (i<4) "Preflight" else if (i==4) "Immediate_Postflight" else "Longterm_Postflight") %>% 
  unlist()

########## AllPost _V_AllPre #######

# Organize sample information
experimental_design <- data.frame(label = colData(i4_D)$subj_time, condition = as.factor(colData(i4_D)$PreVpost), 
                                  replicate = colData(i4_D)$subj_time, timepoint= as.factor(colData(i4_D)$time),
                                  astronaut = as.factor(colData(i4_D)$subject))
# Order data
experimental_design <- experimental_design[order(experimental_design$label),]

# Design matrix
design_matrix <- model.matrix(~ 0 + condition + astronaut, data = experimental_design)

row.names(design_matrix) <- colnames(data_matrix)

arrayw<-arrayWeights(data_matrix, design=design_matrix)


fit <- lmFit(data_matrix, design_matrix,weights=arrayw)
contr <- makeContrasts(PostVPre= conditionPostflight - conditionPreflight, levels = design_matrix)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit)
results_table <- topTable(fit, coef="PostVPre", n=nrow(data_matrix))
write.csv(results_table, "PreVPost.csv")

########## R+1_V_AllPre #######
experimental_design <- data.frame(label = colData(i4_D)$subj_time, condition = as.factor(colData(i4_D)$ImmediateVLongTerm), 
                                  replicate = colData(i4_D)$subj_time, timepoint= as.factor(colData(i4_D)$time),
                                  astronaut = as.factor(colData(i4_D)$subject))

design_matrix <- model.matrix(~ 0 + condition+astronaut, data = experimental_design)
row.names(design_matrix) <- colnames(data_matrix)

arrayw<-arrayWeights(data_matrix,design=design_matrix)
fit <- lmFit(data_matrix, design_matrix,weights=arrayw)

contr <- makeContrasts(ImmPostVSPre=conditionImmediate_Postflight - conditionPreflight,
                       LongTermPostVSImmPost=conditionImmediate_Postflight - conditionLongterm_Postflight,
                       LongTermPostVSPre=conditionLongterm_Postflight - conditionPreflight,
                       levels = design_matrix)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit)
results <- decideTests(fit)

# Uncomment for only significant hits
results_table_imm_sig <- topTable(fit, coef="ImmPostVSPre", n = nrow(data_matrix)) #%>% dplyr::filter(adj.P.Val < 0.05)
write.csv(results_table_imm_sig,"ImmPostVSPre.csv")
#i4_boxplots(results_table_imm_sig, "ImmPostVSPre")


results_table_post_sig <- topTable(fit, coef="LongTermPostVSImmPost", n=nrow(data_matrix)) #%>% dplyr::filter(adj.P.Val < 0.05)
write.csv(results_table_post_sig,"LongTermPostVSImmPost.csv")
#i4_boxplots(results_table_post_sig, "LongTermPostVSImmPost")

results_table_long_sig <- topTable(fit, coef="LongTermPostVSPre", n=nrow(data_matrix)) #%>% dplyr::filter(adj.P.Val < 0.05)
write.csv(results_table_long_sig,"LongTermPostVSPre.csv")
#i4_boxplots(results_table_long_sig, "LongTermPostVSPre")

# R+1_V_All  (immediate post-flight vs all other time points)
colData(i4_D)$ImmediateVAll <- lapply(colData(i4_D)$time, 
                                      function(i) if (i==4) "Immediate_Postflight" else "Other") %>% 
  unlist()

experimental_design <- data.frame(label = colData(i4_D)$subj_time, condition = as.factor(colData(i4_D)$ImmediateVAll), 
                                  replicate = colData(i4_D)$subj_time, timepoint= as.factor(colData(i4_D)$time),
                                  astronaut = as.factor(colData(i4_D)$subject))

design_matrix <- model.matrix(~ 0 + condition+astronaut, data = experimental_design)
row.names(design_matrix) <- colnames(data_matrix)

arrayw<-arrayWeights(data_matrix,design=design_matrix)
fit <- lmFit(data_matrix, design_matrix,weights=arrayw)

contr <- makeContrasts(ImmPostVSAll = conditionImmediate_Postflight - conditionOther,
                       levels = design_matrix)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit)
results <- decideTests(fit)

# Uncomment for only significant hits
results_table_imm_all_sig <- topTable(fit, coef="ImmPostVSAll", n = nrow(data_matrix)) #%>% dplyr::filter(adj.P.Val < 0.05)
write.csv(results_table_imm_all_sig,"ImmPostVSAll.csv")
#i4_boxplots(results_table_imm_all_sig, "ImmPostVSAll")

