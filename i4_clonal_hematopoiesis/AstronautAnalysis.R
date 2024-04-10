library(tidyverse)
library(data.table)
library(vcfR)

astronaut_data <- "data/variantplex_astronaut.csv"
whereincycle <- "data/astronaut_whereincycle_key.csv"
astronaut_pileup_data <- "data/samples.pileup.vcf.gz"
pon_pileup_data <- "data/pon.pileup.vcf.gz"

astronaut <- fread(astronaut_data)
whereincycle <- fread(whereincycle)

dummy.INFO.line <- "##INFO=<ID=AF>"
astronaut_pileup <- read.vcfR(astronaut_pileup_data)
astronaut_pileup@meta[length(astronaut_pileup@meta) + 1] <- dummy.INFO.line
astronaut_pileup <- vcfR2tidy(astronaut_pileup, single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat

pon_pileup <- read.vcfR(pon_pileup_data)
pon_pileup@meta[length(pon_pileup@meta) + 1] <- dummy.INFO.line
pon_pileup <- vcfR2tidy(pon_pileup, single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat

pon_pileup <- pon_pileup %>% 
  mutate(
    key = paste0(CHROM, " ", POS, " ", REF, ">", ALT),
    PON_RefDepth = as.numeric(PON_RefDepth),
    PON_AltDepth = as.numeric(PON_AltDepth),
  ) %>%
  distinct(key, PON_RefDepth, PON_AltDepth)

# Calculate p-value for all 
astronaut_pileup <- astronaut_pileup %>%
  mutate(key = paste0(CHROM, " ", POS, " ", REF, ">", ALT)) %>%
  left_join(
    pon_pileup
  )

astronaut_pileup$PON_Pvalue <- apply(astronaut_pileup %>% select(PON_RefDepth, PON_AltDepth, gt_RD, gt_AD), 1, function(x) {
  if ((x[1]+x[2]==0) | (x[3]+x[4]==0)){
    return(0)
  } else if (x[2]==0 & x[1]!=0) {
    return(0)
  } else if ((x[1]==0 & x[2]!=0) & (x[3]==0 & x[4]!=0)) {
    return(1)
  } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+x[4])) {
    return(1)
  } else {
    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
  }
})

# Merge existing calls into dataframe
astronaut_pileup <- astronaut_pileup %>%
  mutate(
    sample_key = paste0(Indiv, " ", key)
  )

astronaut_pileup <- astronaut_pileup %>%
  rename(
    PILEUP_ALT_DEPTH = gt_AD,
    PILEUP_REF_DEPTH = gt_RD,
    PILEUP_TOTAL_DEPTH = gt_DP,
    PILEUP_VAF = gt_VF
  ) %>%
  select(-c(gt_DPP, gt_DPN, gt_RDP, gt_RDN, gt_ADP, gt_ADN, gt_DPF, gt_RDF, gt_ADF, gt_GT_alleles)) %>%
  left_join(
    astronaut %>% select(sample_key, putative_driver, pd_reason, 
                         FILTER_Mutect, FILTER_Lofreq, FILTER_Vardict, 
                         FP_Filter, FP_Filter_PASS_XGB, FP_Filter_RLD25_XGB, FP_Filter_DETP20_XGB, FP_Filter_MVC4_XGB, FP_Filter_SB1_XGB, FP_Filter_MMQSD50_XGB, FP_Filter_NRC_XGB, FP_Filter_PB10_XGB, FP_Filter_MMQS100_XGB,
                         average_AF, average_AD,
                         pon_FP_pass_XGB, low_AF_pass_XGB, long_indel_pass_XGB, bcbio_pass_XGB, zscore_pass_XGB, all_fp_pass_XGB, di_tri_vard_pass_XGB, long100_indel_pass_XGB, PASS_BY_1, PASS_BY_2, PASS_BY_3
                         )
  ) %>%
  mutate(
    fisher_passed = ifelse(PON_Pvalue <= 2.114164905e-6, TRUE, FALSE),
    was_called = ifelse(is.na(PASS_BY_1), FALSE, TRUE),
    has_pileup = ifelse(PILEUP_ALT_DEPTH != 0, TRUE, FALSE),
    putative_driver = ifelse(is.na(putative_driver), 0, putative_driver)
  )

astronaut_pileup <- astronaut_pileup %>%
  left_join(
    astronaut %>% select(key, Gene) %>% distinct
  ) %>%
  mutate(Gene = ifelse(key == 'chr2 25244297 C>T', "DNMT3A", Gene))

astronaut_pileup <- astronaut_pileup %>%
  separate(sample_key, into = c("SAMPLE"), sep = " ", extra = "drop") %>%
  left_join(
    whereincycle
  )

astronaut <- astronaut %>%
  separate(sample_key, into = c("SAMPLE"), sep = " ", extra = "drop") %>%
  left_join(
    whereincycle
  )

source("utils/timepoints_function.R")

# Define Lenient Pass
astronaut_pileup <- astronaut_pileup %>% 
  mutate(lenient_pass = case_when(
    as.logical(pon_FP_pass_XGB) & 
      as.logical(long100_indel_pass_XGB) & 
      as.logical(long_indel_pass_XGB) & 
      as.logical(di_tri_vard_pass_XGB) & 
      as.logical(bcbio_pass_XGB) & 
      as.logical(zscore_pass_XGB) &
      as.logical(PASS_BY_1) & 
      FP_Filter_DETP20_XGB == 0 & 
      FP_Filter_MMQS100_XGB == 0 &
      FP_Filter_MMQSD50_XGB == 0 &
      FP_Filter_NRC_XGB == 0 &
      FP_Filter_PB10_XGB == 0 &
      FP_Filter_RLD25_XGB == 0 ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(lenient_pass = ifelse(FP_Filter_MVC4_XGB == 1 & PILEUP_ALT_DEPTH > 5, FALSE, lenient_pass)) %>%
  mutate(lenient_pass = ifelse(FP_Filter_SB1_XGB == 1 & PILEUP_ALT_DEPTH > 10, FALSE, lenient_pass)) %>%
  mutate(lenient_pass = ifelse(is.na(lenient_pass), FALSE, lenient_pass))

D_ <- astronaut_pileup %>% group_by(Individual, key, Gene) %>% summarise(result = time_points(variant = across(everything()), ignore_noise = TRUE, rescue = FALSE), in_cohort = ifelse(sum(as.integer(as.logical(putative_driver))) == 0, FALSE, TRUE))
D_ <- D_ %>% separate(result, c("all_time_points", "March", "August", "September_Pre_Flight", "September_Post_Flight", "November", "December", "March_ignore_noise", "August_ignore_noise", "September_Pre_Flight_ignore_noise", "September_Post_Flight_ignore_noise", "November_ignore_noise", "December_ignore_noise", "rescued", "missing", "called", "at_limit"), sep = " ", extra = "merge", fill = "right")
D_ <- D_ %>% mutate_if(is.character, ~na_if(., "NA"))
D_ <- D_ %>% mutate(March = ifelse(is.na(March), March_ignore_noise, March),
                    August = ifelse(is.na(August), August_ignore_noise, August),
                    September_Pre_Flight = ifelse(is.na(September_Pre_Flight), September_Pre_Flight_ignore_noise, September_Pre_Flight),
                    September_Post_Flight = ifelse(is.na(September_Post_Flight), September_Post_Flight_ignore_noise, September_Post_Flight),
                    November = ifelse(is.na(November), November_ignore_noise, November),
                    December = ifelse(is.na(December), December_ignore_noise, December)) %>%
  mutate(
    March = as.numeric(March),
    August = as.numeric(August),
    September_Pre_Flight = as.numeric(September_Pre_Flight),
    September_Post_Flight = as.numeric(September_Post_Flight),
    November = as.numeric(November),
    December = as.numeric(December)
  )

# List of Variants in M
astronaut_variants <- unique(astronaut %>% filter(putative_driver == 1, !is.na(key)) %>% dplyr::select(Individual, key)) %>% mutate(ind_key = paste(Individual, key, sep = "-")) %>% dplyr::select(ind_key)
# List of Variants in Pileup Results
astronaut_pileup_variants <- D_ %>% ungroup() %>% mutate(ind_key = paste(Individual, key, sep = " ")) %>% dplyr::select(ind_key)

# List of Variants that are shared between M and the PR
M_intersect_PR <- dplyr::intersect(unique(astronaut %>% filter(!is.na(key)) %>% dplyr::select(Individual, key)) %>% mutate(ind_key = paste(Individual, key, sep = " ")) %>% dplyr::select(ind_key), 
                                   D_ %>% ungroup() %>% mutate(ind_key = paste(Individual, key, sep = " ")) %>% dplyr::select(ind_key))

# Final list of D we should be using, because these would be "REAL" variants
D <- left_join(M_intersect_PR, D_ %>% mutate(ind_key = paste(Individual, key, sep = " "))) %>% filter(in_cohort == TRUE) %>%
  rbind(D_ %>% mutate(ind_key = paste(Individual, key, sep = " ")) %>% filter(ind_key == '3 chr2 25244297 C>T')) %>%
  mutate(in_cohort = TRUE)

controls <- fread("data/NoTreatment_Timepoints.csv")

setdiff(colnames(D), colnames(controls))

# The specific PD variants we are interested in exploring
ind_key_pd <- c("4 chr2 25234347 G>C", "4 chr2 25235778 C>G", "3 chr2 25241591 C>A", "4 chr2 25246732 GTCGTGGCACACCGGGAACAGCTTCCCCGC>G", "4 chr2 25247647 G>A", "3 chr2 25244297 C>T")

to_plot <- rbind(
  D %>% 
    mutate(
      preVAF = August,
      postVAF = November,
      change_in_VAF = postVAF - preVAF,
      change_in_days = 92,        # There are approximately 92 days between August and November
      fromwhere = "Astronauts",
      change_in_VAF_months = change_in_VAF / 4
    ) %>%
    filter(ind_key %in% ind_key_pd) %>%
    select(change_in_VAF, change_in_days, change_in_VAF_months, fromwhere),
  controls %>%
    select(change_in_VAF, change_in_days, change_in_VAF_months, fromwhere)
)

spider_VAF_months <- to_plot %>%
  filter(
    !is.na(change_in_VAF_months),
    change_in_days < 150
  ) %>%
  ggplot(
    aes(
      x = change_in_days,
      y = change_in_VAF,
      fill = fromwhere
    )
  ) +
  geom_segment(
    aes(
      x = change_in_days,
      y = change_in_VAF,
      xend = 0,
      yend = 0,
      color = fromwhere
    )
  ) +
  geom_point(
    aes(
      x = change_in_days,
      y = change_in_VAF,
      color = fromwhere
    )
  ) +
  labs(
    x = "Days",
    y = "Change in VAF"
  ) +
  theme_bw() +
  theme(
    strip.text = element_blank()
  ) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(-0.025, 0.025)) +
  scale_fill_manual(values = c("#0072B2", "#CC0000")) +
  scale_color_manual(values = c("#0072B2", "#CC0000")) +
  facet_wrap(~fromwhere)
spider_VAF_months
