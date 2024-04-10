require(tidyverse)

individual_timepoints <- function(timepoint) {
    # If the variant was called at this time point, then we should use the VAF that the caller spits out
    if (timepoint["was_called"]) {
        return("called")
        # If it wasn't called, then we have to check if the pileup found it for this time point
    } else if (timepoint["has_pileup"]) {
        # If the variant passes the fisher's exact test, then we can rescue it
        if (as.logical(timepoint["fisher_passed"])) {
            return("rescued")
            # If the variant does not pass. The pileup still found it, but we cannot assume the pileup is right, we have to assume the limit of detection
        } else {
            return("at_limit")
        }
        # If it wasn't called and pileup did not detect it, then this variant for this time point is missing
    } else {
        return("missing")
    }
}

set_vaf <- function(df) {
    # df <- df %>% tidyr::separate(information, into = c("lenient_pass", "average_AF", "was_called", "has_pileup", "fisher_passed", "pileup_vaf"), sep = "\\|")
    vaf <- ifelse(df["lenient_pass"], df["average_AF"], NA)
    vaf_ignore_noise <- ifelse(df["lenient_pass"], df["average_AF"], NA)
    timepoint <- ifelse(is.na(vaf), individual_timepoints(df), "lenient_pass")
    vaf <- ifelse(timepoint == "called", df["average_AF"], vaf)
    vaf_ignore_noise <- ifelse(timepoint == "called", df["average_AF"], vaf_ignore_noise)
    vaf_ignore_noise <- ifelse(!is.na(df["PILEUP_TOTAL_DEPTH"]) & (timepoint == "at_limit" | timepoint == "missing"), 0.5 / as.numeric(df["PILEUP_TOTAL_DEPTH"]), vaf_ignore_noise)
    return(c(vaf, vaf_ignore_noise, timepoint))
}

time_points <- function(variant, ignore_noise = TRUE, rescue = TRUE) {
    variant <- variant %>%
        mutate(
            lenient_pass = case_when(
                as.logical(pon_FP_pass_XGB) &
                    as.logical(long100_indel_pass_XGB) & as.logical(long_indel_pass_XGB) & as.logical(di_tri_vard_pass_XGB) & as.logical(bcbio_pass_XGB) & as.logical(zscore_pass_XGB) &
                    as.logical(PASS_BY_1) &
                    FP_Filter_DETP20_XGB == 0 & FP_Filter_MMQS100_XGB == 0 & FP_Filter_MMQSD50_XGB == 0 & FP_Filter_NRC_XGB == 0 & FP_Filter_PB10_XGB == 0 & FP_Filter_RLD25_XGB == 0 ~ TRUE,
                TRUE ~ FALSE
            )
        ) %>%
        mutate(lenient_pass = ifelse(FP_Filter_MVC4_XGB == 1 & PILEUP_ALT_DEPTH > 5, FALSE, lenient_pass)) %>%
        mutate(lenient_pass = ifelse(FP_Filter_SB1_XGB == 1 & PILEUP_ALT_DEPTH > 10, FALSE, lenient_pass)) %>%
        mutate(lenient_pass = ifelse(is.na(lenient_pass), FALSE, lenient_pass))

    # Check to see if the new definition creates more "PASS" all time points
    # If the variant was CALLED + PASS in ALL time points
    if (sum(as.integer(as.logical(variant$putative_driver))) == length(variant$putative_driver)) {
        all_time_points <- "true"
        # If the variant was CALLED + Lenient PASS in ALL time points
        # This will handle variants that only failed because of LowVAF or SB1/MVC4
    } else if (sum(as.integer(variant$lenient_pass)) == length(variant$lenient_pass)) {
        all_time_points <- "true"
    } else {
        all_time_points <- "false"
    }
    variant <- variant %>% select(lenient_pass, average_AF, was_called, has_pileup, fisher_passed, PILEUP_VAF, PILEUP_TOTAL_DEPTH, whereincycle)
    
    variant <- data.frame(whereincycle = 1:6) %>%
        left_join(variant, by = "whereincycle") %>%
        mutate(
            lenient_pass = ifelse(is.na(lenient_pass), FALSE, lenient_pass),
            was_called = ifelse(is.na(was_called), FALSE, was_called),
            has_pileup = ifelse(is.na(has_pileup), FALSE, has_pileup),
            fisher_passed = ifelse(is.na(fisher_passed), FALSE, fisher_passed)
        )

    res <- as.data.frame(t(apply(variant, 1, set_vaf)))
    names(res) <- c("vaf", "vaf_ignore_noise", "timepoint")
    variant <- cbind(variant, res)

    if ("rescued" %in% variant$timepoint) {
        # If one of the variants was rescued, we should use the pileup IF it has pileup
        variant <- variant %>% 
            mutate(
                vaf_rescue = ifelse(has_pileup, PILEUP_VAF, vaf),
                vaf_ignore_noise_rescue = ifelse(has_pileup, PILEUP_VAF, vaf_ignore_noise)
            )
        rescued <- "true"
    } else {
        variant <- variant %>% 
            mutate(
                vaf_rescue = vaf,
                vaf_ignore_noise_rescue = vaf_ignore_noise
            )
        rescued <- "false"
    }

    missing <- ifelse("missing" %in% variant$timepoint, "true", "false")
    at_limit <- ifelse("at_limit" %in% variant$timepoint, "true", "false")
    called <- ifelse("called" %in% variant$timepoint, "true", "false")

    res <- as.data.frame(t(variant %>% arrange(whereincycle) %>% select(vaf, vaf_ignore_noise, vaf_rescue, vaf_ignore_noise_rescue)))
    return_string <- paste(res[1, ], collapse = " ")
    if (ignore_noise) {
        return_string <- paste(return_string, paste(res[2, ], collapse = " "), sep = " ")
    }
    if (rescue) {
        return_string <- paste(return_string, paste(res[3, ], collapse = " "), sep = " ")
    }
    res <- paste(all_time_points, return_string, rescued, missing, called, at_limit, sep = " ")

    return(res)
}
