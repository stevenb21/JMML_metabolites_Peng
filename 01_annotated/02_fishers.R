# ==============================================================================
# Script Name: fishers.R          (refactor 2025-06-20)
# Author:     Steven Brooks
# Purpose:    Fisher’s exact test on compounds with 20–65 % missingness
# ==============================================================================

# ---- Setup -------------------------------------------------------------------
library(tidyverse)
library(writexl)

# ---- IO helpers --------------------------------------------------------------
run_fisher <- function(data,
                       dataset_name,
                       min_miss = 0.20,
                       max_miss = 0.65,
                       out_path = "../res/annotated") {
  
  message("» Running Fisher tests for ", dataset_name)
  
  # ensure output directory exists
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }
  
  ## 1. keep JMML & Control and tag missingness -------------------------------
  miss_tbl <- data %>%
    filter(sample_group %in% c("JMML", "Control")) %>%
    mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))
  
  ## 2. filter compounds in given missingness window --------------------------
  keep_vec <- miss_tbl %>%
    group_by(compound_name) %>%
    summarise(prop_missing = mean(missing == "NA"), .groups = "drop") %>%
    filter(between(prop_missing, min_miss, max_miss)) %>%
    pull(compound_name)
  
  miss_tbl <- miss_tbl %>% filter(compound_name %in% keep_vec)
  
  ## 3. count → wide (JMML_NA etc.) -------------------------------------------
  fisher_df <- miss_tbl %>%
    count(compound_name, sample_group, missing) %>%
    pivot_wider(
      names_from  = c(sample_group, missing),
      values_from = n,
      values_fill = 0,
      names_sep   = "_") %>%
    # ensure columns exist even if all 0
    rowwise() %>%
    mutate(
      ## 4. define a, b, c, d ---------------------------------------------------
      a_Control_NA = Control_NA,
      b_Control_Present = Control_Present,
      c_JMML_NA = JMML_NA,
      d_JMML_Present = JMML_Present,
      
      ## 5. Fisher test --------------------------------------------------------
      fisher_result = list(fisher.test(matrix(c(a_Control_NA,
                                                b_Control_Present,
                                                c_JMML_NA,
                                                d_JMML_Present),
                                              nrow = 2, byrow = TRUE,
                                              dimnames = list(
                                                Group  = c("Control","JMML"),
                                                Status = c("Missing","Present")
                                              )))),
      odds_ratio    = fisher_result$estimate[[1]],
      ci_lower      = fisher_result$conf.int[1],
      ci_upper      = fisher_result$conf.int[2],
      fisher_pval   = fisher_result$p.value,
      log2_OR       = log2(odds_ratio)
    ) %>%
    ungroup() %>%
    select(compound_name, a_Control_NA,
           b_Control_Present,
           c_JMML_NA,
           d_JMML_Present,
           odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval)
  
  ## 6. save -------------------------------------------------------------------
  out_file <- file.path(out_path,
                        paste0(dataset_name, "_fisher.xlsx"))
  writexl::write_xlsx(fisher_df, out_file)
  message("✓ Results written to: ", out_file)
  
  invisible(fisher_df)
}

# ---- Load data ---------------------------------------------------------------
ZHP   <- readRDS("../data/R_objects/ZHP_annotated_filtered.rds")
RPNPF <- readRDS("../data/R_objects/RPNPF_annotated_filtered.rds")

# ---- Run ---------------------------------------------------------------------
fisher_results_ZHP   <- run_fisher(ZHP,   "ZHP")
fisher_results_RPNPF <- run_fisher(RPNPF, "RPNPF")


