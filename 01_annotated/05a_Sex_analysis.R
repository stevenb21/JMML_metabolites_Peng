# ==== Load Libraries ====
library(tidyverse)  # Core data manipulation and visualization
library(readxl)     # For reading Excel files
library(ggrepel)    # For non-overlapping text labels in ggplot2

if (!dir.exists("../res/annotated/sex")) dir.create("../res/annotated/sex", recursive = TRUE)

# ==== Load Data ====

# Load sample metadata and standardize sample ID format
# JMML info
sample_info <- read.csv("../data/Guthrie_card_additional_data.csv") %>%
  janitor::clean_names() %>%
  transmute(sampleid = tolower(sample_id),
            sex = sex)  # Ensures case-insensitive joins

sample_info$sampleid[sample_info$sampleid == "guth_1"] <- "guth1"
sample_info$sampleid[sample_info$sampleid == "guth_15"] <- "guth15"

control_info <- readxl::read_xlsx("../data/Control_Main_Info.xlsx") %>% 
  janitor::clean_names() %>%
  transmute(sampleid = tolower(id),
            sex = control_gender) # match other dataset

sample_info <- bind_rows(sample_info, control_info)

# ==== ZHP Analysis ====

# Not-imputed dataset
ZHP_filtered <- readRDS("../data/R_objects/ZHP_annotated_filtered.rds") %>% 
  filter(sample_group %in% c("JMML", "Control"))

# Join sex into the filtered ZHP data
ZHP_filtered <- ZHP_filtered %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(sample_info %>% select(sampleid, sex), by = "sampleid")

# Load long-format imputed ZHP data (M2, missing < 20%)
ZHP_imp_long <- readRDS("../data/R_objects/ZHP_imputed.rds")

ZHP_imp_long <- ZHP_imp_long %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(sample_info %>% select(sampleid, sex), by = "sampleid") %>%
  filter(sample_group %in% c("JMML", "Control"))

# ==== Define Missingness Categories (Pre-Imputation) ====

# Compute missingness (% NA or zero) per compound
compound_missingness_ZHP <- ZHP_filtered %>%
  group_by(compound_name) %>%
  summarize(
    n_total = n(),
    n_missing = sum(is.na(intensity) | intensity == 0),
    missing_pct = n_missing / n_total * 100,
    .groups = "drop"
  ) %>%
  mutate(
    missing_category = case_when(
      missing_pct >= 20 & missing_pct <= 65 ~ "M1",
      missing_pct < 20                      ~ "M2",
      TRUE                                  ~ "Excluded"
    )
  )


# Attach to imputed long-format data
ZHP_filtered <- ZHP_filtered %>%
  left_join(compound_missingness_ZHP %>% select(compound_name, missing_category),
            by = "compound_name")



# Sanity check. identical(compounds in M2, compounds in ZHP_imputed)
# M2_compounds <- ZHP_filtered %>% 
#   filter(missing_category == "M2") %>% 
#   distinct(compound_name) %>%
#   pull(compound_name)
# 
# imp_compounds <- ZHP_imp_long %>%
#   distinct(compound_name) %>%
#   pull(compound_name)
# # 
# # identical(M2_compounds, imp_compounds) #TRUE
# 
# not_imp <- ZHP_filtered %>% 
#      filter(missing_category != "M2") %>% 
#   distinct(compound_name) %>%
#   pull(compound_name)
# 
# intersect(imp_compounds, not_imp) #character(0)

# ==== Fisher's Exact Tests (M1) ===============================================

ZHP_M1 <- ZHP_filtered %>%
  filter(missing_category == "M1")


# ==== Fisher's Exact Tests (M1) All samples (JMML + Control) ==================

miss_tbl_full_ZHP <- ZHP_M1 %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_full_ZHP <- miss_tbl_full_ZHP %>%
  count(compound_name, sex, missing) %>%
  pivot_wider(
    names_from  = c(sex, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  # ensure columns exist even if all 0
  replace_na(list(
    F_NA = 0, F_Present = 0,
    M_NA = 0, M_Present = 0)) %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_F_NA = F_NA,
    b_F_Present = F_Present,
    c_M_NA = M_NA,
    d_M_Present = M_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_F_NA,
                                              b_F_Present,
                                              c_M_NA,
                                              d_M_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("F","M"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_F_NA,
         b_F_Present,
         c_M_NA,
         d_M_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)

# ==== Fisher's Exact Tests (M1) Control-only ==================================

miss_tbl_control_ZHP <- ZHP_M1 %>%
  filter(sample_group %in% c("Control")) %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_control_ZHP <- miss_tbl_control_ZHP %>%
  count(compound_name, sex, missing) %>%
  pivot_wider(
    names_from  = c(sex, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  # ensure columns exist even if all 0
  replace_na(list(
    F_NA = 0, F_Present = 0,
    M_NA = 0, M_Present = 0)) %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_F_NA = F_NA,
    b_F_Present = F_Present,
    c_M_NA = M_NA,
    d_M_Present = M_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_F_NA,
                                              b_F_Present,
                                              c_M_NA,
                                              d_M_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("F","M"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_F_NA,
         b_F_Present,
         c_M_NA,
         d_M_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)


# ==== Fisher's Exact Tests (M1) JMML-only =====================================

miss_tbl_JMML_ZHP <- ZHP_M1 %>%
  filter(sample_group %in% c("JMML")) %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_JMML_ZHP <- miss_tbl_JMML_ZHP %>%
  count(compound_name, sex, missing) %>%
  pivot_wider(
    names_from  = c(sex, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  # ensure columns exist even if all 0
  replace_na(list(
    F_NA = 0, F_Present = 0,
    M_NA = 0, M_Present = 0)) %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_F_NA = F_NA,
    b_F_Present = F_Present,
    c_M_NA = M_NA,
    d_M_Present = M_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_F_NA,
                                              b_F_Present,
                                              c_M_NA,
                                              d_M_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("F","M"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_F_NA,
         b_F_Present,
         c_M_NA,
         d_M_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)





# ==== T-tests (M2) ============================================================
# ==== T-tests All samples (JMML + Control) (M2) ===============================
#ZHP_imp_long 

## Group means -----------------------------
group_means_sex_full_ZHP <- ZHP_imp_long %>%
  group_by(compound_name, sex) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_intensity,
    names_prefix = "mean_"
  )

## t-test  -----------------------------
ttest_results_sex_full_ZHP <- ZHP_imp_long %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ sex)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_sex_full_ZHP, by = "compound_name") %>%
  mutate(
    log2_fc = mean_F - mean_M,                        # Effect size: Female - Male
  ) %>% mutate(raw_p = p_value)


# ==== T-tests Control-only (M2) ===============================================

ZHP_imp_long_control <- ZHP_imp_long %>% filter(sample_group == "Control")


group_means_sex_control_ZHP <- ZHP_imp_long_control %>%
  group_by(compound_name, sex) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_intensity,
    names_prefix = "mean_"
  ) 

## t-test  -----------------------------
ttest_results_sex_control_ZHP <- ZHP_imp_long_control %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ sex)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_sex_control_ZHP, by = "compound_name") %>%
  mutate(
    log2_fc = mean_F - mean_M,                        # Effect size: Female - Male
  ) %>% mutate(raw_p = p_value)


# ==== T-tests JMML-only (M2) ==================================================
ZHP_imp_long_JMML_ZHP <- ZHP_imp_long %>% filter(sample_group == "JMML")

group_means_sex_JMML_ZHP <- ZHP_imp_long_JMML_ZHP %>%
  group_by(compound_name, sex) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_intensity,
    names_prefix = "mean_"
  )

## t-test  -----------------------------
ttest_results_sex_JMML_ZHP <- ZHP_imp_long_JMML_ZHP %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ sex)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_sex_JMML_ZHP, by = "compound_name") %>%
  mutate(
    log2_fc = mean_F - mean_M,                        # Effect size: Female - Male
  ) %>% mutate(raw_p = p_value)


# ==== RPNPF Analysis ====

# Not-imputed dataset
RPNPF_filtered <- readRDS("../data/R_objects/RPNPF_annotated_filtered.rds") %>% 
  filter(sample_group %in% c("JMML", "Control"))

# Join sex into the filtered RPNPF data
RPNPF_filtered <- RPNPF_filtered %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(sample_info %>% select(sampleid, sex), by = "sampleid")

# Load long-format imputed RPNPF data (M2, missing < 20%)
RPNPF_imp_long <- readRDS("../data/R_objects/RPNPF_imputed.rds")

RPNPF_imp_long <- RPNPF_imp_long %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(sample_info %>% select(sampleid, sex), by = "sampleid") %>%
  filter(sample_group %in% c("JMML", "Control"))

# ==== Define Missingness Categories (Pre-Imputation) ====

# Compute missingness (% NA or zero) per compound
compound_missingness_RPNPF <- RPNPF_filtered %>%
  group_by(compound_name) %>%
  summarize(
    n_total = n(),
    n_missing = sum(is.na(intensity) | intensity == 0),
    missing_pct = n_missing / n_total * 100,
    .groups = "drop"
  ) %>%
  mutate(
    missing_category = case_when(
      missing_pct >= 20 & missing_pct <= 65 ~ "M1",
      missing_pct < 20                      ~ "M2",
      TRUE                                  ~ "Excluded"
    )
  )


# Attach to imputed long-format data
RPNPF_filtered <- RPNPF_filtered %>%
  left_join(compound_missingness_RPNPF %>% select(compound_name, missing_category),
            by = "compound_name")



# Sanity check. identical(compounds in M2, compounds in RPNPF_imputed)
# M2_compounds <- RPNPF_filtered %>% 
#   filter(missing_category == "M2") %>% 
#   distinct(compound_name) %>%
#   pull(compound_name)
# 
# imp_compounds <- RPNPF_imp_long %>%
#   distinct(compound_name) %>%
#   pull(compound_name)
# # 
# # identical(M2_compounds, imp_compounds) #TRUE
# 
# not_imp <- RPNPF_filtered %>% 
#      filter(missing_category != "M2") %>% 
#   distinct(compound_name) %>%
#   pull(compound_name)
# 
# intersect(imp_compounds, not_imp) #character(0)

# ==== Fisher's Exact Tests (M1) ===============================================

RPNPF_M1 <- RPNPF_filtered %>%
  filter(missing_category == "M1")


# ==== Fisher's Exact Tests (M1) All samples (JMML + Control) ==================

miss_tbl_full_RPNPF <- RPNPF_M1 %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_full_RPNPF <- miss_tbl_full_RPNPF %>%
  count(compound_name, sex, missing) %>%
  pivot_wider(
    names_from  = c(sex, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  # ensure columns exist even if all 0
  replace_na(list(
    F_NA = 0, F_Present = 0,
    M_NA = 0, M_Present = 0)) %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_F_NA = F_NA,
    b_F_Present = F_Present,
    c_M_NA = M_NA,
    d_M_Present = M_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_F_NA,
                                              b_F_Present,
                                              c_M_NA,
                                              d_M_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("F","M"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_F_NA,
         b_F_Present,
         c_M_NA,
         d_M_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)

# ==== Fisher's Exact Tests (M1) Control-only ==================================

miss_tbl_control_RPNPF <- RPNPF_M1 %>%
  filter(sample_group %in% c("Control")) %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_control_RPNPF <- miss_tbl_control_RPNPF %>%
  count(compound_name, sex, missing) %>%
  pivot_wider(
    names_from  = c(sex, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  # ensure columns exist even if all 0
  replace_na(list(
    F_NA = 0, F_Present = 0,
    M_NA = 0, M_Present = 0)) %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_F_NA = F_NA,
    b_F_Present = F_Present,
    c_M_NA = M_NA,
    d_M_Present = M_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_F_NA,
                                              b_F_Present,
                                              c_M_NA,
                                              d_M_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("F","M"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_F_NA,
         b_F_Present,
         c_M_NA,
         d_M_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)


# ==== Fisher's Exact Tests (M1) JMML-only =====================================

miss_tbl_JMML_RPNPF <- RPNPF_M1 %>%
  filter(sample_group %in% c("JMML")) %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_JMML_RPNPF <- miss_tbl_JMML_RPNPF %>%
  count(compound_name, sex, missing) %>%
  pivot_wider(
    names_from  = c(sex, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  # ensure columns exist even if all 0
  replace_na(list(
    F_NA = 0, F_Present = 0,
    M_NA = 0, M_Present = 0)) %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_F_NA = F_NA,
    b_F_Present = F_Present,
    c_M_NA = M_NA,
    d_M_Present = M_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_F_NA,
                                              b_F_Present,
                                              c_M_NA,
                                              d_M_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("F","M"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_F_NA,
         b_F_Present,
         c_M_NA,
         d_M_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)





# ==== T-tests (M2) ============================================================
# ==== T-tests All samples (JMML + Control) (M2) ===============================
#RPNPF_imp_long 

## Group means -----------------------------
group_means_sex_full_RPNPF <- RPNPF_imp_long %>%
  group_by(compound_name, sex) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_intensity,
    names_prefix = "mean_"
  )

## t-test  -----------------------------
ttest_results_sex_full_RPNPF <- RPNPF_imp_long %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ sex)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_sex_full_RPNPF, by = "compound_name") %>%
  mutate(
    log2_fc = mean_F - mean_M,                        # Effect size: Female - Male
  ) %>% mutate(raw_p = p_value)


# ==== T-tests Control-only (M2) ===============================================

RPNPF_imp_long_control <- RPNPF_imp_long %>% filter(sample_group == "Control")


group_means_sex_control_RPNPF <- RPNPF_imp_long_control %>%
  group_by(compound_name, sex) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_intensity,
    names_prefix = "mean_"
  ) 

## t-test  -----------------------------
ttest_results_sex_control_RPNPF <- RPNPF_imp_long_control %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ sex)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_sex_control_RPNPF, by = "compound_name") %>%
  mutate(
    log2_fc = mean_F - mean_M,                        # Effect size: Female - Male
  ) %>% mutate(raw_p = p_value)


# ==== T-tests JMML-only (M2) ==================================================
RPNPF_imp_long_JMML_RPNPF <- RPNPF_imp_long %>% filter(sample_group == "JMML")

group_means_sex_JMML_RPNPF <- RPNPF_imp_long_JMML_RPNPF %>%
  group_by(compound_name, sex) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = sex,
    values_from = mean_intensity,
    names_prefix = "mean_"
  )

## t-test  -----------------------------
ttest_results_sex_JMML_RPNPF <- RPNPF_imp_long_JMML_RPNPF %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ sex)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_sex_JMML_RPNPF, by = "compound_name") %>%
  mutate(
    log2_fc = mean_F - mean_M,                        # Effect size: Female - Male
  ) %>% mutate(raw_p = p_value)




# ==== Combine above  ==========================================================

# Add platform and test labels
fisher_full_RPNPF <- fisher_full_RPNPF %>% 
  mutate(platform = "RPNPF",
         test_type = "Fisher",
         cohort = "Combined_JMML_Control")


fisher_control_RPNPF <- fisher_control_RPNPF %>% 
  mutate(platform = "RPNPF",
         test_type = "Fisher",
         cohort = "Control_only")

fisher_JMML_RPNPF <- fisher_JMML_RPNPF %>% 
  mutate(platform = "RPNPF",
         test_type = "Fisher",
         cohort = "JMML_only")


ttest_results_sex_JMML_RPNPF <- ttest_results_sex_JMML_RPNPF %>% 
  mutate(platform = "RPNPF",
         test_type = "TTest",
         cohort = "JMML_only")

ttest_results_sex_control_RPNPF <- ttest_results_sex_control_RPNPF %>% 
  mutate(platform = "RPNPF",
         test_type = "TTest",
         cohort = "Control_only")

ttest_results_sex_full_RPNPF <- ttest_results_sex_full_RPNPF %>% 
  mutate(platform = "RPNPF",
         test_type = "TTest",
         cohort = "Combined_JMML_Control")


# Add platform and test labels
fisher_full_ZHP <- fisher_full_ZHP %>% 
  mutate(platform = "ZHP",
         test_type = "Fisher",
         cohort = "Combined_JMML_Control")


fisher_control_ZHP <- fisher_control_ZHP %>% 
  mutate(platform = "ZHP",
         test_type = "Fisher",
         cohort = "Control_only")

fisher_JMML_ZHP <- fisher_JMML_ZHP %>% 
  mutate(platform = "ZHP",
         test_type = "Fisher",
         cohort = "JMML_only")


ttest_results_sex_JMML_ZHP <- ttest_results_sex_JMML_ZHP %>% 
  mutate(platform = "ZHP",
         test_type = "TTest",
         cohort = "JMML_only")

ttest_results_sex_control_ZHP <- ttest_results_sex_control_ZHP %>% 
  mutate(platform = "ZHP",
         test_type = "TTest",
         cohort = "Control_only")

ttest_results_sex_full_ZHP <- ttest_results_sex_full_ZHP %>% 
  mutate(platform = "ZHP",
         test_type = "TTest",
         cohort = "Combined_JMML_Control")



# Combine all
combined_results <- bind_rows(
  fisher_full_ZHP,
  fisher_control_ZHP, 
  fisher_JMML_ZHP, 
  ttest_results_sex_JMML_ZHP,
  ttest_results_sex_control_ZHP,
  ttest_results_sex_full_ZHP,
  fisher_full_RPNPF,
  fisher_control_RPNPF, 
  fisher_JMML_RPNPF, 
  ttest_results_sex_JMML_RPNPF,
  ttest_results_sex_control_RPNPF,
  ttest_results_sex_full_RPNPF
)

# Global FDR adjustment
combined_results <- combined_results %>%
  mutate(global_fdr = p.adjust(raw_p, method = "BH"),
         analysis = "Sex")



# ==== Save to results to file =================================================
writexl::write_xlsx(combined_results, "../res/annotated/sex/sex_analysis.xlsx")



