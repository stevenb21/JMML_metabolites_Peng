# ==== Load Libraries ====
library(tidyverse)  # Core data manipulation and visualization
library(readxl)     # For reading Excel files
library(ggrepel)    # For non-overlapping text labels in ggplot2
library(janitor)

if (!dir.exists("../res/annotated/somatic")) dir.create("../res/annotated/somatic", recursive = TRUE)

# ====  Load and clean genomics metadata ==== 
genomics <- read_excel("../data/JMML_DBS_Genomics.xlsx") %>%
  clean_names() %>% #View()
  rename(
    project_id     = project_id,
    upn_id         = upn_stieglitz_lab,
    somatic_status = somatic_mutation_present_in_dbs,
    sex_raw        = sex
  ) %>%
  mutate(
    project_id_lower = str_to_lower(project_id),
    somatic_present = case_when(
      tolower(somatic_status) == "yes" ~ "Somatic_Pos",
      tolower(somatic_status) == "no"  ~ "Somatic_Neg",
      TRUE                                             ~ NA_character_
    ),
    somatic_present = factor(somatic_present, levels = c("Somatic_Neg", "Somatic_Pos"))
  )




# ==== ZHP Analysis ====

# Not-imputed dataset
ZHP_filtered <- readRDS("../data/R_objects/ZHP_annotated_filtered.rds") %>% 
  filter(sample_group %in% c("JMML", "Control"))

# Load long-format imputed ZHP data (M2, missing < 20%)
ZHP_imp_long <- readRDS("../data/R_objects/ZHP_imputed.rds") %>% 
  filter(sample_group %in% c("JMML"))

ZHP_imp_long <- ZHP_imp_long %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(genomics %>% select(project_id_lower, somatic_present), 
            by = c("sampleid" = "project_id_lower"))  %>% 
  filter(somatic_present %in% c("Somatic_Neg", "Somatic_Pos"))

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


# Join Somatic into the filtered ZHP data
ZHP_filtered <- ZHP_filtered %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(genomics %>% select(project_id_lower, somatic_present),
            by = c("sampleid" = "project_id_lower")) %>% 
  filter(somatic_present %in% c("Somatic_Neg", "Somatic_Pos"))

# Attach to imputed long-format data
ZHP_filtered <- ZHP_filtered %>%
  left_join(compound_missingness_ZHP %>% select(compound_name, missing_category),
            by = "compound_name")


# ==== Fisher's Exact Tests (M1) ===============================================

ZHP_M1 <- ZHP_filtered %>%
  filter(missing_category == "M1")


# ==== Fisher's Exact Tests (M1) JMML  =========================================

miss_tbl_full_ZHP <- ZHP_M1 %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_full_ZHP <- miss_tbl_full_ZHP %>%
  count(compound_name, somatic_present, missing) %>%
  pivot_wider(
    names_from  = c(somatic_present, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_Somatic_Neg_NA = Somatic_Neg_NA,
    b_Somatic_Neg_Present = Somatic_Neg_Present,
    c_Somatic_Pos_NA = Somatic_Pos_NA,
    d_Somatic_Pos_Present = Somatic_Pos_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_Somatic_Neg_NA,
                                              b_Somatic_Neg_Present,
                                              c_Somatic_Pos_NA,
                                              d_Somatic_Pos_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("Somatic_Neg","Somatic_Pos"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_Somatic_Neg_NA,
         b_Somatic_Neg_Present,
         c_Somatic_Pos_NA,
         d_Somatic_Pos_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)






# ==== T-tests (M2) ============================================================
# ==== T-tests JMML (M2) =======================================================
#ZHP_imp_long 

## Group means -----------------------------
group_means_somatic_full_ZHP <- ZHP_imp_long %>%
  group_by(compound_name, somatic_present) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = somatic_present,
    values_from = mean_intensity,
    names_prefix = "mean_"
  )

## t-test  -----------------------------
ttest_results_somatic_full_ZHP <- ZHP_imp_long %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ somatic_present)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_somatic_full_ZHP, by = "compound_name") %>%
  mutate(
    log2_fc = mean_Somatic_Pos - mean_Somatic_Neg,                        # Effect size: Somatic_pos - somatic_neg
  ) %>% mutate(raw_p = p_value)




# ==== RPNPF Analysis ====

# Not-imputed dataset
RPNPF_filtered <- readRDS("../data/R_objects/RPNPF_annotated_filtered.rds") %>% 
  filter(sample_group %in% c("JMML"))

# Join sex into the filtered RPNPF data
RPNPF_filtered <- RPNPF_filtered %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(genomics %>% select(project_id_lower, somatic_present),
            by = c("sampleid" = "project_id_lower")) %>% 
  filter(somatic_present %in% c("Somatic_Neg", "Somatic_Pos"))

# Load long-format imputed RPNPF data (M2, missing < 20%)
RPNPF_imp_long <- readRDS("../data/R_objects/RPNPF_imputed.rds") %>% 
  filter(sample_group %in% c("JMML"))

RPNPF_imp_long <- RPNPF_imp_long %>%
  mutate(sampleid = tolower(SampleID)) %>%
  left_join(genomics %>% select(project_id_lower, somatic_present), 
            by = c("sampleid" = "project_id_lower"))  %>% 
  filter(somatic_present %in% c("Somatic_Neg", "Somatic_Pos"))

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


# ==== Fisher's Exact Tests (M1) ===============================================

RPNPF_M1 <- RPNPF_filtered %>%
  filter(missing_category == "M1")


# ==== Fisher's Exact Tests (M1) JMML  =========================================

miss_tbl_full_RPNPF <- RPNPF_M1 %>%
  mutate(missing = if_else(is.na(intensity) | intensity == 0, "NA", "Present"))


##  count → wide (JMML_NA etc.) -------------------------------------------
fisher_full_RPNPF <- miss_tbl_full_RPNPF %>%
  count(compound_name, somatic_present, missing) %>%
  pivot_wider(
    names_from  = c(somatic_present, missing),
    values_from = n,
    values_fill = 0,
    names_sep   = "_") %>%
  rowwise() %>%
  mutate(
    ##  define a, b, c, d ---------------------------------------------------
    a_Somatic_Neg_NA = Somatic_Neg_NA,
    b_Somatic_Neg_Present = Somatic_Neg_Present,
    c_Somatic_Pos_NA = Somatic_Pos_NA,
    d_Somatic_Pos_Present = Somatic_Pos_Present,
    
    ##  Fisher test --------------------------------------------------------
    fisher_result = list(fisher.test(matrix(c(a_Somatic_Neg_NA,
                                              b_Somatic_Neg_Present,
                                              c_Somatic_Pos_NA,
                                              d_Somatic_Pos_Present),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(
                                              Group  = c("Somatic_Neg","Somatic_Pos"),
                                              Status = c("Missing","Present")
                                            )))),
    odds_ratio    = fisher_result$estimate[[1]],
    ci_lower      = fisher_result$conf.int[1],
    ci_upper      = fisher_result$conf.int[2],
    fisher_pval   = fisher_result$p.value,
    log2_OR       = log2(odds_ratio)
  ) %>%
  ungroup() %>%
  select(compound_name, a_Somatic_Neg_NA,
         b_Somatic_Neg_Present,
         c_Somatic_Pos_NA,
         d_Somatic_Pos_Present,
         odds_ratio, ci_lower, ci_upper, log2_OR, fisher_pval) %>% 
  mutate(raw_p = fisher_pval)






# ==== T-tests (M2) ============================================================
# ==== T-tests JMML (M2) =======================================================
#RPNPF_imp_long 

## Group means -----------------------------
group_means_somatic_full_RPNPF <- RPNPF_imp_long %>%
  group_by(compound_name, somatic_present) %>%
  summarize(mean_intensity = mean(intensity_log2, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = somatic_present,
    values_from = mean_intensity,
    names_prefix = "mean_"
  )

## t-test  -----------------------------
ttest_results_somatic_full_RPNPF <- RPNPF_imp_long %>%
  group_by(compound_name) %>%
  summarise(
    t_test = list(t.test(intensity_log2 ~ somatic_present)),  # Run unpaired two-group t-test
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),  # Extract t statistic
    p_value     = map_dbl(t_test, ~ .x$p.value)     # Extract raw p-value
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means_somatic_full_RPNPF, by = "compound_name") %>%
  mutate(
    log2_fc = mean_Somatic_Pos - mean_Somatic_Neg,                        # Effect size: Somatic_pos - somatic_neg
  ) %>% mutate(raw_p = p_value)

# # ==== Combine above  ==========================================================
# 
# # Add platform and test labels
fisher_full_ZHP <- fisher_full_ZHP %>%
  mutate(platform = "ZHP",
         test_type = "Fisher")

ttest_results_somatic_full_ZHP <- ttest_results_somatic_full_ZHP %>% 
  mutate(platform = "ZHP",
         test_type = "TTest")

fisher_full_RPNPF <- fisher_full_RPNPF %>%
  mutate(platform = "RPNPF",
         test_type = "Fisher")

ttest_results_somatic_full_RPNPF <- ttest_results_somatic_full_RPNPF %>% 
  mutate(platform = "RPNPF",
         test_type = "TTest")


 
# Combine all
combined_results <- bind_rows(
  fisher_full_ZHP,
  ttest_results_somatic_full_ZHP,
  fisher_full_RPNPF,
  ttest_results_somatic_full_RPNPF
)

# Global FDR adjustment
combined_results <- combined_results %>%
  mutate(global_fdr = p.adjust(raw_p, method = "BH"),
         analysis = "Somatic")

# # ==== Save to results to file =================================================
writexl::write_xlsx(combined_results, "../res/annotated/somatic/somatic_analysis.xlsx")
