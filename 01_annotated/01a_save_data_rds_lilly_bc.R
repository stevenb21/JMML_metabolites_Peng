# ==============================================================================
# Script Name: save_data_rds.R
# Author: Steven Brooks
# Date: 2025-05-13
# Description: Load the files in, filter, and save to .rds files.
# ==============================================================================

# ==== Setup ===================================================================
library(janitor)
library(tidyverse)
library(readxl)
library(purrr)

if (!dir.exists("../res/annotated")) dir.create("../res/annotated", recursive = TRUE)
if (!dir.exists("../data/R_objects")) dir.create("../res/R_objects", recursive = TRUE)


# ==== Load data ===============================================================

# Annotated
raw_RPNPF <- readxl::read_xlsx("../data/Annotated_RPNPF_JMML_Meyer.xlsx")
raw_ZHP <- readxl::read_xlsx("../data/Annotated_ZHP_JMML_Meyer.xlsx")
raw_meta_info <- readxl::read_xlsx("../data/MetaFileJMML_annotated.xlsx")


#Clean column names
raw_RPNPF <- raw_RPNPF %>% janitor::clean_names()
raw_ZHP <- raw_ZHP %>% janitor::clean_names()

#For  Ectoine: please use row 202
#For 5-Aminopentanoic acid: please use row 222
raw_ZHP <- raw_ZHP[-c(199, 238),]



# ==== ZHP Long format =========================================================

#Annotated
ZHP <- raw_ZHP %>% dplyr::select(-c(molecular_formula,mass,rt_min,rt_sec,adduct,mz,cas_id))

# Create sample group column
raw_meta_info <- raw_meta_info %>%
  mutate(
    sample_group = case_when(
      SampleType == "sample" ~ case_when(
        Sample_Type == "JMML" ~ "JMML",
        Sample_Type == "Control" ~ "Control",
        TRUE ~ "Blank"
      ),
      TRUE ~ SampleType
    )
  )


# Pivot longer to tidy format
ZHP_long <- ZHP %>%
  pivot_longer(
    cols = -compound_name,  # Keep compound_name as identifier
    names_to = "SampleID",  # Column to hold sample names
    values_to = "intensity" # Column to hold the measurement values
  )



# Ensure matching case for SampleID
raw_meta_info <- raw_meta_info %>%
  mutate(SampleID_lower = tolower(SampleID))

ZHP_long <- ZHP_long %>%
  mutate(SampleID_lower = tolower(SampleID)) %>%
  left_join(
    raw_meta_info %>% select(SampleID_lower, sample_group),
    by = "SampleID_lower"
  ) %>%
  select(-SampleID_lower) %>%
  replace_na(list(sample_group = "blank"))


# ==== ZHP Filter 65% Missingness ==============================================


# Get compound names passing the missingness filter (only JMML & Control considered)
valid_compounds <- ZHP_long %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  group_by(compound_name) %>%
  summarise(
    prop_na_or_zero = mean(is.na(intensity) | intensity == 0),
    .groups = "drop"
  ) %>%
  filter(prop_na_or_zero <= 0.65) %>%
  pull(compound_name)

# Step 4: Filter entire dataset to keep only those compounds
ZHP_NA_filtered <- ZHP_long %>%
  filter(compound_name %in% valid_compounds)



# ==== ZHP Filter Noise-to-Signal ==============================================

ZHP_signal_filtered <- ZHP_NA_filtered %>%
  group_by(compound_name) %>%
  summarise(
    median_matrix = median(intensity[sample_group == "matrix"], na.rm = TRUE),
    median_pooledQC = median(intensity[sample_group == "PooledQC"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    noise_to_signal = median_matrix / median_pooledQC
  ) %>%
  filter(
    # Keep if median_pooledQC is NOT NA
    !is.na(median_pooledQC),
    # Keep if noise_to_signal is ≤ 0.3 OR median_matrix is NA
    noise_to_signal <= 0.3 | is.na(median_matrix)
  ) %>%
  select(compound_name)


# Filter original data to keep only good compounds
ZHP_NtS_NA_filtered <- ZHP_NA_filtered %>%
  semi_join(ZHP_signal_filtered, by = "compound_name")


# ==== Hemoglobin Correlation ==================================================

# Load hemoglobin
hemoglobin <- raw_meta_info %>% 
  select(SampleID, sample_group, `Hemoglobin (mg/mL)`) %>% 
  mutate(SampleID = tolower(SampleID)) %>% 
  drop_na()



# Merge hemoglobin with filtered metabolite data
ZHP_with_hgb <- ZHP_NtS_NA_filtered %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  left_join(
    hemoglobin %>% rename(hemoglobin_mg_per_mL = `Hemoglobin (mg/mL)`),
    by = c("SampleID", "sample_group")
  ) %>%
  filter(!is.na(hemoglobin_mg_per_mL), !is.na(intensity))

# Spearman's cor on raw intensity
spearman_hgb_cor_ZHP <- ZHP_with_hgb %>%
  group_by(compound_name) %>%
  summarise(
    n = n(),
    rho = cor(intensity, hemoglobin_mg_per_mL, method = "spearman"),
    .groups = "drop"
  )


#Histogram of rho values
p_hist_rho_ZHP <- ggplot(spearman_hgb_cor_ZHP, aes(x = rho)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Histogram of Spearman's ρ Values\nZHP",
       x = "Spearman's ρ",
       y = "Count") +
  theme_bw()

# Save the histogram
ggsave(
  filename = "../res/annotated/histogram_rho_ZHP.png",
  plot = p_hist_rho_ZHP,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)


# ==== ZHP Power Transformation ================================================


lambdas_QC <- ZHP_NtS_NA_filtered %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  group_by(compound_name) %>%
  summarise(
    lambda = {
      model <- lm(intensity ~ 1, data = cur_data())
      bc <- MASS::boxcox(model, lambda = seq(-2, 2, 0.1), plotit = FALSE)
      lambda_opt <- bc$x[which.max(bc$y)]
      lambda_opt
    },
    .groups = "drop"
  )

# 2. Join λ and apply transformation only when intensity > 0
ZHP_transformed <- ZHP_NtS_NA_filtered %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  left_join(lambdas_QC, by = "compound_name") %>%
  mutate(
    intensity_bc = case_when(
      is.na(intensity) | intensity <= 0 | is.na(lambda) ~ NA_real_,
      abs(lambda) < 1e-8 ~ log(intensity),
      TRUE ~ (intensity^lambda - 1) / lambda
    )
  )


head(ZHP_transformed)

# ==== ZHP Filter Coefficient of Variation =====================================

cv_tbl_ZHP <- ZHP_transformed %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  group_by(compound_name, sample_group) %>%
  summarise(
    mean_intensity_bc = mean(intensity_bc, na.rm = TRUE),
    sd_intensity_bc   = sd(intensity_bc,   na.rm = TRUE),
    cv                = sd_intensity_bc / mean_intensity_bc,
    .groups           = "drop"
  )

# Pivot all three summary columns to wide format
cv_wide_ZHP <- cv_tbl_ZHP %>%
  pivot_wider(
    names_from  = sample_group,
    values_from = c(mean_intensity_bc, sd_intensity_bc, cv),
    names_sep = "_"
  )  %>%
  filter(
    !is.na(cv_Control) & cv_Control < 0.3,
    !is.na(cv_JMML) & cv_JMML < 0.3
  )


# Filter based on conditions
cv_filtered_compounds_ZHP <- cv_wide_ZHP %>%
  pull(compound_name)




ZHP_final_filtered <- ZHP_NtS_NA_filtered %>%
  semi_join(tibble(compound_name = cv_filtered_compounds_ZHP), by = "compound_name")

# ==== ZHP Filtering Summary Table =============================================

ZHP_filter_summary <- tibble(
  step = c(
    "Initial",
    "After ≤65% Missingness",
    "After Pooled Noise-to-Signal ≤ 0.3",
    "After CV(JMML)  ≤ 0.3 & CV(Control) ≤ 0.3"
  ),
  n_metabolites = c(
    n_distinct(ZHP_long$compound_name),
    n_distinct(ZHP_NA_filtered$compound_name),
    n_distinct(ZHP_NtS_NA_filtered$compound_name),
    n_distinct(ZHP_final_filtered$compound_name)
  )
)





# ==== RPNPF Long format =========================================================

#Annotated
RPNPF <- raw_RPNPF %>% select(-c(formula,mass,rt_min,rt_sec,adduct,m_z,cas_id))



# Create sample group column
raw_meta_info <- raw_meta_info %>%
  mutate(
    sample_group = case_when(
      SampleType == "sample" ~ case_when(
        Sample_Type == "JMML" ~ "JMML",
        Sample_Type == "Control" ~ "Control",
        TRUE ~ "Blank"
      ),
      TRUE ~ SampleType
    )
  )


# Pivot longer to tidy format
RPNPF_long <- RPNPF %>%
  pivot_longer(
    cols = -compound_name,  # Keep compound_name as identifier
    names_to = "SampleID",  # Column to hold sample names
    values_to = "intensity" # Column to hold the measurement values
  )

# Ensure matching case for SampleID
raw_meta_info <- raw_meta_info %>%
  mutate(SampleID_lower = tolower(SampleID))

RPNPF_long <- RPNPF_long %>%
  mutate(SampleID_lower = tolower(SampleID)) %>%
  left_join(
    raw_meta_info %>% select(SampleID_lower, sample_group),
    by = "SampleID_lower"
  ) %>%
  select(-SampleID_lower) %>%
  replace_na(list(sample_group = "blank"))


# ==== RPNPF Filter 65% Missingness ==============================================

# Get compound names passing the missingness filter (only JMML & Control considered)
valid_compounds <- RPNPF_long %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  group_by(compound_name) %>%
  summarise(
    prop_na_or_zero = mean(is.na(intensity) | intensity == 0),
    .groups = "drop"
  ) %>%
  filter(prop_na_or_zero <= 0.65) %>%
  pull(compound_name)

# Step 4: Filter entire dataset to keep only those compounds
RPNPF_NA_filtered <- RPNPF_long %>%
  filter(compound_name %in% valid_compounds)





# ==== RPNPF Filter Noise-to-Signal ==============================================

RPNPF_signal_filtered <- RPNPF_NA_filtered %>%
  group_by(compound_name) %>%
  summarise(
    median_matrix = median(intensity[sample_group == "matrix"], na.rm = TRUE),
    median_pooledQC = median(intensity[sample_group == "PooledQC"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    noise_to_signal = median_matrix / median_pooledQC
  ) %>%
  filter(
    # Keep if median_pooledQC is NOT NA
    !is.na(median_pooledQC),
    # Keep if noise_to_signal is ≤ 0.3 OR median_matrix is NA
    noise_to_signal <= 0.3 | is.na(median_matrix)
  ) %>%
  select(compound_name)


# Filter original data to keep only good compounds
RPNPF_NtS_NA_filtered <- RPNPF_NA_filtered %>%
  semi_join(RPNPF_signal_filtered, by = "compound_name")

# ==== Hemoglobin Correlation ==================================================

# Load hemoglobin
hemoglobin <- raw_meta_info %>% 
  select(SampleID, sample_group, `Hemoglobin (mg/mL)`) %>% 
  mutate(SampleID = tolower(SampleID)) %>% 
  drop_na()



# Merge hemoglobin with filtered metabolite data
RPNPF_with_hgb <- RPNPF_NtS_NA_filtered %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  left_join(
    hemoglobin %>% rename(hemoglobin_mg_per_mL = `Hemoglobin (mg/mL)`),
    by = c("SampleID", "sample_group")
  ) %>%
  filter(!is.na(hemoglobin_mg_per_mL), !is.na(intensity))

spearman_hgb_cor_RPNPF <- RPNPF_with_hgb %>%
  group_by(compound_name) %>%
  summarise(
    n = n(),
    rho = cor(intensity, hemoglobin_mg_per_mL, method = "spearman"),
    .groups = "drop"
  ) 


p_hist_rho_RPNPF <- ggplot(spearman_hgb_cor_RPNPF, aes(x = rho)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(title = "Histogram of Spearman's ρ Values\nRPNPF",
       x = "Spearman's ρ",
       y = "Count") +
  theme_bw()

# Save the histogram
ggsave(
  filename = "../res/annotated/histogram_rho_RPNPF.png",
  plot = p_hist_rho_RPNPF,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

# ==== RPNPF Power Transformation ================================================


lambdas_QC <- RPNPF_NtS_NA_filtered %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  group_by(compound_name) %>%
  summarise(
    lambda = {
      model <- lm(intensity ~ 1, data = cur_data())
      bc <- MASS::boxcox(model, lambda = seq(-2, 2, 0.1), plotit = FALSE)
      lambda_opt <- bc$x[which.max(bc$y)]
      lambda_opt
    },
    .groups = "drop"
  )

# 2. Join λ and apply transformation only when intensity > 0
RPNPF_transformed <- RPNPF_NtS_NA_filtered %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  left_join(lambdas_QC, by = "compound_name") %>%
  mutate(
    intensity_bc = case_when(
      is.na(intensity) | intensity <= 0 | is.na(lambda) ~ NA_real_,
      abs(lambda) < 1e-8 ~ log(intensity),
      TRUE ~ (intensity^lambda - 1) / lambda
    )
  )


head(RPNPF_transformed)

# ==== RPNPF Filter Coefficient of Variation =====================================

# Now compute CV for matrix and PooledQC:
cv_tbl_RPNPF <- RPNPF_transformed %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  group_by(compound_name, sample_group) %>%
  summarise(
    mean_intensity_bc = mean(intensity_bc, na.rm = TRUE),
    sd_intensity_bc   = sd(intensity_bc,   na.rm = TRUE),
    cv                = sd_intensity_bc / mean_intensity_bc,
    .groups           = "drop"
  )

# Pivot all three summary columns to wide format
cv_wide_RPNPF <- cv_tbl_RPNPF %>%
  pivot_wider(
    names_from  = sample_group,
    values_from = c(mean_intensity_bc, sd_intensity_bc, cv),
    names_sep = "_"
  ) %>%
  filter(
    !is.na(cv_Control) & cv_Control < 0.3,
    !is.na(cv_JMML) & cv_JMML < 0.3
  )


# Filter based on conditions
cv_filtered_compounds_RPNPF <- cv_wide_RPNPF %>%
  pull(compound_name)



RPNPF_final_filtered <- RPNPF_NtS_NA_filtered %>%
  semi_join(tibble(compound_name = cv_filtered_compounds_RPNPF), by = "compound_name")

# ==== RPNPF Filtering Summary Table =============================================

RPNPF_filter_summary <- tibble(
  step = c(
    "Initial",
    "After ≤65% Missingness",
    "After Pooled Noise-to-Signal ≤ 0.3",
    "After CV(JMML)  ≤ 0.3 & CV(Control) ≤ 0.3"
  ),
  n_metabolites = c(
    n_distinct(RPNPF_long$compound_name),
    n_distinct(RPNPF_NA_filtered$compound_name),
    n_distinct(RPNPF_NtS_NA_filtered$compound_name),
    n_distinct(RPNPF_final_filtered$compound_name)
  )
)





# ==== CV for duplicate metabolites  =============================================

#Are there duplicate metabolites? 
shared_compounds <- intersect(
  unique(RPNPF_final_filtered$compound_name),
  unique(ZHP_final_filtered$compound_name)
)

length(shared_compounds)  # number of overlapping compounds
shared_compounds          # view the compound names



#Compare CVs for ZHP and RPNPF for duplicate metabolites

cv_wide_ZHP <- cv_wide_ZHP %>% 
  mutate(platform = "ZHP")

cv_wide_RPNPF <- cv_wide_RPNPF %>% 
  mutate(platform = "RPNPF")

cv_combined <- bind_rows(cv_wide_ZHP, cv_wide_RPNPF)

# Keep only shared compounds
shared_compounds <- intersect(cv_wide_RPNPF$compound_name, cv_wide_ZHP$compound_name)

cv_combined_shared <- cv_combined %>%
  filter(compound_name %in% shared_compounds)


# pick the platform with the lower max CV across JMML and Control
cv_decision <- cv_combined_shared %>%
  mutate(cv_max = pmax(cv_JMML, cv_Control, na.rm = TRUE)) %>%
  group_by(compound_name) %>%
  slice_min(cv_max, with_ties = FALSE) %>%
  ungroup()

writexl::write_xlsx(cv_decision, "../res/annotated/cv_compare_selected.xlsx")

# Get compounds assigned to ZHP (to remove from RPNPF)
keep_ZHP <- cv_decision %>%
  filter(platform == "ZHP") %>%
  pull(compound_name)

# Remove those compounds from RPNPF
RPNPF_filtered_nodup <- RPNPF_final_filtered %>%
  filter(!(compound_name %in% keep_ZHP))


### Same idea as above
keep_RPNPF <- cv_decision %>%
  filter(platform == "RPNPF") %>%
  pull(compound_name)

# Remove those compounds from RPNPF
ZHP_filtered_nodup <- ZHP_final_filtered %>%
  filter(!(compound_name %in% keep_RPNPF))

# ==== Hemoglobin Table  =======================================================
# Add platform labels
spearman_hgb_cor_RPNPF <- spearman_hgb_cor_RPNPF %>%
  mutate(platform = "RPNPF")

spearman_hgb_cor_ZHP <- spearman_hgb_cor_ZHP %>%
  mutate(platform = "ZHP")

combined_hgb_cor <- bind_rows(spearman_hgb_cor_ZHP, spearman_hgb_cor_RPNPF)

writexl::write_xlsx(combined_hgb_cor, "../res/annotated/hemoglobin_correlated_metabolites.xlsx")



# ==== Save data  ==============================================================
saveRDS(ZHP_filtered_nodup, file = "../data/R_objects/ZHP_annotated_filtered.rds")
saveRDS(RPNPF_filtered_nodup, file = "../data/R_objects/RPNPF_annotated_filtered.rds")


filter_summary <- ZHP_filter_summary %>%
  rename(ZHP_metabolites = n_metabolites) %>%
  left_join(RPNPF_filter_summary %>% rename(RPNPF_metabolites = n_metabolites),
            by = "step")


filter_summary <- bind_rows(
  filter_summary,
  tibble(
    step = "Deduplicated shared metabolites by lowest max CV (JMML, Control)",
    ZHP_metabolites = n_distinct(ZHP_filtered_nodup$compound_name),
    RPNPF_metabolites = n_distinct(RPNPF_filtered_nodup$compound_name)
  )
)

writexl::write_xlsx(filter_summary, "../res/annotated/filter_summary.xlsx")
print(filter_summary)


