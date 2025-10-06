library(tidyverse)
library(readxl)
library(janitor)

# ==== Load metadata ====
# We would like to attach sex and case status as well to the dataframe.

# Load sample metadata and standardize sample ID format
# JMML info
##–– Sample info ––----
sample_info <- read.csv("../data/Guthrie_card_additional_data.csv") %>%
  janitor::clean_names() %>% 
  transmute(
    sampleid = tolower(sample_id),
    sex,
    age_at_colctn,
    race_ethnicity,
    bd_blood_collection_time_x = as.character(bd_blood_collection_time_x)
  ) %>%
  mutate(
    sampleid = recode(sampleid,
                      "guth_1"  = "guth1",
                      "guth_15" = "guth15"
    )
  )


##–– Control info ––----
control_info <- read_xlsx("../data/Control_Main_Info.xlsx") %>% 
  janitor::clean_names() %>%
  mutate(
    sampleid = tolower(id),
    sex      = control_gender
  ) %>% 
  rename(
    age_at_colctn              = control_age_at_colctn,
    race_ethnicity             = control_race_ethncty_dscr,
    bd_blood_collection_time_x = control_bd_blood_collection_time_x
  ) %>% 
  mutate(
    bd_blood_collection_time_x = as.character(bd_blood_collection_time_x)
  ) %>% 
  select(sampleid, sex, age_at_colctn, race_ethnicity, bd_blood_collection_time_x)

combined_info <- bind_rows(sample_info, control_info) %>%
  rename_with(~ paste0("meta_", .))




##–– Meta info ––----
raw_meta_info <- readxl::read_xlsx("../data/MetaFileJMML_annotated.xlsx") %>% 
  clean_names() %>% 
  transmute(meta_sampleid = tolower(id),
            meta_hemoglobin = hemoglobin_mg_m_l) %>% 
  drop_na()


combined_info <- combined_info %>% 
  left_join(raw_meta_info, by = "meta_sampleid")


### Clean up metadata:

combined_info$meta_sex <- as.factor(combined_info$meta_sex)


combined_info <- combined_info %>%
  mutate(meta_race_ethnicity = tolower(meta_race_ethnicity),
         meta_race_ethnicity = case_when(
           str_detect(meta_race_ethnicity, "white")        ~ "White",
           str_detect(meta_race_ethnicity, "black")        ~ "Black",
           str_detect(meta_race_ethnicity, "asian|chinese|filipino|east indian|vietnamese") ~ "Asian",
           str_detect(meta_race_ethnicity, "hispanic|latino") ~ "Hispanic",
           str_detect(meta_race_ethnicity, "more than")    ~ "Multiracial",
           str_detect(meta_race_ethnicity, "unknown")      ~ "Unknown",
           TRUE ~ "Other"
         ),
         meta_race_ethnicity = factor(meta_race_ethnicity))



# Missing age at collection for guth_1 and guth_15:
median_age <- median(combined_info$meta_age_at_colctn, na.rm = TRUE)

combined_info <- combined_info %>%
  mutate(meta_age_at_colctn = ifelse(is.na(meta_age_at_colctn), median_age, meta_age_at_colctn))


#===== ZHP Load ================================================================

# For each row we want samples, and for each column we want metabolite intensities and sample metadata.

## ==== M1 Data (Missingness <20%) =============================================

# We want the imputed box-cox 'M1 data' (Missingness <20%).
ZHP_imp_long <- readRDS("../data/R_objects/ZHP_imputed.rds")


## ==== 'M2 data' (20% - 65% Missing) ====

# We want to make this categorical.

## Not-imputed dataset
ZHP_filtered <- readRDS("../data/R_objects/ZHP_annotated_filtered.rds") %>% 
  filter(sample_group %in% c("JMML", "Control"))


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
      missing_pct >= 20 & missing_pct <= 65 ~ "M2",
      missing_pct < 20                      ~ "M1",
      TRUE                                  ~ "Excluded"
    )
  )


# Attach to imputed long-format data
ZHP_filtered <- ZHP_filtered %>%
  left_join(compound_missingness_ZHP %>% select(compound_name, missing_category),
            by = "compound_name") %>% 
  filter(missing_category == "M2") %>% 
  mutate(M2_Detected = ifelse(!is.na(intensity) & intensity != 0, 1, 0))

## ==== Load metadata ====


## Combine the continuous and categorical features =====


## ----------  Binary features  ----------
bin_wide_ZHP <- ZHP_filtered %>% 
  select(SampleID, compound_name, M2_Detected, sample_group) %>% 
  pivot_wider(names_from  = compound_name,
              values_from = M2_Detected,
              names_prefix = "bin_")  

## ----------  Box-Cox intensities  ----------
bc_wide_ZHP <- ZHP_imp_long %>% 
  select(SampleID, compound_name, intensity_bc, sample_group) %>% 
  pivot_wider(names_from  = compound_name,
              values_from = intensity_bc,
              names_prefix = "bc_")


full_ZHP <- bin_wide_ZHP %>% 
  full_join(bc_wide_ZHP, by = c("SampleID", "sample_group")) %>% 
  mutate(SampleID = tolower(SampleID)) %>% 
  left_join(combined_info, by = c("SampleID" = "meta_sampleid")) %>% 
  mutate(meta_sample_group = sample_group,
         meta_SampleID = SampleID) %>% 
  select(-sample_group, -SampleID)


#===== RPNPF Load ================================================================

# For each row we want samples, and for each column we want metabolite intensities and sample metadata.

## ==== M1 Data (Missingness <20%) =============================================

# We want the imputed box-cox 'M1 data' (Missingness <20%).
RPNPF_imp_long <- readRDS("../data/R_objects/RPNPF_imputed.rds")


## ==== 'M2 data' (20% - 65% Missing) ====

# We want to make this categorical.

## Not-imputed dataset
RPNPF_filtered <- readRDS("../data/R_objects/RPNPF_annotated_filtered.rds") %>% 
  filter(sample_group %in% c("JMML", "Control"))


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
      missing_pct >= 20 & missing_pct <= 65 ~ "M2",
      missing_pct < 20                      ~ "M1",
      TRUE                                  ~ "Excluded"
    )
  )


# Attach to imputed long-format data
RPNPF_filtered <- RPNPF_filtered %>%
  left_join(compound_missingness_RPNPF %>% select(compound_name, missing_category),
            by = "compound_name") %>% 
  filter(missing_category == "M2") %>% 
  mutate(M2_Detected = ifelse(!is.na(intensity) & intensity != 0, 1, 0))



## Combine the continuous and categorical features =====


## ----------  Binary features  ----------
bin_wide_RPNPF <- RPNPF_filtered %>% 
  select(SampleID, compound_name, M2_Detected, sample_group) %>% 
  pivot_wider(names_from  = compound_name,
              values_from = M2_Detected,
              names_prefix = "bin_")  

## ----------  Box-Cox intensities  ----------
bc_wide_RPNPF <- RPNPF_imp_long %>% 
  select(SampleID, compound_name, intensity_bc, sample_group) %>% 
  pivot_wider(names_from  = compound_name,
              values_from = intensity_bc,
              names_prefix = "bc_")


full_RPNPF <- bin_wide_RPNPF %>% 
  full_join(bc_wide_RPNPF, by = c("SampleID", "sample_group")) %>% 
  mutate(SampleID = tolower(SampleID)) %>% 
  left_join(combined_info, by = c("SampleID" = "meta_sampleid")) %>% 
  mutate(meta_sample_group = sample_group,
         meta_SampleID = SampleID) %>% 
  select(-sample_group, -SampleID)

#saveRDS(full_RPNPF, file = "../datamine_results/full_RPNPF_matrix.rds")

# ======== Combine ZHP and RPNPF ========

# — helper to build each platform’s full matrix with platform-specific prefixes —  
make_full_matrix <- function(imp_long, filtered, platform) {
  # define your prefixes
  bin_pref <- paste0("bin_", tolower(platform), "_")
  bc_pref  <- paste0("bc_",  tolower(platform), "_")
  
  # binary (M2) features  
  bin_wide <- filtered %>% 
    select(SampleID, compound_name, M2_Detected, sample_group) %>% 
    mutate(SampleID = tolower(SampleID)) %>% 
    pivot_wider(names_from  = compound_name,
                values_from = M2_Detected,
                names_prefix = bin_pref)
  
  # continuous (M1) features  
  bc_wide <- imp_long %>% 
    select(SampleID, compound_name, intensity_bc, sample_group) %>% 
    mutate(SampleID = tolower(SampleID)) %>% 
    pivot_wider(names_from  = compound_name,
                values_from = intensity_bc,
                names_prefix = bc_pref)
  
  # join features + attach metadata  
  full_mat <- bin_wide %>% 
    full_join(bc_wide, by = c("SampleID", "sample_group")) %>% 
    left_join(combined_info,
              by = c("SampleID" = "meta_sampleid")) %>% 
    mutate(
      meta_SampleID      = SampleID,
      meta_sample_group  = sample_group
    ) %>% 
    select(starts_with("meta_"), everything(), -SampleID, -sample_group)
  
  return(full_mat)
}

# — build each —  
full_ZHP    <- make_full_matrix(ZHP_imp_long,    ZHP_filtered,    "ZHP")
full_RPNPF  <- make_full_matrix(RPNPF_imp_long, RPNPF_filtered, "RPNPF")




# — now merge into one wide matrix —  
# join on all of your meta_ columns; features have distinct prefixes so won’t collide  
full_combined <- full_ZHP %>% 
  left_join(full_RPNPF, 
            by = intersect(names(full_ZHP), names(full_RPNPF)))

# Summary table for EDA:
race_table <- full_combined %>% 
  select(meta_race_ethnicity, meta_sample_group)

# - Simplify race/ethnicity for one-hot encoding downstream
full_combined <- full_combined %>%  mutate(
  meta_race_ethnicity = case_when(
    meta_race_ethnicity == "White"    ~ "White",
    meta_race_ethnicity == "Hispanic" ~ "Hispanic",
    TRUE                        ~ "Other"
  )
)

# — save out —  
saveRDS(full_combined, "../res/datamine/full_combined_matrix.rds")