# ==============================================================================
# Script Name: QRILC_imputation.R
# Author: Steven Brooks
# Date: 2025-06-04
  # Description: Perform QRILC on the least 20% missing values.
# The idea is to:
# 1) log2(x+1) the raw intensities.
# 2) Perform QRLIC on the log-transformed intensities.
# 3) Invert the log transformation on the imputed intensities (2^QRLIC(x))
# 4) Box-cox the 'raw' imputed intensities.


### Also, this script assumes both RPNPF and ZHP do not have zero values.
# We have to change the logic for the unannotated analysis.
# ==============================================================================

# ==== Setup ===================================================================
library(tidyverse)
library(readxl)
library(imputeLCMD)
library(writexl)


# ==== Load data ===============================================================

ZHP <- readRDS(file = "../data/R_objects/ZHP_annotated_filtered.rds")
RPNPF <- readRDS(file = "../data/R_objects/RPNPF_annotated_filtered.rds")

# ==== Process ZHP =============================================================

## ==== QRILC imputation on compounds with <20% missingness =====================

ZHP_filtered <- ZHP %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  mutate(sample_group = factor(sample_group, levels = c("Control", "JMML")))


compound_missingness <- ZHP_filtered %>%
  group_by(compound_name) %>%
  summarize(
    total = n(),
    zeros = sum(intensity == 0, na.rm = TRUE), 
    NAs = sum(is.na(intensity)),
    missing = sum(is.na(intensity) | intensity == 0),
    missing_pct = missing / total
  ) %>%
  filter(missing_pct < 0.2)


ZHP_clean <- ZHP_filtered %>%
  filter(compound_name %in% compound_missingness$compound_name)

ZHP_wide <- ZHP_clean %>%
  select(SampleID, compound_name, intensity) %>%
  mutate(
    intensity_log2 = if_else(
      is.na(intensity) | intensity == 0,
      NA_real_,
      log2(intensity + 1)
    )
  ) %>% 
  select(SampleID, compound_name, intensity_log2) %>%
  pivot_wider(names_from = compound_name, values_from = intensity_log2)


# Drop SampleID for imputation
intensity_matrix <- ZHP_wide %>% select(-SampleID) %>% as.matrix()


# Apply QRLIC imputation
set.seed(123)
imputed_matrix <- impute.QRILC(intensity_matrix)[[1]]

# Set any negative imputed values to zero
imputed_matrix[imputed_matrix < 0] <- 0


# Combine back with SampleID
ZHP_imputed <- bind_cols(SampleID = ZHP_wide$SampleID, as.data.frame(imputed_matrix))

inv_ZHP_imputed <- ZHP_imputed %>%
  mutate(across(-SampleID, ~ 2^.x - 1))


## ====  Box-Cox Transform ZHP  =================================================

# Prepare sample group info
group_info <- ZHP_filtered %>%
  filter(SampleID %in% ZHP_imputed$SampleID) %>%
  distinct(SampleID, sample_group)

# Join group info to imputed data
imputed_long_ZHP <- inv_ZHP_imputed %>%
  pivot_longer(-SampleID, names_to = "compound_name", values_to = "intensity") %>%
  left_join(group_info, by = "SampleID") %>%
  filter(sample_group %in% c("JMML", "Control"))  %>% 
  mutate(sample_group = factor(sample_group, levels = c("Control", "JMML")))
  

# Box-cox transform before t.test:
lambdas <- imputed_long_ZHP %>%
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
ZHP_transformed <- imputed_long_ZHP %>%
  left_join(lambdas, by = "compound_name") %>%
  mutate(
    intensity_bc = case_when(
      is.na(intensity) | intensity <= 0 | is.na(lambda) ~ NA_real_,
      abs(lambda) < 1e-8 ~ log(intensity),
      TRUE ~ (intensity^lambda - 1) / lambda
    )
  )


#Save the imputed raw intensities, log2() intensities, and box-cox intensities
ZHP_transformed <- ZHP_transformed %>% 
  mutate(intensity_log2 = log2(intensity + 1))

## ====  t.test ZHP  ============================================================

# Calculate means by group
group_means <- ZHP_transformed %>%
  group_by(compound_name, sample_group) %>%
  summarize(mean_intensity = mean(intensity_log2), .groups = "drop") %>%
  pivot_wider(names_from = sample_group, values_from = mean_intensity, names_prefix = "mean_")

# Run t-tests
ttest_results_ZHP <- ZHP_transformed %>%
  group_by(compound_name) %>%
  summarize(
    t_test = list(t.test(intensity_bc ~ sample_group)),
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),
    p_value = map_dbl(t_test, ~ .x$p.value)
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means, by = "compound_name") %>% 
  mutate(log2_fc = mean_JMML - mean_Control)

## ------- ZHP PCA Before imputation --------------------------------------------
# PCA Validation of Normalization and Imputation Procedures
#
# PCA is performed on both unimputed and imputed datasets to assess the impact
# of data preprocessing on overall structure and variance.
#
# Step 1: Normalize the unimputed data
# - Normalization is applied to raw intensity values (e.g., log2, Box-Cox).
# - PCA is run with samples as rows and metabolites as columns.
# - This serves as a baseline before any imputation is applied.

# Box-cox transform before t.test:
lambdas_unimp <- ZHP_clean %>%
  group_by(compound_name) %>%
  summarise(
    lambda = {
      model    <- lm(intensity ~ 1, data = cur_data())
      bc       <- MASS::boxcox(model, lambda = seq(-2, 2, 0.1), plotit = FALSE)
      bc$x[which.max(bc$y)]
    },
    .groups = "drop"
  )

# 2. Join λ and apply Box–Cox only when intensity >  0
ZHP_unimp_transformed <- ZHP_clean %>%
  left_join(lambdas_unimp, by = "compound_name") %>%
  mutate(
    intensity_bc = case_when(
      # missing or non-positive intensities → NA
      is.na(intensity) | intensity <= 0 | is.na(lambda)    ~ NA_real_,
      # λ ≈ 0 → log transform
      abs(lambda) < 1e-8                                    ~ log(intensity),
      # otherwise Box–Cox
      TRUE                                                  ~ (intensity^lambda - 1) / lambda
    )
  )


# Step 2: Impute missing values
# - Missing values are imputed using [insert method, e.g., KNN or QRILC].
# Step 3: Normalize the imputed data (if required by method)
# - Normalization is re-applied if the imputation method alters data scale.

# See above code: QRILC + Box-Cox
head(ZHP_transformed)


# Step 4: PCA on the imputed dataset
# - PCA is rerun on the normalized, imputed data.


# 1. Pivot to wide, then make SampleID the rownames before dropping metadata
ZHP_unimp_pca_mat <- ZHP_unimp_transformed %>%
  select(SampleID, sample_group, compound_name, intensity_bc) %>%
  pivot_wider(
    names_from  = compound_name,
    values_from = intensity_bc
  ) %>%
  column_to_rownames(var = "SampleID") %>%   # ← move SampleID into rownames
  select(-sample_group) %>%                  # ← drop sample_group once rownames set
  as.matrix()

# 2. Clean columns with any NA or infinite
keep_cols <- colSums(is.na(ZHP_unimp_pca_mat) | is.infinite(ZHP_unimp_pca_mat)) == 0
pca_unimp_matrix_clean <- ZHP_unimp_pca_mat[, keep_cols]

# 3. Run PCA
pca_unimp_res <- prcomp(pca_unimp_matrix_clean, center = TRUE, scale. = TRUE)

# 4. Extract scores (they’ll carry rownames = your SampleIDs)
pca_scores_unimp <- as.data.frame(pca_unimp_res$x) %>%
  rownames_to_column(var = "SampleID")

# 5. Re-attach sample_group (if you need it for plotting)
pca_scores_unimp <- pca_scores_unimp %>%
  left_join(
    ZHP_unimp_transformed %>% 
      distinct(SampleID, sample_group),
    by = "SampleID"
  )


p_before_ZHP <- ggplot(pca_scores_unimp, aes(PC1, PC2, color = sample_group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", round(summary(pca_unimp_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_unimp_res)$importance[2,2]*100,1), "%)"),
    color = "Sample Group",
    title = "PCA Before QRILC Imputation\non Box-Cox Transformed Intensities (ZHP)"
  ) +
  theme_minimal(base_size = 12)

print(p_before_ZHP)


if (!dir.exists("../res/annotated/QC")) dir.create("../res/annotated/QC", recursive = TRUE)

# 5. Save plot to PNG
ggsave(
  filename = "../res/annotated/QC/PCA_boxcox_before_imputation_ZHP.png",
  plot     = p_before_ZHP,
  width    = 6,
  height   = 5,
  dpi      = 300
)

# Purpose:
# - Compare clustering and variance explained (PC1, PC2, etc.) before and after imputation.
# - Assess whether group structure (e.g., JMML vs. Control) is preserved or clarified.
# - Ensure imputation does not introduce artificial patterns or mask biological signal.

## ------- ZHP PCA After imputation ---------------------------------------------


#  Pivot to wide, then make SampleID the rownames before dropping metadata
ZHP_imp_pca_mat <- ZHP_transformed %>%
  select(SampleID, sample_group, compound_name, intensity_bc) %>%
  pivot_wider(
    names_from  = compound_name,
    values_from = intensity_bc
  ) %>%
  column_to_rownames(var = "SampleID") %>%   # ← move SampleID into rownames
  select(-sample_group) %>%                  # ← drop sample_group once rownames set
  as.matrix()

#  Clean columns with any NA or infinite
keep_cols <- colSums(is.na(ZHP_imp_pca_mat) | is.infinite(ZHP_imp_pca_mat)) == 0
pca_imp_matrix_clean <- ZHP_imp_pca_mat[, keep_cols]

#  Run PCA
pca_imp_res <- prcomp(pca_imp_matrix_clean, center = TRUE, scale. = TRUE)

# 4. Extract scores (they’ll carry rownames = your SampleIDs)
pca_scores_imp <- as.data.frame(pca_imp_res$x) %>%
  rownames_to_column(var = "SampleID")

# 5. Re-attach sample_group (if you need it for plotting)
pca_scores_imp <- pca_scores_imp %>%
  left_join(
    ZHP_transformed %>% 
      distinct(SampleID, sample_group),
    by = "SampleID"
  )


p_after_ZHP <- ggplot(pca_scores_imp, aes(PC1, PC2, color = sample_group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", round(summary(pca_imp_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_imp_res)$importance[2,2]*100,1), "%)"),
    color = "Sample Group",
    title = "PCA After QRILC Imputation\nand Box-Cox Transformation Intensities (ZHP)"
  ) +
  theme_minimal(base_size = 12)

print(p_after_ZHP)

# 5. Save plot to PNG
ggsave(
  filename = "../res/annotated/QC/PCA_boxcox_after_imputation_ZHP.png",
  plot     = p_after_ZHP,
  width    = 6,
  height   = 5,
  dpi      = 300
)



# ==== Process RPNPF =============================================================

## ==== QRILC imputation on compounds with <20% missingness =====================

RPNPF_filtered <- RPNPF %>%
  filter(sample_group %in% c("JMML", "Control")) %>% 
  mutate(sample_group = factor(sample_group, levels = c("Control", "JMML")))


compound_missingness <- RPNPF_filtered %>%
  group_by(compound_name) %>%
  summarize(
    total = n(),
    zeros = sum(intensity == 0, na.rm = TRUE), 
    NAs = sum(is.na(intensity)),
    missing = sum(is.na(intensity) | intensity == 0),
    missing_pct = missing / total
  ) %>%
  filter(missing_pct < 0.2)


RPNPF_clean <- RPNPF_filtered %>%
  filter(compound_name %in% compound_missingness$compound_name)

RPNPF_wide <- RPNPF_clean %>%
  select(SampleID, compound_name, intensity) %>%
  mutate(
    intensity_log2 = if_else(
      is.na(intensity) | intensity == 0,
      NA_real_,
      log2(intensity + 1)
    )
  ) %>% 
  select(SampleID, compound_name, intensity_log2) %>%
  pivot_wider(names_from = compound_name, values_from = intensity_log2)


# Drop SampleID for imputation
intensity_matrix <- RPNPF_wide %>% select(-SampleID) %>% as.matrix()


# Apply QRLIC imputation
set.seed(123)
imputed_matrix <- impute.QRILC(intensity_matrix)[[1]]

# Set any negative imputed values to zero
imputed_matrix[imputed_matrix < 0] <- 0


# Combine back with SampleID
RPNPF_imputed <- bind_cols(SampleID = RPNPF_wide$SampleID, as.data.frame(imputed_matrix))

inv_RPNPF_imputed <- RPNPF_imputed %>%
  mutate(across(-SampleID, ~ 2^.x - 1))


## ====  Box-Cox Transform RPNPF  =================================================

# Prepare sample group info
group_info <- RPNPF_filtered %>%
  filter(SampleID %in% RPNPF_imputed$SampleID) %>%
  distinct(SampleID, sample_group)

# Join group info to imputed data
imputed_long_RPNPF <- inv_RPNPF_imputed %>%
  pivot_longer(-SampleID, names_to = "compound_name", values_to = "intensity") %>%
  left_join(group_info, by = "SampleID") %>%
  filter(sample_group %in% c("JMML", "Control"))  %>% 
  mutate(sample_group = factor(sample_group, levels = c("Control", "JMML")))


# Box-cox transform before t.test:
lambdas <- imputed_long_RPNPF %>%
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
RPNPF_transformed <- imputed_long_RPNPF %>%
  left_join(lambdas, by = "compound_name") %>%
  mutate(
    intensity_bc = case_when(
      is.na(intensity) | intensity <= 0 | is.na(lambda) ~ NA_real_,
      abs(lambda) < 1e-8 ~ log(intensity),
      TRUE ~ (intensity^lambda - 1) / lambda
    )
  )


#Save the imputed raw intensities, log2() intensities, and box-cox intensities
RPNPF_transformed <- RPNPF_transformed %>% 
  mutate(intensity_log2 = log2(intensity + 1))

## ====  t.test RPNPF  ============================================================

# Calculate means by group
group_means <- RPNPF_transformed %>%
  group_by(compound_name, sample_group) %>%
  summarize(mean_intensity = mean(intensity_log2), .groups = "drop") %>%
  pivot_wider(names_from = sample_group, values_from = mean_intensity, names_prefix = "mean_")

# Run t-tests
ttest_results_RPNPF <- RPNPF_transformed %>%
  group_by(compound_name) %>%
  summarize(
    t_test = list(t.test(intensity_bc ~ sample_group)),
    .groups = "drop"
  ) %>%
  mutate(
    t_statistic = map_dbl(t_test, ~ .x$statistic),
    p_value = map_dbl(t_test, ~ .x$p.value)
  ) %>%
  select(compound_name, t_statistic, p_value) %>%
  left_join(group_means, by = "compound_name") %>% 
  mutate(log2_fc = mean_JMML - mean_Control)

## ------- RPNPF PCA Before imputation --------------------------------------------
# PCA Validation of Normalization and Imputation Procedures
#
# PCA is performed on both unimputed and imputed datasets to assess the impact
# of data preprocessing on overall structure and variance.
#
# Step 1: Normalize the unimputed data
# - Normalization is applied to raw intensity values (e.g., log2, Box-Cox).
# - PCA is run with samples as rows and metabolites as columns.
# - This serves as a baseline before any imputation is applied.

# Box-cox transform before t.test:
lambdas_unimp <- RPNPF_clean %>%
  group_by(compound_name) %>%
  summarise(
    lambda = {
      model    <- lm(intensity ~ 1, data = cur_data())
      bc       <- MASS::boxcox(model, lambda = seq(-2, 2, 0.1), plotit = FALSE)
      bc$x[which.max(bc$y)]
    },
    .groups = "drop"
  )

# 2. Join λ and apply Box–Cox only when intensity >  0
RPNPF_unimp_transformed <- RPNPF_clean %>%
  left_join(lambdas_unimp, by = "compound_name") %>%
  mutate(
    intensity_bc = case_when(
      # missing or non-positive intensities → NA
      is.na(intensity) | intensity <= 0 | is.na(lambda)    ~ NA_real_,
      # λ ≈ 0 → log transform
      abs(lambda) < 1e-8                                    ~ log(intensity),
      # otherwise Box–Cox
      TRUE                                                  ~ (intensity^lambda - 1) / lambda
    )
  )


# Step 2: Impute missing values
# - Missing values are imputed using [insert method, e.g., KNN or QRILC].
# Step 3: Normalize the imputed data (if required by method)
# - Normalization is re-applied if the imputation method alters data scale.

# See above code: QRILC + Box-Cox
head(RPNPF_transformed)


# Step 4: PCA on the imputed dataset
# - PCA is rerun on the normalized, imputed data.


# 1. Pivot to wide, then make SampleID the rownames before dropping metadata
RPNPF_unimp_pca_mat <- RPNPF_unimp_transformed %>%
  select(SampleID, sample_group, compound_name, intensity_bc) %>%
  pivot_wider(
    names_from  = compound_name,
    values_from = intensity_bc
  ) %>%
  column_to_rownames(var = "SampleID") %>%   # ← move SampleID into rownames
  select(-sample_group) %>%                  # ← drop sample_group once rownames set
  as.matrix()

# 2. Clean columns with any NA or infinite
keep_cols <- colSums(is.na(RPNPF_unimp_pca_mat) | is.infinite(RPNPF_unimp_pca_mat)) == 0
pca_unimp_matrix_clean <- RPNPF_unimp_pca_mat[, keep_cols]

# 3. Run PCA
pca_unimp_res <- prcomp(pca_unimp_matrix_clean, center = TRUE, scale. = TRUE)

# 4. Extract scores (they’ll carry rownames = your SampleIDs)
pca_scores_unimp <- as.data.frame(pca_unimp_res$x) %>%
  rownames_to_column(var = "SampleID")

# 5. Re-attach sample_group (if you need it for plotting)
pca_scores_unimp <- pca_scores_unimp %>%
  left_join(
    RPNPF_unimp_transformed %>% 
      distinct(SampleID, sample_group),
    by = "SampleID"
  )


p_before_RPNPF <- ggplot(pca_scores_unimp, aes(PC1, PC2, color = sample_group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", round(summary(pca_unimp_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_unimp_res)$importance[2,2]*100,1), "%)"),
    color = "Sample Group",
    title = "PCA Before QRILC Imputation\non Box-Cox Transformed Intensities (RPNPF)"
  ) +
  theme_minimal(base_size = 12)

print(p_before_RPNPF)

# 5. Save plot to PNG
ggsave(
  filename = "../res/annotated/QC/PCA_boxcox_before_imputation_RPNPF.png",
  plot     = p_before_RPNPF,
  width    = 6,
  height   = 5,
  dpi      = 300
)

# Purpose:
# - Compare clustering and variance explained (PC1, PC2, etc.) before and after imputation.
# - Assess whether group structure (e.g., JMML vs. Control) is preserved or clarified.
# - Ensure imputation does not introduce artificial patterns or mask biological signal.

## ------- RPNPF PCA After imputation ---------------------------------------------


#  Pivot to wide, then make SampleID the rownames before dropping metadata
RPNPF_imp_pca_mat <- RPNPF_transformed %>%
  select(SampleID, sample_group, compound_name, intensity_bc) %>%
  pivot_wider(
    names_from  = compound_name,
    values_from = intensity_bc
  ) %>%
  column_to_rownames(var = "SampleID") %>%   # ← move SampleID into rownames
  select(-sample_group) %>%                  # ← drop sample_group once rownames set
  as.matrix()

#  Clean columns with any NA or infinite
keep_cols <- colSums(is.na(RPNPF_imp_pca_mat) | is.infinite(RPNPF_imp_pca_mat)) == 0
pca_imp_matrix_clean <- RPNPF_imp_pca_mat[, keep_cols]

#  Run PCA
pca_imp_res <- prcomp(pca_imp_matrix_clean, center = TRUE, scale. = TRUE)

# 4. Extract scores (they’ll carry rownames = your SampleIDs)
pca_scores_imp <- as.data.frame(pca_imp_res$x) %>%
  rownames_to_column(var = "SampleID")

# 5. Re-attach sample_group (if you need it for plotting)
pca_scores_imp <- pca_scores_imp %>%
  left_join(
    RPNPF_transformed %>% 
      distinct(SampleID, sample_group),
    by = "SampleID"
  )


p_after_RPNPF <- ggplot(pca_scores_imp, aes(PC1, PC2, color = sample_group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", round(summary(pca_imp_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_imp_res)$importance[2,2]*100,1), "%)"),
    color = "Sample Group",
    title = "PCA After QRILC Imputation\nand Box-Cox Transformation Intensities (RPNPF)"
  ) +
  theme_minimal(base_size = 12)

print(p_after_RPNPF)

# 5. Save plot to PNG
ggsave(
  filename = "../res/annotated/QC/PCA_boxcox_after_imputation_RPNPF.png",
  plot     = p_after_RPNPF,
  width    = 6,
  height   = 5,
  dpi      = 300
)





# ====  Combining Above  =======================================================

#Load fisher:
#source("02_fishers.R")


# Add platform and test labels
fisher_results_ZHP <- fisher_results_ZHP %>%
  mutate(platform = "ZHP", test_type = "Fisher")

fisher_results_RPNPF <- fisher_results_RPNPF %>%
  mutate(platform = "RPNPF", test_type = "Fisher")

ttest_results_ZHP <- ttest_results_ZHP %>%
  mutate(platform = "ZHP", test_type = "TTest")

ttest_results_RPNPF <- ttest_results_RPNPF %>%
  mutate(platform = "RPNPF", test_type = "TTest")


# Rename t-test p-values to match Fisher column names
ttest_results_ZHP <- ttest_results_ZHP %>%
  rename(effect = t_statistic, raw_p = p_value)

ttest_results_RPNPF <- ttest_results_RPNPF %>%
  rename(effect = t_statistic, raw_p = p_value)

fisher_results_ZHP <- fisher_results_ZHP %>%
  rename(effect = odds_ratio, raw_p = fisher_pval)

fisher_results_RPNPF <- fisher_results_RPNPF %>%
  rename(effect = odds_ratio, raw_p = fisher_pval)


# Combine all
combined_results <- bind_rows(
  ttest_results_ZHP,
  ttest_results_RPNPF,
  fisher_results_ZHP,
  fisher_results_RPNPF
)

# Global FDR adjustment
combined_results <- combined_results %>%
  mutate(global_fdr = p.adjust(raw_p, method = "BH"))

print(combined_results %>% filter(global_fdr < .05))

write_xlsx(combined_results, "../res/annotated/final_results_global_fdr.xlsx")

saveRDS(ZHP_transformed, "../data/R_objects/ZHP_imputed.rds")
saveRDS(RPNPF_transformed, "../data/R_objects/RPNPF_imputed.rds")

