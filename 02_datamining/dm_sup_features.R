# Feature Selection / Dimensionality reduction
# ==== Load Libraries ====
library(tidyverse)
library(limma)
library(vip)
library(randomForest)
library(glmnet)

# ==== combined Analysis ====
full_combined <- readRDS(file = "../res/datamine/full_combined_matrix.rds")

# For each platform we have a matrix (full_combined, full_RPNPF) with samples as rows.
# 
# We have three sets of features in the column space: 
#   
#   - The imputed boxcox intensities (bc_ prefix) for the metabolites with less than 20% missingness.
# 
# - The binary detected or not detected metabolites (prefix bin_) in the 20% - 65% missingness. 
# 
# - The sample metadata (prefix meta_), which includes: sex, race/ethnicity, age, and blood collection time.


## ==== 3.1 Rank metabolites by class‐separation power ====

### ==== Identify Metabolite Features ====
metab_features_combined <-  grep("^bc_", colnames(full_combined), value = TRUE)

### ==== Extract Class Labels ====
group_combined <- factor(full_combined$meta_sample_group, levels = c("Control", "JMML"))
design_combined <- model.matrix(~ 0 + group_combined)
colnames(design_combined) <- levels(group_combined)

### ==== Build Expression Matrix (features × samples) ====
expr_matrix_combined <- t(as.matrix(full_combined[, metab_features_combined]))

### ==== Fit Linear Model ====
fit_combined <- lmFit(expr_matrix_combined, design_combined)

### ==== Define 2-Class Contrast ====
contrast_matrix_combined <- makeContrasts(Diff = JMML - Control, levels = design_combined)
fit2_combined <- contrasts.fit(fit_combined, contrast_matrix_combined)
fit2_combined <- eBayes(fit2_combined)

### ==== Extract and Rank Results ====
results_combined <- topTable(fit2_combined, number = Inf, sort.by = "P") %>%
  rownames_to_column("metabolite") %>%
  rename(
    pval = P.Value,
    fdr  = adj.P.Val
  )

### ==== Save Ranked Metabolites ====
write_csv(results_combined, "../res/datamine/ranked_metabolites_combined.csv")


## ==== 3.2 Multivariate importance on full dataset ====


# 1) identify features and response
metab_features <- full_combined %>%
  select(-c(meta_SampleID, meta_sample_group, meta_bd_blood_collection_time_x)) %>% 
  colnames(.)


X <- full_combined[, metab_features]

y <- factor(full_combined$meta_sample_group, levels = c("Control", "JMML"))

set.seed(123)

# 2) fit RF with importance turned on
rf_model <- randomForest(
  x          = X,
  y          = y,
  importance = TRUE
)

# === OPTION A: extract built-in OOB permutation importance ===
imp_oob <- importance(rf_model, type = 1)        # type=1 is MeanDecreaseAccuracy
imp_df <- tibble(
  Feature    = rownames(imp_oob),
  Importance = imp_oob[, 1]
) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 30)

p1 <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "",
    x     = NULL,
    y     = "Mean Decrease in Accuracy"
  )

ggsave(
  filename = "../res/datamine/combined_rf_oob_importance.png",
  plot     = p1,
  width    = 7,
  height   = 5,
  dpi      = 300
)

