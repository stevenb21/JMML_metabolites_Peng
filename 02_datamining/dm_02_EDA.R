# ==== Load Libraries & Data ====
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)


# ==== Box-cox + Binary Analysis ====
full_boxcox_bin <- readRDS(file ="../res/datamine/full_combined_matrix.rds")

## ==== Prepare Data ====
# # Drop metadata columns, retain intensity (no prefix) + binary (bin_) features
metab_boxcox_bin <- full_boxcox_bin %>% 
  column_to_rownames("meta_SampleID") %>% 
  select(-starts_with("meta_"))    # drop all the meta_ columns
# that leaves the intensities (bc_) + bin_ columns

# run PCA, centering and scaling everything
pca_boxcox_bin <- metab_boxcox_bin %>% 
  prcomp(center = TRUE, scale. = TRUE)


# extract scores and join back on the metadata  sample_group
scores_boxcox_bin <- as_tibble(pca_boxcox_bin$x, rownames = "meta_SampleID") %>% 
  left_join(
    full_boxcox_bin %>% select(starts_with("meta_")), 
    by = "meta_SampleID"
  )


## ==== PCA  ====
### ==== PCA of Case status ====

p <- ggplot(scores_boxcox_bin, aes(PC1, PC2, color = meta_sample_group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("JMML" = "#E74C3C", "Control" = "#3498DB")) +
  labs(x = paste0("PC1 (", round(summary(pca_boxcox_bin)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_boxcox_bin)$importance[2,2] * 100, 1), "%)"),
       color = "") +
  theme_minimal()

# save to .png 
ggsave("../res/datamine/boxcox_bin_pca_case.png", p, width = 6, height = 5, dpi = 300)



### ==== PCA of Sex ====

# 4. plot PC1 vs PC2, coloring by sex, shaping by cohort, etc.
p <- ggplot(scores_boxcox_bin, aes(PC1, PC2, color = meta_sex)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("F" = "#E67E22", "M" = "#1ABC9C")) +
  labs(x = paste0("PC1 (", round(summary(pca_boxcox_bin)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_boxcox_bin)$importance[2,2] * 100, 1), "%)"),
       color = "Sex") +
  theme_minimal()


# save to .png:
ggsave("../res/datamine/boxcox_bin_pca_sex.png", p, width = 6, height = 5, dpi = 300)


## ======= Heatmaps of top 50 variable metabolites ====

# Compute variance per feature
feature_variance_boxcox_bin <- metab_boxcox_bin %>% 
  summarise(across(everything(), var, na.rm = TRUE)) %>% 
  pivot_longer(cols = everything(), names_to = "feature", values_to = "variance")

# Select top 50 most variable features
top50_features_boxcox_bin <- feature_variance_boxcox_bin %>% 
  arrange(desc(variance)) %>% 
  slice_head(n = 50) %>% 
  pull(feature)

# Subset and scale the matrix for heatmap
heatmap_matrix_boxcox_bin <- metab_boxcox_bin %>% 
  select(all_of(top50_features_boxcox_bin)) %>% 
  scale()

# Create annotation dataframe for side colors
annotation_boxcox_bin <- full_boxcox_bin %>% 
  select(meta_SampleID, meta_sample_group, meta_sex) %>% 
  column_to_rownames("meta_SampleID")

# Ensure annotation rows match matrix
annotation_boxcox_bin <- annotation_boxcox_bin[rownames(heatmap_matrix_boxcox_bin), , drop = FALSE]

# Define color palettes
ann_colors <- list(
  meta_sample_group = c(JMML = "#E74C3C", Control = "#3498DB"),
  meta_sex = c(F = "#E67E22", M = "#1ABC9C")
)

# Generate and save heatmap
pheatmap(
  mat = heatmap_matrix_boxcox_bin,
  annotation_row = annotation_boxcox_bin,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
  show_rownames = FALSE,
  fontsize_col = 8,
  filename = "../res/datamine/boxcox_bin_heatmap_top50.png",
  width = 8,
  height = 10
)

