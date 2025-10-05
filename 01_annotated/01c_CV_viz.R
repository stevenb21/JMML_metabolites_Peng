# 05_cv_distribution_viz.R
# ─────────────────────────
# Visualize CV of Box–Cox–transformed intensities, showing kept vs filtered-out compounds

# assumes this just ran:
#source("01a_save_data_rds_lilly_bc.R")

# ==== ZHP  =========================================

# ==== Plot Coefficient of Variation =========================================

# Create a histogram to visualize the distribution of CV values
p_cv_hist_JMML_ZHP  <- ggplot(cv_tbl_ZHP %>% filter(sample_group == "JMML"), aes(x = cv)) +
  geom_histogram(binwidth = 0.01, fill =  "#E74C3C", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "blue", size = 1) +  # Vertical line at CV = 0.3
  labs(title = "Histogram of Coefficient of Variation (CV) Before Filtering",
       subtitle = "JMML Samples | Box-Cox Transformed | ZHP",
       x = "Coefficient of Variation (CV)",
       y = "Count") +
  theme_bw()


# Save the histogram
ggsave(
  filename = "../res/annotated/QC/histogram_cv_before_filtering_JMML_ZHP.png",
  plot = p_cv_hist_JMML_ZHP,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)


# ==== Plot Coefficient of Variation =========================================

# Create a histogram to visualize the distribution of CV values
p_cv_hist_Control_ZHP  <- ggplot(cv_tbl_ZHP %>% filter(sample_group == "Control"), aes(x = cv)) +
  geom_histogram(binwidth = 0.01, fill =  "#3498DB", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "blue", size = 1) +  # Vertical line at CV = 0.3
  labs(title = "Histogram of Coefficient of Variation (CV) Before Filtering",
       subtitle = "Control Samples | Box-Cox Transformed | ZHP",
       x = "Coefficient of Variation (CV)",
       y = "Count") +
  theme_bw()


# Save the histogram
ggsave(
  filename = "../res/annotated/QC/histogram_cv_before_filtering_Control_ZHP.png",
  plot = p_cv_hist_Control_ZHP,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)



## ==== RPNPF  =========================================

# ==== Plot Coefficient of Variation =========================================

# Create a histogram to visualize the distribution of CV values
p_cv_hist_JMML_RPNPF  <- ggplot(cv_tbl_RPNPF %>% filter(sample_group == "JMML"), aes(x = cv)) +
  geom_histogram(binwidth = 0.01, fill =  "#E74C3C", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "blue", size = 1) +  # Vertical line at CV = 0.3
  labs(title = "Histogram of Coefficient of Variation (CV) Before Filtering",
       subtitle = "JMML Samples | Box-Cox Transformed | RPNPF",
       x = "Coefficient of Variation (CV)",
       y = "Count") +
  theme_bw()


# Save the histogram
ggsave(
  filename = "../res/annotated/QC/histogram_cv_before_filtering_JMML_RPNPF.png",
  plot = p_cv_hist_JMML_RPNPF,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)


# ==== Plot Coefficient of Variation =========================================

# Create a histogram to visualize the distribution of CV values
p_cv_hist_Control_RPNPF  <- ggplot(cv_tbl_RPNPF %>% filter(sample_group == "Control"), aes(x = cv)) +
  geom_histogram(binwidth = 0.01, fill =  "#3498DB", color = "black", alpha = 0.8) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "blue", size = 1) +  # Vertical line at CV = 0.3
  labs(title = "Histogram of Coefficient of Variation (CV) Before Filtering",
       subtitle = "Control Samples | Box-Cox Transformed | RPNPF",
       x = "Coefficient of Variation (CV)",
       y = "Count") +
  theme_bw()


# Save the histogram
ggsave(
  filename = "../res/annotated/QC/histogram_cv_before_filtering_Control_RPNPF.png",
  plot = p_cv_hist_Control_RPNPF,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)
