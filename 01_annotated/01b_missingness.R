# Missingness ZHP
library(naniar)

# assumes this just ran:
#source("01a_save_data_rds_lilly_bc.R")

# ==== ZHP Missingness =========================================================

# Filter and reshape
ZHP_wide <- ZHP_long %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  pivot_wider(
    names_from = compound_name,
    values_from = intensity
  )

# Reorder columns to fix display (optional)
ZHP_wide <- ZHP_wide %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

#### ==== Missingness Heatmap ==============================================
p_miss_clean_ZHP <- vis_miss(ZHP_wide, sort_miss = FALSE) +
  labs(
    title = "ZHP Metabolite Missingness (JMML + Control)",
    subtitle = "Each row = sample, each column = metabolite",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  )


# Step 4: Save
ggsave(
  filename = "../res/annotated/QC/ZHP_missingness_heatmap.png",
  plot = p_miss_clean_ZHP,
  width = 10,
  height = 6,
  dpi = 300
)



# Recompute missingness
ZHP_missingness_summary <- ZHP_long %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  group_by(compound_name) %>%
  summarise(
    prop_missing_or_zero = mean(is.na(intensity) | intensity == 0),
    .groups = "drop"
  )

# Threshold range
thresholds <- seq(-0.01, 1, by = 0.01)

# Compute number removed at each threshold
ZHP_removed_curve <- tibble(
  threshold = thresholds,
  n_removed = map_int(thresholds, function(t) {
    sum(ZHP_missingness_summary$prop_missing_or_zero > t)
  })
)


M3_tick <- .65 + (1 - .65)/2
M2_tick <- .2 + (.65 - .2)/2
M1_tick <- .1


a <- 0.01/2
# Plot
p_removed_ZHP <- ggplot(ZHP_removed_curve, aes(x = threshold, y = n_removed)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(0.20, 0.65), linetype = "dashed", color = "red") +
  
  # Region annotations
  annotate("text", x = 0.10, y = max(ZHP_removed_curve$n_removed) * 0.75,
           label = "T-Test on\nImputed Values", size = 4, color = "gray20") +
  annotate("text", x = 0.425, y = max(ZHP_removed_curve$n_removed) * 0.5,
           label = "Fisher's Exact Test\non Missingness", size = 4, color = "gray20") +
  annotate("text", x = 0.825, y = max(ZHP_removed_curve$n_removed) * 0.5,
           label = "Removed from Analysis", size = 4, color = "gray20") +
  
  
  # Custom x-axis ticks
  scale_x_continuous(
    name = "Missingness Threshold",
    breaks = c(0, 0.1, 0.2, 0.425, 0.65, 0.825, 1.0),
    labels = c("0", "M1", "0.2", "M2", "0.65", "M3", "1.0")
  ) +
  
  labs(
    title = "Number of Metabolites Classified by Missingness Threshold - ZHP",
    y = "Number of Metabolites Removed"
  ) +
  theme_bw(base_size = 12)




# Save to PNG
ggsave(
  "../res/annotated/QC/ZHP_missingness_threshold_excluded.png",
  plot = p_removed_ZHP,
  width = 9,
  height = 5,
  dpi = 300
)




a <- 0.01/2
# Plot
p_removed_ZHP_ttest <- ggplot(ZHP_removed_curve, aes(x = threshold, y = n_removed)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(0.20, 0.65), linetype = "dashed", color = "red") +
  
  # Region annotations
  annotate("text", x = 0.10, y = max(ZHP_removed_curve$n_removed) * 0.75,
           label = "T-Test on\nImputed Values", size = 4, color = "gray20") +
  annotate("text", x = 0.425, y = max(ZHP_removed_curve$n_removed) * 0.5,
           label = "Fisher's Exact Test\non Missingness", size = 4, color = "gray20") +
  annotate("text", x = 0.825, y = max(ZHP_removed_curve$n_removed) * 0.5,
           label = "Removed from Analysis", size = 4, color = "gray20") +
  
  #Color the part you want
  geom_rect(aes(xmin = 0, xmax = 0.20, ymin = -Inf, ymax = Inf),
              fill = "dodgerblue", alpha = a) +
  # geom_rect(aes(xmin = 0.20, xmax = 0.65, ymin = -Inf, ymax = Inf),
  #           fill = "darkorange", alpha = a) +
  # geom_rect(aes(xmin = 0.65, xmax = 1.0, ymin = -Inf, ymax = Inf),
  #           fill = "gray60", alpha = a) +
  
  # Custom x-axis ticks
  scale_x_continuous(
    name = "Missingness Threshold",
    breaks = c(0, 0.1, 0.2, 0.425, 0.65, 0.825, 1.0),
    labels = c("0", "M1", "0.2", "M2", "0.65", "M3", "1.0")
  ) +
  
  labs(
    title = "Number of Metabolites Classified by Missingness Threshold - ZHP",
    y = "Number of Metabolites Removed"
  ) +
  theme_bw(base_size = 12)


ggsave("../res/annotated/QC/missingness_ZHP_ttest.png", plot = p_removed_ZHP_ttest, width = 9)


# Plot
p_removed_ZHP_fisher <- ggplot(ZHP_removed_curve, aes(x = threshold, y = n_removed)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(0.20, 0.65), linetype = "dashed", color = "red") +
  
  # Region annotations
  annotate("text", x = 0.10, y = max(ZHP_removed_curve$n_removed) * 0.75,
           label = "T-Test on\nImputed Values", size = 4, color = "gray20") +
  annotate("text", x = 0.425, y = max(ZHP_removed_curve$n_removed) * 0.5,
           label = "Fisher's Exact Test\non Missingness", size = 4, color = "gray20") +
  annotate("text", x = 0.825, y = max(ZHP_removed_curve$n_removed) * 0.5,
           label = "Removed from Analysis", size = 4, color = "gray20") +
  
  #Color the part you want
  # geom_rect(aes(xmin = 0, xmax = 0.20, ymin = -Inf, ymax = Inf),
  #           fill = "dodgerblue", alpha = a) +
  geom_rect(aes(xmin = 0.20, xmax = 0.65, ymin = -Inf, ymax = Inf),
            fill = "dodgerblue", alpha = a) +
  # geom_rect(aes(xmin = 0.65, xmax = 1.0, ymin = -Inf, ymax = Inf),
  #           fill = "gray60", alpha = a) +
  
  # Custom x-axis ticks
  scale_x_continuous(
    name = "Missingness Threshold",
    breaks = c(0, 0.1, 0.2, 0.425, 0.65, 0.825, 1.0),
    labels = c("0", "M1", "0.2", "M2", "0.65", "M3", "1.0")
  ) +
  
  labs(
    title = "Number of Metabolites Classified by Missingness Threshold - ZHP",
    y = "Number of Metabolites Removed"
  ) +
  theme_bw(base_size = 12)

ggsave("../res/annotated/QC/missingness_ZHP_fisher.png", plot = p_removed_ZHP_fisher, width = 9)


# ==== RPNPF Missingness =========================================================

  
  # Filter and reshape
  RPNPF_wide <- RPNPF_long %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  pivot_wider(
    names_from = compound_name,
    values_from = intensity
  )

# Reorder columns to fix display (optional)
RPNPF_wide <- RPNPF_wide %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

#### ==== Missingness Heatmap ==============================================
p_miss_clean_RPNPF <- vis_miss(RPNPF_wide, sort_miss = FALSE) +
  labs(
    title = "RPNPF Metabolite Missingness (JMML + Control)",
    subtitle = "Each row = sample, each column = metabolite",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  )


# Save
ggsave(
  filename = "../res/annotated/QC/RPNPF_missingness_heatmap.png",
  plot = p_miss_clean_RPNPF,
  width = 10,
  height = 6,
  dpi = 300
)



# Recompute missingness
RPNPF_missingness_summary <- RPNPF_long %>%
  filter(sample_group %in% c("JMML", "Control")) %>%
  group_by(compound_name) %>%
  summarise(
    prop_missing_or_zero = mean(is.na(intensity) | intensity == 0),
    .groups = "drop"
  )

# Threshold range
thresholds <- seq(-0.01, 1, by = 0.01)

# Compute number removed at each threshold
RPNPF_removed_curve <- tibble(
  threshold = thresholds,
  n_removed = map_int(thresholds, function(t) {
    sum(RPNPF_missingness_summary$prop_missing_or_zero > t)
  })
)


M3_tick <- .65 + (1 - .65)/2
M2_tick <- .2 + (.65 - .2)/2
M1_tick <- .1


a <- 0.01/2
# Plot
p_removed_RPNPF <- ggplot(RPNPF_removed_curve, aes(x = threshold, y = n_removed)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(0.20, 0.65), linetype = "dashed", color = "red") +
  
  # Region annotations
  annotate("text", x = 0.10, y = max(RPNPF_removed_curve$n_removed) * 0.75,
           label = "T-Test on\nImputed Values", size = 4, color = "gray20") +
  annotate("text", x = 0.425, y = max(RPNPF_removed_curve$n_removed) * 0.5,
           label = "Fisher's Exact Test\non Missingness", size = 4, color = "gray20") +
  annotate("text", x = 0.825, y = max(RPNPF_removed_curve$n_removed) * 0.5,
           label = "Removed from Analysis", size = 4, color = "gray20") +

  
  # Custom x-axis ticks
  scale_x_continuous(
    name = "Missingness Threshold",
    breaks = c(0, 0.1, 0.2, 0.425, 0.65, 0.825, 1.0),
    labels = c("0", "M1", "0.2", "M2", "0.65", "M3", "1.0")
  ) +
  
  labs(
    title = "Number of Metabolites Classified by Missingness Threshold - RPNPF",
    y = "Number of Metabolites Removed"
  ) +
  theme_bw(base_size = 12)



# Save to PNG
ggsave(
  "../res/annotated/QC/RPNPF_missingness_threshold_excluded.png",
  plot = p_removed_RPNPF,
  width = 9,
  height = 5,
  dpi = 300
)




a <- 0.01/2
# Plot
p_removed_RPNPF_ttest <- ggplot(RPNPF_removed_curve, aes(x = threshold, y = n_removed)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(0.20, 0.65), linetype = "dashed", color = "red") +
  
  # Region annotations
  annotate("text", x = 0.10, y = max(RPNPF_removed_curve$n_removed) * 0.75,
           label = "T-Test on\nImputed Values", size = 4, color = "gray20") +
  annotate("text", x = 0.425, y = max(RPNPF_removed_curve$n_removed) * 0.5,
           label = "Fisher's Exact Test\non Missingness", size = 4, color = "gray20") +
  annotate("text", x = 0.825, y = max(RPNPF_removed_curve$n_removed) * 0.5,
           label = "Removed from Analysis", size = 4, color = "gray20") +
  
  #Color the part you want
  geom_rect(aes(xmin = 0, xmax = 0.20, ymin = -Inf, ymax = Inf),
            fill = "dodgerblue", alpha = a) +
  
  # Custom x-axis ticks
  scale_x_continuous(
    name = "Missingness Threshold",
    breaks = c(0, 0.1, 0.2, 0.425, 0.65, 0.825, 1.0),
    labels = c("0", "M1", "0.2", "M2", "0.65", "M3", "1.0")
  ) +
  
  labs(
    title = "Number of Metabolites Classified by Missingness Threshold - RPNPF",
    y = "Number of Metabolites Removed"
  ) +
  theme_bw(base_size = 12)



# Plot
p_removed_RPNPF_fisher <- ggplot(RPNPF_removed_curve, aes(x = threshold, y = n_removed)) +
  geom_line(color = "black") +
  geom_vline(xintercept = c(0.20, 0.65), linetype = "dashed", color = "red") +
  
  # Region annotations
  annotate("text", x = 0.10, y = max(RPNPF_removed_curve$n_removed) * 0.75,
           label = "T-Test on\nImputed Values", size = 4, color = "gray20") +
  annotate("text", x = 0.425, y = max(RPNPF_removed_curve$n_removed) * 0.5,
           label = "Fisher's Exact Test\non Missingness", size = 4, color = "gray20") +
  annotate("text", x = 0.825, y = max(RPNPF_removed_curve$n_removed) * 0.5,
           label = "Removed from Analysis", size = 4, color = "gray20") +
  
  #Color the part you want
  geom_rect(aes(xmin = 0.20, xmax = 0.65, ymin = -Inf, ymax = Inf),
            fill = "dodgerblue", alpha = a) +
  
  # Custom x-axis ticks
  scale_x_continuous(
    name = "Missingness Threshold",
    breaks = c(0, 0.1, 0.2, 0.425, 0.65, 0.825, 1.0),
    labels = c("0", "M1", "0.2", "M2", "0.65", "M3", "1.0")
  ) +
  
  labs(
    title = "Number of Metabolites Classified by Missingness Threshold - RPNPF",
    y = "Number of Metabolites Removed"
  ) +
  theme_bw(base_size = 12)



