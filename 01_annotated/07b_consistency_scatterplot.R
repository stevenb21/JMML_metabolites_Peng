#source("05b_Sex_analysis_viz.R")

# ==== ZHP Visualization  ======================================================

# ==== TTest  ==================================================================
# ==== Scatterplot: JMML vs Control Log Fold Change ====
# Create a scatterplot of sex log2 fold changes comparing JMML-only to Control-only

# Join JMML-only and Control-only t-test results
ZHP_ttest_scatter_df <- ttest_results_sex_JMML_ZHP %>%  # defined in Sex_analysis_viz.R
  select(compound_name, log2_fc, raw_p) %>%
  rename(log2_fc_JMML = log2_fc) %>%
  inner_join(
    ttest_results_sex_control_ZHP %>%
      select(compound_name, log2_fc) %>%
      rename(log2_fc_Control = log2_fc),
    by = "compound_name"
  ) %>%
  # Determine significance categories
  mutate(
    sig_JMML    = raw_p < 0.05,
    sig_Control = ttest_results_sex_control_ZHP$raw_p < 0.05,
    sig = case_when(
      sig_JMML & sig_Control  ~ "Significant in both",
      sig_JMML                ~ "Significant in JMML-only",
      sig_Control             ~ "Significant in Control-only",
      TRUE                    ~ "Not Significant"
    )
  )

# Define shared color scale
sig_colors <- c(
  "Significant in both"     = "purple",  
  "Significant in JMML-only" = "#E74C3C", 
  "Significant in Control-only" = "#3498DB", 
  "Not Significant"         = "gray70"
)


# Plot
p_scatter_ttest_ZHP <- ggplot(ZHP_ttest_scatter_df, aes(x = log2_fc_JMML, y = log2_fc_Control, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(
    data = subset(ZHP_ttest_scatter_df, sig != "Not Significant"),
    aes(label = compound_name), 
    size = 3,
    max.overlaps = 15
  ) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Sex Log2 Fold Change:\nJMML-only vs Control-only",
    subtitle = "TTest using raw-p values (ZHP)",
    x = "log2(Female / Male) in JMML-only",
    y = "log2(Female / Male) in Control-only",
    color = "Significance"
  ) +
  theme_bw(base_size = 14)

p_scatter_ttest_ZHP

# Save plot
ggsave(
  filename = "../res/annotated/sex/ZHP_ttest_scatter_logFC_JMML_vs_Control.png",
  plot    = p_scatter_ttest_ZHP,
  width   = 7,
  height  = 5,
  dpi     = 300
)

# ==== Fisher  =================================================================
# ==== Scatterplot: JMML vs Control Log Fold Change ====
# Create a scatterplot of sex log2 fold changes comparing JMML-only to Control-only

ZHP_fish_scatter_df <- fish_plot_ZHP_sex_JMML %>%
  select(compound_name, log2_OR, raw_p) %>%
  rename(log2_OR_JMML = log2_OR,
         p_JMML       = raw_p) %>%
  inner_join(
    fish_plot_ZHP_sex_Control %>%
      select(compound_name, log2_OR, raw_p) %>%
      rename(log2_OR_Control = log2_OR,
             p_Control       = raw_p),
    by = "compound_name"
  ) %>%
  mutate(
    sig_JMML    = p_JMML   < 0.05,
    sig_Control = p_Control < 0.05,
    sig = case_when(
      sig_JMML & sig_Control   ~ "Significant in both",
      sig_JMML                 ~ "Significant in JMML-only",
      sig_Control              ~ "Significant in Control-only",
      TRUE                     ~ "Not Significant"
    )
  )


# Plot
p_scatter_fish_ZHP <- ggplot(ZHP_fish_scatter_df, aes(x = log2_OR_JMML, y = log2_OR_Control, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(
    data = subset(ZHP_fish_scatter_df, sig != "Not Significant"),
    aes(label = compound_name), 
    size = 3,
    max.overlaps = 15
  ) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Sex Log2 OR:\nJMML-only vs Control-only",
    subtitle = "Fisher's exact using raw-p values (ZHP)",
    x = "log2(OR: Female / Male) in JMML-only",
    y = "log2(OR Female / Male) in Control-only",
    color = "Significance"
  ) +
  theme_bw(base_size = 14)

p_scatter_fish_ZHP

# Save plot
ggsave(
  filename = "../res/annotated/sex/ZHP_fish_scatter_logOR_JMML_vs_Control.png",
  plot    = p_scatter_fish_ZHP,
  width   = 7,
  height  = 5,
  dpi     = 300
)



# ==== RPNPF Visualization  ====================================================

# ==== TTest  ==================================================================
# ==== Scatterplot: JMML vs Control Log Fold Change ====
# Create a scatterplot of sex log2 fold changes comparing JMML-only to Control-only

# Join JMML-only and Control-only t-test results
RPNPF_ttest_scatter_df <- ttest_results_sex_JMML_RPNPF %>%  # defined in Sex_analysis_viz.R
  select(compound_name, log2_fc, raw_p) %>%
  rename(log2_fc_JMML = log2_fc) %>%
  inner_join(
    ttest_results_sex_control_RPNPF %>%
      select(compound_name, log2_fc) %>%
      rename(log2_fc_Control = log2_fc),
    by = "compound_name"
  ) %>%
  # Determine significance categories
  mutate(
    sig_JMML    = raw_p < 0.05,
    sig_Control = ttest_results_sex_control_RPNPF$raw_p < 0.05,
    sig = case_when(
      sig_JMML & sig_Control  ~ "Significant in both",
      sig_JMML                ~ "Significant in JMML-only",
      sig_Control             ~ "Significant in Control-only",
      TRUE                    ~ "Not Significant"
    )
  )

# Define shared color scale
sig_colors <- c(
  "Significant in both"     = "purple",  
  "Significant in JMML-only" = "#E74C3C", 
  "Significant in Control-only" = "#3498DB", 
  "Not Significant"         = "gray70"
)


# Plot
p_scatter_ttest_RPNPF <- ggplot(RPNPF_ttest_scatter_df, aes(x = log2_fc_JMML, y = log2_fc_Control, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(
    data = subset(RPNPF_ttest_scatter_df, sig != "Not Significant"),
    aes(label = compound_name), 
    size = 3,
    max.overlaps = 15
  ) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Sex Log2 Fold Change:\nJMML-only vs Control-only",
    subtitle = "TTest using raw-p values (RPNPF)",
    x = "log2(Female / Male) in JMML-only",
    y = "log2(Female / Male) in Control-only",
    color = "Significance"
  ) +
  theme_bw(base_size = 14)

p_scatter_ttest_RPNPF

# Save plot
ggsave(
  filename = "../res/annotated/sex/RPNPF_ttest_scatter_logFC_JMML_vs_Control.png",
  plot    = p_scatter_ttest_RPNPF,
  width   = 7,
  height  = 5,
  dpi     = 300
)

# ==== Fisher  =================================================================
# ==== Scatterplot: JMML vs Control Log Fold Change ====
# Create a scatterplot of sex log2 fold changes comparing JMML-only to Control-only

RPNPF_fish_scatter_df <- fish_plot_RPNPF_sex_JMML %>%
  select(compound_name, log2_OR, raw_p) %>%
  rename(log2_OR_JMML = log2_OR,
         p_JMML       = raw_p) %>%
  inner_join(
    fish_plot_RPNPF_sex_Control %>%
      select(compound_name, log2_OR, raw_p) %>%
      rename(log2_OR_Control = log2_OR,
             p_Control       = raw_p),
    by = "compound_name"
  ) %>%
  mutate(
    sig_JMML    = p_JMML   < 0.05,
    sig_Control = p_Control < 0.05,
    sig = case_when(
      sig_JMML & sig_Control   ~ "Significant in both",
      sig_JMML                 ~ "Significant in JMML-only",
      sig_Control              ~ "Significant in Control-only",
      TRUE                     ~ "Not Significant"
    )
  )


# Plot
p_scatter_fish_RPNPF <- ggplot(RPNPF_fish_scatter_df, aes(x = log2_OR_JMML, y = log2_OR_Control, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(
    data = subset(RPNPF_fish_scatter_df, sig != "Not Significant"),
    aes(label = compound_name), 
    size = 3,
    max.overlaps = 15
  ) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Sex Log2 OR:\nJMML-only vs Control-only",
    subtitle = "Fisher's exact using raw-p values (RPNPF)",
    x = "log2(OR: Female / Male) in JMML-only",
    y = "log2(OR Female / Male) in Control-only",
    color = "Significance"
  ) +
  theme_bw(base_size = 14)

p_scatter_fish_RPNPF

# Save plot
ggsave(
  filename = "../res/annotated/sex/RPNPF_fish_scatter_logOR_JMML_vs_Control.png",
  plot    = p_scatter_fish_RPNPF,
  width   = 7,
  height  = 5,
  dpi     = 300
)

