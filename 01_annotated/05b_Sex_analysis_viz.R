library(tidyverse)
library(ggrepel)


#source("05a_Sex_analysis.R")
#Make sure you run the above line if you are debugging!

combined_results <- readxl::read_xlsx("../res/annotated/sex/sex_analysis.xlsx")

# -------------------------------------------------------------------------


#Flip effect size of fisher to match the T-Test effect direction. 
# set it to: (log(OR) female positive)
combined_results <- combined_results %>%
  mutate(
    log2_OR = case_when(
      test_type == "Fisher" ~ -log2_OR,
      TRUE                  ~ log2_OR
    )
  )


combined_results_ZHP <- combined_results %>% filter(platform == "ZHP")
combined_results_RPNPF <- combined_results %>% filter(platform == "RPNPF")

# ==== ZHP  data  ==============================================================

# ==== ZHP - Combined data (JMML + Control) ====================================
ttest_Combined_ZHP <- combined_results_ZHP %>% 
  filter(test_type == "TTest", cohort == "Combined_JMML_Control")

fish_Combined_ZHP <- combined_results_ZHP %>%
  filter(cohort == "Combined_JMML_Control") %>% 
  filter(test_type == "Fisher")

# ==== ZHP - TTest Visualizations ==============================================

# ==== Volcano Plot of Sex Differences ====

ttest_results_Combined_ZHP <- ttest_Combined_ZHP %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Females",
      log2_fc < 0           ~ "Higher in Males",
      TRUE                    ~ NA_character_
    ),
    sig = ifelse(p_value < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

p_sex_volcano_Combined_ZHP <- ggplot(
  ttest_results_Combined_ZHP,
  aes(x = log2_fc, y = -log10(p_value))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_Combined_ZHP %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5, max.overlaps = 15, box.padding = 0.3,
  ) +
  scale_color_manual(values = c(
    "Higher in Females" = "#E67E22",
    "Higher in Males" = "#1ABC9C",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (ZHP Platform)\nCombined JMML + Control",
    caption = paste0("n = ", nrow(ttest_results_Combined_ZHP), " metabolites tested via two-sided t-test."),
    x = "log2(Female / Male)\n↑ Higher in females, ↓ Higher in males",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_Combined_ZHP

# Save volcano plot
ggsave("../res/annotated/sex/ZHP_Combined_ttest_volcano.png",
       p_sex_volcano_Combined_ZHP, width = 10, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_Combined_ZHP <- ttest_results_Combined_ZHP %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_ZHP_Combined <- ggplot(
  ZHP_imp_long %>% filter(compound_name %in% top5_Combined_ZHP),
  aes(x = sex, y = intensity_log2, fill = sex)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#E67E22", "M" = "#1ABC9C")) +
  labs(
    title = "Top 5 Sex-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (ZHP Platform)\nCombined JMML + Control",
    x = "Biological Sex",
    y = "log2(Intensity)",
    fill = "Sex"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_ZHP_Combined

# Save violin plot
ggsave(
  filename = "../res/annotated/sex/ZHP_Combined_top5_sex_ttest.png",
  plot = p_violin_ZHP_Combined,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish Combined JMML + Control ========================================


# Add enrichment direction labels
fish_plot_ZHP_sex_Combined <- fish_Combined_ZHP %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Females",
      log2_OR <= -1              ~ "Detected more in Males",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))




#Males = "#1ABC9C"
#Females = "#E67E22" 

fish_p_ZHP_sex_Combined_JMML_Control <- ggplot(fish_plot_ZHP_sex_Combined,
                                               aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_ZHP_sex_Combined %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Males" = "#1ABC9C",
    "Detected more in Females" = "#E67E22",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (ZHP Platform) | Combined JMML + Control",
    caption = paste0("n = ", nrow(fish_plot_ZHP_sex_Combined), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Female / Male)\n↑ More detected in females, ↓ More detected in males",
    y = "-log10(p-value)",
    color = "Relative Detection\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_ZHP_sex_Combined_JMML_Control

# top_fish_Control_ZHP  <- fish_plot_ZHP_sex_Control_JMML_Control %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_Control_ZHP)



ggsave(filename = "../res/annotated/sex/ZHP_Combined_fisher_volcano_sex.png",
       plot = fish_p_ZHP_sex_Combined_JMML_Control, width = 10)


# ==== ZHP - Control ===========================================================

ttest_Control_ZHP <- combined_results_ZHP %>% 
  filter(cohort == "Control_only") %>% 
  filter(test_type == "TTest")

fish_Control_ZHP <- combined_results_ZHP %>%
  filter(cohort == "Control_only") %>% 
  filter(test_type == "Fisher")

# ==== Viz Fish Control ========================================================

# ==== ZHP - TTest Visualizations ==============================================

# ==== Volcano Plot of Sex Differences ====

ttest_results_Control_ZHP <- ttest_Control_ZHP %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Females",
      log2_fc < 0           ~ "Higher more in Males",
      TRUE                    ~ NA_character_
    ),
    sig = ifelse(p_value < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

p_sex_volcano_Control_ZHP <- ggplot(
  ttest_results_Control_ZHP,
  aes(x = log2_fc, y = -log10(p_value))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_Control_ZHP %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5, max.overlaps = 15, box.padding = 0.3,
  ) +
  scale_color_manual(values = c(
    "Higher in Females" = "#E67E22",
    "Higher more in Males" = "#1ABC9C",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (ZHP Platform) | Control",
    caption = paste0("n = ", nrow(ttest_results_Control_ZHP), " metabolites tested via two-sided t-test."),
    x = "log2(Female / Male)\n↑ Higher in females, ↓ Higher in males",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_Control_ZHP

# Save volcano plot
ggsave("../res/annotated/sex/ZHP_Control_ttest_volcano.png",
       p_sex_volcano_Control_ZHP, width = 8, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_Control_ZHP <- ttest_results_Control_ZHP %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_ZHP_Control <- ggplot(
  ZHP_imp_long %>% filter(compound_name %in% top5_Control_ZHP),
  aes(x = sex, y = intensity_log2, fill = sex)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#E67E22", "M" = "#1ABC9C")) +
  labs(
    title = "Top 5 Sex-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (ZHP Platform) | Control",
    x = "Biological Sex",
    y = "log2(Intensity)",
    fill = "Sex"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_ZHP_Control

# Save violin plot
ggsave(
  filename = "../res/annotated/sex/ZHP_Control_top5_sex_ttest.png",
  plot = p_violin_ZHP_Control,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish Control JMML + Control ========================================


# Add enrichment direction labels
fish_plot_ZHP_sex_Control <- fish_Control_ZHP %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Females",
      log2_OR <= -1              ~ "Detected more in Males",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))




#Males = "#1ABC9C"
#Females = "#E67E22" 

fish_p_ZHP_sex_Control <- ggplot(fish_plot_ZHP_sex_Control,
                                               aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_ZHP_sex_Control %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Males" = "#1ABC9C",
    "Detected more in Females" = "#E67E22",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (ZHP Platform) | Control",
    caption = paste0("n = ", nrow(fish_plot_ZHP_sex_Control), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Female / Male)\n↑ More detected in females, ↓ More detected in males",
    y = "-log10(p-value)",
    color = "Relative Detection\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_ZHP_sex_Control

# top_fish_Control_ZHP  <- fish_plot_ZHP_sex_Control_JMML_Control %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_Control_ZHP)



ggsave(filename = "../res/annotated/sex/ZHP_Control_fisher_volcano_sex.png",
       plot = fish_p_ZHP_sex_Control, width = 10)

# ==== ZHP - JMML ===========================================================

ttest_JMML_ZHP <- combined_results_ZHP %>% 
  filter(cohort == "JMML_only") %>% 
  filter(test_type == "TTest")

fish_JMML_ZHP <- combined_results_ZHP %>%
  filter(cohort == "JMML_only") %>% 
  filter(test_type == "Fisher")


# ==== ZHP - TTest Visualizations ==============================================

# ==== Volcano Plot of Sex Differences ====

ttest_results_JMML_ZHP <- ttest_JMML_ZHP %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Females",
      log2_fc < 0           ~ "Higher in Males",
      TRUE                    ~ NA_character_
    ),
    sig = ifelse(p_value < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

p_sex_volcano_JMML_ZHP <- ggplot(
  ttest_results_JMML_ZHP,
  aes(x = log2_fc, y = -log10(p_value))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_JMML_ZHP %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5, max.overlaps = 15, box.padding = 0.3,
  ) +
  scale_color_manual(values = c(
    "Higher in Females" = "#E67E22",
    "Higher in Males" = "#1ABC9C",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (ZHP Platform) | JMML",
    caption = paste0("n = ", nrow(ttest_results_JMML_ZHP), " metabolites tested via two-sided t-test."),
    x = "log2(Female / Male)\n↑ Higher in females, ↓ Higher in males",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_JMML_ZHP

# Save volcano plot
ggsave("../res/annotated/sex/ZHP_JMML_ttest_volcano.png",
       p_sex_volcano_JMML_ZHP, width = 8, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_JMML_ZHP <- ttest_results_JMML_ZHP %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_ZHP_JMML <- ggplot(
  ZHP_imp_long %>% filter(compound_name %in% top5_JMML_ZHP),
  aes(x = sex, y = intensity_log2, fill = sex)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#E67E22", "M" = "#1ABC9C")) +
  labs(
    title = "Top 5 Sex-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (ZHP Platform) | JMML",
    x = "Biological Sex",
    y = "log2(Intensity)",
    fill = "Sex"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_ZHP_JMML

# Save violin plot
ggsave(
  filename = "../res/annotated/sex/ZHP_JMML_top5_sex_ttest.png",
  plot = p_violin_ZHP_JMML,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish JMML  ==========================================================


# Add enrichment direction labels
fish_plot_ZHP_sex_JMML <- fish_JMML_ZHP %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Females",
      log2_OR <= -1              ~ "Detected more in Males",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))




#Males = "#1ABC9C"
#Females = "#E67E22" 

fish_p_ZHP_sex_JMML <- ggplot(fish_plot_ZHP_sex_JMML,
                                 aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_ZHP_sex_JMML %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Males" = "#1ABC9C",
    "Detected more in Females" = "#E67E22",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (ZHP Platform) | JMML",
    caption = paste0("n = ", nrow(fish_plot_ZHP_sex_JMML), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Female / Male)\n↑ More detected in females, ↓ More detected in males",
    y = "-log10(p-value)",
    color = "Relative Detection\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_ZHP_sex_JMML

# top_fish_JMML_ZHP  <- fish_plot_ZHP_sex_JMML_JMML_JMML %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_JMML_ZHP)



ggsave(filename = "../res/annotated/sex/ZHP_JMML_fisher_volcano_sex.png",
       plot = fish_p_ZHP_sex_JMML, width = 10)

# ==== RPNPF  data  ==============================================================

# ==== RPNPF - Combined data (JMML + Control) ====================================
ttest_Combined_RPNPF <- combined_results_RPNPF %>% 
  filter(test_type == "TTest", cohort == "Combined_JMML_Control")

fish_Combined_RPNPF <- combined_results_RPNPF |> 
  filter(cohort == "Combined_JMML_Control") %>% 
  filter(test_type == "Fisher")

# ==== RPNPF - TTest Visualizations ==============================================

# ==== Volcano Plot of Sex Differences ====

ttest_results_Combined_RPNPF <- ttest_Combined_RPNPF %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Females",
      log2_fc < 0           ~ "Higher in Males",
      TRUE                    ~ NA_character_
    ),
    sig = ifelse(p_value < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

p_sex_volcano_Combined_RPNPF <- ggplot(
  ttest_results_Combined_RPNPF,
  aes(x = log2_fc, y = -log10(p_value))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_Combined_RPNPF %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5, max.overlaps = 15, box.padding = 0.3,
  ) +
  scale_color_manual(values = c(
    "Higher in Females" = "#E67E22",
    "Higher in Males" = "#1ABC9C",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform)\nCombined JMML + Control",
    caption = paste0("n = ", nrow(ttest_results_Combined_RPNPF), " metabolites tested via two-sided t-test."),
    x = "log2(Female / Male)\n↑ Higher in females, ↓ Higher in males",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_Combined_RPNPF

# Save volcano plot
ggsave("../res/annotated/sex/RPNPF_Combined_ttest_volcano.png",
       p_sex_volcano_Combined_RPNPF, width = 10, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_Combined_RPNPF <- ttest_results_Combined_RPNPF %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_RPNPF_Combined <- ggplot(
  RPNPF_imp_long %>% filter(compound_name %in% top5_Combined_RPNPF),
  aes(x = sex, y = intensity_log2, fill = sex)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#E67E22", "M" = "#1ABC9C")) +
  labs(
    title = "Top 5 Sex-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform)\nCombined JMML + Control",
    x = "Biological Sex",
    y = "log2(Intensity)",
    fill = "Sex"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_RPNPF_Combined

# Save violin plot
ggsave(
  filename = "../res/annotated/sex/RPNPF_Combined_top5_sex_ttest.png",
  plot = p_violin_RPNPF_Combined,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish Combined JMML + Control ========================================


# Add enrichment direction labels
fish_plot_RPNPF_sex_Combined <- fish_Combined_RPNPF %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Females",
      log2_OR <= -1              ~ "Detected more in Males",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))




#Males = "#1ABC9C"
#Females = "#E67E22" 

fish_p_RPNPF_sex_Combined_JMML_Control <- ggplot(fish_plot_RPNPF_sex_Combined,
                                               aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_RPNPF_sex_Combined %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Males" = "#1ABC9C",
    "Detected more in Females" = "#E67E22",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (RPNPF Platform) | Combined JMML + Control",
    caption = paste0("n = ", nrow(fish_plot_RPNPF_sex_Combined), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Female / Male)\n↑ More detected in females, ↓ More detected in males",
    y = "-log10(p-value)",
    color = "Relative Detection\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_RPNPF_sex_Combined_JMML_Control

# top_fish_Control_RPNPF  <- fish_plot_RPNPF_sex_Control_JMML_Control %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_Control_RPNPF)



ggsave(filename = "../res/annotated/sex/RPNPF_Combined_fisher_volcano_sex.png",
       plot = fish_p_RPNPF_sex_Combined_JMML_Control, width = 10)


# ==== RPNPF - Control ===========================================================

ttest_Control_RPNPF <- combined_results_RPNPF %>% 
  filter(cohort == "Control_only") %>% 
  filter(test_type == "TTest")

fish_Control_RPNPF <- combined_results_RPNPF %>%
  filter(cohort == "Control_only") %>% 
  filter(test_type == "Fisher")

# ==== Viz Fish Control ========================================================

# ==== RPNPF - TTest Visualizations ==============================================

# ==== Volcano Plot of Sex Differences ====

ttest_results_Control_RPNPF <- ttest_Control_RPNPF %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Females",
      log2_fc < 0           ~ "Higher more in Males",
      TRUE                    ~ NA_character_
    ),
    sig = ifelse(p_value < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

p_sex_volcano_Control_RPNPF <- ggplot(
  ttest_results_Control_RPNPF,
  aes(x = log2_fc, y = -log10(p_value))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_Control_RPNPF %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5, max.overlaps = 15, box.padding = 0.3,
  ) +
  scale_color_manual(values = c(
    "Higher in Females" = "#E67E22",
    "Higher more in Males" = "#1ABC9C",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform) | Control",
    caption = paste0("n = ", nrow(ttest_results_Control_RPNPF), " metabolites tested via two-sided t-test."),
    x = "log2(Female / Male)\n↑ Higher in females, ↓ Higher in males",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_Control_RPNPF

# Save volcano plot
ggsave("../res/annotated/sex/RPNPF_Control_ttest_volcano.png",
       p_sex_volcano_Control_RPNPF, width = 8, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_Control_RPNPF <- ttest_results_Control_RPNPF %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_RPNPF_Control <- ggplot(
  RPNPF_imp_long %>% filter(compound_name %in% top5_Control_RPNPF),
  aes(x = sex, y = intensity_log2, fill = sex)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#E67E22", "M" = "#1ABC9C")) +
  labs(
    title = "Top 5 Sex-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform) | Control",
    x = "Biological Sex",
    y = "log2(Intensity)",
    fill = "Sex"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_RPNPF_Control

# Save violin plot
ggsave(
  filename = "../res/annotated/sex/RPNPF_Control_top5_sex_ttest.png",
  plot = p_violin_RPNPF_Control,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish Control JMML + Control ========================================


# Add enrichment direction labels
fish_plot_RPNPF_sex_Control <- fish_Control_RPNPF %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Females",
      log2_OR <= -1              ~ "Detected more in Males",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))




#Males = "#1ABC9C"
#Females = "#E67E22" 

fish_p_RPNPF_sex_Control <- ggplot(fish_plot_RPNPF_sex_Control,
                                 aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_RPNPF_sex_Control %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Males" = "#1ABC9C",
    "Detected more in Females" = "#E67E22",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (RPNPF Platform) | Control",
    caption = paste0("n = ", nrow(fish_plot_RPNPF_sex_Control), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Female / Male)\n↑ More detected in females, ↓ More detected in males",
    y = "-log10(p-value)",
    color = "Relative Detection\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_RPNPF_sex_Control

# top_fish_Control_RPNPF  <- fish_plot_RPNPF_sex_Control_JMML_Control %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_Control_RPNPF)



ggsave(filename = "../res/annotated/sex/RPNPF_Control_fisher_volcano_sex.png",
       plot = fish_p_RPNPF_sex_Control, width = 10)

# ==== RPNPF - JMML ===========================================================

ttest_JMML_RPNPF <- combined_results_RPNPF %>% 
  filter(cohort == "JMML_only") %>% 
  filter(test_type == "TTest")

fish_JMML_RPNPF <- combined_results_RPNPF %>%
  filter(cohort == "JMML_only") %>% 
  filter(test_type == "Fisher")


# ==== RPNPF - TTest Visualizations ==============================================

# ==== Volcano Plot of Sex Differences ====

ttest_results_JMML_RPNPF <- ttest_JMML_RPNPF %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Females",
      log2_fc < 0           ~ "Higher in Males",
      TRUE                    ~ NA_character_
    ),
    sig = ifelse(p_value < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

p_sex_volcano_JMML_RPNPF <- ggplot(
  ttest_results_JMML_RPNPF,
  aes(x = log2_fc, y = -log10(p_value))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_JMML_RPNPF %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5, max.overlaps = 15, box.padding = 0.3,
  ) +
  scale_color_manual(values = c(
    "Higher in Females" = "#E67E22",
    "Higher in Males" = "#1ABC9C",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform) | JMML",
    caption = paste0("n = ", nrow(ttest_results_JMML_RPNPF), " metabolites tested via two-sided t-test."),
    x = "log2(Female / Male)\n↑ Higher in females, ↓ Higher in males",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_JMML_RPNPF

# Save volcano plot
ggsave("../res/annotated/sex/RPNPF_JMML_ttest_volcano.png",
       p_sex_volcano_JMML_RPNPF, width = 8, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_JMML_RPNPF <- ttest_results_JMML_RPNPF %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_RPNPF_JMML <- ggplot(
  RPNPF_imp_long %>% filter(compound_name %in% top5_JMML_RPNPF),
  aes(x = sex, y = intensity_log2, fill = sex)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#E67E22", "M" = "#1ABC9C")) +
  labs(
    title = "Top 5 Sex-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform) | JMML",
    x = "Biological Sex",
    y = "log2(Intensity)",
    fill = "Sex"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_RPNPF_JMML

# Save violin plot
ggsave(
  filename = "../res/annotated/sex/RPNPF_JMML_top5_sex_ttest.png",
  plot = p_violin_RPNPF_JMML,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish JMML  ==========================================================


# Add enrichment direction labels
fish_plot_RPNPF_sex_JMML <- fish_JMML_RPNPF %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Females",
      log2_OR <= -1              ~ "Detected more in Males",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))




#Males = "#1ABC9C"
#Females = "#E67E22" 

fish_p_RPNPF_sex_JMML <- ggplot(fish_plot_RPNPF_sex_JMML,
                              aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_RPNPF_sex_JMML %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Males" = "#1ABC9C",
    "Detected more in Females" = "#E67E22",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Sex-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (RPNPF Platform) | JMML",
    caption = paste0("n = ", nrow(fish_plot_RPNPF_sex_JMML), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Female / Male)\n↑ More detected in females, ↓ More detected in males",
    y = "-log10(p-value)",
    color = "Relative Detection\n(F vs. M)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_RPNPF_sex_JMML

# top_fish_JMML_RPNPF  <- fish_plot_RPNPF_sex_JMML_JMML_JMML %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_JMML_RPNPF)



ggsave(filename = "../res/annotated/sex/RPNPF_JMML_fisher_volcano_sex.png",
       plot = fish_p_RPNPF_sex_JMML, width = 10)

