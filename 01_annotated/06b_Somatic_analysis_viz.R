library(tidyverse)
library(ggrepel)


#source("06a_Somatic_analysis.R")
combined_results <- readxl::read_xlsx("../res/annotated/somatic/somatic_analysis.xlsx")



combined_results_ZHP <- combined_results %>% filter(platform == "ZHP")
combined_results_RPNPF <- combined_results %>% filter(platform == "RPNPF")

# ==== ZHP  data  ==============================================================

ttest_Combined_ZHP <- combined_results_ZHP %>% 
  filter(test_type == "TTest")

fish_Combined_ZHP <- combined_results_ZHP %>%
  filter(test_type == "Fisher")

# ==== ZHP - TTest Visualizations ==============================================

# ==== Volcano Plot of Somatic Differences ====

#"Somatic" = "#C0392B",
#"Non-somatic" = "#F1948A"

ttest_results_Combined_ZHP <- ttest_Combined_ZHP %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Somatic+",
      log2_fc < 0           ~ "Higher in Somatic-",
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
    "Higher in Somatic+" = "#C0392B",
    "Higher in Somatic-" = "#F1948A",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Somatic-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (ZHP Platform)",
    caption = paste0("n = ", nrow(ttest_results_Combined_ZHP), " metabolites tested via two-sided t-test."),
    x = "log2(Somatic+ / Somatic-)\n↑ Higher in Somatic+, ↓ Higher in Somatic-",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(Somatic+ vs. Somatic-)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_Combined_ZHP

# Save volcano plot
ggsave("../res/annotated/somatic/ZHP_Combined_ttest_volcano.png",
       p_sex_volcano_Combined_ZHP, width = 10, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_Combined_ZHP <- ttest_results_Combined_ZHP %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_ZHP_Combined <- ggplot(
  ZHP_imp_long %>% filter(compound_name %in% top5_Combined_ZHP),
  aes(x = somatic_present, y = intensity_log2, fill = somatic_present)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("Somatic_Neg" = "#F1948A", "Somatic_Pos" = "#C0392B")) +
  labs(
    title = "Top 5 Somatic-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (ZHP Platform)",
    x = "Somatic Status",
    y = "log2(Intensity)",
    fill = "Somatic Status"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_ZHP_Combined

# Save violin plot
ggsave(
  filename = "../res/annotated/somatic/ZHP_Combined_top5_somatic_ttest.png",
  plot = p_violin_ZHP_Combined,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish JMML  ========================================


# Add enrichment direction labels
fish_plot_ZHP_somatic_Combined <- fish_Combined_ZHP %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Somatic+",
      log2_OR <= -1              ~ "Detected more in Somatic-",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))





#"Somatic" = "#C0392B",
#"Non-somatic" = "#F1948A"

fish_p_ZHP_somatic_Combined_JMML_Control <- ggplot(fish_plot_ZHP_somatic_Combined,
                                               aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_ZHP_somatic_Combined %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Somatic+" = "#C0392B",
    "Detected more in Somatic-" = "#F1948A",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Somatic-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (ZHP Platform)",
    caption = paste0("n = ", nrow(fish_plot_ZHP_somatic_Combined), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Somatic+ / Somatic-)\n↑ More detected in Somatic+, ↓ More detected in Somatic-",
    y = "-log10(p-value)",
    color = "Relative Detection\n(Somatic+ vs. Somatic-)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_ZHP_somatic_Combined_JMML_Control

# top_fish_Control_ZHP  <- fish_plot_ZHP_somatic_Control_JMML_Control %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_Control_ZHP)



ggsave(filename = "../res/annotated/somatic/ZHP_Combined_fisher_volcano_somatic.png",
       plot = fish_p_ZHP_somatic_Combined_JMML_Control, width = 10)


# ==== RPNPF  data  ==============================================================

ttest_Combined_RPNPF <- combined_results_RPNPF %>% 
  filter(test_type == "TTest")

fish_Combined_RPNPF <- combined_results_RPNPF %>%
  filter(test_type == "Fisher")

# ==== RPNPF - TTest Visualizations ==============================================

# ==== Volcano Plot of Somatic Differences ====

#"Somatic" = "#C0392B",
#"Non-somatic" = "#F1948A"

ttest_results_Combined_RPNPF <- ttest_Combined_RPNPF %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in Somatic+",
      log2_fc < 0           ~ "Higher in Somatic-",
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
    "Higher in Somatic+" = "#C0392B",
    "Higher in Somatic-" = "#F1948A",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Somatic-based Metabolite Differences",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform)",
    caption = paste0("n = ", nrow(ttest_results_Combined_RPNPF), " metabolites tested via two-sided t-test."),
    x = "log2(Somatic+ / Somatic-)\n↑ Higher in Somatic+, ↓ Higher in Somatic-",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity\n(Somatic+ vs. Somatic-)",
    shape = "Statistically Significant\n(p < 0.05)"
  ) +
  theme_bw(base_size = 14)


p_sex_volcano_Combined_RPNPF

# Save volcano plot
ggsave("../res/annotated/somatic/RPNPF_Combined_ttest_volcano.png",
       p_sex_volcano_Combined_RPNPF, width = 10, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====

# Identify top 5 metabolites by raw p-value
top5_Combined_RPNPF <- ttest_results_Combined_RPNPF %>%
  arrange(p_value) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_RPNPF_Combined <- ggplot(
  RPNPF_imp_long %>% filter(compound_name %in% top5_Combined_RPNPF),
  aes(x = somatic_present, y = intensity_log2, fill = somatic_present)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("Somatic_Neg" = "#F1948A", "Somatic_Pos" = "#C0392B")) +
  labs(
    title = "Top 5 Somatic-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform)",
    x = "Somatic Status",
    y = "log2(Intensity)",
    fill = "Somatic Status"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )
p_violin_RPNPF_Combined

# Save violin plot
ggsave(
  filename = "../res/annotated/somatic/RPNPF_Combined_top5_somatic_ttest.png",
  plot = p_violin_RPNPF_Combined,
  width = 10, height = 6, units = "in", dpi = 300
)

# ==== Viz Fish Combined JMML + Control ========================================


# Add enrichment direction labels
fish_plot_RPNPF_somatic_Combined <- fish_Combined_RPNPF %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in Somatic+",
      log2_OR <= -1              ~ "Detected more in Somatic-",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")),
         neglog10_p = -log10(raw_p))





#"Somatic" = "#C0392B",
#"Non-somatic" = "#F1948A"

fish_p_RPNPF_somatic_Combined_JMML_Control <- ggplot(fish_plot_RPNPF_somatic_Combined,
                                               aes(x = log2_OR, y = neglog10_p)) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +       # p = 0.05 threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(data = fish_plot_RPNPF_somatic_Combined %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in Somatic+" = "#C0392B",
    "Detected more in Somatic-" = "#F1948A",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "Somatic-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (RPNPF Platform)",
    caption = paste0("n = ", nrow(fish_plot_RPNPF_somatic_Combined), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: Somatic+ / Somatic-)\n↑ More detected in Somatic+, ↓ More detected in Somatic-",
    y = "-log10(p-value)",
    color = "Relative Detection\n(Somatic+ vs. Somatic-)",
    shape = "Statistically Significant\n(p < 0.05)"
  )+
  theme_bw(base_size = 14)

fish_p_RPNPF_somatic_Combined_JMML_Control

# top_fish_Control_RPNPF  <- fish_plot_RPNPF_somatic_Control_JMML_Control %>%
#   arrange(desc(neglog10_p)) %>%
#   slice_head(n = 1) %>%
#   transmute(Compound = compound_name,
#             `Missing in Female Samples` = a_F_NA,
#             `Detected in Female Samples` = b_F_Present,
#             `Missing in Male Samples` = c_M_NA,
#             `Detected in Male Samples` = d_M_Present)
# 
# print(top_fish_Control_RPNPF)



ggsave(filename = "../res/annotated/somatic/RPNPF_Combined_fisher_volcano_somatic.png",
       plot = fish_p_RPNPF_somatic_Combined_JMML_Control, width = 10)
