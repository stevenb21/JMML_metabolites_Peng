library(tidyverse)
library(readxl)
library(ggrepel)



# ==== Load data ===============================================================

#source("03_QLIC_imputation.R")
#assumes the objects in 03_QLIC are in the environment

final_results <- read_xlsx("../res/annotated/final_results_global_fdr.xlsx")


# ==== ZHP Visualizations ======================================================

# ==== FDR Visualizations ======================================================

# ==== Viz Fish FDR ============================================================

fish_ZHP <- final_results %>%
  filter(test_type == "Fisher") %>%
  filter(platform == "ZHP")

# Add enrichment direction labels
fish_plot_ZHP_fdr <- fish_ZHP %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in JMML",
      log2_OR <= -1              ~ "Detected more in Control",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(global_fdr < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))


fish_p_ZHP_fdr <- ggplot(fish_plot_ZHP_fdr, aes(x = log2_OR, y = neglog10_fdr)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray") +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7) +
  geom_text_repel(data = fish_plot_ZHP_fdr %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in JMML" = "#E74C3C",
    "Detected more in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (ZHP Platform)",
    caption = paste0("n = ", nrow(fish_plot_ZHP_fdr), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: JMML vs. Control)\n↑ More detected in JMML, ↓ More detected in Control",
    y = "-log10(FDR)",
    color = "Relative Detection",
    shape = "Statistically Significant\n(FDR < 0.05)"
  )+
  theme_bw(base_size = 14)

print(fish_p_ZHP_fdr)

ggsave(filename = "../res/annotated/fisher_volcano_ZHP_FDR.png", plot = fish_p_ZHP_fdr)



# ==== ZHP - TTest Visualizations FDR ==========================================

ttest_ZHP <- final_results %>%
  filter(test_type == "TTest") %>%
  filter(platform == "ZHP")

# ==== Volcano Plot of Differences ====

ttest_results_Combined_ZHP_fdr <- ttest_ZHP %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in JMML",
      log2_fc < 0           ~ "Higher in Control",
      TRUE                    ~ "Ambiguous"
    ),
    sig = ifelse(global_fdr < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

ZHP_ttest_volcano_FDR <- ggplot(
  ttest_results_Combined_ZHP_fdr,
  aes(x = log2_fc, y = -log10(global_fdr))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_Combined_ZHP_fdr %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5,
    max.overlaps = 15,
    box.padding = 0.3
  ) +
  scale_color_manual(values = c(
    "Higher in JMML" = "#E74C3C",
    "Higher in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "T-test on Imputed Box-Cox Transformed Intensities\n(ZHP Platform)",
    caption = paste0("n = ", nrow(ttest_results_Combined_ZHP_fdr), " metabolites tested via two-sided t-test."),
    x = "log2FC(JMML / Control)\n↑ Higher in JMML, ↓ Higher in Control",
    y = "-log10(FDR)\nfrom t-test",
    color = "Relative Intensity",
    shape = "Statistically Significant\n(FDR < 0.05)"
  ) +
  theme_bw(base_size = 14)


print(ZHP_ttest_volcano_FDR)

# Save volcano plot
ggsave("../res/annotated/ZHP_ttest_volcano_FDR.png",
       ZHP_ttest_volcano_FDR, width = 10, height = 6, dpi = 300)




# ==== Top 5 Violin Plots of Most Significant Sex Differences ====


# Identify top 5 metabolites by raw p-value
top5_Combined_ZHP <- ttest_results_Combined_ZHP_fdr %>%
  arrange(raw_p) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_ZHP <- ggplot(
  ZHP_transformed %>% filter(compound_name %in% top5_Combined_ZHP),
  aes(x = sample_group, y = intensity_log2, fill = sample_group)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("JMML" = "#E74C3C", "Control" = "#3498DB")) +
  labs(
    title = "Top 5 JMML-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (ZHP Platform)",
    x = "Case Status",
    y = "log2(Intensity)",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

p_violin_ZHP

# Save violin plot
ggsave(
  filename = "../res/annotated/ZHP_top5_ttest.png",
  plot = p_violin_ZHP,
  width = 10, height = 6, units = "in", dpi = 300
)






# ==== Raw-p value Visualizations ==============================================


# ==== Viz Fish ZHP P-value  ===================================================

fish_ZHP <- final_results %>%
  filter(test_type == "Fisher") %>%
  filter(platform == "ZHP")

# Add enrichment direction labels
fish_plot_ZHP_p <- fish_ZHP %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in JMML",
      log2_OR <= -1              ~ "Detected more in Control",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_p = -log10(raw_p),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))


# Drop columns that are entirely NA
fish_plot_ZHP_p_clean <- fish_plot_ZHP_p %>%
  select(where(~ !all(is.na(.))))  %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write.csv(fish_plot_ZHP_p_clean, file = "../res/annotated/clean_fisher_ZHP_table.csv")

fish_p_ZHP <- ggplot(fish_plot_ZHP_p,  aes(x = log2_OR, y = neglog10_p)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray") +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7) +
  # geom_text_repel(data = fish_plot_ZHP_p %>% filter(sig == "Significant"), 
  #                 aes(label = compound_name), 
  #                 size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in JMML" = "#E74C3C",
    "Detected more in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (ZHP Platform)",
    caption = paste0("n = ", nrow(fish_plot_ZHP_p ), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: JMML vs. Control)\n↑ More detected in JMML, ↓ More detected in Control",
    y = "-log10(P-value)",
    color = "Relative Detection",
    shape = "Statistically Significant\n(P-value < 0.05)"
  )+
  theme_bw(base_size = 14)


ggsave(filename = "../res/annotated/fisher_volcano_ZHP_pval.png", plot = fish_p_ZHP)




# ==== ZHP - TTest Visualizations FDR ==========================================

ttest_ZHP <- final_results %>%
  filter(test_type == "TTest") %>%
  filter(platform == "ZHP")

# ==== Volcano Plot of Sex Differences ====

ttest_results_Combined_ZHP_p <- ttest_ZHP %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in JMML",
      log2_fc < 0           ~ "Higher in Control",
      TRUE                    ~ "Ambiguous"
    ),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))


# Drop columns that are entirely NA
ttest_results_Combined_ZHP_clean <- ttest_results_Combined_ZHP_p %>%
  select(where(~ !all(is.na(.))))  %>%
  mutate(across(where(is.numeric), ~ round(., 2)))


write.csv(ttest_results_Combined_ZHP_clean, 
          file = "../res/annotated/ttest_results_Combined_ZHP_clean.csv")


ZHP_ttest_volcano_pval <- ggplot(
  ttest_results_Combined_ZHP_p ,
  aes(x = log2_fc, y = -log10(raw_p))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  scale_color_manual(values = c(
    "Higher in JMML" = "#E74C3C",
    "Higher in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "T-test on Imputed Box-Cox Transformed Intensities\n(ZHP Platform)",
    caption = paste0("n = ", nrow(ttest_results_Combined_ZHP_p), " metabolites tested via two-sided t-test."),
    x = "log2FC(JMML / Control)\n↑ Higher in JMML, ↓ Higher in Control",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity",
    shape = "Statistically Significant\n(p-value < 0.05)"
  ) +
  theme_bw(base_size = 14)


print(ZHP_ttest_volcano_pval)

# Save volcano plot
ggsave("../res/annotated/ZHP_ttest_volcano_pval.png",
       ZHP_ttest_volcano_pval, width = 10, height = 6, dpi = 300)





# ==== RPNPF Visualizations ======================================================

# ==== FDR Visualizations ======================================================

# ==== Viz Fish FDR ============================================================

fish_RPNPF <- final_results %>%
  filter(test_type == "Fisher") %>%
  filter(platform == "RPNPF")

# Add enrichment direction labels
fish_plot_RPNPF_fdr <- fish_RPNPF %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in JMML",
      log2_OR <= -1              ~ "Detected more in Control",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_fdr = -log10(global_fdr),
    sig = ifelse(global_fdr < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))


fish_p_RPNPF_fdr <- ggplot(fish_plot_RPNPF_fdr, aes(x = log2_OR, y = neglog10_fdr)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray") +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7) +
  geom_text_repel(data = fish_plot_RPNPF_fdr %>% filter(sig == "Significant"), 
                  aes(label = compound_name), 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in JMML" = "#E74C3C",
    "Detected more in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (RPNPF Platform)",
    caption = paste0("n = ", nrow(fish_plot_RPNPF_fdr), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: JMML vs. Control)\n↑ More detected in JMML, ↓ More detected in Control",
    y = "-log10(FDR)",
    color = "Relative Detection",
    shape = "Statistically Significant\n(FDR < 0.05)"
  )+
  theme_bw(base_size = 14)

print(fish_p_RPNPF_fdr)

ggsave(filename = "../res/annotated/fisher_volcano_RPNPF_FDR.png", plot = fish_p_RPNPF_fdr)

# ==== RPNPF - TTest Visualizations FDR ==========================================

ttest_RPNPF <- final_results %>%
  filter(test_type == "TTest") %>%
  filter(platform == "RPNPF")

# ==== Volcano Plot of Sex Differences ====

ttest_results_Combined_RPNPF_fdr <- ttest_RPNPF %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in JMML",
      log2_fc < 0           ~ "Higher in Control",
      TRUE                    ~ "Ambiguous"
    ),
    sig = ifelse(global_fdr < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

RPNPF_ttest_volcano_FDR <- ggplot(
  ttest_results_Combined_RPNPF_fdr,
  aes(x = log2_fc, y = -log10(global_fdr))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  geom_text_repel(
    data = ttest_results_Combined_RPNPF_fdr %>% filter(sig == "Significant"),
    aes(label = compound_name),
    size = 2.5,
    max.overlaps = 15,
    box.padding = 0.3
  ) +
  scale_color_manual(values = c(
    "Higher in JMML" = "#E74C3C",
    "Higher in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "T-test on Imputed Box-Cox Transformed Intensities\n(RPNPF Platform)",
    caption = paste0("n = ", nrow(ttest_results_Combined_RPNPF_fdr), " metabolites tested via two-sided t-test."),
    x = "log2FC(JMML / Control)\n↑ Higher in JMML, ↓ Higher in Control",
    y = "-log10(FDR)\nfrom t-test",
    color = "Relative Intensity",
    shape = "Statistically Significant\n(FDR < 0.05)"
  ) +
  theme_bw(base_size = 14)


print(RPNPF_ttest_volcano_FDR)

# Save volcano plot
ggsave("../res/annotated/RPNPF_ttest_volcano_FDR.png",
       RPNPF_ttest_volcano_FDR, width = 10, height = 6, dpi = 300)

# ==== Top 5 Violin Plots of Most Significant Sex Differences ====


# Identify top 5 metabolites by raw p-value
top5_Combined_RPNPF <- ttest_results_Combined_RPNPF_fdr %>%
  arrange(raw_p) %>%
  slice_head(n = 5) %>%
  pull(compound_name)

p_violin_RPNPF <- ggplot(
  RPNPF_transformed %>% filter(compound_name %in% top5_Combined_RPNPF),
  aes(x = sample_group, y = intensity_log2, fill = sample_group)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("JMML" = "#E74C3C", "Control" = "#3498DB")) +
  labs(
    title = "Top 5 JMML-associated Metabolites",
    subtitle = "T-test on Imputed Intensities (RPNPF Platform)",
    x = "Case Status",
    y = "log2(Intensity)",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

p_violin_RPNPF

# Save violin plot
ggsave(
  filename = "../res/annotated/RPNPF_top5_ttest.png",
  plot = p_violin_RPNPF,
  width = 10, height = 6, units = "in", dpi = 300
)






# ==== Raw-p value Visualizations ==============================================


# ==== Viz Fish RPNPF P-value  ===================================================

fish_RPNPF <- final_results %>%
  filter(test_type == "Fisher") %>%
  filter(platform == "RPNPF")

# Add enrichment direction labels
fish_plot_RPNPF_p <- fish_RPNPF %>%
  mutate(
    direction = case_when(
      is.na(log2_OR)             ~ "Ambiguous",
      log2_OR >= 1               ~ "Detected more in JMML",
      log2_OR <= -1              ~ "Detected more in Control",
      abs(log2_OR) < 1           ~ "Ambiguous"
    ),
    neglog10_p = -log10(raw_p),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")) %>% 
  mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

# Drop columns that are entirely NA
fish_plot_RPNPF_p_clean <- fish_plot_RPNPF_p %>%
  select(where(~ !all(is.na(.))))  %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write.csv(fish_plot_RPNPF_p_clean, file = "../res/annotated/clean_fisher_RPNPF_table.csv")

fish_p_RPNPF <- ggplot(fish_plot_RPNPF_p,  aes(x = log2_OR, y = neglog10_p)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray") +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7) +
  # geom_text_repel(data = fish_plot_RPNPF_p %>% filter(sig == "Significant"), 
  #                 aes(label = compound_name), 
  #                 size = 3, max.overlaps = 15) +
  scale_color_manual(values = c(
    "Detected more in JMML" = "#E74C3C",
    "Detected more in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "Fisher’s Exact Test on Missingness (RPNPF Platform)",
    caption = paste0("n = ", nrow(fish_plot_RPNPF_p ), 
                     " metabolites tested using Fisher’s Exact Test"),
    x = "log2(Odds Ratio: JMML vs. Control)\n↑ More detected in JMML, ↓ More detected in Control",
    y = "-log10(P-value)",
    color = "Relative Detection",
    shape = "Statistically Significant\n(P-value < 0.05)"
  )+
  theme_bw(base_size = 14)

print(fish_p_RPNPF)

ggsave(filename = "../res/annotated/fisher_volcano_RPNPF_pval.png", plot = fish_p_RPNPF)




# ==== RPNPF - TTest Visualizations FDR ==========================================

ttest_RPNPF <- final_results %>%
  filter(test_type == "TTest") %>%
  filter(platform == "RPNPF")

# ==== Volcano Plot of Sex Differences ====

ttest_results_Combined_RPNPF_p <- ttest_RPNPF %>%
  mutate(
    direction = case_when(
      log2_fc > 0            ~ "Higher in JMML",
      log2_fc < 0           ~ "Higher in Control",
      TRUE                    ~ "Ambiguous"
    ),
    sig = ifelse(raw_p < 0.05, "Significant", "Not Significant")
  ) %>% mutate(sig = factor(sig, levels = c("Not Significant", "Significant")))

RPNPF_ttest_volcano_pval <- ggplot(
  ttest_results_Combined_RPNPF_p ,
  aes(x = log2_fc, y = -log10(raw_p))
) +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
  scale_color_manual(values = c(
    "Higher in JMML" = "#E74C3C",
    "Higher in Control" = "#3498DB",
    "Ambiguous" = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) + # ● = not sig, ▲ = sig
  labs(
    title = "JMML-based Metabolite Differences",
    subtitle = "T-test on Imputed Box-Cox Transformed Intensities\n(RPNPF Platform)",
    caption = paste0("n = ", nrow(ttest_results_Combined_RPNPF_p), " metabolites tested via two-sided t-test."),
    x = "log2FC(JMML / Control)\n↑ Higher in JMML, ↓ Higher in Control",
    y = "-log10(p-value)\nfrom t-test",
    color = "Relative Intensity",
    shape = "Statistically Significant\n(p-value < 0.05)"
  ) +
  theme_bw(base_size = 14)


print(RPNPF_ttest_volcano_pval)

# Save volcano plot
ggsave("../res/annotated/RPNPF_ttest_volcano_pval.png",
       RPNPF_ttest_volcano_pval, width = 10, height = 6, dpi = 300)



# ==== Significant Hits ========================================================


# ==== Viz ttest sig hit FDR 2-hydro ===========================================

t_res <- final_results %>% filter(test_type == "TTest") %>% filter(global_fdr < .05)


ZHP_transformed <- readRDS("../data/R_objects/ZHP_imputed.rds") 



ttest_sig <- ggplot(ZHP_transformed %>% filter(compound_name == "2-hydroxychrysene"),
                    aes(x = sample_group, y = intensity_log2, fill = sample_group)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, shape = 21, stroke = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "white") +
  scale_fill_manual(values = c("Control" = "#3498DB", "JMML" = "#E74C3C")) +
  labs(
    title = "",
    x = "",
    y = "log2(Intensity)",
    fill = ""
  ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold")
    )
  print(ttest_sig)
  
#Distribution of log2(intensity) of 2-hydroxychrysene



ggsave(filename = "../res/annotated/ttest_violin_2-hydroxychrysene.png", plot = ttest_sig)


# ==== Viz violin sig hit Fisher ===========================================

fish_res_top <- final_results %>%
  filter(test_type == "Fisher") %>%
  filter(global_fdr < .05) %>% 
  pull(compound_name)


# Load the non-imputed intensities (M2) 
ZHP_filtered_fish <- readRDS("../data/R_objects/ZHP_annotated_filtered.rds") %>% 
  filter(sample_group %in% c("JMML", "Control")) %>% 
  filter(compound_name %in% fish_res_top) %>% 
  mutate(intensity_log2 = log2(intensity + 1))

writexl::write_xlsx(ZHP_filtered_fish, path = "../res/annotated/Fisher_FDR_hits.xlsx")

fish_res_violin_ZHP <- ggplot(
  ZHP_filtered_fish,
  aes(x = sample_group, y = intensity_log2, fill = sample_group)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = c("JMML" = "#E74C3C", "Control" = "#3498DB")) +
  labs(
    title = "Top 4 JMML-associated Metabolites ",
    subtitle = "Fisher's Exact Test on Missingness (ZHP Platform)",
    x = "Case Status",
    y = "log2(Intensity)",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

fish_res_violin_ZHP

# Save violin plot
ggsave(
  filename = "../res/annotated/ZHP_top4violin_fish_res.png",
  plot = fish_res_violin_ZHP,
  width = 10, height = 6, units = "in", dpi = 300
)


# ==== Lilly Poster ==========================================
set.seed(42)
# build the plot
fish_p_ZHP_fdr_pub <- ggplot(fish_plot_ZHP_fdr, aes(x = log2_OR, y = neglog10_fdr)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray80") +
  geom_vline(xintercept = c(-1, 1),   linetype = "dotted", color = "gray80") +
  geom_point(aes(color = direction, shape = sig), alpha = 0.7) +
  geom_text_repel(
    data          = filter(fish_plot_ZHP_fdr, sig == "Significant"),
    aes(label     = compound_name),
    size          = 3,
    box.padding   = 0.6,
    point.padding = 0.4,
    force         = 2,
    max.overlaps  = Inf
  ) +
  scale_color_manual(values = c(
    "Detected more in JMML"    = "#E74C3C",
    "Detected more in Control" = "#3498DB",
    "Ambiguous"                = "gray50"
  )) +
  scale_shape_manual(values = c("Not Significant" = 16, "Significant" = 17)) +
  labs(
    x = expression(Log[2](OR)~"(JMML/Control)"),
    y     = expression(-Log[10](FDR)),
    color = "Relative Detection",
    shape = "Statistically Significant\n(FDR < 0.05)"
  ) +
  theme_classic(base_size = 14, base_family = "Times") +
  theme(
    axis.title.x     = element_text(face   = "bold",
                                    size   = 14,
                                    margin = ggplot2::margin(t = 8)),
    axis.title.y     = element_text(face   = "bold",
                                    size   = 14,
                                    margin = ggplot2::margin(r = 8)),
    axis.text        = element_text(size = 12),
    legend.key.size  = unit(0.8, "lines"),
    legend.position = "bottom",
    legend.box      = "vertical",
    legend.title    = element_text(face = "bold"),
    legend.text     = element_text(size = 11)
  ) +
  guides(
    color = guide_legend(title.position = "top", nrow = 1),
    shape = guide_legend(title.position = "top", nrow = 1)
  )


print(fish_p_ZHP_fdr_pub)

ggsave(filename = "../res/annotated/fisher_volcano_ZHP_FDR_pub.png", plot = fish_p_ZHP_fdr_pub)
