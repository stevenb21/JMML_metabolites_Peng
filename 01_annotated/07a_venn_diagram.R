# ==== Load Libraries ====
library(tidyverse)
library(readxl)
library(ggVennDiagram)

# ==== Read and Process Sex Analysis Results ====
sex_res <- read_xlsx("../res/annotated/sex/sex_analysis.xlsx")
sex_hits_all <- sex_res %>%
  filter(cohort == "Combined_JMML_Control", raw_p < 0.05) %>%
  mutate(
    direction = case_when(
      test_type == "TTest" & log2_fc >  0 ~ "↑ Females",
      test_type == "TTest" & log2_fc <  0 ~ "↓ Females",
      test_type == "Fisher" & log2_OR > 0 ~ "↑ Females",
      test_type == "Fisher" & log2_OR < 0 ~ "↓ Females",
      TRUE                                  ~ ""
    )
  ) %>%
  transmute(hit = paste(compound_name, direction)) %>%
  distinct(hit) %>%
  pull(hit)

sex_hits_JMML <- sex_res %>%
  filter(cohort == "JMML_only", raw_p < 0.05) %>%
  mutate(
    direction = case_when(
      test_type == "TTest" & log2_fc >  0 ~ "↑ Females",
      test_type == "TTest" & log2_fc <  0 ~ "↓ Females",
      test_type == "Fisher" & log2_OR > 0 ~ "↑ Females",
      test_type == "Fisher" & log2_OR < 0 ~ "↓ Females",
      TRUE                                  ~ ""
    )
  ) %>%
  transmute(hit = paste(compound_name, direction)) %>%
  distinct(hit) %>%
  pull(hit)



sex_hits_Control <- sex_res %>%
  filter(cohort == "Control_only", raw_p < 0.05) %>%
  mutate(
    direction = case_when(
      test_type == "TTest" & log2_fc >  0 ~ "↑ Females",
      test_type == "TTest" & log2_fc <  0 ~ "↓ Females",
      test_type == "Fisher" & log2_OR > 0 ~ "↑ Females",
      test_type == "Fisher" & log2_OR < 0 ~ "↓ Females",
      TRUE                                  ~ ""
    )
  ) %>%
  transmute(hit = paste(compound_name, direction)) %>%
  distinct(hit) %>%
  pull(hit)


# ==== Read and Process Somatic Mutation Analysis Results ====
somatic_res <- read_xlsx("../res/annotated/somatic/somatic_analysis.xlsx")
somatic_hits <- somatic_res %>%
  filter(raw_p < 0.05) %>%
  mutate(
    direction = case_when(
      test_type == "TTest" & log2_fc >  0 ~ "↑ Somatic+",
      test_type == "TTest" & log2_fc <  0 ~ "↓ Somatic+",
      test_type == "Fisher" & log2_OR > 0 ~ "↑ Somatic+",
      test_type == "Fisher" & log2_OR < 0 ~ "↓ Somatic+",
      TRUE                                  ~ ""
    )
  ) %>%
  transmute(hit = paste(compound_name, direction)) %>%
  distinct(hit) %>%
  pull(hit)

# ==== Read and Process JMML vs Control Analysis Results ====
jmml_res <- read_xlsx("../res/annotated/final_results_global_fdr.xlsx")

jmml_hits <- jmml_res %>%
  filter(raw_p < 0.05) %>%
  mutate(
    direction = case_when(
      test_type == "TTest" & log2_fc >  0 ~ "↑ JMML",
      test_type == "TTest" & log2_fc <  0 ~ "↓ JMML",
      test_type == "Fisher" & log2_OR > 0 ~ "↑ JMML",
      test_type == "Fisher" & log2_OR < 0 ~ "↓ JMML",
      TRUE                                  ~ ""
    )
  ) %>%  
  transmute(hit = paste(compound_name, direction)) %>%
  distinct(hit) %>%
  pull(hit)

# ==== Create Venn List ====
venn_list <- list(
  Sex_Combined             = sex_hits_all,
#  Sex_JMML             = sex_hits_JMML,
#  Sex_Control = sex_hits_Control,
  Somatic         = somatic_hits,
  JMML_vs_Control = jmml_hits
)

# ==== Plot Venn Diagram ====
p_venn <- ggVennDiagram(
  venn_list,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Overlap of Significant Compounds\nAcross Analyses (raw p-value)",
    x = NULL, y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 70)
  ) + scale_x_continuous(expand = expansion(mult = .2))


print(p_venn)

# ==== Save Venn Diagram ====
ggsave(
  filename = "../res/annotated/venn_significant_hits_raw_p.png",
  plot     = p_venn,
  width    = 6,
  height   = 6,
  dpi      = 300
)

# ==== Create Venn List ====
venn_list_sex <- list(
  Sex_Combined             = sex_hits_all,
    Sex_JMML             = sex_hits_JMML,
   Sex_Control = sex_hits_Control
#  Somatic         = somatic_hits,
#  JMML_vs_Control = jmml_hits
)

# ==== Plot Venn Diagram ====
p_sex <- ggVennDiagram(
  venn_list_sex,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Overlap of Significant Compounds\nAcross Sex (raw p-value)",
    x = NULL, y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 70)
  ) + scale_x_continuous(expand = expansion(mult = .2))


print(p_sex)

# ==== Save Venn Diagram ====
ggsave(
  filename = "../res/annotated/venn_significant_hits_raw_p_sex.png",
  plot     = p_sex,
  width    = 6,
  height   = 6,
  dpi      = 300
)

