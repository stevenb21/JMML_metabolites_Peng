# Load libraries
library(tidyverse)
library(readxl)
library(writexl)

# Load differential results
final_results <- read_xlsx("../res/annotated/final_results_global_fdr.xlsx")

# Filter raw p < 0.05 and remove "Bisphenol AP"
sig_compounds <- final_results %>%
  filter(raw_p < 0.05) %>%
  distinct(compound_name) %>%
  pull(compound_name)

# Turn your character vector into a tibble before using mutate()
pathway_template <- tibble(compound_name = sig_compounds)

# Save to CSV for manual editing
write_csv(pathway_template, "../res/annotated/pathway_table_template.csv")
