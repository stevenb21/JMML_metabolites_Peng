# JMML Newborn DBS Metabolomics (Research Workflow)

This repo is a **research analysis workflow (not a packaged library)** for biomarker discovery and classification modeling using **untargeted LC–MS metabolomics** from **newborn dried blood spots (DBS)**. Dataset: **58 JMML cases + 58 controls** with two acquisition modes (**RPNPF** and **ZHP**).

![Workflow overview](Workflow.png)

**Poster:** [JMML_poster_StevenBrooks_2025.pdf](JMML_poster_StevenBrooks_2025.pdf) :contentReference[oaicite:1]{index=1}

## What’s inside
- `01_annotated/` — QC + statistical testing on annotated metabolites (`main.R` runs the numbered pipeline). :contentReference[oaicite:2]{index=2}
- `02_datamining/` — feature engineering + classification experiments (JMML vs control)
- `analysis/` — `Annotated_Metabolites.Rmd` and `Data_mining.Rmd` reports

## How to run
Run the pipeline:
- source("main.R")
- Knit: `analysis/Annotated_Metabolites.Rmd` and `analysis/Data_mining.Rmd`
  
Raw data are not included in this public repo; results/intermediates are generated locally.
