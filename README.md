# Adiposity-associated serum lipid signatures correlate with efficacy and toxicity of CAR-T cell therapy in lymphoid malignancies

This repository contains the R code used for the manuscript:

**Adiposity-associated serum lipid signatures correlate with efficacy and toxicity of CAR-T cell therapy in lymphoid malignancies**

The analyses integrate longitudinal serum lipidomics, clinical outcomes, CT-derived body-composition measures, validation-cohort analyses, human CAR-T cell lipid supplementation experiments, and supernatant metabolomics from CAR-T co-culture assays.

## Repository contents

| File | Description |
|---|---|
| `0_packages.R` | Loads required R packages and analysis-wide settings. Run this first. |
| `1_surv.R` | Discovery-cohort baseline analyses, adiposity-associated lipid selection, PFS/OS modeling, Kaplan-Meier plots, and lipid-class outcome summaries. |
| `2_crs.R` | CRS-associated metabolite analyses, including early post-infusion lipid associations and lipid-class summary analyses. |
| `3_time.R` | Longitudinal serum lipid dynamics across baseline, early post-infusion, and day 14 time points. |
| `4_validation.R` | Validation-cohort analyses, including application of discovery-derived lipid signatures and association testing. |
| `5_invitro.R` | Human CAR-T lipid supplementation experiments, including killing assays and flow-cytometry activation, redox, and metabolic marker analyses. |
| `6_supernatant.R` | Supernatant metabolomics analyses from CAR-T co-culture experiments. |
| `7_baseline_char.R` | Baseline characteristics and summary tables for the validation cohort. |
| `.gitignore` | Specifies local files and output folders that should not be tracked by Git. |
| `CART_metabolomics.Rproj` | RStudio project file. |
| `README.md` | Repository overview and usage instructions. |

## Study overview

This study evaluates whether circulating lipidomic profiles associated with host adiposity correlate with clinical outcomes after CAR-T cell therapy in lymphoid malignancies.

Serum lipidomics was analyzed longitudinally at baseline, early post-infusion, and day 14. Lipid features were integrated with CT-derived visceral adipose tissue, subcutaneous adipose tissue, and total adipose tissue measurements, as well as clinical endpoints including response, progression-free survival, overall survival, cytokine release syndrome, and immune effector cell-associated neurotoxicity syndrome where available.

Revision-added experiments include human CAR-T cell lipid supplementation assays using BCMA CAR-T cells and NCI-H929 target cells, as well as supernatant metabolomics from CAR-T co-culture experiments.

## Analysis modules

### 1. Baseline discovery and survival analyses

Implemented in:

    source("1_surv.R")

This script performs the main discovery-cohort baseline analyses, including:

- clinical and metabolomics data loading;
- QC-based metabolite filtering;
- adiposity-associated metabolite selection;
- correlation with VAT, SAT, and TAT compartments;
- Cox proportional-hazards modeling for PFS and OS;
- Kaplan-Meier visualization;
- lipid-class-level outcome summaries.

### 2. CRS-associated metabolite analyses

Implemented in:

    source("2_crs.R")

This script evaluates lipid features associated with cytokine release syndrome, including:

- early post-infusion metabolite analyses;
- CRS-grade association testing;
- adiposity-stratified metabolite comparisons;
- lipid-class summary analyses;
- CRS-associated figure generation.

### 3. Longitudinal lipid dynamics

Implemented in:

    source("3_time.R")

This script evaluates changes in lipid species across CAR-T therapy time points:

- baseline / Day 0;
- early post-infusion / Day 3–5;
- intermediate post-infusion / Day 14.

The analysis generates longitudinal plots stratified by adiposity group and lipid class where applicable.

### 4. Validation-cohort analyses

Implemented in:

    source("4_validation.R")

This script applies discovery-derived features and signatures to an independent validation cohort, including:

- validation-cohort clinical annotation;
- application of discovery-derived lipid signatures and cutoffs;
- CRS and survival association analyses;
- validation-specific plots and tables.

### 5. Human CAR-T lipid supplementation experiments

Implemented in:

    source("5_invitro.R")

This script analyzes the human in vitro experiments added during manuscript revision, including:

- BCMA CAR-T cells with 4-1BB and CD28ζ costimulatory domains;
- NCI-H929 target cells;
- untransduced T-cell controls where applicable;
- acylcarnitine supplementation;
- physiologic and supraphysiologic lipid-supplementation conditions;
- CAR-T killing score;
- activation markers including CD69, CD25, and PD-1;
- redox and metabolic markers including cROS, GSH, MitoRed, MitoGreen, and mSOX;
- construct-specific metabolic remodeling.

### 6. Supernatant metabolomics from co-culture experiments

Implemented in:

    source("6_supernatant.R")

This script analyzes metabolomics data generated from the supernatants of the CAR-T co-culture experiments added during revision. These analyses complement the flow-cytometry and killing-assay readouts by evaluating extracellular metabolite changes after lipid supplementation and CAR-T/target-cell co-culture.

### 7. Baseline characteristics

Implemented in:

    source("7_baseline_char.R")

This script generates baseline characteristic summaries for the validation cohort.

## How to run the analyses

Open the repository in RStudio using:

    CART_metabolomics.Rproj

Then run the package setup script:

    source("0_packages.R")

After that, run the analysis scripts as needed.

For example, to run the main baseline survival analysis:

    source("1_surv.R")

To run the revision-added in vitro and supernatant analyses:

    source("5_invitro.R")
    source("6_supernatant.R")

Scripts should be run from the repository root directory so that relative file paths resolve correctly.

## Recommended run order

For a full analysis run, use:

    source("0_packages.R")
    source("1_surv.R")
    source("2_crs.R")
    source("3_time.R")
    source("4_validation.R")
    source("5_invitro.R")
    source("6_supernatant.R")
    source("7_baseline_char.R")

Some scripts may require intermediate files generated by earlier scripts or local input files that are not included in this public repository.

## Requirements

The analyses were developed in R.

Recommended:

- R >= 4.3.0
- RStudio
- Local access to the required clinical, body-composition, metabolomics, and in vitro input files

Required packages are loaded in:

    source("0_packages.R")

Major package groups include:

- data handling: `readxl`, `openxlsx`, `readr`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`, `forcats`;
- visualization: `ggplot2`, `ggpubr`, `ggrepel`, `ggprism`, `ComplexHeatmap`, `ggvenn`, `ggVennDiagram`;
- survival and clinical outcomes: `survival`, `survminer`, `ggsurvfit`, `tidycmprsk`;
- statistics and modeling: `rstatix`, `broom`, `emmeans`, `pROC`, `cutpointr`, `metafor`, `caret`, `randomForest`;
- metabolomics utilities: `MetaboAnalystR`, `pls`, `igraph`.

## Input data

The scripts expect local analysis-ready input files used during manuscript development. These include:

- clinical annotation files;
- CT-derived body-composition data;
- processed serum metabolomics data;
- validation-cohort data;
- in vitro CAR-T lipid supplementation readouts;
- supernatant metabolomics data.

Patient-level clinical data and raw metabolomics files will be deposited for public access once the study is published.


## Output files

The scripts generate manuscript and supplementary outputs, including:

- Kaplan-Meier plots;
- Cox model summaries;
- CRS-associated metabolite plots;
- longitudinal lipid dynamics plots;
- lipid-class meta-analysis plots;
- validation-cohort figures and tables;
- in vitro CAR-T killing and flow-marker plots;
- supernatant metabolomics outputs.

Generated figures, tables, and intermediate output folders are not intended to be version-controlled unless specifically needed for manuscript review.

## Metabolomics preprocessing

Serum metabolomics data were processed using semi-quantitative LC-MS workflows. Internal standards were used to monitor extraction, injection, instrument stability, and drift. Lamivudine-normalized relative abundance values were used for downstream analyses where indicated.

Metabolite features were filtered according to quality-control criteria before statistical analysis. Lipid-class annotations and metabolite labels were harmonized before figure generation.


## Citation

Please cite the associated manuscript once published.

## Contact

Email: davidm_cordasdossantos@dfci.harvard.edu

For questions about the analysis code or data access, please contact the corresponding study authors.
