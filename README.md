# Longitudinal depressive symptoms and the number of unhealthy lifestyle factors in adolescents

This repository provides analysis code for the study: **“Longitudinal depressive symptoms and the number of unhealthy lifestyle factors in adolescents”**.

## Repository contents

- **01_mplus_RI-CLPMs_codes.inp**  
  Mplus input scripts for Random Intercept Cross-Lagged Panel Models (RI-CLPMs), including:
  - Model A: baseline model (no covariates)
  - Model B: time-invariant covariates (T1)
  - Model C: time-invariant + time-varying covariates (T1–T3)
  - Model D: autoregressive paths constrained equal across time
  - Model E: cross-lagged paths constrained equal across time
  - Model F: autoregressive + cross-lagged paths constrained equal across time

- **02_rcs_CESD22_to_unlifestyle24.R**  
  Restricted cubic spline (RCS) analysis evaluating the (non)linear association between baseline depressive symptoms (CESD_22) and follow-up unhealthy lifestyle count (unlifestyle_24), adjusted for covariates.

- **03_lcmm_CESD_trajectories_fit-indices.R**  
  Latent class mixed model (LCMM) scripts for depressive symptom trajectories, including model fitting (1–4 classes) and extraction of fit indices and posterior probabilities.

- **04_mixed-models_GLMM_LMM_pairwise_diagnostics.R**  
  Mixed-effects models (GLMM/LMM) and pairwise comparisons (via emmeans), including diagnostics (e.g., DHARMa) and sensitivity models (e.g., NB/ZINB).

## How to run (overview)

> **Note:** The original data are not publicly available (see *Data availability* below).  
> With approved access to the data, the recommended running order is:

1. **01 (Mplus RI-CLPMs):** run the `.inp` script(s) in Mplus to obtain RI-CLPM estimates.
2. **02 (RCS):** run the spline model and generate the effect curve (with 95% CI).
3. **03 (LCMM):** fit 1–4 class trajectory models and export posterior probabilities / fit summaries.
4. **04 (Mixed-effects models):** run GLMM/LMM and pairwise comparisons; export result tables.

## Data availability

The data used in this study were obtained from the provincial surveillance system on student health and contain sensitive information on adolescents. In accordance with ethical requirements and data protection regulations, the datasets are not publicly available. De-identified data may be accessed upon reasonable request for research purposes. Interested researchers may contact the corresponding author, Dr. Lianlong Yu (Email: lianlong00a@163.com), at the Shandong Center for Disease Control and Prevention. Access to the data requires on-site review and analysis at the Shandong CDC in compliance with institutional and provincial ethical guidelines.

## Notes for reproducibility

- The scripts assume the analysis datasets contain variables consistent with the file names (e.g., `CESD_22`, `CESD_24`, `unlifestyle_22`, `unlifestyle_24`, and covariates used in each model).
- File paths in scripts should be set according to your local environment (avoid using `setwd("Desktop")` if you are adapting the code).

## License

This project is licensed under the MIT License (see `LICENSE`).
