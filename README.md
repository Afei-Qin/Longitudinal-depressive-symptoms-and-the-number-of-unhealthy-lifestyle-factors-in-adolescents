# Longitudinal depressive symptoms and the number of unhealthy lifestyle factors in adolescents

This repository provides analysis code for the study:  
“Longitudinal depressive symptoms and the number of unhealthy lifestyle factors in adolescents”.

## Software and environment
- **R**: 4.5.1  
- **Mplus**: 8.3  

## Repository contents

### 01_mplus_RI-CLPMs_codes.inp
Mplus input scripts for Random Intercept Cross-Lagged Panel Models (RI-CLPMs), including:
- **Model A**: baseline RI-CLPM (no covariates)
- **Model B**: RI-CLPM with time-invariant covariates at T1 (age, sex, ethnicity, living status, region, campus)
- **Model C**: RI-CLPM with time-invariant covariates at T1 + time-varying covariates at T1–T3 (ACE and BMI)
- **Model D**: Model C with autoregressive paths constrained equal across time
- **Model E**: Model C with cross-lagged paths constrained equal across time
- **Model F**: Model C with both autoregressive and cross-lagged paths constrained equal across time

**Expected input:** analysis dataset referenced in the script (e.g., `fulldataRICLPM.csv`, not included).  
**Expected output:** Mplus model estimates and fit statistics (saved via Mplus output files).

### 02_rcs_CESD22_to_unlifestyle24.R
Restricted cubic spline (RCS) analysis evaluating the (non)linear association between baseline depressive symptoms (**CESD_22**) and follow-up unhealthy lifestyle count (**unlifestyle_24**), adjusted for covariates.

**Expected input:** R analysis dataset object (e.g., `be`) containing `CESD_22`, `unlifestyle_24`, and covariates used in the script.  
**Expected output:** spline curve plot (and model comparison output in console).

### 03_lcmm_CESD_trajectories_fit-indices.R
Latent class mixed model (LCMM) scripts for depressive symptom trajectories, including model fitting (1–4 classes) and extraction of fit indices and posterior probabilities.

**Expected input:** R analysis dataset (as defined in the script) with repeated CESD measures and required covariates/IDs.  
**Expected output:** fit indices and posterior probabilities (as written/exported in the script).

### 04_mixed-models_GLMM_LMM_pairwise_diagnostics.R
Mixed-effects models (GLMM/LMM) and pairwise comparisons (via **emmeans**), including diagnostics (e.g., **DHARMa**) and sensitivity models (e.g., NB/ZINB).

**Expected input:** R analysis dataset object (e.g., `da`) with clustering variable (e.g., `school_id`), outcomes, and covariates used in the models.  
**Expected output:** model summaries in console and exported result tables (CSV files, as specified in the script).

## How to run (overview)

> Note: The original data are not publicly available (see Data availability below).

With approved access to the data, the recommended running order is:
1. **01 (Mplus RI-CLPMs)**: run the `.inp` script(s) in Mplus 8.3 to obtain RI-CLPM estimates.
2. **02 (RCS)**: run the spline model and generate the effect curve (with 95% CI).
3. **03 (LCMM)**: fit 1–4 class trajectory models and export posterior probabilities / fit summaries.
4. **04 (Mixed-effects models)**: run GLMM/LMM and pairwise comparisons; export result tables.

## Data availability

The data used in this study were obtained from the provincial surveillance system on student health and contain sensitive information on adolescents. In accordance with ethical requirements and data protection regulations, the datasets are not publicly available. De-identified data may be accessed upon reasonable request for research purposes. Interested researchers may contact the corresponding author, **Dr. Lianlong Yu** (Email: **lianlong00a@163.com**), at the Shandong Center for Disease Control and Prevention. Access to the data requires on-site review and analysis at the Shandong CDC in compliance with institutional and provincial ethical guidelines.  
**The code for this study is publicly available in this repository.**

## Notes for reproducibility

- The scripts assume analysis datasets contain variables consistent with the file names and model specifications (e.g., `CESD_22`, `CESD_24`, `unlifestyle_22`, `unlifestyle_24`, and covariates used in each model).
- Update file paths/object names in the scripts to match your local environment and dataset naming (avoid hard-coded paths when adapting the code).
- If you encounter package-version differences, consider recording your session info (`sessionInfo()`) when reproducing results.

## License

This project is licensed under the MIT License (see `LICENSE`).
