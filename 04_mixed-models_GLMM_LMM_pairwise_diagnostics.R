## 04_mixed-models_GLMM_LMM_pairwise_diagnostics.R
## Generalized linear mixed-effects models (GLMM), LMM, pairwise comparisons, and diagnostics
## Minimal necessary fixes for public reproducibility:
## 1) create results/ folder and write outputs there (no scattered files)
## 2) make unlifestyle_cat a factor (0/1) for logistic GLMM
## 3) remove View() (non-reproducible in non-interactive runs)
## 4) keep everything else as close as possible to your original script

## ---- Packages ----
library(lme4)
library(lmerTest)
library(performance)
library(dplyr)
library(broom.mixed)
library(emmeans)
library(DHARMa)
library(pscl)
library(glmmTMB)

## ---- Output folder (minimal fix) ----
dir.create("results", showWarnings = FALSE, recursive = TRUE)

## ---- Data prep (same as your original) ----
## 'da' is expected to exist in the environment
da <- da %>%
  mutate(
    school_id  = factor(school_id),
    male       = factor(male),
    minzu      = factor(minzu),
    region     = factor(region),
    acedemic   = factor(acedemic_22),
    living     = factor(living_22),
    Livecampus = factor(Livecampus_22),
    class      = factor(class)
  )

## Minimal fix: make binary outcome a factor for logistic GLMM
da$unlifestyle_cat <- factor(ifelse(da$unlifestyle_24 >= 4, 1, 0), levels = c(0, 1))

# primary poisson GLMM #
model_poisson <- glmer(
  unlifestyle_24 ~ age_z + male + minzu + region +
    acedemic + bmi_z + living + Livecampus + unlifestyle_22 + class +
    (1 | school_id),
  family = poisson,
  data = da
)
summary(model_poisson)

emm_model_poisson <- emmeans(
  model_poisson, ~ class,
  nuisance = c("age_z", "male", "minzu", "region", "acedemic", "bmi_z", "living",
               "Livecampus", "unlifestyle_22")
)

contrast_list <- list(
  "2-1" = c(-1, 1, 0),
  "3-1" = c(-1, 0, 1),
  "3-2" = c(0, -1, 1)
)

contr_model_poisson <- contrast(emm_model_poisson, contrast_list) %>%
  summary(infer = c(TRUE, TRUE))

results_model_poisson <- as.data.frame(contr_model_poisson) %>%
  mutate(
    variable = "unlifestyle_24",
    p_adj_BH = p.adjust(p.value, method = "BH"),
    beta = estimate,
    beta_CI_lower = estimate - 1.96 * SE,
    beta_CI_upper = estimate + 1.96 * SE
  ) %>%
  select(contrast, beta, beta_CI_lower, beta_CI_upper, p.value, p_adj_BH, variable)

summary_model_poisson <- summary(model_poisson)
est_poisson <- summary_model_poisson$coefficients[, "Estimate"]
se_poisson  <- summary_model_poisson$coefficients[, "Std. Error"]
pval_poisson <- summary_model_poisson$coefficients[, "Pr(>|z|)"]
lower_poisson <- est_poisson - 1.96 * se_poisson
upper_poisson <- est_poisson + 1.96 * se_poisson
result_table_poisson <- data.frame(
  beta = est_poisson,
  Std.Error = se_poisson,
  CI.Lower = lower_poisson,
  CI.Upper = upper_poisson,
  P.value = pval_poisson
)

contr_model_poissonOR <- contrast(emm_model_poisson, contrast_list, type = "response") %>%
  summary(infer = c(TRUE, TRUE))

results_model_poissonOR <- as.data.frame(contr_model_poissonOR) %>%
  mutate(
    variable = "unlifestyle_24",
    OR = ratio,
    OR_CI_lower = asymp.LCL,
    OR_CI_upper = asymp.UCL,
    p_adj_BH = p.adjust(p.value, method = "BH")
  ) %>%
  select(contrast, OR, OR_CI_lower, OR_CI_upper, p.value, p_adj_BH, variable)

summary_model_poisson <- summary(model_poisson)
est_poisson <- summary_model_poisson$coefficients[, "Estimate"]
se_poisson  <- summary_model_poisson$coefficients[, "Std. Error"]
pval_poisson <- summary_model_poisson$coefficients[, "Pr(>|z|)"]
OR_poisson <- exp(est_poisson)
OR_CI_lower_poisson <- exp(est_poisson - 1.96 * se_poisson)
OR_CI_upper_poisson <- exp(est_poisson + 1.96 * se_poisson)
result_table_poissonOR <- data.frame(
  Variable = rownames(summary_model_poisson$coefficients),
  OR = OR_poisson,
  OR_CI_Lower = OR_CI_lower_poisson,
  OR_CI_Upper = OR_CI_upper_poisson,
  P_value = pval_poisson
)

## Minimal fix: write outputs into results/
write.csv(result_table_poissonOR, file.path("results", "result_table_poissonOR.csv"), row.names = FALSE)
write.csv(results_model_poissonOR, file.path("results", "results_model_poissonOR.csv"), row.names = FALSE)
###########

# 1）Model S1: baseline CESD and unlifestyle_24
da$CESD_22_z <- scale(da$CESD_22)
model_poisson_CESD <- glmer(
  unlifestyle_24 ~ age_z + male + minzu + region +
    acedemic + bmi_z + living + Livecampus + unlifestyle_22 + CESD_22_z +
    (1 | school_id),
  family = poisson,
  data = da
)
summary(model_poisson_CESD)

summary_model_poisson_CESD <- summary(model_poisson_CESD)
est_poisson_CESD <- summary_model_poisson_CESD$coefficients[, "Estimate"]
se_poisson_CESD  <- summary_model_poisson_CESD$coefficients[, "Std. Error"]
pval_poisson_CESD <- summary_model_poisson_CESD$coefficients[, "Pr(>|z|)"]
OR_poisson_CESD <- exp(est_poisson_CESD)
OR_CI_lower_poisson_CESD <- exp(est_poisson_CESD - 1.96 * se_poisson_CESD)
OR_CI_upper_poisson_CESD <- exp(est_poisson_CESD + 1.96 * se_poisson_CESD)
result_table_poisson_CESD <- data.frame(
  Variable = rownames(summary_model_poisson_CESD$coefficients),
  OR = OR_poisson_CESD,
  OR_CI_Lower = OR_CI_lower_poisson_CESD,
  OR_CI_Upper = OR_CI_upper_poisson_CESD,
  P_value = pval_poisson_CESD
)
write.csv(result_table_poisson_CESD, file.path("results", "result_table_poisson_CESD.csv"), row.names = FALSE)

# 2）Model S2: LMM
model_normal <- lmer(
  unlifestyle_24 ~ age_z + male + minzu + region +
    acedemic + bmi_z + living + Livecampus + unlifestyle_22 + class +
    (1 | school_id),
  data = da
)
summary(model_normal)

emm_model_normal <- emmeans(
  model_normal, ~ class,
  nuisance = c("age_z", "male", "minzu", "region", "acedemic", "bmi_z", "living",
               "Livecampus", "unlifestyle_22")
)

contrast_list <- list(
  "2-1" = c(-1, 1, 0),
  "3-1" = c(-1, 0, 1),
  "3-2" = c(0, -1, 1)
)

contr_model_normal <- contrast(emm_model_normal, contrast_list) %>%
  summary(infer = c(TRUE, TRUE))

results_model_normal <- as.data.frame(contr_model_normal) %>%
  mutate(
    variable = "unlifestyle_24",
    p_adj_BH = p.adjust(p.value, method = "BH"),
    beta = estimate,
    beta_CI_lower = estimate - 1.96 * SE,
    beta_CI_upper = estimate + 1.96 * SE
  ) %>%
  select(contrast, beta, beta_CI_lower, beta_CI_upper, p.value, p_adj_BH, variable)

summary_model_normal <- summary(model_normal)
est_normal <- summary_model_normal$coefficients[, "Estimate"]
se_normal  <- summary_model_normal$coefficients[, "Std. Error"]
pval_normal <- summary_model_normal$coefficients[, "Pr(>|t|)"]
lower_normal <- est_normal - 1.96 * se_normal
upper_normal <- est_normal + 1.96 * se_normal
result_table_normal <- data.frame(
  beta = est_normal,
  Std.Error = se_normal,
  CI.Lower = lower_normal,
  CI.Upper = upper_normal,
  P.value = pval_normal
)

# 3）Model S3: GLMM-logit
model_cate <- glmer(
  unlifestyle_cat ~ age_z + male + minzu + region +
    acedemic + bmi_z + living + Livecampus + unlifestyle_22 + class +
    (1 | school_id),
  family = binomial(link = "logit"),
  data = da,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
summary(model_cate)

emm_model_cate <- emmeans(
  model_cate, ~ class,
  nuisance = c("age_z", "male", "minzu", "region", "acedemic", "bmi_z", "living",
               "Livecampus", "unlifestyle_22")
)

contrast_list <- list(
  "2-1" = c(-1, 1, 0),
  "3-1" = c(-1, 0, 1),
  "3-2" = c(0, -1, 1)
)

contr_model_cate <- contrast(emm_model_cate, contrast_list, type = "response") %>%
  summary(infer = c(TRUE, TRUE))

results_model_cate <- as.data.frame(contr_model_cate) %>%
  mutate(
    variable = "unlifestyle_cat",
    OR = odds.ratio,
    OR_CI_lower = asymp.LCL,
    OR_CI_upper = asymp.UCL,
    p_adj_BH = p.adjust(p.value, method = "BH")
  ) %>%
  select(contrast, OR, OR_CI_lower, OR_CI_upper, p.value, p_adj_BH, variable)

summary_model_cate <- summary(model_cate)
est_cate <- summary_model_cate$coefficients[, "Estimate"]
se_cate  <- summary_model_cate$coefficients[, "Std. Error"]
pval_cate <- summary_model_cate$coefficients[, "Pr(>|z|)"]
OR_cate <- exp(est_cate)
OR_CI_lower_cate <- exp(est_cate - 1.96 * se_cate)
OR_CI_upper_cate <- exp(est_cate + 1.96 * se_cate)
result_table_cate <- data.frame(
  Variable = rownames(summary_model_cate$coefficients),
  OR = OR_cate,
  OR_CI_Lower = OR_CI_lower_cate,
  OR_CI_Upper = OR_CI_upper_cate,
  P_value = pval_cate
)

all_results0 <- bind_rows(result_table_poisson, result_table_poissonOR, result_table_normal, result_table_cate)

## Minimal fix: write outputs into results/
write.csv(results_model_poisson, file.path("results", "results_model_poisson.csv"), row.names = FALSE)
write.csv(results_model_normal,  file.path("results", "results_model_normal.csv"),  row.names = FALSE)
write.csv(results_model_cate,    file.path("results", "results_model_cate.csv"),    row.names = FALSE)
write.csv(all_results0,          file.path("results", "all_results0.csv"),          row.names = FALSE)

# Checking for Over-Diversification Issues #
sim_poisson <- simulateResiduals(fittedModel = model_poisson, n = 1000)
plot(sim_poisson)
testOverdispersion(sim_poisson)

## The baseline lifestyle quantity and the unordered grouping of follow-up depression and trajectory
## 1) baseline number of unhealthy lifestyle factors and T3 depressive symptoms
model_poisson24 <- glmer(
  CESD_24 ~ age_z + male + minzu + region +
    acedemic + bmi_z + living + Livecampus + CESD_22 + unlifestyle_22 +
    (1 | school_id),
  family = poisson,
  data = da
)
summary(model_poisson24)

sim_poisson24 <- simulateResiduals(fittedModel = model_poisson24, n = 1000)
plot(sim_poisson24)
testOverdispersion(sim_poisson24)

model_nb24 <- glmer.nb(
  CESD_24 ~ age_z + male + minzu + region +
    acedemic + bmi_z + living + Livecampus + CESD_22 + unlifestyle_22 +
    (1 | school_id),
  data = da,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5),
    check.conv.grad = .makeCC("warning", tol = 0.02),
    check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)
  )
)
summary(model_nb24)

model_zinb24 <- glmmTMB(
  CESD_24 ~ age_z + male + minzu + region +
    acedemic + bmi_z + living + Livecampus + CESD_22 + unlifestyle_22 +
    (1 | school_id),
  ziformula = ~1,
  family = nbinom2,
  data = da
)
summary(model_zinb24)
summary_model_zinb24 <- summary(model_zinb24)

cond_coef <- summary_model_zinb24$coefficients$cond
est_cond <- cond_coef[, "Estimate"]
se_cond <- cond_coef[, "Std. Error"]
p_cond <- cond_coef[, "Pr(>|z|)"]

zi_coef <- summary_model_zinb24$coefficients$zi
est_zi <- zi_coef[, "Estimate"]
se_zi <- zi_coef[, "Std. Error"]
p_zi <- zi_coef[, "Pr(>|z|)"]

IRR_cond <- exp(est_cond)
IRR_cond_CI_lower <- exp(est_cond - 1.96 * se_cond)
IRR_cond_CI_upper <- exp(est_cond + 1.96 * se_cond)

OR_zi <- exp(est_zi)
OR_zi_CI_lower <- exp(est_zi - 1.96 * se_zi)
OR_zi_CI_upper <- exp(est_zi + 1.96 * se_zi)

results_cond <- data.frame(
  Variable = rownames(cond_coef),
  Estimate = est_cond,
  IRR = IRR_cond,
  IRR_CI_lower = IRR_cond_CI_lower,
  IRR_CI_upper = IRR_cond_CI_upper,
  p_value = p_cond
)
results_cond

results_zi <- data.frame(
  Variable = rownames(zi_coef),
  Estimate = est_zi,
  OR = OR_zi,
  OR_CI_lower = OR_zi_CI_lower,
  OR_CI_upper = OR_zi_CI_upper,
  p_value = p_zi
)
results_zi

sim_nb24 <- simulateResiduals(fittedModel = model_nb24, n = 1000)
sim_zinb24 <- simulateResiduals(fittedModel = model_zinb24, n = 1000)

compareLL <- function(model1, model2) {
  ll1 <- logLik(model1)
  ll2 <- logLik(model2)
  test_stat <- -2 * (as.numeric(ll1) - as.numeric(ll2))
  df <- attr(ll2, "df") - attr(ll1, "df")
  p_value <- pchisq(test_stat, df, lower.tail = FALSE)
  return(data.frame(Test_Statistic = test_stat, DF = df, P_Value = p_value))
}

lr_test <- compareLL(model_nb24, model_zinb24)
cat("Likelihood Ratio Test (NB vs ZINB):\n")
print(lr_test)

## 2) baseline number of unhealthy lifestyle factors and depressive symptom trajectories
da <- da %>%
  mutate(
    school_id  = factor(school_id),
    male       = factor(male),
    minzu      = factor(minzu),
    region     = factor(region),
    acedemic   = factor(acedemic_22),
    living     = factor(living_22),
    Livecampus = factor(Livecampus_22)
  )

run_glmer_pairwise_corrected <- function(data, response_var, vars_of_interest, nuisance_vars = NULL, adjust_method = "BH") {
  levels_resp <- sort(unique(data[[response_var]]))
  if (length(levels_resp) < 2) stop("There are only two levels of the dependent variable, and thus no comparison can be made.")
  comb_resp <- t(combn(levels_resp, 2))
  all_results <- data.frame()
  models_list <- list()

  for (i in seq_len(nrow(comb_resp))) {
    g_small <- comb_resp[i, 1]
    g_big   <- comb_resp[i, 2]
    cat("The comparison is being processed:", g_big, "vs", g_small, "\n")
    subdat <- subset(data, data[[response_var]] %in% c(g_small, g_big))
    subdat[[response_var]] <- factor(subdat[[response_var]], levels = c(g_small, g_big))
    if (nrow(subdat) < 10) {
      cat("Due to insufficient data, the comparison has been skipped:", g_big, "vs", g_small, "\n")
      next
    }
    all_vars <- c(vars_of_interest, nuisance_vars)
    fml <- as.formula(
      paste0(response_var, " ~ ", paste(all_vars, collapse = " + "), " + (1 | school_id)")
    )
    tryCatch({
      fit <- glmer(
        fml, data = subdat, family = binomial(link = "logit"),
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )
      models_list[[paste0(g_big, "_vs_", g_small)]] <- fit

      for (v in vars_of_interest) {
        emm <- emmeans(fit, specs = v, type = "response")
        contr <- pairs(emm) %>%
          summary(infer = TRUE) %>%
          as.data.frame()
        contr_df <- contr %>%
          mutate(
            variable = v,
            contrast_resp = paste0(g_big, "_vs_", g_small),
            comparison = paste(contrast)
          ) %>%
          rename(
            estimate_em = estimate,
            SE_em = SE,
            p_value_em = p.value
          )
        all_results <- bind_rows(all_results, contr_df)
      }
    }, error = function(e) {
      cat("Model fitting failed:", g_big, "vs", g_small, "-", e$message, "\n")
    })
  }

  if (nrow(all_results) > 0) {
    all_results <- all_results %>%
      group_by(variable) %>%
      mutate(p_adj_BH = p.adjust(p_value_em, method = adjust_method)) %>%
      ungroup()
  }
  return(list(
    models = models_list,
    emmeans_results = all_results
  ))
}

extract_fixed_effects <- function(models_list) {
  all_coef <- data.frame()
  for (nm in names(models_list)) {
    fit <- models_list[[nm]]
    tryCatch({
      coef_df <- broom.mixed::tidy(fit, effects = "fixed") %>%
        mutate(
          OR = exp(estimate),
          CI_lower = exp(estimate - 1.96 * std.error),
          CI_upper = exp(estimate + 1.96 * std.error),
          contrast_resp = nm
        )
      all_coef <- bind_rows(all_coef, coef_df)
    }, error = function(e) {
      cat("Failed to extract the coefficient.:", nm, "-", e$message, "\n")
    })
  }
  return(all_coef)
}

result_corrected <- run_glmer_pairwise_corrected(
  da,
  response_var = "class",
  vars_of_interest = c("unlifestyle_22"),
  nuisance_vars = c("age_z", "male", "minzu", "region", "acedemic", "bmi_z", "living", "Livecampus", "CESD_22")
)

coef_results <- extract_fixed_effects(result_corrected$models)

## Minimal fix: remove View() for non-interactive reproducibility
print(coef_results)
write.csv(coef_results, file.path("results", "coef_results_pairwise_class.csv"), row.names = FALSE)
write.csv(result_corrected$emmeans_results, file.path("results", "emmeans_pairwise_class.csv"), row.names = FALSE)
