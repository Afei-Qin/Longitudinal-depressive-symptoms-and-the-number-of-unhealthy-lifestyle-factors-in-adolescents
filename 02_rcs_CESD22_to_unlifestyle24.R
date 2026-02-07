## Restricted cubic spline (RCS): CESD_22 -> unlifestyle_24 ##
# Packages
library(ggplot2)
library(rms)

# Data: be (analysis dataset)
# NOTE: The raw data are not publicly available due to ethics/privacy restrictions
# This script assumes an analysis dataset named `be` has been prepared with the variables used below
# Make categorical covariates factors
be$male_22        <- factor(be$male_22)
be$minzu_22       <- factor(be$minzu_22)
be$acedemic_22    <- factor(be$acedemic_22)
be$living_22      <- factor(be$living_22)
be$region_22      <- factor(be$region_22)
be$Livecampus_22  <- factor(be$Livecampus_22)

# rms needs datadist()
dd <- datadist(be)
options(datadist = "dd")

# Fit models with different numbers of knots
fit3 <- ols(unlifestyle_24 ~ rcs(CESD_22, 3) + z_age_22 + male_22 + minzu_22 +
              acedemic_22 + living_22 + region_22 + Livecampus_22 + z_bmi_22 +
              unlifestyle_22,
            data = be, x = TRUE, y = TRUE)

fit4 <- ols(unlifestyle_24 ~ rcs(CESD_22, 4) + z_age_22 + male_22 + minzu_22 +
              acedemic_22 + living_22 + region_22 + Livecampus_22 + z_bmi_22 +
              unlifestyle_22,
            data = be, x = TRUE, y = TRUE)

fit5 <- ols(unlifestyle_24 ~ rcs(CESD_22, 5) + z_age_22 + male_22 + minzu_22 +
              acedemic_22 + living_22 + region_22 + Livecampus_22 + z_bmi_22 +
              unlifestyle_22,
            data = be, x = TRUE, y = TRUE)

# Compare AIC
AIC(fit3, fit4, fit5)

# Choose the final model (here: 3 knots)
fit <- fit3
an  <- anova(fit)
print(an)

# rms default plot (optional)
plot(Predict(fit, CESD_22), anova = an, pval = TRUE)

# ggplot curve with 95% CI
beta <- Predict(fit, CESD_22, ref.zero = TRUE)
beta_df <- as.data.frame(beta)

# ---- minimal, robust p extraction ----
rn <- trimws(rownames(an))                           # 去掉行名前后空格
pcol <- grep("^P", colnames(an), value = TRUE)[1]    # 找到P列（通常就是"P"）
idx_nl <- which(tolower(rn) == "nonlinear")          # 找 Nonlinear 行（忽略空格差异）

p_nl <- if (length(idx_nl) == 1) an[idx_nl, pcol] else NA
lab  <- if (is.na(p_nl)) "p for non-linearity: NA" else
  if (p_nl < 0.001) "p for non-linearity < 0.001"
  else paste0("p for non-linearity = ", format.pval(p_nl, digits = 3, eps = 0.001))
# -------------------------------------

ggplot(beta_df, aes(x = CESD_22, y = yhat)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  theme_classic() +
  labs(x = "Baseline depressive symptoms (CESD_22)",
       y = "Beta (95% CI)") +
  annotate("text", x = 12, y = 0.5, label = lab, size = 4)
