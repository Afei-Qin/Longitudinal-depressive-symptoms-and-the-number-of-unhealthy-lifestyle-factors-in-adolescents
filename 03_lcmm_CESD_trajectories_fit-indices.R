## LCMM ##
library(lcmm)
library(dplyr)

set.seed(123)

## ---- Minimal fix 1: avoid hard-coded Desktop path (use repo-relative output folders) ----
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("models",  showWarnings = FALSE, recursive = TRUE)
dir.create("logs",    showWarnings = FALSE, recursive = TRUE)

result_csv <- file.path("results", "CESD.csv")
log_file   <- file.path("logs", "CESD.txt")

## ---- Data types (same as your original) ----
df_long_CESD$ID <- as.integer(df_long_CESD$ID)
df_long_CESD$Wave <- as.factor(df_long_CESD$Wave)
df_long_CESD$Sex <- as.factor(df_long_CESD$Sex)
df_long_CESD$Age <- as.numeric(df_long_CESD$Age)
df_long_CESD$CESD <- as.numeric(df_long_CESD$CESD)

model_list <- list()
model_BIC <- data.frame(
  Variable = character(),
  Model_Type = character(),
  Classes = integer(),
  BIC = numeric(),
  LogLik = numeric(),
  Entropy = numeric(),
  stringsAsFactors = FALSE
)

if (!file.exists(result_csv)) {
  write.csv(model_BIC, result_csv, row.names = FALSE)
}

## ---- Minimal fix 2: make APPA/OCC/Entropy robust for ng>2 (no hard-coded prob1) ----
calculate_fit_indices <- function(model) {
  if (model$ng == 1) {
    return(list(
      APPA = NA, OCC = NA, entropy = NA,
      class_size = 1, BIC = model$BIC, loglik = model$loglik
    ))
  }

  postprobs <- model$pprob
  prob_cols <- grep("^prob", names(postprobs), value = TRUE)
  if (length(prob_cols) == 0) stop("No probability column found")

  prob_matrix <- as.matrix(postprobs[, prob_cols, drop = FALSE])
  assigned_class <- apply(prob_matrix, 1, which.max)

  # APPA: mean posterior probability among those assigned to class k
  APPA <- numeric(model$ng)
  for (k in 1:model$ng) {
    idx <- which(assigned_class == k)
    APPA[k] <- if (length(idx) > 0) mean(prob_matrix[idx, k]) else NA_real_
  }

  # class size (proportion)
  class_size <- table(assigned_class) / length(assigned_class)

  # OCC
  OCC <- sapply(1:model$ng, function(k) {
    pi_k <- if (as.character(k) %in% names(class_size)) as.numeric(class_size[as.character(k)]) else NA_real_
    appa_k <- APPA[k]
    if (is.na(appa_k) || is.na(pi_k) || pi_k %in% c(0, 1) || appa_k %in% c(0, 1)) return(NA_real_)
    (appa_k / (1 - appa_k)) / (pi_k / (1 - pi_k))
  })

  # normalized entropy based on full posterior matrix (0-1; higher is better)
  eps <- 1e-10
  entropy <- 1 - sum(prob_matrix * log(pmax(prob_matrix, eps))) / (nrow(prob_matrix) * log(model$ng))

  return(list(
    APPA = APPA,
    OCC = OCC,
    entropy = entropy,
    class_size = class_size,
    BIC = model$BIC,
    loglik = model$loglik
  ))
}

var_use <- c("CESD")
var_use_name <- c("CESD")

for (i in 1:length(var_use)) {
  data_model <- df_long_CESD[, c("ID", "Wave", var_use[i], "Age")]
  colnames(data_model) <- c("ID", "Wave", "IU", "Age")
  data_model <- na.omit(data_model)

  cat("🔹 Start fitting the variables:", var_use_name[i], "\n")

  ng_list <- 1:4
  total_models <- length(ng_list)
  start_time_all <- Sys.time()
  elapsed_vec <- c()

  for (idx in seq_along(ng_list)) {
    g <- ng_list[idx]
    model_name <- paste0(var_use_name[i], "_model1_", g)

    ## ---- Minimal fix 3: save models into /models folder ----
    model_file <- file.path("models", paste0(model_name, ".rds"))

    # If it exists, then skip it.
    if (file.exists(model_file)) {
      cat("⏩ Skip the existing model:", model_file, "\n")
      model_list[[model_name]] <- readRDS(model_file)
      next
    }

    # Progress bar
    pct <- round((idx - 1) / total_models * 100, 1)
    cat(sprintf("\n📊 Progress: [%-20s] %3.0f%% | Class number %d/%d | %s\n",
                paste(rep("=", idx - 1), collapse = ""),
                pct, g, total_models, format(Sys.time(), "%H:%M:%S")))

    # Model fitting
    cat("▶️Start fitting the model:", model_name, "\n")
    start_time <- Sys.time()

    model <- tryCatch({
      if (g == 1) {
        hlme(
          fixed = IU ~ Age + I(Age^2),
          subject = "ID",
          random = ~ Age,
          ng = 1,
          data = data_model,
          maxiter = 500
        )
      } else {
        hlme(
          fixed = IU ~ Age + I(Age^2),
          subject = "ID",
          random = ~ Age,
          mixture = ~ Age + I(Age^2),
          ng = g,
          nwg = TRUE,
          data = data_model,
          B = model_list[[paste0(var_use_name[i], "_model1_1")]],
          maxiter = 500
        )
      }
    },
    error = function(e) {
      cat("❌ Model failure:", e$message, "\n")
      write(paste(Sys.time(), "Model failure:", model_name, e$message), file = log_file, append = TRUE)
      return(NULL)
    })

    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
    elapsed_vec <- c(elapsed_vec, elapsed)
    cat(sprintf("✅ Completed %s | Took %.1f minutes| %s\n", model_name, elapsed, format(end_time, "%H:%M:%S")))

    # Estimated time
    avg_time <- mean(elapsed_vec)
    remain <- avg_time * (total_models - idx)
    cat(sprintf("⏱️ On average, each model takes %.1f minutes. The remaining time is expected to be %.1f minutes\n", avg_time, remain))

    if (!is.null(model)) {
      saveRDS(model, model_file)
      model_list[[model_name]] <- model
      fit_indices <- calculate_fit_indices(model)

      sum_row <- data.frame(
        Variable = var_use_name[i],
        Model_Type = "Random_Quadratic_NoCovariate",
        Classes = g,
        BIC = fit_indices$BIC,
        LogLik = fit_indices$loglik,
        Entropy = fit_indices$entropy,
        stringsAsFactors = FALSE
      )

      write.table(sum_row, result_csv, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
    }
  }

  total_elapsed <- as.numeric(difftime(Sys.time(), start_time_all, units = "mins"))
  cat(sprintf("\n🎉 All models have been fitted for %s | Total time consumed: %.1f minutes\n", var_use_name[i], total_elapsed))
}

## ---- Minimal fix 4: read models from /models and write postprobs into /results ----
model_2 <- readRDS(file.path("models", "CESD_model1_2.rds"))
postprobs_2 <- model_2$pprob
summary_2 <- postprobs_2 %>%
  group_by(class) %>%
  summarise(
    n = n(),
    Percent = round(n() / nrow(postprobs_2) * 100, 1),
    APPA = round(mean(get(paste0("prob", first(class))), na.rm = TRUE), 3)
  )
print(summary_2)
write.csv(postprobs_2, file.path("results", "CESD_model1_2.csv"), row.names = TRUE)

model_3 <- readRDS(file.path("models", "CESD_model1_3.rds"))
postprobs_3 <- model_3$pprob
summary_3 <- postprobs_3 %>%
  group_by(class) %>%
  summarise(
    n = n(),
    Percent = round(n() / nrow(postprobs_3) * 100, 1),
    APPA = round(mean(get(paste0("prob", first(class))), na.rm = TRUE), 3)
  )
print(summary_3)
write.csv(postprobs_3, file.path("results", "CESD_model1_3.csv"), row.names = TRUE)

model_4 <- readRDS(file.path("models", "CESD_model1_4.rds"))
postprobs_4 <- model_4$pprob
summary_4 <- postprobs_4 %>%
  group_by(class) %>%
  summarise(
    n = n(),
    Percent = round(n() / nrow(postprobs_4) * 100, 1),
    APPA = round(mean(get(paste0("prob", first(class))), na.rm = TRUE), 3)
  )
print(summary_4)
write.csv(postprobs_4, file.path("results", "CESD_model1_4.csv"), row.names = TRUE)

calculate_comprehensive_fit_indices <- function(model) {
  if (is.null(model)) {
    return(list(
      APPA = NA, OCC = NA, entropy = NA, relative_entropy = NA,
      class_size = NA, BIC = NA, loglik = NA, n_classes = NA
    ))
  }

  if (model$ng == 1) {
    return(list(
      APPA = NA, OCC = NA, entropy = NA, relative_entropy = NA,
      class_size = 1, BIC = model$BIC, loglik = model$loglik, n_classes = 1
    ))
  }

  # Obtain the posterior probability
  postprobs <- model$pprob

  # Debugging information: View the posterior probability structure
  cat("Posterior probability data structure:\n")
  print(head(postprobs))
  cat("Column name:", names(postprobs), "\n")
  cat("Number of Categories:", model$ng, "\n")

  # Correctly extract the probability columns
  prob_cols <- grep("^prob", names(postprobs), value = TRUE)
  cat("Probability columns:", prob_cols, "\n")

  if (length(prob_cols) == 0) {
    stop("No probability column found")
  }

  # Calculate the maximum posterior probability and the corresponding category for each individual
  prob_matrix <- as.matrix(postprobs[, prob_cols, drop = FALSE])
  max_prob <- apply(prob_matrix, 1, max)
  assigned_class <- apply(prob_matrix, 1, which.max)

  # Calculate APPA - Correction Method
  APPA <- numeric(model$ng)
  for (k in 1:model$ng) {
    class_members <- which(assigned_class == k)
    if (length(class_members) > 0) {
      APPA[k] <- mean(prob_matrix[class_members, k])
    } else {
      APPA[k] <- NA
    }
  }

  # Calculate the proportion of categories class_size
  class_size <- table(assigned_class) / nrow(postprobs)

  # Calculate OCC - Add stricter checks
  OCC <- numeric(model$ng)
  for (k in 1:model$ng) {
    if (as.character(k) %in% names(class_size)) {
      pi_k <- class_size[as.character(k)]
      appa_k <- APPA[k]

      # Add boundary checks
      if (is.na(appa_k) || is.na(pi_k) || pi_k == 0 || pi_k == 1 || appa_k == 1 || appa_k == 0) {
        OCC[k] <- NA
      } else {
        OCC[k] <- (appa_k / (1 - appa_k)) / (pi_k / (1 - pi_k))
      }
    } else {
      OCC[k] <- NA
    }
  }

  # Calculate entropy - using the maximum a posteriori probability
  entropy <- -sum(max_prob * log(pmax(max_prob, 1e-10))) / nrow(postprobs)

  # Calculate Relative Entropy
  relative_entropy <- 1 - (entropy / log(model$ng))

  return(list(
    APPA = APPA,
    OCC = OCC,
    entropy = entropy,
    relative_entropy = relative_entropy,
    class_size = as.numeric(class_size),
    BIC = model$BIC,
    loglik = model$loglik,
    n_classes = model$ng,
    assigned_class = assigned_class,
    max_prob = max_prob
  ))
}

for (i in 1:length(var_use)) {
  for (g in 1:4) {
    model_name <- paste0(var_use_name[i], "_model1_", g)
    model_file <- file.path("models", paste0(model_name, ".rds"))

    if (file.exists(model_file)) {
      cat("\n📊 Model:", model_name, "\n")
      cat("--------------------------------------------\n")

      model <- readRDS(model_file)

      # Add error handling
      result <- tryCatch({
        fit_indices <- calculate_comprehensive_fit_indices(model)

        if (!is.null(model)) {
          # Output basic information
          cat(sprintf("Number of categories: %d\n", fit_indices$n_classes))
          cat(sprintf("BIC: %.1f\n", fit_indices$BIC))
          cat(sprintf("Log-likelihood value: %.1f\n", fit_indices$loglik))

          if (fit_indices$n_classes > 1) {
            cat(sprintf("Entropy: %.3f\n", fit_indices$entropy))
            cat(sprintf("Relative Entropy: %.3f\n", fit_indices$relative_entropy))
            cat("\nDetailed information for each category:\n")

            for (class_idx in 1:fit_indices$n_classes) {
              cat(sprintf("  Category %d:\n", class_idx))
              cat(sprintf("    - Category proportion: %.3f\n", fit_indices$class_size[class_idx]))
              cat(sprintf("    - APPA: %.3f\n", fit_indices$APPA[class_idx]))
              cat(sprintf("    - OCC: %.3f\n", fit_indices$OCC[class_idx]))
            }

            # Output Summary Indicators
            valid_APPA <- fit_indices$APPA[!is.na(fit_indices$APPA)]
            valid_OCC <- fit_indices$OCC[!is.na(fit_indices$OCC)]

            if (length(valid_APPA) > 0) {
              cat(sprintf("\nSummary Indicators:\n"))
              cat(sprintf("  mean APPA: %.3f\n", mean(valid_APPA)))
              cat(sprintf("  min APPA: %.3f\n", min(valid_APPA)))
            }

            if (length(valid_OCC) > 0) {
              cat(sprintf("  mean OCC: %.3f\n", mean(valid_OCC, na.rm = TRUE)))
              cat(sprintf("  min OCC: %.3f\n", min(valid_OCC, na.rm = TRUE)))
            }

            # Quality Assessment
            cat(sprintf("\nClassification Quality Assessment:\n"))
            if (fit_indices$relative_entropy > 0.8) {
              cat("  ✓ Relative Entropy > 0.8: Excellent classification quality\n")
            } else if (fit_indices$relative_entropy > 0.6) {
              cat("  ~ Relative Entropy 0.6 - 0.8: Moderate classification quality\n")
            } else {
              cat("  ⚠ Relative entropy < 0.6: Poor classification quality\n")
            }

            if (length(valid_APPA) > 0 && min(valid_APPA) > 0.7) {
              cat("  ✓ All APPA values > 0.7: Good category distinction\n")
            } else if (length(valid_APPA) > 0 && min(valid_APPA) > 0.5) {
              cat("  ~ Minimum APPA value between 0.5 and 0.7: Moderate category distinction\n")
            } else {
              cat("  ⚠ Minimum APPA < 0.5: Poor category discrimination\n")
            }
          } else {
            cat("Single-class model, without classification quality indicators\n")
          }
        }
      }, error = function(e) {
        cat("❌ Calculation error:", e$message, "\n")

        if (!is.null(model)) {
          cat(sprintf("Number of categories: %d\n", model$ng))
          cat(sprintf("BIC: %.1f\n", model$BIC))
          cat(sprintf("对数似然值: %.1f\n", model$loglik))
        }
      })

      cat("--------------------------------------------\n")
    }
  }
}
