# ------------------------------------------------------------
# File: utility_functions.R
# Description: This file contains a collection of utility functions
#              that are commonly used across various scripts.
#
# Author: Thiago Gonçalves-Souza
# Date: September-30-2024
# R version 4.4.1
#
# Functions included:
#   - calc_partial_sd: Compute Partial Standard Deviations Considering Correlation
#   - standardize_estimates: Function to standardize estimates
#   - collect_results: Function to collect results for one diversity method
#   - add_estimates: function to add the calculated estimates to a data.frame
#   - model_averaging: Model Averaging Using AIC Weights
#   - model_average_scaled_importance: Function to model-average scaled importance
#   - extract_specific_stats: function to extract relevant stats from a glmmTMB object
#   - lm_results: function to process multiple lm models and extract the relevant terms
#   - glmm_results: function to process multiple glmmm models and extract the relevant terms
#   - extract_model_results: Function to extract results from meta-analysis models 
#   - overall_rma_results: Function to extract relevant results from the model summary of the rma function
#   - extract_glmmTMB_summary: extract the glmmTMB results from multiple models # it's a bit different from the previous one
#   - extract_confint: Extract confidence intervals from lm models in the moderators test
#
# Reference: Cade, B.S. (2015). Ecology 96, https://doi.org/10.1890/14-1639.1	
#
# ------------------------------------------------------------


### Compute Partial Standard Deviations Considering Correlation

calc_partial_sd <- function(model, predictors) {
  vcov_mat <- summary(model)$vcov$cond
  partial_sds <- sapply(predictors, function(predictor) {
    if (predictor %in% rownames(vcov_mat)) {
      sqrt(vcov_mat[predictor, predictor])
    } else {
      NA
    }
  })
  return(partial_sds)
}

# Function to standardize estimates
standardize_estimates <- function(model, predictors, partial_sds) {
  coefs <- summary(model)$coefficients$cond
  standardized_estimates <- sapply(1:length(predictors), function(i) {
    predictor <- predictors[i]
    if (!is.na(partial_sds[i]) && predictor %in% rownames(coefs)) {
      estimate <- coefs[predictor, "Estimate"]
      return(estimate * partial_sds[i])
    } else {
      return(NA)
    }
  })
  return(standardized_estimates)
}

# Function to collect results for one diversity method
collect_results <- function(mod1, mod2, mod3, predictors, partial_sds1, partial_sds2, partial_sds3, method_name) {
  results <- data.frame(
    Model = character(),
    Predictor = character(),
    Estimate = numeric(),
    Std_Error = numeric(),
    Method = character()
  )
  
  add_estimates <- function(model, predictors, partial_sds, model_name) {
    coefs <- summary(model)$coefficients$cond
    for (i in 1:length(predictors)) {
      predictor <- predictors[i]
      if (!is.na(partial_sds[i]) && predictor %in% rownames(coefs)) {
        standardized_estimate <- coefs[predictor, "Estimate"] * partial_sds[i]
        std_error <- coefs[predictor, "Std. Error"]
        results <<- rbind(
          results,
          data.frame(
            Model = model_name,
            Predictor = predictor,
            Estimate = standardized_estimate,
            Std_Error = std_error,
            Method = method_name
          )
        )
      }
    }
  }

  # Collect estimates for each model
  add_estimates(mod1, predictors, partial_sds1, "mod1")
  add_estimates(mod2, predictors, partial_sds2, "mod2")
  add_estimates(mod3, predictors, partial_sds3, "mod3")
  
  return(results)
}



# Function to model-average scaled importance
model_average_scaled_importance <- function(scaled_importance, weights) {
  weighted_avg <- sum(scaled_importance * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
  return(weighted_avg)
}


# Calculate t-Statistics and Scale Ratios

calculate_t_stats <- function(results) {
  # Calculate t-statistics
  results$t_stat <- with(results, Estimate / Std_Error)
  
  # Calculate the ratio of absolute values of t-statistics for each predictor
  results$abs_t_stat <- abs(results$t_stat)
  
  # Scale ratios relative to the maximum within each model
  results$scaled_importance <- ave(results$abs_t_stat, results$Model, FUN = function(x) x / max(x, na.rm = TRUE))
  
  return(results)
}



# Model Averaging Using AIC Weights

model_averaging <- function(results, model_list) {
  # Calculate AIC weights
  aic_weights <- aictab(model_list, modnames = c("mod1", "mod2", "mod3"))
  
  # Extract weights
  weights <- aic_weights$AICcWt
  
  # Adjust dataframe to match weights
  results$model_weight <- rep(NA, nrow(results))
  for (i in 1:nrow(results)) {
    results$model_weight[i] <- weights[which(aic_weights$Modnames == results$Model[i])]
  }
  
  # Model-average scaled importance
  model_average_scaled_importance <- function(scaled_importance, weights) {
    weighted_avg <- sum(scaled_importance * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
    return(weighted_avg)
  }
  
  # Weighted average of the importance ratios across models
  weighted_importance_patch_type <- model_average_scaled_importance(results$scaled_importance[results$Predictor == "patch_typefragmented"], results$model_weight)
  weighted_importance_habitat_amount <- model_average_scaled_importance(results$scaled_importance[results$Predictor == "habitat_amount_cent"], results$model_weight)
  
  # Combine and view the weighted importance ratios
  weighted_importance <- data.frame(
    Predictor = c("patch_typefragmented", "habitat_amount_cent"),
    Scaled_Importance = c(weighted_importance_patch_type, weighted_importance_habitat_amount)
  )
  
  return(weighted_importance)
}



# Function to extract relevant stats from a specific row (fixed effect) in a glmmTMB object

extract_specific_stats <- function(model, included_moderator, row_index = 2) {
  summary_model <- summary(model)
  
  # Extract the fixed effects coefficients
  coef_table <- summary_model$coefficients$cond
  
  # Select the specified row (e.g., the second fixed effect)
  coef_stats <- coef_table[row_index, ]
  
  # Create a data frame with the results
  df <- data.frame(
    Included_moderator = included_moderator,
    Estimate = coef_stats["Estimate"],
    `Std. Error` = coef_stats["Std. Error"],
    `z value` = coef_stats["z value"],
    `Pr(>|z|)` = coef_stats["Pr(>|z|)"]
  )
  
  # Ensure the column names are readable
  colnames(df) <- c("Included_moderator", "Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  return(df)
}



# Function to process multiple models
glmm_results <- function(..., row_index = 2) {
  models <- list(...)
  
  # Extract the diversity types (assuming the naming convention of mod1, mod2, etc.)
  diversity_types <- sapply(substitute(list(...))[-1], deparse)
  
  # Create an empty list to store results
  results_list <- list()
  
  # Loop through each model and extract relevant statistics
  for (i in seq_along(models)) {
    results_list[[i]] <- extract_specific_stats(models[[i]], diversity_types[i], row_index)
  }
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}

  # Function to extract the glmmTMB results from multiple models 

extract_glmmTMB_summary <- function(..., model_names) {
  models <- list(...)
  if (length(models) != length(model_names)) {
    stop("The number of models and model names must be the same.")
  }
  
  results <- list()
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    model_name <- model_names[i]
    
    # Extract summary
    model_summary <- summary(model)
    
    # Extract intercept and coefficients
    coefs <- model_summary$coefficients$cond
    
    # Convert to data frame
    coefs_df <- as.data.frame(coefs)
    
    # Create a new column for the model name
    coefs_df$model <- model_name
    
    # Create a new column for the term names
    coefs_df$term <- rownames(coefs_df)
    rownames(coefs_df) <- NULL
    
    # Create a sequence column
    coefs_df$id <- seq_len(nrow(coefs_df))
    
    # Reorder columns to have 'id', 'model', 'term' as the first columns
    coefs_df <- coefs_df[, c("id", "model", "term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    
    results[[i]] <- coefs_df
  }
  
  # Combine all results into a single data frame
  output_df <- do.call(rbind, results)
  
  return(output_df)
}


  # Function to extract results from meta-analysis models


extract_model_results <- function(model, diversity_index, diversity_type) {
  # Check if the object is a class of rma from metafor package
  if (!inherits(model, "rma")) {
    stop("The provided model is not an 'rma' object from the metafor package.")
  }
  
  # Extracting the summary of the model
  model_summary <- summary(model)
  
  # Extract the coefficients (variables)
  variable_names <- rownames(model_summary$beta)
  
  # Convert the coefficients matrix to a data frame
  coefficients_df <- data.frame(
    variable = variable_names,
    estimate = model_summary$beta,
    se = model_summary$se,
    zval = model_summary$zval,
    pval = model_summary$pval,
    ci.lb = model_summary$ci.lb,
    ci.ub = model_summary$ci.ub
  )
  
  # Adding the columns for diversity index and type at the beginning
  results_df <- data.frame(
    diversity_index = diversity_index,
    diversity_type = diversity_type,
    coefficients_df,
    row.names = NULL
  )
  
  return(results_df)
}


# Function to extract relevant results from the model summary of the rma function

overall_rma_results <- function(model) {
  # Check if model is a valid rma object
  if (!inherits(model, "rma")) {
    stop("The model is not a valid 'rma' object.")
  }
  
  # Extract the summary statistics
  summary_model <- summary(model)
  
  # Extract coefficients
  coefficients <- summary_model$b
  
  # Extract standard errors
  se <- summary_model$se
  
  # Extract z-values
  zval <- summary_model$zval
  
  # Extract p-values
  pval <- summary_model$pval

  # Combine into a data frame
  results <- data.frame(
    Estimate = coefficients,
    SE = se,
    Z.Value = zval,
    P.Value = pval
  )

  # Rename the row name "intrcpt" to "overall_model"
  rownames(results) <- "overall_model"
  
  return(results)
}


# Function to extract relevant results from the model summary of the lm function using multiple models


lm_results <- function(..., row_index = 2) {
  models <- list(...)
  
  # Extract the diversity types (assuming the naming convention of mod1, mod2, etc.)
  included_moderators <- sapply(substitute(list(...))[-1], deparse)
  
  # Create an empty list to store results
  results_list <- list()
  
  # Function to extract specific statistics from lm model
  extract_specific_stats <- function(model, included_moderator, row_index) {
    summary_model <- summary(model)
    
    # Extract specific results
    estimate <- summary_model$coefficients[row_index, "Estimate"]
    std_error <- summary_model$coefficients[row_index, "Std. Error"]
    t_value <- summary_model$coefficients[row_index, "t value"]
    p_value <- summary_model$coefficients[row_index, "Pr(>|t|)"]
    
    # Combine into a data frame
    data.frame(
      Included_moderator = included_moderator,
      Estimate = estimate,
      Std_Error = std_error,
      T_value = t_value,
      P_value = p_value
    )
  }
  
  # Loop through each model and extract relevant statistics
  for (i in seq_along(models)) {
    results_list[[i]] <- extract_specific_stats(models[[i]], included_moderators[i], row_index)
  }
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}

### Extract confidence intervals from lm models in the moderators test

extract_confint <- function(..., row_index = 2) {
  models <- list(...)
  
  # Extract the diversity types (assuming the naming convention of mod1, mod2, etc.)
  included_moderators <- sapply(substitute(list(...))[-1], deparse)
  
  # Create an empty list to store results
  results_list <- list()
  
  # Function to extract confidence intervals from lm model
  extract_confidence_intervals <- function(model, included_moderator, row_index) {
    conf_intervals <- confint(model)
    
    # Extract the specific confidence intervals, excluding intercept
    lower_ci <- conf_intervals[row_index, 1]
    upper_ci <- conf_intervals[row_index, 2]
    
    # Combine into a data frame
    data.frame(
      Included_moderator = included_moderator,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci
    )
  }
  
  # Loop through each model and extract relevant confidence intervals
  for (i in seq_along(models)) {
    results_list[[i]] <- extract_confidence_intervals(models[[i]], included_moderators[i], row_index)
  }
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}


