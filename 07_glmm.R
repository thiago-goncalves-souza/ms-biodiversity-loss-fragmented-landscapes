# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: GLMMs to test the effects of landscape type on species diversity
# Author: Thiago Gonçalves Souza
# Date: [October 10, 2024]
# Notes: R version 4.4.1 [2024-06-14 ucrt]
# ==============================================================================


# Using renv for package management
renv::restore() # it requires installing the package renv

# Packages ---------------------------------------------------------------------

library(glmmTMB)
library(cowplot)
library(DHARMa)
library(tidyverse)
library(ggplot2)
library(lme4)
library(svglite)
library(gtools)
source("utility_functions.R")

### Data -----------------------------------------------------------------------

# continents 

study_continents <- read.csv("data/landfrag_37sset_continents.csv")

# time since fragmentation 

time_since <- read.csv("data/landfrag_37sset_time_since_frag.csv") %>% 
  dplyr::arrange(refshort)

# Diversity of 2 (pair = all pairs)
frag_cont_div2 <- read.csv("processed_data/diversity_of_2.csv") 

# Diversity of 2 (pair = closest fragments)

frag_cont_div2_cl <- read.csv("processed_data/diversity_of_2_close.csv") 

# Diversity of 2 (pair = all pairs) rarefied estimates (q = 0, q = 2)

frag_cont_rar_div2 <- read.csv("processed_data/frag_cont_raref_2.csv") 

# Diversity of 2 (pair = closest fragments) rarefied estimates (q = 0, q = 2)

frag_cont_rar_div2_cl <- read.csv("processed_data/frag_cont_raref_close_2.csv") 


### Data preparation -----------------------------------------------------------


# habitat amount class by study
frag_cont_div2 %>% 
  dplyr::filter(diversity_index == "alpha") %>% 
  dplyr::select(refshort, habitat_amount) %>% 
  dplyr::group_by(refshort) %>% 
  dplyr::summarise(aver_ha = mean(habitat_amount)) %>% 
  dplyr::mutate(ha_class = quantcut(aver_ha, 2)) %>% 
  dplyr::select(refshort, ha_class) -> land_hab_amount_class


### (summary tables) averages per diversity and landscape type  ----------------


frag_cont_div2 %>% 
  dplyr::select(patch_type, diversity_index, diversity_value) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean = mean(diversity_value),
                   .groups = "keep") %>% 
  pivot_wider(names_from = patch_type, values_from = mean) %>%
  dplyr::mutate(percentage_diff = ((fragmented - continuous) / continuous) * 100,
         diversity_type = "All pairs") %>%
  dplyr::select(diversity_type, diversity_index, percentage_diff) -> all_pairs_perc_diff

frag_cont_rar_div2 %>% 
  dplyr::select(patch_type, diversity_index, diversity_value) %>% 
  dplyr::group_by(patch_type, diversity_index) %>% 
  summarise(mean = mean(diversity_value),
            .groups = "keep") %>% 
  pivot_wider(names_from = patch_type, values_from = mean) %>%
  dplyr::mutate(percentage_diff = ((fragmented - continuous) / continuous) * 100,
         diversity_type = "Nearest pairs") %>%
  dplyr::select(diversity_type, diversity_index, percentage_diff) -> near_pairs_perc_diff

frag_cont_rar_div2 %>% 
  dplyr::filter(q_order == "q = 0") %>% 
  dplyr::select(patch_type, diversity_index, diversity_value) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean = mean(diversity_value),
                   .groups = "keep") %>% 
  pivot_wider(names_from = patch_type, values_from = mean) %>%
  dplyr::mutate(percentage_diff = ((fragmented - continuous) / continuous) * 100,
         diversity_type = "All pairs (q = 0)") %>%
  dplyr::select(diversity_type, diversity_index, percentage_diff) -> all_pairs_q0_perc_diff

frag_cont_rar_div2 %>% 
  dplyr::filter(q_order == "q = 2") %>% 
  dplyr::select(patch_type, diversity_index, diversity_value) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean = mean(diversity_value),
                   .groups = "keep") %>% 
  pivot_wider(names_from = patch_type, values_from = mean) %>%
  dplyr::mutate(percentage_diff = ((fragmented - continuous) / continuous) * 100,
         diversity_type = "All pairs (q = 2)") %>%
  dplyr::select(diversity_type, diversity_index, percentage_diff) -> all_pairs_q2_perc_diff

frag_cont_rar_div2_cl %>% 
  dplyr::filter(q_order == "q = 0") %>% 
  dplyr::select(patch_type, diversity_index, diversity_value) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean = mean(diversity_value),
                   .groups = "keep") %>% 
  pivot_wider(names_from = patch_type, values_from = mean) %>%
  dplyr::mutate(percentage_diff = ((fragmented - continuous) / continuous) * 100,
         diversity_type = "Nearest pairs (q = 0)") %>%
  dplyr::select(diversity_type, diversity_index, percentage_diff) -> near_pairs_q0_perc_diff


frag_cont_rar_div2_cl %>% 
  dplyr::filter(q_order == "q = 2") %>% 
  dplyr::select(patch_type, diversity_index, diversity_value) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean = mean(diversity_value),
                   .groups = "keep") %>% 
  pivot_wider(names_from = patch_type, values_from = mean) %>%
  dplyr::mutate(percentage_diff = ((fragmented - continuous) / continuous) * 100,
         diversity_type = "Nearest pairs (q = 2)") %>%
  dplyr::select(diversity_type, diversity_index, percentage_diff) -> near_pairs_q2_perc_diff


perc_table_results <- rbind.data.frame(
  all_pairs_perc_diff,
  all_pairs_q0_perc_diff,
  all_pairs_q2_perc_diff,
  near_pairs_perc_diff,
  near_pairs_q0_perc_diff,
  near_pairs_q2_perc_diff
)

perc_table_results %>% 
  dplyr::select(diversity_index, percentage_diff) %>% 
  group_by(diversity_index) %>% 
  dplyr::summarise(min = min(percentage_diff), # note they are reversed because some are negative!
                   max = max(percentage_diff),
                   mean_diff = mean(percentage_diff),
                   .groups = "keep")


### GLMMs ----------------------------------------------------------------------
### 1. Fitting models using all plot pairs ###


## Effects on Alpha-diversity -----------------------------------------------------


## Data Preparation ------------------------------------------------------------

# Filter Alpha-diversity data

frag_cont_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "alpha") -> frag_cont_div2_alpha_glmm


# Prepare moderators ofAlpha-diversity 

frag_cont_div2_alpha_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_alpha

## Model Fitting ---------------------------------------------------------------

mod1_div2_alpha <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                           data = frag_cont_div2_alpha_glmm)
summary(mod1_div2_alpha)

## Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_alpha_final_res <- simulateResiduals(fittedModel = mod1_div2_alpha, plot = F)
plot(mod_div2_alpha_final_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_alpha_final_res) # no overdispersion

## Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response

mod1_div2_alpha_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                               data = frag_cont_div2_alpha_glmm)

# Check the transformed model

mod_div2_alpha_log_res <- simulateResiduals(fittedModel = mod1_div2_alpha_log, plot = F)
plot(mod_div2_alpha_log_res) # residuals normally distributed
testDispersion(mod_div2_alpha_final_res) # no overdispersion

## Final Model Fitting ---------------------------------------------------------
# Refit final model with REML

mod1_div2_alpha_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                     REML = TRUE,
                                     data = frag_cont_div2_alpha_glmm)

summary(mod1_div2_alpha_log_final) # report results


# Fit models with habitat amount as a random variable

mod1_div2_alpha_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                      control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                     REML = TRUE,
                                     data = frag_cont_div2_alpha_glmm)

summary(mod1_div2_alpha_ha_final) # report results

## Extracting Slopes -----------------------------------------------------------

mod1_div2_alpha_2_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                data = frag_cont_div2_alpha_glmm)

ha_slope_div2_alpha <- ranef(mod1_div2_alpha_2_ha)$cond$refshort$habitat_amount


## Effects of moderators -------------------------------------------------------

# Habitat amount class 

mod_hac_div2_alpha <- lm(ha_slope_div2_alpha ~ ha_class, data = moders_div2_alpha)
mod_hac_div2_alpha_res <- simulateResiduals(fittedModel = mod_hac_div2_alpha, plot = F)
plot(mod_hac_div2_alpha_res) # residuals normally distributed
testDispersion(mod_hac_div2_alpha_res) # no overdispersion

summary(mod_hac_div2_alpha) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_alpha$effort_diff_cent <- scale(log(moders_div2_alpha$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_alpha <- lm(ha_slope_div2_alpha ~ effort_diff_cent, data = moders_div2_alpha)

mod_eff_div2_alpha_res <- simulateResiduals(fittedModel = mod_eff_div2_alpha, plot = F)
plot(mod_eff_div2_alpha_res) # residuals normally distributedm with heterogeneous variances specially in large values (huge sampling effort in one study)
testDispersion(mod_eff_div2_alpha_res) # no overdispersion


summary(mod_eff_div2_alpha) # report results # 


# time since fragmentation

mod_time_div2_alpha <- lm(ha_slope_div2_alpha ~ time_since_fragmentation, data = moders_div2_alpha)
mod_time_div2_alpha_res <- simulateResiduals(fittedModel = mod_time_div2_alpha, plot = F)
plot(mod_time_div2_alpha_res) # residuals normally distributed
testDispersion(mod_hac_div2_alpha_res) # no overdispersion

summary(mod_time_div2_alpha) # report results

# continent 

mod_cont_div2_alpha <- lm(ha_slope_div2_alpha ~ continent_cat, data = moders_div2_alpha)
mod_cont_div2_alpha_res <- simulateResiduals(fittedModel = mod_cont_div2_alpha, plot = F)
plot(mod_cont_div2_alpha_res) # residuals normally distributed
testDispersion(mod_cont_div2_alpha_res) # no overdispersion

summary(mod_cont_div2_alpha)  # report results


## save the results to export later

lm_results(mod_hac_div2_alpha, mod_eff_div2_alpha, mod_time_div2_alpha, mod_cont_div2_alpha) %>% 
  dplyr::mutate(Diversity_type = "Alpha (all pairs)", .before = "Included_moderator") -> moder_results_div2_alpha


## Effects on Beta-diversity ---------------------------------------------------

## Data Preparation ------------------------------------------------------------

frag_cont_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "beta") -> frag_cont_div2_beta_glmm


## Prepare moderators for Beta-diversity

frag_cont_div2_beta_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_beta


## Model fitting ---------------------------------------------------------------
# Fit initial model

mod1_div2_beta <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                          data = frag_cont_div2_beta_glmm)
summary(mod1_div2_beta)

# Model checking
# Simulate residuals and check assumptions

mod_div2_beta_final_res <- simulateResiduals(fittedModel = mod1_div2_beta, plot = F)
plot(mod_div2_beta_final_res) # residuals normally distributed 
testDispersion(mod_div2_beta_final_res) # no overdispersion


# Final Model Fitting
# Refit final model with REML

mod1_div2_beta_final <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                REML = TRUE,
                                data = frag_cont_div2_beta_glmm)
summary(mod1_div2_beta_final) # report results


# Fit models with habitat amount as a random variable

mod1_div2_beta_ha_final <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                             control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                             REML = TRUE,
                             data = frag_cont_div2_beta_glmm)

summary(mod1_div2_beta_ha_final) # report results

## Extracting slopes ----------------------------------------------------------

mod1_div2_beta_ha <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                             control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                             data = frag_cont_div2_beta_glmm)

ha_slope_div2_beta <- ranef(mod1_div2_beta_ha)$cond$refshort$habitat_amount

## Effects of Moderators -------------------------------------------------------

# Habitat amount class 

mod_hac_div2_beta <- lm(ha_slope_div2_beta ~ ha_class, data = moders_div2_beta)
mod_hac_div2_beta_res <- simulateResiduals(fittedModel = mod_hac_div2_beta, plot = F)
plot(mod_hac_div2_beta_res) # residuals normally distributed
testDispersion(mod_hac_div2_beta_res) # no overdispersion

summary(mod_hac_div2_beta) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_beta$effort_diff_cent <- scale(log(moders_div2_beta$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_beta <- lm(ha_slope_div2_beta ~ effort_diff_cent, data = moders_div2_beta)

mod_eff_div2_beta_res <- simulateResiduals(fittedModel = mod_eff_div2_beta, plot = F)
plot(mod_eff_div2_beta_res) # residuals normally distributed
testDispersion(mod_eff_div2_beta_res) # no overdispersion


summary(mod_eff_div2_beta) # report results # 


# time since fragmentation

mod_time_div2_beta <- lm(ha_slope_div2_beta ~ time_since_fragmentation, data = moders_div2_beta)
mod_time_div2_beta_res <- simulateResiduals(fittedModel = mod_time_div2_beta, plot = F)
plot(mod_time_div2_beta_res) # residuals normally distributed
testDispersion(mod_hac_div2_beta_res) # no overdispersion

summary(mod_time_div2_beta) # report results

# continent 

mod_cont_div2_beta <- lm(ha_slope_div2_beta ~ continent_cat, data = moders_div2_beta)
mod_cont_div2_beta_res <- simulateResiduals(fittedModel = mod_cont_div2_beta, plot = F)
plot(mod_cont_div2_beta_res) # residuals normally distributed
testDispersion(mod_cont_div2_beta_res) # no overdispersion

summary(mod_cont_div2_beta)  # report results


## save the results to export 

lm_results(mod_hac_div2_beta, mod_eff_div2_beta, mod_time_div2_beta, mod_cont_div2_beta) %>% 
  dplyr::mutate(Diversity_type = "Beta (all pairs)", .before = "Included_moderator") -> moder_results_div2_beta


## Effects on Gamma-Diversity  ----------------------------------------------------

## Data Preparation ------------------------------------------------------------
# Filter Gamma-diversity data

frag_cont_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "gamma") -> frag_cont_div2_gamma_glmm


# Prepare moderators for Gamma-diversity

frag_cont_div2_gamma_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_gamma

## Model Fitting ---------------------------------------------------------------
# Fit initial model

mod1_div2_gamma <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                           data = frag_cont_div2_gamma_glmm)
summary(mod1_div2_gamma)

## Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_gamma_final_res <- simulateResiduals(fittedModel = mod1_div2_gamma, plot = F)
plot(mod_div2_gamma_final_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_gamma_final_res) # no overdispersion

## Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response

mod1_div2_gamma_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                               data = frag_cont_div2_gamma_glmm)
summary(mod1_div2_gamma_log)

# Check the transformed model
mod_div2_gamma_log_res <- simulateResiduals(fittedModel = mod1_div2_gamma_log, plot = F)
plot(mod_div2_gamma_log_res) # residuals normally distributed
testDispersion(mod_div2_gamma_log_res) # no overdispersion

## Final Model Fitting ---------------------------------------------------------
# Refit final model with REML

mod1_div2_gamma_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                     REML = TRUE,
                                     data = frag_cont_div2_gamma_glmm)
summary(mod1_div2_gamma_log_final) # report results

# Fit models with habitat amount as a random variable

mod1_div2_gamma_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                              control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                              REML = TRUE,
                              data = frag_cont_div2_gamma_glmm)

summary(mod1_div2_gamma_ha_final) # report results
 
## Extracting Slopes -----------------------------------------------------------

mod1_div2_gamma_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                              control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                              data = frag_cont_div2_gamma_glmm)
ha_slope_div2_gamma <- ranef(mod1_div2_gamma_ha)$cond$refshort$habitat_amount

## Effects of Moderators -------------------------------------------------------

# Habitat amount class 

mod_hac_div2_gamma <- lm(ha_slope_div2_gamma ~ ha_class, data = moders_div2_gamma)
mod_hac_div2_gamma_res <- simulateResiduals(fittedModel = mod_hac_div2_gamma, plot = F)
plot(mod_hac_div2_gamma_res) # residuals normally distributed
testDispersion(mod_hac_div2_gamma_res) # no overdispersion

summary(mod_hac_div2_gamma) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_gamma$effort_diff_cent <- scale(log(moders_div2_gamma$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_gamma <- lm(ha_slope_div2_gamma ~ effort_diff_cent, data = moders_div2_gamma)

mod_eff_div2_gamma_res <- simulateResiduals(fittedModel = mod_eff_div2_gamma, plot = F)
plot(mod_eff_div2_gamma_res) # residuals normally distributed
testDispersion(mod_eff_div2_gamma_res) # no overdispersion


summary(mod_eff_div2_gamma) # report results


# time since fragmentation

mod_time_div2_gamma <- lm(ha_slope_div2_gamma ~ time_since_fragmentation, data = moders_div2_gamma)
mod_time_div2_gamma_res <- simulateResiduals(fittedModel = mod_time_div2_gamma, plot = F)
plot(mod_time_div2_gamma_res) # residuals normally distributed
testDispersion(mod_hac_div2_gamma_res) # no overdispersion

summary(mod_time_div2_gamma) # report results

# continent 

mod_cont_div2_gamma <- lm(ha_slope_div2_gamma ~ continent_cat, data = moders_div2_gamma)
mod_cont_div2_gamma_res <- simulateResiduals(fittedModel = mod_cont_div2_gamma, plot = F)
plot(mod_cont_div2_gamma_res) # residuals normally distributed
testDispersion(mod_cont_div2_gamma_res) # no overdispersion

summary(mod_cont_div2_gamma)  # report results


## save the results to export later

lm_results(mod_hac_div2_gamma, mod_eff_div2_gamma, mod_time_div2_gamma, mod_cont_div2_gamma) %>% 
  dplyr::mutate(Diversity_type = "Gamma (all pairs)", .before = "Included_moderator") -> moder_results_div2_gamma


### GLMMs ----------------------------------------------------------------------
### 2. Fitting models using only the nearest pairs

## Effects on Alpha-diversity

## Data Preparation ------------------------------------------------------------
# Filter Alpha-diversity data

frag_cont_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "alpha") -> div_cl2_alpha_glmm

# Prepare moderators for Alpha-diversity

div_cl2_alpha_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_alpha


## Model Fitting ---------------------------------------------------------------
# Fit initial model

mod1_div2_cl_alpha_2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                data = div_cl2_alpha_glmm)
summary(mod1_div2_cl_alpha_2)

## Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_cl_alpha_final_res <- simulateResiduals(fittedModel = mod1_div2_cl_alpha_2, plot = F)
plot(mod_div2_cl_alpha_final_res) # residuals not normally distributed - log it
testDispersion(mod_div2_cl_alpha_final_res) # no overdispersion

## Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response

mod1_div2_cl_alpha_2_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                    data = div_cl2_alpha_glmm)

# Check the transformed model

mod_div2_cl_alpha_log_res <- simulateResiduals(fittedModel = mod1_div2_cl_alpha_2_log, plot = F)
plot(mod_div2_cl_alpha_log_res) # residuals normally distributed
testDispersion(mod_div2_cl_alpha_final_res) # no overdispersion

## Final Model Fitting ---------------------------------------------------------
# Refit final model with REML

mod1_div2_cl_alpha_2_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                          REML = TRUE,
                                          data = div_cl2_alpha_glmm)
summary(mod1_div2_cl_alpha_2_log_final) # report results


# Fit models with habitat amount as a random variable
mod1_div2_cl_alpha_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                   REML = TRUE,
                                   data = div_cl2_alpha_glmm)

summary(mod1_div2_cl_alpha_ha_final) # report results


## Extracting Slopes -----------------------------------------------------------

mod1_div2_cl_alpha_2_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                   data = div_cl2_alpha_glmm)

ha_slope_div2_cl_alpha <- ranef(mod1_div2_cl_alpha_2_ha)$cond$refshort$habitat_amount

## Effects of Moderators -------------------------------------------------------

# Habitat amount class 

mod_hac_div2_cl_alpha <- lm(ha_slope_div2_cl_alpha ~ ha_class, data = moders_div2_cl_alpha)
mod_hac_div2_cl_alpha_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_alpha, plot = F)
plot(mod_hac_div2_cl_alpha_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_alpha_res) # no overdispersion

summary(mod_hac_div2_cl_alpha) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_alpha$effort_diff_cent <- scale(log(moders_div2_cl_alpha$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_alpha <- lm(ha_slope_div2_cl_alpha ~ effort_diff_cent, data = moders_div2_cl_alpha)

mod_eff_div2_cl_alpha_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_alpha, plot = F)
plot(mod_eff_div2_cl_alpha_res) # residuals normally distributed
testDispersion(mod_eff_div2_cl_alpha_res) # no overdispersion


summary(mod_eff_div2_cl_alpha) # report results # 

# time since fragmentation

mod_time_div2_cl_alpha <- lm(ha_slope_div2_cl_alpha ~ time_since_fragmentation, data = moders_div2_cl_alpha)
mod_time_div2_cl_alpha_res <- simulateResiduals(fittedModel = mod_time_div2_cl_alpha, plot = F)
plot(mod_time_div2_cl_alpha_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_alpha_res) # no overdispersion

summary(mod_time_div2_cl_alpha) # report results

# continent 

mod_cont_div2_cl_alpha <- lm(ha_slope_div2_cl_alpha ~ continent_cat, data = moders_div2_cl_alpha)
mod_cont_div2_cl_alpha_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_alpha, plot = F)
plot(mod_cont_div2_cl_alpha_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_alpha_res) # no overdispersion

summary(mod_cont_div2_cl_alpha)  # report results


## save the results to export later

lm_results(mod_hac_div2_cl_alpha, mod_eff_div2_cl_alpha, mod_time_div2_cl_alpha, mod_cont_div2_cl_alpha) %>% 
  dplyr::mutate(Diversity_type = "Alpha (nearest pairs)", .before = "Included_moderator") -> moder_results_div2_cl_alpha


## Effects on Beta-Diversity --------------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter beta-diversity data

frag_cont_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "beta") -> div_cl2_beta_glmm

# Prepare moderators for beta-diversity
div_cl2_beta_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_beta

### Model Fitting ---------------------------------------------------------------
# Fit initial model

mod1_div2_cl_beta_2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = div_cl2_beta_glmm)
summary(mod1_div2_cl_beta_2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_cl_beta_res <- simulateResiduals(fittedModel = mod1_div2_cl_beta_2, plot = F)
plot(mod_div2_cl_beta_res) # residuals normally distributed
testDispersion(mod_div2_cl_beta_res) # no overdispersion 

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML

mod1_div2_cl_beta_final <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                   data = div_cl2_beta_glmm,
                                   REML = TRUE)
summary(mod1_div2_cl_beta_final) # report results


# Fit models with habitat amount as a random variable

mod1_div2_cl_beta_ha_final <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                  REML = TRUE,
                                  data = div_cl2_beta_glmm)

summary(mod1_div2_cl_beta_ha_final) # report results 

### Extracting Slopes -----------------------------------------------------------

mod1_div2_cl_beta_2_ha <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                  data = div_cl2_beta_glmm)


ha_slope_div2_cl_beta <- ranef(mod1_div2_cl_beta_2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_beta <- lm(ha_slope_div2_cl_beta ~ ha_class, data = moders_div2_cl_beta)
mod_hac_div2_cl_beta_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_beta, plot = F)
plot(mod_hac_div2_cl_beta_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_beta_res) # no overdispersion

summary(mod_hac_div2_cl_beta) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_beta$effort_diff_cent <- scale(log(moders_div2_cl_beta$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_beta <- lm(ha_slope_div2_cl_beta ~ effort_diff_cent, data = moders_div2_cl_beta)

mod_eff_div2_cl_beta_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_beta, plot = F)
plot(mod_eff_div2_cl_beta_res) # residuals normally distributed
testDispersion(mod_eff_div2_cl_beta_res) # no overdispersion


summary(mod_eff_div2_cl_beta) # report results # 

# time since fragmentation

mod_time_div2_cl_beta <- lm(ha_slope_div2_cl_beta ~ time_since_fragmentation, data = moders_div2_cl_beta)
mod_time_div2_cl_beta_res <- simulateResiduals(fittedModel = mod_time_div2_cl_beta, plot = F)
plot(mod_time_div2_cl_beta_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_beta_res) # no overdispersion

summary(mod_time_div2_cl_beta) # report results

# continent 

mod_cont_div2_cl_beta <- lm(ha_slope_div2_cl_beta ~ continent_cat, data = moders_div2_cl_beta)
mod_cont_div2_cl_beta_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_beta, plot = F)
plot(mod_cont_div2_cl_beta_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_beta_res) # no overdispersion

summary(mod_cont_div2_cl_beta)  # report results

## save the results to export later

lm_results(mod_hac_div2_cl_beta, mod_eff_div2_cl_beta, mod_time_div2_cl_beta, mod_cont_div2_cl_beta) %>% 
  dplyr::mutate(Diversity_type = "Beta (nearest pairs)", .before = "Included_moderator") -> moder_results_div2_cl_beta

## Effects on Gamma-Diversity --------------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter gamma-diversity data

frag_cont_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "gamma") -> div_cl2_gamma_glmm

# Prepare moderators for gamma-diversity

div_cl2_gamma_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_gamma

### Model Fitting ---------------------------------------------------------------
# Fit initial model

mod1_div2_cl_gamma_2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = div_cl2_gamma_glmm)
summary(mod1_div2_cl_gamma_2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_cl_gamma_res <- simulateResiduals(fittedModel = mod1_div2_cl_gamma_2, plot = F)
plot(mod_div2_cl_gamma_res) # residuals normally distributed
testDispersion(mod_div2_cl_gamma_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML

mod1_div2_cl_gamma_final <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                    data = div_cl2_gamma_glmm,
                                    REML = TRUE)

summary(mod1_div2_cl_gamma_final) # report results

# Fit models with habitat amount as a random variable

mod1_div2_cl_gamma_ha_final <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                         control=glmmTMBControl(optimizer=optim,
                                                                optArgs=list(method="L-BFGS-B")),
                                         REML = TRUE, data = div_cl2_gamma_glmm)

summary(mod1_div2_cl_gamma_ha_final)


### Extracting Slopes -----------------------------------------------------------

mod1_div2_cl_gamma_2_ha <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                   data = div_cl2_gamma_glmm)


ha_slope_div2_cl_gamma <- ranef(mod1_div2_cl_gamma_2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_gamma <- lm(ha_slope_div2_cl_gamma ~ ha_class, data = moders_div2_cl_gamma)
mod_hac_div2_cl_gamma_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_gamma, plot = F)
plot(mod_hac_div2_cl_gamma_res) # residuals not normally distributed and variances heterogeneous
testDispersion(mod_hac_div2_cl_gamma_res) # no overdispersion

summary(mod_hac_div2_cl_gamma) # report results (note: this model does not have normal distribution)


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_gamma$effort_diff_cent <- scale(log(moders_div2_cl_gamma$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_gamma <- lm(ha_slope_div2_cl_gamma ~ effort_diff_cent, data = moders_div2_cl_gamma)

mod_eff_div2_cl_gamma_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_gamma, plot = F)
plot(mod_eff_div2_cl_gamma_res) # residuals normally distributed and heterogeneous variance 
testDispersion(mod_eff_div2_cl_gamma_res) # no overdispersion


summary(mod_eff_div2_cl_gamma) # report results (note: this model does not have homogeneous variance)

# time since fragmentation

mod_time_div2_cl_gamma <- lm(ha_slope_div2_cl_gamma ~ time_since_fragmentation, data = moders_div2_cl_gamma)
mod_time_div2_cl_gamma_res <- simulateResiduals(fittedModel = mod_time_div2_cl_gamma, plot = F)
plot(mod_time_div2_cl_gamma_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_gamma_res) # no overdispersion

summary(mod_time_div2_cl_gamma) # report results

# continent 

mod_cont_div2_cl_gamma <- lm(ha_slope_div2_cl_gamma ~ continent_cat, data = moders_div2_cl_gamma)
mod_cont_div2_cl_gamma_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_gamma, plot = F)
plot(mod_cont_div2_cl_gamma_res) # residuals not normally distributed
testDispersion(mod_cont_div2_cl_gamma_res) # no overdispersion

summary(mod_cont_div2_cl_gamma)  # report results  (note: this model does not have normal distribution)

## save the results to export later

lm_results(mod_hac_div2_cl_gamma, mod_eff_div2_cl_gamma, mod_time_div2_cl_gamma, mod_cont_div2_cl_gamma) %>% 
  dplyr::mutate(Diversity_type = "Gamma (nearest pairs)", .before = "Included_moderator") -> moder_results_div2_cl_gamma


### GLMMs ----------------------------------------------------------------------
### 3. Fitting models using q = 0 and q = 2 using all pairs ### 

## Effects on Alpha-Diversity (q = 0) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter alpha-diversity data 

frag_cont_rar_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "alpha" & q_order == "q = 0") -> raref_alpha_q0_glmm

# Prepare moderators for alpha-diversity
raref_alpha_q0_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_alpha_r0

### Model Fitting ---------------------------------------------------------------
# Fit initial model

mod1_div2_alpha_r0 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_alpha_q0_glmm)
summary(mod1_div2_alpha_r0)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_alpha_r0_res <- simulateResiduals(fittedModel = mod1_div2_alpha_r0, plot = F)
plot(mod_div2_alpha_r0_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_alpha_r0_res) # no overdispersion


### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response

mod1_div2_alpha_r0_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_alpha_q0_glmm)

# Check the transformed model
mod_div2_alpha_r0_log_res <- simulateResiduals(fittedModel = mod1_div2_alpha_r0_log, plot = F)
plot(mod_div2_alpha_r0_log_res) # residuals normally distributed
testDispersion(mod_div2_alpha_r0_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_alpha_r0_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                        REML = TRUE,
                                        data = raref_alpha_q0_glmm)
summary(mod1_div2_alpha_r0_log_final) # report results


# Fit models with habitat amount as a random variable
mod1_div2_alpha_r0_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 REML = TRUE,
                                 data = raref_alpha_q0_glmm)

summary(mod1_div2_alpha_r0_ha_final)

### Extracting Slopes -----------------------------------------------------------
mod1_div2_alpha_r0_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 data = raref_alpha_q0_glmm)


ha_slope_div2_alpha_r0 <- ranef(mod1_div2_alpha_r0_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_alpha_r0 <- lm(ha_slope_div2_alpha_r0 ~ ha_class, data = moders_div2_alpha_r0)
mod_hac_div2_alpha_r0_res <- simulateResiduals(fittedModel = mod_hac_div2_alpha_r0, plot = F)
plot(mod_hac_div2_alpha_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_alpha_r0_res) # no overdispersion

summary(mod_hac_div2_alpha_r0) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_alpha_r0$effort_diff_cent <- scale(log(moders_div2_alpha_r0$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_alpha_r0 <- lm(ha_slope_div2_alpha_r0 ~ effort_diff_cent, data = moders_div2_alpha_r0)

mod_eff_div2_alpha_r0_res <- simulateResiduals(fittedModel = mod_eff_div2_alpha_r0, plot = F)
plot(mod_eff_div2_alpha_r0_res) # residuals normally distributed
testDispersion(mod_eff_div2_alpha_r0_res) # no overdispersion


summary(mod_eff_div2_alpha_r0) # report results # 


# time since fragmentation

mod_time_div2_alpha_r0 <- lm(ha_slope_div2_alpha_r0 ~ time_since_fragmentation, data = moders_div2_alpha_r0)
mod_time_div2_alpha_r0_res <- simulateResiduals(fittedModel = mod_time_div2_alpha_r0, plot = F)
plot(mod_time_div2_alpha_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_alpha_r0_res) # no overdispersion

summary(mod_time_div2_alpha_r0) # report results

# continent 

mod_cont_div2_alpha_r0 <- lm(ha_slope_div2_alpha_r0 ~ continent_cat, data = moders_div2_alpha_r0)
mod_cont_div2_alpha_r0_res <- simulateResiduals(fittedModel = mod_cont_div2_alpha_r0, plot = F)
plot(mod_cont_div2_alpha_r0_res) # residuals normally distributed
testDispersion(mod_cont_div2_alpha_r0_res) # no overdispersion

summary(mod_cont_div2_alpha_r0)  # report results


## save the results to export later
lm_results(mod_hac_div2_alpha_r0, mod_eff_div2_alpha_r0, mod_time_div2_alpha_r0, mod_cont_div2_alpha_r0) %>% 
  dplyr::mutate(Diversity_type = "Alpha (all pairs) - q = 0", .before = "Included_moderator") -> moder_results_div2_alpha_r0

## Effects on Beta-Diversity (q = 0) -------------------------------------------

### Data Preparation (q = 0) ---------------------------------------------------
# Filter beta-diversity data
frag_cont_rar_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "beta" & q_order == "q = 0") -> raref_beta_q0_glmm

# Prepare moderators for beta-diversity
raref_beta_q0_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_beta_r0

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_beta_r0 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_beta_q0_glmm)
summary(mod1_div2_beta_r0)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_beta_r0_res <- simulateResiduals(fittedModel = mod1_div2_beta_r0, plot = F)
plot(mod_div2_beta_r0_res) # residuals normally distributed
testDispersion(mod_div2_beta_r0_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_beta_r0_final <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                   REML = TRUE,
                                   data = raref_beta_q0_glmm)
summary(mod1_div2_beta_r0_final) # report results

# Fit models with habitat amount as a random variable
mod1_div2_beta_r0_ha_final <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                REML = TRUE,
                                data = raref_beta_q0_glmm)

summary(mod1_div2_beta_r0_ha_final) # report results

### Extracting Slopes -----------------------------------------------------------
mod1_div2_beta_r0_ha <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                                data = raref_beta_q0_glmm)


ha_slope_div2_beta_r0 <- ranef(mod1_div2_beta_r0_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_beta_r0 <- lm(ha_slope_div2_beta_r0 ~ ha_class, data = moders_div2_beta_r0)
mod_hac_div2_beta_r0_res <- simulateResiduals(fittedModel = mod_hac_div2_beta_r0, plot = F)
plot(mod_hac_div2_beta_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_beta_r0_res) # no overdispersion

summary(mod_hac_div2_beta_r0) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_beta_r0$effort_diff_cent <- scale(log(moders_div2_beta_r0$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_beta_r0 <- lm(ha_slope_div2_beta_r0 ~ effort_diff_cent, data = moders_div2_beta_r0)

mod_eff_div2_beta_r0_res <- simulateResiduals(fittedModel = mod_eff_div2_beta_r0, plot = F)
plot(mod_eff_div2_beta_r0_res) # residuals normally distributed
testDispersion(mod_eff_div2_beta_r0_res) # no overdispersion

summary(mod_eff_div2_beta_r0) # report results # 


# time since fragmentation

mod_time_div2_beta_r0 <- lm(ha_slope_div2_beta_r0 ~ time_since_fragmentation, data = moders_div2_beta_r0)
mod_time_div2_beta_r0_res <- simulateResiduals(fittedModel = mod_time_div2_beta_r0, plot = F)
plot(mod_time_div2_beta_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_beta_r0_res) # no overdispersion

summary(mod_time_div2_beta_r0) # report results

# continent 

mod_cont_div2_beta_r0 <- lm(ha_slope_div2_beta_r0 ~ continent_cat, data = moders_div2_beta_r0)
mod_cont_div2_beta_r0_res <- simulateResiduals(fittedModel = mod_cont_div2_beta_r0, plot = F)
plot(mod_cont_div2_beta_r0_res) # residuals normally distributed
testDispersion(mod_cont_div2_beta_r0_res) # no overdispersion

summary(mod_cont_div2_beta_r0)  # report results


## save the results to export later

lm_results(mod_hac_div2_beta_r0, mod_eff_div2_beta_r0, mod_time_div2_beta_r0, mod_cont_div2_beta_r0) %>% 
  dplyr::mutate(Diversity_type = "Beta (all pairs) - q = 0", .before = "Included_moderator") -> moder_results_div2_beta_r0

## Effects on Gamma-Diversity (q = 0) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter gamma-diversity data

frag_cont_rar_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "gamma" & q_order == "q = 0") -> raref_gamma_q0_glmm

# Prepare moderators for gamma-diversity
raref_gamma_q0_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_gamma_r0

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_gamma_r0 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_gamma_q0_glmm)
summary(mod1_div2_gamma_r0)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_gamma_r0_res <- simulateResiduals(fittedModel = mod1_div2_gamma_r0, plot = F)
plot(mod_div2_gamma_r0_res) # # residuals not normally distributed and variance heterogeneous - we need to log transform it
testDispersion(mod_div2_gamma_r0_res) # no overdispersion

### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response
mod1_div2_gamma_r0_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_gamma_q0_glmm)
summary(mod1_div2_gamma_r0_log)

# Check the transformed model
mod_div2_gamma_r0_log_res <- simulateResiduals(fittedModel = mod1_div2_gamma_r0_log, plot = F)
plot(mod_div2_gamma_r0_log_res) # residuals normally distributed
testDispersion(mod_div2_gamma_r0_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_gamma_r0_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                        REML = TRUE,
                                        data = raref_gamma_q0_glmm)
summary(mod1_div2_gamma_r0_log_final) # report results

# Fit models with habitat amount as a random variable
mod1_div2_gamma_r0_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                 REML = TRUE,
                                 data = raref_gamma_q0_glmm)

summary(mod1_div2_gamma_r0_ha_final) # report results


### Extracting Slopes -----------------------------------------------------------
mod1_div2_gamma_r0_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                 data = raref_gamma_q0_glmm)


ha_slope_div2_gamma_r0 <- ranef(mod1_div2_gamma_r0_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_gamma_r0 <- lm(ha_slope_div2_gamma_r0 ~ ha_class, data = moders_div2_gamma_r0)
mod_hac_div2_gamma_r0_res <- simulateResiduals(fittedModel = mod_hac_div2_gamma_r0, plot = F)
plot(mod_hac_div2_gamma_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_gamma_r0_res) # no overdispersion

summary(mod_hac_div2_gamma_r0) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_gamma_r0$effort_diff_cent <- scale(log(moders_div2_gamma_r0$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_gamma_r0 <- lm(ha_slope_div2_gamma_r0 ~ effort_diff_cent, data = moders_div2_gamma_r0)

mod_eff_div2_gamma_r0_res <- simulateResiduals(fittedModel = mod_eff_div2_gamma_r0, plot = F)
plot(mod_eff_div2_gamma_r0_res) # residuals normally distributed
testDispersion(mod_eff_div2_gamma_r0_res) # no overdispersion


summary(mod_eff_div2_gamma_r0) # report results # 


# time since fragmentation

mod_time_div2_gamma_r0 <- lm(ha_slope_div2_gamma_r0 ~ time_since_fragmentation, data = moders_div2_gamma_r0)
mod_time_div2_gamma_r0_res <- simulateResiduals(fittedModel = mod_time_div2_gamma_r0, plot = F)
plot(mod_time_div2_gamma_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_gamma_r0_res) # no overdispersion

summary(mod_time_div2_gamma_r0) # report results

# continent 

mod_cont_div2_gamma_r0 <- lm(ha_slope_div2_gamma_r0 ~ continent_cat, data = moders_div2_gamma_r0)
mod_cont_div2_gamma_r0_res <- simulateResiduals(fittedModel = mod_cont_div2_gamma_r0, plot = F)
plot(mod_cont_div2_gamma_r0_res) # residuals normally distributed
testDispersion(mod_cont_div2_gamma_r0_res) # no overdispersion

summary(mod_cont_div2_gamma_r0)  # report results

## save the results to export later

lm_results(mod_hac_div2_gamma_r0, mod_eff_div2_gamma_r0, mod_time_div2_gamma_r0, mod_cont_div2_gamma_r0) %>% 
  dplyr::mutate(Diversity_type = "Gamma (all pairs) - q = 0", .before = "Included_moderator") -> moder_results_div2_gamma_r0


## Effects on Alpha-Diversity (q = 2) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter alpha-diversity data 
frag_cont_rar_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "alpha" & q_order == "q = 2") -> raref_alpha_q2_glmm

# Prepare moderators for alpha-diversity
raref_alpha_q2_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_alpha_r2

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_alpha_r2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_alpha_q2_glmm)
summary(mod1_div2_alpha_r2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_alpha_r2_res <- simulateResiduals(fittedModel = mod1_div2_alpha_r2, plot = F)
plot(mod_div2_alpha_r2_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_alpha_r2_res) # no overdispersion

### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response
mod1_div2_alpha_r2_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_alpha_q2_glmm)

# Check the transformed model
mod_div2_alpha_r2_log_res <- simulateResiduals(fittedModel = mod1_div2_alpha_r2_log, plot = F)
plot(mod_div2_alpha_r2_log_res) # residuals normally distributed
testDispersion(mod_div2_alpha_r2_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_alpha_r2_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                        REML = TRUE,
                                        data = raref_alpha_q2_glmm)
summary(mod1_div2_alpha_r2_log_final) # report results

# Fit models with habitat amount as a random variable

mod1_div2_alpha_r2_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 control=glmmTMBControl(optimizer=optim,
                                                        optArgs=list(method="L-BFGS-B")),
                                 REML = TRUE,
                                 data = raref_alpha_q2_glmm)

summary(mod1_div2_alpha_r2_ha_final) # report results 

### Extracting Slopes -----------------------------------------------------------

mod1_div2_alpha_r2_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 control=glmmTMBControl(optimizer=optim,
                                                        optArgs=list(method="L-BFGS-B")),
                                 data = raref_alpha_q2_glmm)

ha_slope_div2_alpha_r2 <- ranef(mod1_div2_alpha_r2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_alpha_r2 <- lm(ha_slope_div2_alpha_r2 ~ ha_class, data = moders_div2_alpha_r2)
mod_hac_div2_alpha_r2_res <- simulateResiduals(fittedModel = mod_hac_div2_alpha_r2, plot = F)
plot(mod_hac_div2_alpha_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_alpha_r2_res) # no overdispersion

summary(mod_hac_div2_alpha_r2) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_alpha_r2$effort_diff_cent <- scale(log(moders_div2_alpha_r2$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_alpha_r2 <- lm(ha_slope_div2_alpha_r2 ~ effort_diff_cent, data = moders_div2_alpha_r2)

mod_eff_div2_alpha_r2_res <- simulateResiduals(fittedModel = mod_eff_div2_alpha_r2, plot = F)
plot(mod_eff_div2_alpha_r2_res) # residuals normally distributed
testDispersion(mod_eff_div2_alpha_r2_res) # no overdispersion


summary(mod_eff_div2_alpha_r2) # report results # 


# time since fragmentation

mod_time_div2_alpha_r2 <- lm(ha_slope_div2_alpha_r2 ~ time_since_fragmentation, data = moders_div2_alpha_r2)
mod_time_div2_alpha_r2_res <- simulateResiduals(fittedModel = mod_time_div2_alpha_r2, plot = F)
plot(mod_time_div2_alpha_r2_res) # residuals normally distributed, but heterogeneous variance
testDispersion(mod_hac_div2_alpha_r2_res) # no overdispersion

summary(mod_time_div2_alpha_r2) # report results (note: this model does not have homogeneous variance)

# continent 

mod_cont_div2_alpha_r2 <- lm(ha_slope_div2_alpha_r2 ~ continent_cat, data = moders_div2_alpha_r2)
mod_cont_div2_alpha_r2_res <- simulateResiduals(fittedModel = mod_cont_div2_alpha_r2, plot = F)
plot(mod_cont_div2_alpha_r2_res) # residuals normally distributed
testDispersion(mod_cont_div2_alpha_r2_res) # no overdispersion

summary(mod_cont_div2_alpha_r2)  # report results

## save the results to export later

lm_results(mod_hac_div2_alpha_r2, mod_eff_div2_alpha_r2, mod_time_div2_alpha_r2, mod_cont_div2_alpha_r2) %>% 
  dplyr::mutate(Diversity_type = "Alpha (all pairs) - q = 2", .before = "Included_moderator") -> moder_results_div2_alpha_r2


## Effects on Beta-Diversity (q = 2) -------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter beta-diversity data
frag_cont_rar_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "beta" & q_order == "q = 2") -> raref_beta_q2_glmm

# Prepare moderators for beta-diversity
raref_beta_q2_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_beta_r2

### Model Fitting ---------------------------------------------------------------
# Fit initial model

mod1_div2_beta_r2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_beta_q2_glmm)
summary(mod1_div2_beta_r2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_beta_r2_res <- simulateResiduals(fittedModel = mod1_div2_beta_r2, plot = F)
plot(mod_div2_beta_r2_res) # residuals normally distributed 
testDispersion(mod_div2_beta_r2_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_beta_r2_final <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                   REML = TRUE,
                                   data = raref_beta_q2_glmm)
summary(mod1_div2_beta_r2_final) # report results

# Fit models with habitat amount as a random variable
mod1_div2_beta_r2_ha_final <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                control=glmmTMBControl(optimizer=optim,
                                                       optArgs=list(method="L-BFGS-B")),
                                REML = TRUE,
                                data = raref_beta_q2_glmm)

summary(mod1_div2_beta_r2_ha_final) # report results

### Extracting Slopes -----------------------------------------------------------
mod1_div2_beta_r2_ha <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                control=glmmTMBControl(optimizer=optim,
                                                       optArgs=list(method="L-BFGS-B")),
                                data = raref_beta_q2_glmm)


ha_slope_div2_beta_r2 <- ranef(mod1_div2_beta_r2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_beta_r2 <- lm(ha_slope_div2_beta_r2 ~ ha_class, data = moders_div2_beta_r2)
mod_hac_div2_beta_r2_res <- simulateResiduals(fittedModel = mod_hac_div2_beta_r2, plot = F)
plot(mod_hac_div2_beta_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_beta_r2_res) # no overdispersion

summary(mod_hac_div2_beta_r2) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_beta_r2$effort_diff_cent <- scale(log(moders_div2_beta_r2$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_beta_r2 <- lm(ha_slope_div2_beta_r2 ~ effort_diff_cent, data = moders_div2_beta_r2)

mod_eff_div2_beta_r2_res <- simulateResiduals(fittedModel = mod_eff_div2_beta_r2, plot = F)
plot(mod_eff_div2_beta_r2_res) # residuals normally distributed
testDispersion(mod_eff_div2_beta_r2_res) # no overdispersion


summary(mod_eff_div2_beta_r2) # report results # 


# time since fragmentation

mod_time_div2_beta_r2 <- lm(ha_slope_div2_beta_r2 ~ time_since_fragmentation, data = moders_div2_beta_r2)
mod_time_div2_beta_r2_res <- simulateResiduals(fittedModel = mod_time_div2_beta_r2, plot = F)
plot(mod_time_div2_beta_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_beta_r2_res) # no overdispersion

summary(mod_time_div2_beta_r2) # report results 

# continent 

mod_cont_div2_beta_r2 <- lm(ha_slope_div2_beta_r2 ~ continent_cat, data = moders_div2_beta_r2)
mod_cont_div2_beta_r2_res <- simulateResiduals(fittedModel = mod_cont_div2_beta_r2, plot = F)
plot(mod_cont_div2_beta_r2_res) # residuals normally distributed
testDispersion(mod_cont_div2_beta_r2_res) # no overdispersion

summary(mod_cont_div2_beta_r2)  # report results


## save the results to export later

lm_results(mod_hac_div2_beta_r2, mod_eff_div2_beta_r2, mod_time_div2_beta_r2, mod_cont_div2_beta_r2) %>% 
  dplyr::mutate(Diversity_type = "Beta (all pairs) - q = 2", .before = "Included_moderator") -> moder_results_div2_beta_r2


## Effects on Gamma-Diversity (q = 2) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter gamma-diversity data
frag_cont_rar_div2 %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "gamma" & q_order == "q = 2") -> raref_gamma_q2_glmm

# Prepare moderators for gamma-diversity
raref_gamma_q2_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_gamma_r2

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_gamma_r2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_gamma_q2_glmm)
summary(mod1_div2_gamma_r2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_gamma_r2_res <- simulateResiduals(fittedModel = mod1_div2_gamma_r2, plot = F)
plot(mod_div2_gamma_r2_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_gamma_r2_res) # no overdispersion

### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response
mod1_div2_gamma_r2_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_gamma_q2_glmm)
summary(mod1_div2_gamma_r2_log)

# Check the transformed model
mod_div2_gamma_r2_log_res <- simulateResiduals(fittedModel = mod1_div2_gamma_r2_log, plot = F)
plot(mod_div2_gamma_r2_log_res) # residuals normally distributed  
testDispersion(mod_div2_gamma_r2_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_gamma_r2_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                    REML = TRUE,
                                    data = raref_gamma_q2_glmm)
summary(mod1_div2_gamma_r2_final) # report results

# Fit models with habitat amount as a random variable

mod1_div2_gamma_r2_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 control=glmmTMBControl(optimizer=optim,
                                                        optArgs=list(method="L-BFGS-B")),
                                 REML = TRUE,
                                 data = raref_gamma_q2_glmm)

summary(mod1_div2_gamma_r2_ha_final) # report results

### Extracting Slopes -----------------------------------------------------------
mod1_div2_gamma_r2_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                 control=glmmTMBControl(optimizer=optim,
                                                        optArgs=list(method="L-BFGS-B")),
                                 data = raref_gamma_q2_glmm)

ha_slope_div2_gamma_r2 <- ranef(mod1_div2_gamma_r2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_gamma_r2 <- lm(ha_slope_div2_gamma_r2 ~ ha_class, data = moders_div2_gamma_r2)
mod_hac_div2_gamma_r2_res <- simulateResiduals(fittedModel = mod_hac_div2_gamma_r2, plot = F)
plot(mod_hac_div2_gamma_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_gamma_r2_res) # no overdispersion

summary(mod_hac_div2_gamma_r2) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_gamma_r2$effort_diff_cent <- scale(log(moders_div2_gamma_r2$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_gamma_r2 <- lm(ha_slope_div2_gamma_r2 ~ effort_diff_cent, data = moders_div2_gamma_r2)

mod_eff_div2_gamma_r2_res <- simulateResiduals(fittedModel = mod_eff_div2_gamma_r2, plot = F)
plot(mod_eff_div2_gamma_r2_res) # residuals normally distributed
testDispersion(mod_eff_div2_gamma_r2_res) # no overdispersion


summary(mod_eff_div2_gamma_r2) # report results # 


# time since fragmentation

mod_time_div2_gamma_r2 <- lm(ha_slope_div2_gamma_r2 ~ time_since_fragmentation, data = moders_div2_gamma_r2)
mod_time_div2_gamma_r2_res <- simulateResiduals(fittedModel = mod_time_div2_gamma_r2, plot = F)
plot(mod_time_div2_gamma_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_gamma_r2_res) # no overdispersion

summary(mod_time_div2_gamma_r2) # report results 

# continent 

mod_cont_div2_gamma_r2 <- lm(ha_slope_div2_gamma_r2 ~ continent_cat, data = moders_div2_gamma_r2)
mod_cont_div2_gamma_r2_res <- simulateResiduals(fittedModel = mod_cont_div2_gamma_r2, plot = F)
plot(mod_cont_div2_gamma_r2_res) # residuals normally distributed
testDispersion(mod_cont_div2_gamma_r2_res) # no overdispersion

summary(mod_cont_div2_gamma_r2)  # report results

## save the results to export later 
lm_results(mod_hac_div2_gamma_r2, mod_eff_div2_gamma_r2, mod_time_div2_gamma_r2, mod_cont_div2_gamma_r2) %>% 
  dplyr::mutate(Diversity_type = "Gamma (all pairs) - q = 2", .before = "Included_moderator") -> moder_results_div2_gamma_r2



### GLMMs ----------------------------------------------------------------------
### 4. Fitting models using q = 0 and q = 2 using the nearest plot pairs ### 

## Effects on Alpha-Diversity (q = 0) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter alpha-diversity data
frag_cont_rar_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "alpha" & q_order == "q = 0") -> raref_cl_alpha_q0_glmm

# Prepare moderators for alpha-diversity
raref_cl_alpha_q0_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_alpha_r0

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_cl_alpha_r0 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_cl_alpha_q0_glmm)
summary(mod1_div2_cl_alpha_r0)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_cl_alpha_r0_res <- simulateResiduals(fittedModel = mod1_div2_cl_alpha_r0, plot = F)
plot(mod_div2_cl_alpha_r0_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_cl_alpha_r0_res) # no overdispersion+


### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response
mod1_div2_cl_alpha_r0_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_cl_alpha_q0_glmm)

# Check the transformed model
mod_div2_cl_alpha_r0_log_res <- simulateResiduals(fittedModel = mod1_div2_cl_alpha_r0_log, plot = F)
plot(mod_div2_cl_alpha_r0_log_res) # residuals normally distributed
testDispersion(mod_div2_cl_alpha_r0_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_cl_alpha_r0_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                           REML = TRUE,
                                           data = raref_cl_alpha_q0_glmm)
summary(mod1_div2_cl_alpha_r0_log_final) # report results

# Fit models with habitat amount as a random variable
mod1_div2_cl_alpha_r0_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                          REML = TRUE, data = raref_cl_alpha_q0_glmm)

summary(mod1_div2_cl_alpha_r0_ha_final) # report results

### Extracting Slopes -----------------------------------------------------------
mod1_div2_cl_alpha_r0_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                    data = raref_cl_alpha_q0_glmm)


ha_slope_div2_cl_alpha_r0 <- ranef(mod1_div2_cl_alpha_r0_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_alpha_r0 <- lm(ha_slope_div2_cl_alpha_r0 ~ ha_class, data = moders_div2_cl_alpha_r0)
mod_hac_div2_cl_alpha_r0_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_alpha_r0, plot = F)
plot(mod_hac_div2_cl_alpha_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_alpha_r0_res) # no overdispersion

summary(mod_hac_div2_cl_alpha_r0) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_alpha_r0$effort_diff_cent <- scale(log(moders_div2_cl_alpha_r0$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_alpha_r0 <- lm(ha_slope_div2_cl_alpha_r0 ~ effort_diff_cent, data = moders_div2_cl_alpha_r0)

mod_eff_div2_cl_alpha_r0_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_alpha_r0, plot = F)
plot(mod_eff_div2_cl_alpha_r0_res) # residuals normally distributed
testDispersion(mod_eff_div2_cl_alpha_r0_res) # no overdispersion


summary(mod_eff_div2_cl_alpha_r0) # report results # 


# time since fragmentation

mod_time_div2_cl_alpha_r0 <- lm(ha_slope_div2_cl_alpha_r0 ~ time_since_fragmentation, data = moders_div2_cl_alpha_r0)
mod_time_div2_cl_alpha_r0_res <- simulateResiduals(fittedModel = mod_time_div2_cl_alpha_r0, plot = F)
plot(mod_time_div2_cl_alpha_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_alpha_r0_res) # no overdispersion

summary(mod_time_div2_cl_alpha_r0) # report results

# continent 

mod_cont_div2_cl_alpha_r0 <- lm(ha_slope_div2_cl_alpha_r0 ~ continent_cat, data = moders_div2_cl_alpha_r0)
mod_cont_div2_cl_alpha_r0_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_alpha_r0, plot = F)
plot(mod_cont_div2_cl_alpha_r0_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_alpha_r0_res) # no overdispersion

summary(mod_cont_div2_cl_alpha_r0)  # report results

## save the results to export later

lm_results(mod_hac_div2_cl_alpha_r0, mod_eff_div2_cl_alpha_r0, mod_time_div2_cl_alpha_r0, mod_cont_div2_cl_alpha_r0) %>% 
  dplyr::mutate(Diversity_type = "Alpha (nearest pairs) - q = 0", .before = "Included_moderator") -> moder_results_div2_cl_alpha_r0

## Effects on Beta-Diversity (q = 0) -------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter beta-diversity data

frag_cont_rar_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "beta" & q_order == "q = 0") -> raref_cl_beta_q0_glmm

# Prepare moderators for beta-diversity
raref_cl_beta_q0_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_beta_r0

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_cl_beta_r0 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_cl_beta_q0_glmm)
summary(mod1_div2_cl_beta_r0)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_cl_beta_r0_res <- simulateResiduals(fittedModel = mod1_div2_cl_beta_r0, plot = F)
plot(mod_div2_cl_beta_r0_res) # residuals normally distributed
testDispersion(mod_div2_cl_beta_r0_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_cl_beta_r0_final <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                      REML = TRUE,
                                      data = raref_cl_beta_q0_glmm)
summary(mod1_div2_cl_beta_r0_final) # report results


# Fit models with habitat amount as a random variable

mod1_div2_cl_beta_r0_ha_final <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                   REML = TRUE,
                                   data = raref_cl_beta_q0_glmm) 

summary(mod1_div2_cl_beta_r0_ha_final) # report results (note: Model convergence problem independent of the optimizer)

### Extracting Slopes -----------------------------------------------------------
mod1_div2_cl_beta_r0_ha <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
                                   data = raref_cl_beta_q0_glmm) # Model convergence problem


ha_slope_div2_cl_beta_r0 <- ranef(mod1_div2_cl_beta_r0_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_beta_r0 <- lm(ha_slope_div2_cl_beta_r0 ~ ha_class, data = moders_div2_cl_beta_r0)
mod_hac_div2_cl_beta_r0_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_beta_r0, plot = F)
plot(mod_hac_div2_cl_beta_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_beta_r0_res) # no overdispersion

summary(mod_hac_div2_cl_beta_r0) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_beta_r0$effort_diff_cent <- scale(log(moders_div2_cl_beta_r0$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_beta_r0 <- lm(ha_slope_div2_cl_beta_r0 ~ effort_diff_cent, data = moders_div2_cl_beta_r0)

mod_eff_div2_cl_beta_r0_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_beta_r0, plot = F)
plot(mod_eff_div2_cl_beta_r0_res) # residuals normally distributed, but heterogeneous variance
testDispersion(mod_eff_div2_cl_beta_r0_res) # no overdispersion


summary(mod_eff_div2_cl_beta_r0) # report results (note: model does not have homogeneous variance) 


# time since fragmentation

mod_time_div2_cl_beta_r0 <- lm(ha_slope_div2_cl_beta_r0 ~ time_since_fragmentation, data = moders_div2_cl_beta_r0)
mod_time_div2_cl_beta_r0_res <- simulateResiduals(fittedModel = mod_time_div2_cl_beta_r0, plot = F)
plot(mod_time_div2_cl_beta_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_beta_r0_res) # no overdispersion

summary(mod_time_div2_cl_beta_r0) # report results

# continent 

mod_cont_div2_cl_beta_r0 <- lm(ha_slope_div2_cl_beta_r0 ~ continent_cat, data = moders_div2_cl_beta_r0)
mod_cont_div2_cl_beta_r0_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_beta_r0, plot = F)
plot(mod_cont_div2_cl_beta_r0_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_beta_r0_res) # no overdispersion

summary(mod_cont_div2_cl_beta_r0)  # report results


## save the results to export later

lm_results(mod_hac_div2_cl_beta_r0, mod_eff_div2_cl_beta_r0, mod_time_div2_cl_beta_r0, mod_cont_div2_cl_beta_r0) %>% 
  dplyr::mutate(Diversity_type = "Beta (nearest pairs) - q = 0", .before = "Included_moderator") -> moder_results_div2_cl_beta_r0

## Effects on Gamma-Diversity (q = 0) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter gamma-diversity data
frag_cont_rar_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "gamma" & q_order == "q = 0") -> raref_cl_gamma_q0_glmm

# Prepare moderators for gamma-diversity
raref_cl_gamma_q0_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_gamma_r0

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_cl_gamma_r0 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_cl_gamma_q0_glmm)
summary(mod1_div2_cl_gamma_r0)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_cl_gamma_r0_res <- simulateResiduals(fittedModel = mod1_div2_cl_gamma_r0, plot = F)
plot(mod_div2_cl_gamma_r0_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_cl_gamma_r0_res) # no overdispersion

### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response
mod1_div2_cl_gamma_r0_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_cl_gamma_q0_glmm)
summary(mod1_div2_cl_gamma_r0_log)

# Check the transformed model
mod_div2_cl_gamma_r0_log_res <- simulateResiduals(fittedModel = mod1_div2_cl_gamma_r0_log, plot = F)
plot(mod_div2_cl_gamma_r0_log_res) # residuals normally distributed
testDispersion(mod_div2_cl_gamma_r0_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_cl_gamma_r0_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                           REML = TRUE,
                                           data = raref_cl_gamma_q0_glmm)
summary(mod1_div2_cl_gamma_r0_log_final) # report results

# Fit models with habitat amount as a random variable
mod1_div2_cl_gamma_r0_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                    control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                    REML = TRUE,
                                    data = raref_cl_gamma_q0_glmm)

summary(mod1_div2_cl_gamma_r0_ha_final) # report results 

### Extracting Slopes -----------------------------------------------------------
mod1_div2_cl_gamma_r0_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                    control=glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                                    data = raref_cl_gamma_q0_glmm)


ha_slope_div2_cl_gamma_r0 <- ranef(mod1_div2_cl_gamma_r0_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_gamma_r0 <- lm(ha_slope_div2_cl_gamma_r0 ~ ha_class, data = moders_div2_cl_gamma_r0)
mod_hac_div2_cl_gamma_r0_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_gamma_r0, plot = F)
plot(mod_hac_div2_cl_gamma_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_gamma_r0_res) # no overdispersion

summary(mod_hac_div2_cl_gamma_r0) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_gamma_r0$effort_diff_cent <- scale(log(moders_div2_cl_gamma_r0$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_gamma_r0 <- lm(ha_slope_div2_cl_gamma_r0 ~ effort_diff_cent, data = moders_div2_cl_gamma_r0)

mod_eff_div2_cl_gamma_r0_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_gamma_r0, plot = F)
plot(mod_eff_div2_cl_gamma_r0_res) # residuals normally distributed
testDispersion(mod_eff_div2_cl_gamma_r0_res) # no overdispersion


summary(mod_eff_div2_cl_gamma_r0) # report results # 


# time since fragmentation

mod_time_div2_cl_gamma_r0 <- lm(ha_slope_div2_cl_gamma_r0 ~ time_since_fragmentation, data = moders_div2_cl_gamma_r0)
mod_time_div2_cl_gamma_r0_res <- simulateResiduals(fittedModel = mod_time_div2_cl_gamma_r0, plot = F)
plot(mod_time_div2_cl_gamma_r0_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_gamma_r0_res) # no overdispersion

summary(mod_time_div2_cl_gamma_r0) # report results

# continent 

mod_cont_div2_cl_gamma_r0 <- lm(ha_slope_div2_cl_gamma_r0 ~ continent_cat, data = moders_div2_cl_gamma_r0)
mod_cont_div2_cl_gamma_r0_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_gamma_r0, plot = F)
plot(mod_cont_div2_cl_gamma_r0_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_gamma_r0_res) # no overdispersion

summary(mod_cont_div2_cl_gamma_r0)  # report results


## save the results to export later

lm_results(mod_hac_div2_cl_gamma_r0, mod_eff_div2_cl_gamma_r0, mod_time_div2_cl_gamma_r0, mod_cont_div2_cl_gamma_r0) %>% 
  dplyr::mutate(Diversity_type = "Gamma (nearest pairs) - q = 0", .before = "Included_moderator") -> moder_results_div2_cl_gamma_r0


## Effects on Alpha-Diversity (q = 2) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter alpha-diversity data

frag_cont_rar_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "alpha" & q_order == "q = 2") -> raref_cl_alpha_q2_glmm

# Prepare moderators for alpha-diversity
raref_cl_alpha_q2_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_alpha_r2

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_cl_alpha_r2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_cl_alpha_q2_glmm)
summary(mod1_div2_cl_alpha_r2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_cl_alpha_r2_res <- simulateResiduals(fittedModel = mod1_div2_cl_alpha_r2, plot = F)
plot(mod_div2_cl_alpha_r2_res) # residuals not normally distributed - we need to log transform it
testDispersion(mod_div2_cl_alpha_r2_res) # no overdispersion

### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response
mod1_div2_cl_alpha_r2_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_cl_alpha_q2_glmm)

# Check the transformed model
mod_div2_cl_alpha_r2_log_res <- simulateResiduals(fittedModel = mod1_div2_cl_alpha_r2_log, plot = F)
plot(mod_div2_cl_alpha_r2_log_res) # residuals normally distributed
testDispersion(mod_div2_cl_alpha_r2_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_cl_alpha_r2_log_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                           REML = TRUE,
                                           data = raref_cl_alpha_q2_glmm)
summary(mod1_div2_cl_alpha_r2_log_final) # report results

# Fit models with habitat amount as a random variable
mod1_div2_cl_alpha_r2_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                    control=glmmTMBControl(optimizer=optim,
                                                           optArgs=list(method="L-BFGS-B")),
                                    REML = TRUE,
                                    data = raref_cl_alpha_q2_glmm)

summary(mod1_div2_cl_alpha_r2_ha_final) # report results


### Extracting Slopes -----------------------------------------------------------
mod1_div2_cl_alpha_r2_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                    control=glmmTMBControl(optimizer=optim,
                                                           optArgs=list(method="L-BFGS-B")),
                                    data = raref_cl_alpha_q2_glmm)


ha_slope_div2_cl_alpha_r2 <- ranef(mod1_div2_cl_alpha_r2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_alpha_r2 <- lm(ha_slope_div2_cl_alpha_r2 ~ ha_class, data = moders_div2_cl_alpha_r2)
mod_hac_div2_cl_alpha_r2_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_alpha_r2, plot = F)
plot(mod_hac_div2_cl_alpha_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_alpha_r2_res) # no overdispersion

summary(mod_hac_div2_cl_alpha_r2) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_alpha_r2$effort_diff_cent <- scale(log(moders_div2_cl_alpha_r2$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_alpha_r2 <- lm(ha_slope_div2_cl_alpha_r2 ~ effort_diff_cent, data = moders_div2_cl_alpha_r2)

mod_eff_div2_cl_alpha_r2_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_alpha_r2, plot = F)
plot(mod_eff_div2_cl_alpha_r2_res) # residuals normally distributed
testDispersion(mod_eff_div2_cl_alpha_r2_res) # no overdispersion


summary(mod_eff_div2_cl_alpha_r2) # report results


# time since fragmentation

mod_time_div2_cl_alpha_r2 <- lm(ha_slope_div2_cl_alpha_r2 ~ time_since_fragmentation, data = moders_div2_cl_alpha_r2)
mod_time_div2_cl_alpha_r2_res <- simulateResiduals(fittedModel = mod_time_div2_cl_alpha_r2, plot = F)
plot(mod_time_div2_cl_alpha_r2_res) # residuals normally distributed, but heterogeneous variance
testDispersion(mod_hac_div2_cl_alpha_r2_res) # no overdispersion

summary(mod_time_div2_cl_alpha_r2) # report results (note: model does not have homogeneous variance)

# continent 

mod_cont_div2_cl_alpha_r2 <- lm(ha_slope_div2_cl_alpha_r2 ~ continent_cat, data = moders_div2_cl_alpha_r2)
mod_cont_div2_cl_alpha_r2_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_alpha_r2, plot = F)
plot(mod_cont_div2_cl_alpha_r2_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_alpha_r2_res) # no overdispersion

summary(mod_cont_div2_cl_alpha_r2)  # report results

## save the results to export later

lm_results(mod_hac_div2_cl_alpha_r2, mod_eff_div2_cl_alpha_r2, mod_time_div2_cl_alpha_r2, mod_cont_div2_cl_alpha_r2) %>% 
  dplyr::mutate(Diversity_type = "Alpha (nearest pairs) - q = 2", .before = "Included_moderator") -> moder_results_div2_cl_alpha_r2

## Effects on Beta-Diversity (q = 2) -------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter beta-diversity data
frag_cont_rar_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "beta" & q_order == "q = 2") -> raref_cl_beta_q2_glmm

# Prepare moderators for beta-diversity
raref_cl_beta_q2_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_beta_r2

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_cl_beta_r2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_cl_beta_q2_glmm)
summary(mod1_div2_cl_beta_r2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions
mod_div2_cl_beta_r2_res <- simulateResiduals(fittedModel = mod1_div2_cl_beta_r2, plot = F)
plot(mod_div2_cl_beta_r2_res) # residuals normally distributed 
testDispersion(mod_div2_cl_beta_r2_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_cl_beta_r2_final <- glmmTMB(diversity_value ~ patch_type + (1|refshort), 
                                      REML = TRUE,
                                      data = raref_cl_beta_q2_glmm)
summary(mod1_div2_cl_beta_r2_final) # report results

# Fit models with habitat amount as a random variable
mod1_div2_cl_beta_r2_ha_final <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                   control=glmmTMBControl(optimizer=optim,
                                                          optArgs=list(method="L-BFGS-B")),
                                   REML = TRUE,
                                   data = raref_cl_beta_q2_glmm)

summary(mod1_div2_cl_beta_r2_ha_final) # report results 

### Extracting Slopes -----------------------------------------------------------
mod1_div2_cl_beta_r2_ha <- glmmTMB(diversity_value ~ patch_type + (habitat_amount|refshort), 
                                   control=glmmTMBControl(optimizer=optim,
                                                          optArgs=list(method="L-BFGS-B")),
                                   data = raref_cl_beta_q2_glmm)


ha_slope_div2_cl_beta_r2 <- ranef(mod1_div2_cl_beta_r2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_beta_r2 <- lm(ha_slope_div2_cl_beta_r2 ~ ha_class, data = moders_div2_cl_beta_r2)
mod_hac_div2_cl_beta_r2_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_beta_r2, plot = F)
plot(mod_hac_div2_cl_beta_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_beta_r2_res) # no overdispersion

summary(mod_hac_div2_cl_beta_r2) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_beta_r2$effort_diff_cent <- scale(log(moders_div2_cl_beta_r2$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_beta_r2 <- lm(ha_slope_div2_cl_beta_r2 ~ effort_diff_cent, data = moders_div2_cl_beta_r2)

mod_eff_div2_cl_beta_r2_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_beta_r2, plot = F)
plot(mod_eff_div2_cl_beta_r2_res) # residuals normally distributed
testDispersion(mod_eff_div2_cl_beta_r2_res) # no overdispersion


summary(mod_eff_div2_cl_beta_r2) # report results


# time since fragmentation

mod_time_div2_cl_beta_r2 <- lm(ha_slope_div2_cl_beta_r2 ~ time_since_fragmentation, data = moders_div2_cl_beta_r2)
mod_time_div2_cl_beta_r2_res <- simulateResiduals(fittedModel = mod_time_div2_cl_beta_r2, plot = F)
plot(mod_time_div2_cl_beta_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_beta_r2_res) # no overdispersion

summary(mod_time_div2_cl_beta_r2) # report results

# continent 

mod_cont_div2_cl_beta_r2 <- lm(ha_slope_div2_cl_beta_r2 ~ continent_cat, data = moders_div2_cl_beta_r2)
mod_cont_div2_cl_beta_r2_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_beta_r2, plot = F)
plot(mod_cont_div2_cl_beta_r2_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_beta_r2_res) # no overdispersion

summary(mod_cont_div2_cl_beta_r2)  # report results

## save the results to export later

lm_results(mod_hac_div2_cl_beta_r2, mod_eff_div2_cl_beta_r2, mod_time_div2_cl_beta_r2, mod_cont_div2_cl_beta_r2) %>% 
  dplyr::mutate(Diversity_type = "Beta (nearest pairs) - q = 2", .before = "Included_moderator") -> moder_results_div2_cl_beta_r2


## Effects on Gamma-Diversity (q = 2) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter gamma-diversity data
frag_cont_rar_div2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_index, q_order, diversity_value, habitat_amount, number_patches, n_pairs) %>% 
  filter(diversity_index == "gamma" & q_order == "q = 2") -> raref_cl_gamma_q2_glmm

# Prepare moderators for gamma-diversity
raref_cl_gamma_q2_glmm %>% 
  dplyr::select(refshort, patch_type, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = n_pairs,
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(effort_diff = abs(fragmented_n_pairs-continuous_n_pairs)) %>% 
  dplyr::select(refshort, effort_diff) %>% 
  left_join(., time_since, by = "refshort") %>% 
  left_join(., study_continents, by = "refshort") %>% 
  left_join(., land_hab_amount_class, by = "refshort") %>% 
  dplyr::mutate(continent_cat = ifelse(continent == "South America", "South America", "Other")) -> moders_div2_cl_gamma_r2

### Model Fitting ---------------------------------------------------------------
# Fit initial model
mod1_div2_cl_gamma_r2 <- glmmTMB(diversity_value ~ patch_type + (1|refshort), data = raref_cl_gamma_q2_glmm)
summary(mod1_div2_cl_gamma_r2)

### Model Checking --------------------------------------------------------------
# Simulate residuals and check assumptions

mod_div2_cl_gamma_r2_res <- simulateResiduals(fittedModel = mod1_div2_cl_gamma_r2, plot = F)
plot(mod_div2_cl_gamma_r2_res) # residuals not normally distributed and heterogeneous variance - we need to log transform it
testDispersion(mod_div2_cl_gamma_r2_res) # no overdispersion

### Model Re-fitting ------------------------------------------------------------
# Fit model with log-transformed response
mod1_div2_cl_gamma_r2_log <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), data = raref_cl_gamma_q2_glmm)
summary(mod1_div2_cl_gamma_r2_log)

# Check the transformed model
mod_div2_cl_gamma_r2_log_res <- simulateResiduals(fittedModel = mod1_div2_cl_gamma_r2_log, plot = F)
plot(mod_div2_cl_gamma_r2_log_res) # residuals normally distributed  
testDispersion(mod_div2_cl_gamma_r2_log_res) # no overdispersion

### Final Model Fitting ---------------------------------------------------------
# Refit final model with REML
mod1_div2_cl_gamma_r2_final <- glmmTMB(log(diversity_value) ~ patch_type + (1|refshort), 
                                       REML = TRUE,
                                       data = raref_cl_gamma_q2_glmm)
summary(mod1_div2_cl_gamma_r2_final) # report results


# Fit models with habitat amount as a random variable
mod1_div2_cl_gamma_r2_ha_final <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                    control=glmmTMBControl(optimizer=optim,
                                                           optArgs=list(method="L-BFGS-B")),
                                    REML = TRUE,
                                    data = raref_cl_gamma_q2_glmm)

summary(mod1_div2_cl_gamma_r2_ha_final) # report results

### Extracting Slopes -----------------------------------------------------------
mod1_div2_cl_gamma_r2_ha <- glmmTMB(log(diversity_value) ~ patch_type + (habitat_amount|refshort), 
                                    control=glmmTMBControl(optimizer=optim,
                                                           optArgs=list(method="L-BFGS-B")),
                                    data = raref_cl_gamma_q2_glmm)


ha_slope_div2_cl_gamma_r2 <- ranef(mod1_div2_cl_gamma_r2_ha)$cond$refshort$habitat_amount

### Effects of Moderators -------------------------------------------------------
# Habitat amount class 

mod_hac_div2_cl_gamma_r2 <- lm(ha_slope_div2_cl_gamma_r2 ~ ha_class, data = moders_div2_cl_gamma_r2)
mod_hac_div2_cl_gamma_r2_res <- simulateResiduals(fittedModel = mod_hac_div2_cl_gamma_r2, plot = F)
plot(mod_hac_div2_cl_gamma_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_gamma_r2_res) # no overdispersion

summary(mod_hac_div2_cl_gamma_r2) # report results


# Difference in the sampling effort between fragmented and continuous landscapes

moders_div2_cl_gamma_r2$effort_diff_cent <- scale(log(moders_div2_cl_gamma_r2$effort_diff+1), center = TRUE, scale = FALSE) 

mod_eff_div2_cl_gamma_r2 <- lm(ha_slope_div2_cl_gamma_r2 ~ effort_diff_cent, data = moders_div2_cl_gamma_r2)

mod_eff_div2_cl_gamma_r2_res <- simulateResiduals(fittedModel = mod_eff_div2_cl_gamma_r2, plot = F)
plot(mod_eff_div2_cl_gamma_r2_res) # residuals normally distributed
testDispersion(mod_eff_div2_cl_gamma_r2_res) # no overdispersion


summary(mod_eff_div2_cl_gamma_r2) # report results


# time since fragmentation

mod_time_div2_cl_gamma_r2 <- lm(ha_slope_div2_cl_gamma_r2 ~ time_since_fragmentation, data = moders_div2_cl_gamma_r2)
mod_time_div2_cl_gamma_r2_res <- simulateResiduals(fittedModel = mod_time_div2_cl_gamma_r2, plot = F)
plot(mod_time_div2_cl_gamma_r2_res) # residuals normally distributed
testDispersion(mod_hac_div2_cl_gamma_r2_res) # no overdispersion

summary(mod_time_div2_cl_gamma_r2) # report results 

# continent 

mod_cont_div2_cl_gamma_r2 <- lm(ha_slope_div2_cl_gamma_r2 ~ continent_cat, data = moders_div2_cl_gamma_r2)
mod_cont_div2_cl_gamma_r2_res <- simulateResiduals(fittedModel = mod_cont_div2_cl_gamma_r2, plot = F)
plot(mod_cont_div2_cl_gamma_r2_res) # residuals normally distributed
testDispersion(mod_cont_div2_cl_gamma_r2_res) # no overdispersion

summary(mod_cont_div2_cl_gamma_r2)  # report results


## save the results to export later

lm_results(mod_hac_div2_cl_gamma_r2, mod_eff_div2_cl_gamma_r2, mod_time_div2_cl_gamma_r2, mod_cont_div2_cl_gamma_r2) %>% 
  dplyr::mutate(Diversity_type = "Gamma (nearest pairs) - q = 2", .before = "Included_moderator") -> moder_results_div2_cl_gamma_r2



#### EXPORT ALL RESULTS - main models and habitat amount models with moderators 

## Main models
## study as random variable
extracted_summaries_study <- extract_glmmTMB_summary(mod1_div2_alpha_log_final, mod1_div2_beta_final, mod1_div2_gamma_log_final, 
                                               mod1_div2_alpha_r0_log_final, mod1_div2_beta_r0_final, mod1_div2_gamma_r0_log_final,
                                               mod1_div2_alpha_r2_log_final, mod1_div2_beta_r2_final, mod1_div2_gamma_r2_final,
                                               mod1_div2_cl_alpha_2_log_final, mod1_div2_cl_beta_final, mod1_div2_cl_gamma_final,
                                               mod1_div2_cl_alpha_r0_log_final, mod1_div2_cl_beta_r0_final, mod1_div2_cl_gamma_r0_log_final,
                                               mod1_div2_cl_alpha_r2_log_final, mod1_div2_cl_beta_r2_final, mod1_div2_cl_gamma_r2_final,
                                               model_names = c("Alpha (all pairs)", "Beta (all pairs)", "Gamma (all pairs)",
                                                               "Alpha (all pairs) - q = 0", "Beta (all pairs) - q = 0", "Gamma (all pairs) - q = 0",
                                                               "Alpha (all pairs) - q = 2", "Beta (all pairs) - q = 2", "Gamma (all pairs) - q = 2",
                                                               "Alpha (nearest pairs)", "Beta (nearest pairs)",  "Gamma (nearest pairs)", 
                                                               "Alpha (nearest pairs) - q = 0", "Beta (nearest pairs) - q = 0", "Gamma (nearest pairs) - q = 0",
                                                               "Alpha (nearest pairs) - q = 2", "Beta (nearest pairs) - q = 2",  "Gamma (nearest pairs) - q = 2"))


## habitat amoubt and study as random variable
extracted_summaries_hamount <- extract_glmmTMB_summary(mod1_div2_alpha_ha_final, mod1_div2_beta_ha_final, mod1_div2_gamma_ha_final,
                                                       mod1_div2_alpha_r0_ha_final, mod1_div2_beta_r0_ha_final, mod1_div2_gamma_r0_ha_final,
                                                       mod1_div2_alpha_r2_ha_final, mod1_div2_beta_r2_ha_final, mod1_div2_gamma_r2_ha_final,
                                                       mod1_div2_cl_alpha_ha_final, mod1_div2_cl_beta_ha_final, mod1_div2_cl_gamma_ha_final,
                                                       mod1_div2_cl_alpha_r0_ha_final, mod1_div2_cl_beta_r0_ha_final, mod1_div2_cl_gamma_r0_ha_final,
                                                       mod1_div2_cl_alpha_r2_ha_final, mod1_div2_cl_beta_r2_ha_final, mod1_div2_cl_gamma_r2_ha_final,
                                                     model_names = c("Alpha (all pairs)", "Beta (all pairs)", "Gamma (all pairs)",
                                                                     "Alpha (all pairs) - q = 0", "Beta (all pairs) - q = 0", "Gamma (all pairs) - q = 0",
                                                                     "Alpha (all pairs) - q = 2", "Beta (all pairs) - q = 2", "Gamma (all pairs) - q = 2",
                                                                     "Alpha (nearest pairs)", "Beta (nearest pairs)",  "Gamma (nearest pairs)", 
                                                                     "Alpha (nearest pairs) - q = 0", "Beta (nearest pairs) - q = 0", "Gamma (nearest pairs) - q = 0",
                                                                     "Alpha (nearest pairs) - q = 2", "Beta (nearest pairs) - q = 2",  "Gamma (nearest pairs) - q = 2"))




## habitat amount models with moderators 

moderator_results <- rbind.data.frame(
  moder_results_div2_alpha,
  moder_results_div2_beta,
  moder_results_div2_gamma,
  moder_results_div2_alpha_r0,
  moder_results_div2_beta_r0,
  moder_results_div2_gamma_r0,
  moder_results_div2_alpha_r2,
  moder_results_div2_beta_r2,
  moder_results_div2_gamma_r2,
  moder_results_div2_cl_alpha,
  moder_results_div2_cl_beta,
  moder_results_div2_cl_gamma,
  moder_results_div2_cl_alpha_r0,
  moder_results_div2_cl_beta_r0,
  moder_results_div2_cl_gamma_r0,
  moder_results_div2_cl_alpha_r2,
  moder_results_div2_cl_beta_r2,
  moder_results_div2_cl_gamma_r2
) %>% 
  dplyr::mutate(moderator = ifelse(grepl("eff",Included_moderator), "Number of pairs",
                                   ifelse(grepl("hac",Included_moderator), "Habitat amount",
                                          ifelse(grepl("time",Included_moderator), "Time after fragmentation", "Continent"))), 
                .before=Included_moderator) %>% 
  dplyr::select(-Included_moderator)

rownames(moderator_results) <- seq_len(nrow(moderator_results))

# The two files above are already exported in the folder results 

# write.csv(extracted_summaries_study, "results/extracted_summaries_study.csv")
# write.csv(extracted_summaries_hamount, "results/extracted_summaries_hamount.csv")
# write.csv(moderator_results, "results/moderator_results.csv")



