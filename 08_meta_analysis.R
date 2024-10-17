# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: Meta-analysis to test the effects of landscape type on species diversity
# Author: Thiago Gonçalves Souza
# Date: [October 10, 2024]
# Notes: R version 4.4.1 [2024-06-14 ucrt]
# ==============================================================================


# Using renv for package management
renv::restore() # it requires installing the package renv

# Packages ---------------------------------------------------------------------

library(tidyverse)
library(metafor)
library(ggplot2)
# library(orchaRd) # devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(gtools)
library(ggbeeswarm)
# library(Rmisc)
source("utility_functions.R")

# Data -------------------------------------------------------------------------

# continents 

study_continents <- read.csv("data/landfrag_37sset_continents.csv")

# time since fragmentation

time_since <- read.csv("data/landfrag_37sset_time_since_frag.csv")

## diversity of 2 based on all fragment / plot pairs 

frag_cont_div2 <- read.csv("processed_data/diversity_of_2.csv") 

# diversity of 2 based only on the nearest pairs 

frag_cont_div2_cl <- read.csv("processed_data/diversity_of_2_close.csv") 

# rarefied diversity of 2 - all pairs 

frag_cont_rare_div2 <- read.csv("processed_data/frag_cont_raref_2.csv") 

# rarefied diversity of 2 - nearest pairs 

frag_cont_rare_div2_cl <- read.csv("processed_data/frag_cont_raref_close_2.csv") 

### Data preparation -----------------------------------------------------------

### landscape data to be merged in all diversity matrices

frag_cont_div2 %>% 
  dplyr::filter(diversity_index == "alpha") %>% 
  dplyr::select(refshort, habitat_amount) %>% 
  dplyr::group_by(refshort) %>% 
  dplyr::summarise(aver_ha = mean(habitat_amount)) %>% 
  dplyr::mutate(ha_class = quantcut(aver_ha, 2)) %>% 
  dplyr::select(refshort, ha_class) -> land_hab_amount_class


landscape_info <- frag_cont_div2 %>% 
  dplyr::filter(diversity_index == "alpha") %>% 
  dplyr::select(refshort, patch_type, habitat_amount, number_patches) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(habitat_amount, number_patches),
              names_glue = "{patch_type}_{.value}") %>% 
  dplyr::mutate(ha_diff = log(continuous_habitat_amount) - log(fragmented_habitat_amount),
                np_diff = log(continuous_number_patches) - log(fragmented_number_patches)) %>% 
  dplyr::select(refshort, ha_diff, np_diff)


landscape_info <- left_join(landscape_info, study_continents, by = "refshort")
landscape_info <- left_join(landscape_info, land_hab_amount_class, by = "refshort")
landscape_info <- left_join(landscape_info, time_since, by = "refshort")


### Meta-analysis --------------------------------------------------------------
### 1. models using all plot pairs ###

## Effects on Alpha-Diversity --------------------------------------------------

### Data Preparation -----------------------------------------------------------
# Filter alpha-diversity data

frag_cont_div2 %>% 
  filter(diversity_index == "alpha") -> alpha2

transformed_data_alpha2 <- alpha2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011) 
alpha2_logR <- escalc(measure="ROM", 
                      m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                      sd1i = continuous_sd, sd2i = fragmented_sd,
                      n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                      data = transformed_data_alpha2)


# combine matrices to add habitat amount

alpha2_logR <- left_join(alpha2_logR, landscape_info, by = "refshort")


### Model Fitting --------------------------------------------------------------
# Meta-Analysis via Linear (Mixed-Effects) Models 

mod_overall_alpha2 <- rma.uni(yi, vi,
                              data = alpha2_logR,
                              method = "REML") # with random effect 

summary(mod_overall_alpha2) # reported results 

# tables for plotting later

overall_alpha2_res <- overall_rma_results(mod_overall_alpha2)
overall_alpha2_est <- overall_alpha2_res %>% 
  dplyr::mutate(diversity = "Alpha", 
                type = "All pairs",
                .before = Estimate)
overall_alpha2_dat <- mod_overall_alpha2$data %>% 
  dplyr::mutate(diversity = "Alpha",
                type = "All pairs",
                .before = yi)

### Effects of Moderators -------------------------------------------------------
## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_alpha2 <- rma.uni(yi, vi,
                            data = alpha2_logR,
                            mods = ~ ha_diff,
                            method = "REML") # with random effect 

summary(mod_mods1_alpha2) # reported results 

## Habitat amount (quantiles)

mod_mods2_alpha2 <- rma.uni(yi, vi,
                            data = alpha2_logR,
                            mods = ~ ha_class -1,
                            method = "REML") # with random effect 

summary(mod_mods2_alpha2) # reported results 

## Comparing South America with other continents 

alpha2_logR$continent_bin <- ifelse(alpha2_logR$continent == "South America", "South America", "Others")


mod_mods3_alpha2 <- rma.uni(yi, vi,
                            data = alpha2_logR,
                            mods = ~ continent_bin - 1,
                            method = "REML") # with random effect 

summary(mod_mods3_alpha2) # reported results 

## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_alpha2 <- rma.uni(yi, vi,
                            data = alpha2_logR,
                            mods = ~ time_since_fragmentation - 1,
                            method = "REML") # with random effect 

summary(mod_mods4_alpha2) # reported results 

### Sensitivity Analysis -------------------------------------------------------
# influential study diagnostics 

mod_infl_alpha2 <- influence.rma.uni(mod_overall_alpha2)
plot(mod_infl_alpha2) # study 26 is an outlier

# double check if the results are the same by removing the study 26

mod_overall_alpha_rem26 <- rma.uni(yi, vi,
                              data = alpha2_logR[-26,],
                              method = "REML") 

summary(mod_overall_alpha_rem26) # same result. The estimate is even more strongly positive

# calculate the inverse of effective sample size
alpha2_logR$inv_n <- with(alpha2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
alpha2_logR$sqrt_inv_n <- with(alpha2_logR, sqrt(inv_n))

# Time-lag bias test

alpha2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", alpha2_logR$refshort))
alpha2_logR$year_c <- as.vector(scale(alpha2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_alpha2_check <- rma(yi, vi,
                        mods = ~1 + sqrt_inv_n + year_c,
                        data = alpha2_logR,
                        method = "REML")
summary(mod_alpha2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_alpha2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Alpha", 
         type = "All pairs",
         refshort = alpha2_logR$refshort,
         continent_bin = alpha2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_alpha2
leave1out_mod_alpha2

### Output preparation ---------------------------------------------------------
## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type

results_mods1_alpha2 <- extract_model_results(mod_mods1_alpha2, 
                                              diversity_index = "Alpha", 
                                              diversity_type = "All pairs") 


results_mods2_alpha2 <-extract_model_results(mod_mods2_alpha2, 
                                             diversity_index = "Alpha", 
                                             diversity_type = "All pairs") 


results_mods3_alpha2 <- extract_model_results(mod_mods3_alpha2, 
                                              diversity_index = "Alpha", 
                                              diversity_type = "All pairs") 



results_mods4_alpha2 <- extract_model_results(mod_mods4_alpha2, 
                                              diversity_index = "Alpha", 
                                              diversity_type = "All pairs") 

results_alpha2_check <- extract_model_results(mod_alpha2_check, 
                                              diversity_index = "Alpha", 
                                              diversity_type = "All pairs") 


mod_results_alpha2 <- rbind.data.frame(results_mods1_alpha2, 
                                       results_mods2_alpha2, 
                                       results_mods3_alpha2,
                                       results_mods4_alpha2,
                                       results_alpha2_check)

## Effects on Beta-Diversity --------------------------------------------------
### Beta diversity (all pairs)

### Data Preparation ------------------------------------------------------------
# Filter beta-diversity data

frag_cont_div2 %>% 
  filter(diversity_index == "beta") -> beta2

transformed_data_beta2 <- beta2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011)

beta2_logR <- escalc(measure="ROM", 
                     m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                     sd1i = continuous_sd, sd2i = fragmented_sd,
                     n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                     data = transformed_data_beta2)


# Add habitat amount information

beta2_logR <- left_join(beta2_logR, landscape_info, by = "refshort")

### Model Fitting --------------------------------------------------------------
# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_beta2 <- rma.uni(yi, vi,
                             data = beta2_logR,
                             method = "REML") # with random effect 

summary(mod_overall_beta2) # reported results 

# Tables for plotting later

overall_beta2_res <- overall_rma_results(mod_overall_beta2)
overall_beta2_est <- overall_beta2_res %>% 
  dplyr::mutate(diversity = "Beta",
                type = "All pairs",
                .before = Estimate)
overall_beta2_dat <- mod_overall_beta2$data %>% 
  dplyr::mutate(diversity = "Beta", 
                type = "All pairs",
                .before = yi)

### Effects of Moderators -------------------------------------------------------
## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_beta2 <- rma.uni(yi, vi,
                           data = beta2_logR,
                           mods = ~ ha_diff,
                           method = "REML") # with random effect 

summary(mod_mods1_beta2) # reported results 

## Habitat amount (quantiles)

mod_mods2_beta2 <- rma.uni(yi, vi,
                           data = beta2_logR,
                           mods = ~ ha_class -1,
                           method = "REML") # with random effect 

summary(mod_mods2_beta2) # reported results 

## Comparing South America with other continents 

beta2_logR$continent_bin <- ifelse(beta2_logR$continent == "South America", "South America", "Others")

mod_mods3_beta2 <- rma.uni(yi, vi,
                           data = beta2_logR,
                           mods = ~ continent_bin - 1,
                           method = "REML") # with random effect 

summary(mod_mods3_beta2) # reported results 

## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_beta2 <- rma.uni(yi, vi,
                           data = beta2_logR,
                           mods = ~ time_since_fragmentation - 1,
                           method = "REML") # with random effect 

summary(mod_mods4_beta2) # reported results 

### Sensitivity Analysis -------------------------------------------------------
# influential study diagnostics 
mod_infl_beta2 <- influence.rma.uni(mod_overall_beta2)
plot(mod_infl_beta2) # no outlier

# calculate the inverse of effective sample size
beta2_logR$inv_n <- with(beta2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
beta2_logR$sqrt_inv_n <- with(beta2_logR, sqrt(inv_n))

# Time-lag bias test

beta2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", beta2_logR$refshort))
beta2_logR$year_c <- as.vector(scale(beta2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_beta2_check <- rma(yi, vi,
                       mods = ~1 + sqrt_inv_n + year_c,
                       data = beta2_logR,
                       method = "REML")
summary(mod_beta2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out analysis
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_beta2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Beta", 
         type = "All pairs",
         refshort = beta2_logR$refshort,
         continent_bin = beta2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_beta2
leave1out_mod_beta2

### Output Preparation ---------------------------------------------------------
# Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_beta2 <- extract_model_results(mod_mods1_beta2, 
                                             diversity_index = "Beta", 
                                             diversity_type = "All pairs") 


results_mods2_beta2 <-extract_model_results(mod_mods2_beta2, 
                                            diversity_index = "Beta", 
                                            diversity_type = "All pairs") 


results_mods3_beta2 <- extract_model_results(mod_mods3_beta2, 
                                             diversity_index = "beta", 
                                             diversity_type = "All pairs") 



results_mods4_beta2 <- extract_model_results(mod_mods4_beta2, 
                                             diversity_index = "Beta", 
                                             diversity_type = "All pairs") 

results_beta2_check <- extract_model_results(mod_beta2_check, 
                                             diversity_index = "Beta", 
                                             diversity_type = "All pairs") 


mod_results_beta2 <- rbind.data.frame(results_mods1_beta2, 
                                      results_mods2_beta2, 
                                      results_mods3_beta2,
                                      results_mods4_beta2,
                                      results_beta2_check)


## Effects on Gamma-Diversity --------------------------------------------------


### Data Preparation ------------------------------------------------------------
# Filter gamma-diversity data
frag_cont_div2 %>% 
  filter(diversity_index == "gamma") -> gamma2

transformed_data_gamma2 <- gamma2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011)

gamma2_logR <- escalc(measure="ROM", 
                      m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                      sd1i = continuous_sd, sd2i = fragmented_sd,
                      n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                      data = transformed_data_gamma2)


# Add habitat amount information

gamma2_logR <- left_join(gamma2_logR, landscape_info, by = "refshort")


### Model Fitting --------------------------------------------------------------
# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_gamma2 <- rma.uni(yi, vi,
                              data = gamma2_logR,
                              method = "REML") # with random effect 

summary(mod_overall_gamma2) # reported results 

# tables for plotting later

overall_gamma2_res <- overall_rma_results(mod_overall_gamma2)
overall_gamma2_est <- overall_gamma2_res %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "All pairs",
                .before = Estimate)
overall_gamma2_dat <- mod_overall_gamma2$data %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "All pairs",
                .before = yi)

### Effects of Moderators -------------------------------------------------------
## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_gamma2 <- rma.uni(yi, vi,
                            data = gamma2_logR,
                            mods = ~ ha_diff,
                            method = "REML") # with random effect 

summary(mod_mods1_gamma2) # reported results 

## Habitat amount (quantiles)

mod_mods2_gamma2 <- rma.uni(yi, vi,
                            data = gamma2_logR,
                            mods = ~ ha_class -1,
                            method = "REML") # with random effect 

summary(mod_mods2_gamma2) # reported results 

## Comparing South America with other continents 

gamma2_logR$continent_bin <- ifelse(gamma2_logR$continent == "South America", "South America", "Others")

mod_mods3_gamma2 <- rma.uni(yi, vi,
                            data = gamma2_logR,
                            mods = ~ continent_bin - 1,
                            method = "REML") # with random effect 

summary(mod_mods3_gamma2) # reported results 

## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_gamma2 <- rma.uni(yi, vi,
                            data = gamma2_logR,
                            mods = ~ time_since_fragmentation - 1,
                            method = "REML") # with random effect 

summary(mod_mods4_gamma2) # reported results 

### Sensitivity Analysis -------------------------------------------------------
# influential study diagnostics 

mod_infl_gamma2 <- influence.rma.uni(mod_overall_gamma2)
plot(mod_infl_gamma2) # no outlier

# calculate the inverse of effective sample size
gamma2_logR$inv_n <- with(gamma2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
gamma2_logR$sqrt_inv_n <- with(gamma2_logR, sqrt(inv_n))

# Time-lag bias test

gamma2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", gamma2_logR$refshort))
gamma2_logR$year_c <- as.vector(scale(gamma2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_gamma2_check <- rma(yi, vi,
                        mods = ~1 + sqrt_inv_n + year_c,
                        data = gamma2_logR,
                        method = "REML")
summary(mod_gamma2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_gamma2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Gamma", 
         type = "All pairs",
         refshort = gamma2_logR$refshort,
         continent_bin = gamma2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_gamma2
leave1out_mod_gamma2


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_gamma2 <- extract_model_results(mod_mods1_gamma2, 
                                              diversity_index = "Gamma", 
                                              diversity_type = "All pairs") 


results_mods2_gamma2 <-extract_model_results(mod_mods2_gamma2, 
                                             diversity_index = "Gamma", 
                                             diversity_type = "All pairs") 


results_mods3_gamma2 <- extract_model_results(mod_mods3_gamma2, 
                                              diversity_index = "Gamma", 
                                              diversity_type = "All pairs") 



results_mods4_gamma2 <- extract_model_results(mod_mods4_gamma2, 
                                              diversity_index = "Gamma", 
                                              diversity_type = "All pairs") 

results_gamma2_check <- extract_model_results(mod_gamma2_check, 
                                              diversity_index = "Gamma", 
                                              diversity_type = "All pairs") 


mod_results_gamma2 <- rbind.data.frame(results_mods1_gamma2, 
                                       results_mods2_gamma2, 
                                       results_mods3_gamma2,
                                       results_mods4_gamma2,
                                       results_gamma2_check)

### Meta-analysis --------------------------------------------------------------
### 2. models using the nearest plot pairs ###

## Effects on Alpha-Diversity --------------------------------------------------

frag_cont_div2_cl %>% 
  filter(diversity_index == "alpha") -> alpha2_cl

transformed_data_alpha2_cl <- alpha2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011)
alpha2_cl_logR <- escalc(measure="ROM", 
                         m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                         sd1i = continuous_sd, sd2i = fragmented_sd,
                         n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                         data = transformed_data_alpha2_cl)


# Add habitat amount information
alpha2_cl_logR <- left_join(alpha2_cl_logR, landscape_info, by = "refshort")

### Model Fitting --------------------------------------------------------------
# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_alpha2_cl <- rma.uni(yi, vi,
                                 data = alpha2_cl_logR,
                                 method = "REML") # with random effect 

summary(mod_overall_alpha2_cl) # reported results 

# Tables for plotting later

overall_alpha2_cl_res <- overall_rma_results(mod_overall_alpha2_cl)
overall_alpha2_cl_est <- overall_alpha2_cl_res %>% 
  dplyr::mutate(diversity = "Alpha",
                type = "Nearest pairs",
                .before = Estimate)
overall_alpha2_cl_dat <- mod_overall_alpha2_cl$data %>% 
  dplyr::mutate(diversity = "Alpha",
                type = "Nearest pairs",
                .before = yi)

### Effects of Moderators -------------------------------------------------------
## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_alpha2_cl <- rma.uni(yi, vi,
                               data = alpha2_cl_logR,
                               mods = ~ ha_diff,
                               method = "REML") # with random effect 

summary(mod_mods1_alpha2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_alpha2_cl <- rma.uni(yi, vi,
                               data = alpha2_cl_logR,
                               mods = ~ ha_class -1,
                               method = "REML") # with random effect 

summary(mod_mods2_alpha2_cl) # reported results 

## Comparing South America with other continents 

alpha2_cl_logR$continent_bin <- ifelse(alpha2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_alpha2_cl <- rma.uni(yi, vi,
                               data = alpha2_cl_logR,
                               mods = ~ continent_bin - 1,
                               method = "REML") # with random effect 

summary(mod_mods3_alpha2_cl) # reported results 

## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_alpha2_cl <- rma.uni(yi, vi,
                               data = alpha2_cl_logR,
                               mods = ~ time_since_fragmentation - 1,
                               method = "REML") # with random effect 

summary(mod_mods4_alpha2_cl) # reported results 

### Sensitivity Analysis -------------------------------------------------------

# influential study diagnostics 

mod_infl_alpha2_cl <- influence.rma.uni(mod_overall_alpha2_cl)
plot(mod_infl_alpha2_cl) # no outlier

# calculate the inverse of effective sample size
alpha2_cl_logR$inv_n <- with(alpha2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
alpha2_cl_logR$sqrt_inv_n <- with(alpha2_cl_logR, sqrt(inv_n))

# Time-lag bias test

alpha2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", alpha2_cl_logR$refshort))
alpha2_cl_logR$year_c <- as.vector(scale(alpha2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_alpha2_cl_check <- rma(yi, vi,
                           mods = ~1 + sqrt_inv_n + year_c,
                           data = alpha2_cl_logR,
                           method = "REML")
summary(mod_alpha2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_alpha2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Alpha", 
         type = "Nearest pairs",
         refshort = alpha2_cl_logR$refshort,
         continent_bin = alpha2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_alpha2_cl
leave1out_mod_alpha2_cl

### Output Preparation ---------------------------------------------------------
## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type
results_mods1_alpha2_cl <- extract_model_results(mod_mods1_alpha2_cl, 
                                                 diversity_index = "Alpha", 
                                                 diversity_type = "Nearest pairs") 


results_mods2_alpha2_cl <-extract_model_results(mod_mods2_alpha2_cl, 
                                                diversity_index = "Alpha", 
                                                diversity_type = "Nearest pairs") 


results_mods3_alpha2_cl <- extract_model_results(mod_mods3_alpha2_cl, 
                                                 diversity_index = "Alpha", 
                                                 diversity_type = "Nearest pairs") 



results_mods4_alpha2_cl <- extract_model_results(mod_mods4_alpha2_cl, 
                                                 diversity_index = "Alpha", 
                                                 diversity_type = "Nearest pairs") 

results_alpha2_cl_check <- extract_model_results(mod_alpha2_cl_check, 
                                                 diversity_index = "Alpha", 
                                                 diversity_type = "Nearest pairs") 


mod_results_alpha2_cl <- rbind.data.frame(results_mods1_alpha2_cl, 
                                          results_mods2_alpha2_cl, 
                                          results_mods3_alpha2_cl,
                                          results_mods4_alpha2_cl,
                                          results_alpha2_cl_check)



## Effects on Beta-Diversity --------------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter beta-diversity data

frag_cont_div2_cl %>% 
  filter(diversity_index == "beta") -> beta2_cl

transformed_data_beta2_cl <- beta2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011)
beta2_cl_logR <- escalc(measure="ROM", 
                        m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                        sd1i = continuous_sd, sd2i = fragmented_sd,
                        n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                        data = transformed_data_beta2_cl)


# Add habitat amount information
beta2_cl_logR <- left_join(beta2_cl_logR, landscape_info, by = "refshort")

### Model Fitting --------------------------------------------------------------
# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_beta2_cl <- rma.uni(yi, vi,
                                data = beta2_cl_logR,
                                method = "REML") # with random effect 

summary(mod_overall_beta2_cl) # reported results 

# Tables for plotting later

overall_beta2_cl_res <- overall_rma_results(mod_overall_beta2_cl)
overall_beta2_cl_est <- overall_beta2_cl_res %>% 
  dplyr::mutate(diversity = "Beta", 
                type = "Nearest pairs",
                .before = Estimate)
overall_beta2_cl_dat <- mod_overall_beta2_cl$data %>% 
  dplyr::mutate(diversity = "Beta", 
                type = "Nearest pairs",
                .before = yi)

### Effects of Moderators -------------------------------------------------------
## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_beta2_cl <- rma.uni(yi, vi,
                              data = beta2_cl_logR,
                              mods = ~ ha_diff,
                              method = "REML") # with random effect 

summary(mod_mods1_beta2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_beta2_cl <- rma.uni(yi, vi,
                              data = beta2_cl_logR,
                              mods = ~ ha_class -1,
                              method = "REML") # with random effect 

summary(mod_mods2_beta2_cl) # reported results 

## Comparing South America with other continents 

beta2_cl_logR$continent_bin <- ifelse(beta2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_beta2_cl <- rma.uni(yi, vi,
                              data = beta2_cl_logR,
                              mods = ~ continent_bin - 1,
                              method = "REML") # with random effect 

summary(mod_mods3_beta2_cl) # reported results 

## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_beta2_cl <- rma.uni(yi, vi,
                              data = beta2_cl_logR,
                              mods = ~ time_since_fragmentation - 1,
                              method = "REML") # with random effect 

summary(mod_mods4_beta2_cl) # reported results 

### Sensitivity Analysis -------------------------------------------------------
# influential study diagnostics 

mod_infl_beta2_cl <- influence.rma.uni(mod_overall_beta2_cl)
plot(mod_infl_beta2_cl) # no outlier

# calculate the inverse of effective sample size
beta2_cl_logR$inv_n <- with(beta2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
beta2_cl_logR$sqrt_inv_n <- with(beta2_cl_logR, sqrt(inv_n))

# Time-lag bias test

beta2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", beta2_cl_logR$refshort))
beta2_cl_logR$year_c <- as.vector(scale(beta2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_beta2_cl_check <- rma(yi, vi,
                          mods = ~1 + sqrt_inv_n + year_c,
                          data = beta2_cl_logR,
                          method = "REML")
summary(mod_beta2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_beta2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Beta", 
         type = "Nearest pairs",
         refshort = beta2_cl_logR$refshort,
         continent_bin = beta2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_beta2_cl
leave1out_mod_beta2_cl

### Output Preparation ---------------------------------------------------------
## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_beta2_cl <- extract_model_results(mod_mods1_beta2_cl, 
                                                diversity_index = "Beta", 
                                                diversity_type = "Nearest pairs") 


results_mods2_beta2_cl <-extract_model_results(mod_mods2_beta2_cl, 
                                               diversity_index = "Beta", 
                                               diversity_type = "Nearest pairs") 


results_mods3_beta2_cl <- extract_model_results(mod_mods3_beta2_cl, 
                                                diversity_index = "Beta", 
                                                diversity_type = "Nearest pairs") 



results_mods4_beta2_cl <- extract_model_results(mod_mods4_beta2_cl, 
                                                diversity_index = "Beta", 
                                                diversity_type = "Nearest pairs") 

results_beta2_cl_check <- extract_model_results(mod_beta2_cl_check, 
                                                diversity_index = "Beta", 
                                                diversity_type = "Nearest pairs") 


mod_results_beta2_cl <- rbind.data.frame(results_mods1_beta2_cl, 
                                         results_mods2_beta2_cl, 
                                         results_mods3_beta2_cl,
                                         results_mods4_beta2_cl,
                                         results_beta2_cl_check)

## Effects on Gamma-Diversity --------------------------------------------------


### Data Preparation ------------------------------------------------------------
# Filter gamma-diversity data
frag_cont_div2_cl %>% 
  filter(diversity_index == "gamma") -> gamma2_cl

transformed_data_gamma2_cl <- gamma2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011)
gamma2_cl_logR <- escalc(measure="ROM", 
                         m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                         sd1i = continuous_sd, sd2i = fragmented_sd,
                         n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                         data = transformed_data_gamma2_cl)

# Add habitat amount information
gamma2_cl_logR <- left_join(gamma2_cl_logR, landscape_info, by = "refshort")

### Model Fitting --------------------------------------------------------------
# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_gamma2_cl <- rma.uni(yi, vi,
                                 data = gamma2_cl_logR,
                                 method = "REML") # with random effect 

summary(mod_overall_gamma2_cl) # reported results 


# Tables for plotting later

overall_gamma2_cl_res <- overall_rma_results(mod_overall_gamma2_cl)
overall_gamma2_cl_est <- overall_gamma2_cl_res %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "Nearest pairs",
                .before = Estimate)
overall_gamma2_cl_dat <- mod_overall_gamma2_cl$data %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "Nearest pairs",
                .before = yi)

### Effects of Moderators -------------------------------------------------------
## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_gamma2_cl <- rma.uni(yi, vi,
                               data = gamma2_cl_logR,
                               mods = ~ ha_diff,
                               method = "REML") # with random effect 

summary(mod_mods1_gamma2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_gamma2_cl <- rma.uni(yi, vi,
                               data = gamma2_cl_logR,
                               mods = ~ ha_class -1,
                               method = "REML") # with random effect 

summary(mod_mods2_gamma2_cl) # reported results 

## Comparing South America with other continents 

gamma2_cl_logR$continent_bin <- ifelse(gamma2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_gamma2_cl <- rma.uni(yi, vi,
                               data = gamma2_cl_logR,
                               mods = ~ continent_bin - 1,
                               method = "REML") # with random effect 

summary(mod_mods3_gamma2_cl) # reported results 

## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_gamma2_cl <- rma.uni(yi, vi,
                               data = gamma2_cl_logR,
                               mods = ~ time_since_fragmentation - 1,
                               method = "REML") # with random effect 

summary(mod_mods4_gamma2_cl) # reported results 


### Sensitivity Analysis -------------------------------------------------------
# influential study diagnostics 

mod_infl_gamma2_cl <- influence.rma.uni(mod_overall_gamma2_cl)
plot(mod_infl_gamma2_cl) # no outlier

# calculate the inverse of effective sample size
gamma2_cl_logR$inv_n <- with(gamma2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
gamma2_cl_logR$sqrt_inv_n <- with(gamma2_cl_logR, sqrt(inv_n))

# Time-lag bias test

gamma2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", gamma2_cl_logR$refshort))
gamma2_cl_logR$year_c <- as.vector(scale(gamma2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_gamma2_cl_check <- rma(yi, vi,
                           mods = ~1 + sqrt_inv_n + year_c,
                           data = gamma2_cl_logR,
                           method = "REML")
summary(mod_gamma2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_gamma2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Gamma", 
         type = "Nearest pairs",
         refshort = gamma2_cl_logR$refshort,
         continent_bin = gamma2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_gamma2_cl
leave1out_mod_gamma2_cl

### Output Preparation ---------------------------------------------------------
## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type
results_mods1_gamma2_cl <- extract_model_results(mod_mods1_gamma2_cl, 
                                                 diversity_index = "Gamma", 
                                                 diversity_type = "Nearest pairs") 


results_mods2_gamma2_cl <-extract_model_results(mod_mods2_gamma2_cl, 
                                                diversity_index = "Gamma", 
                                                diversity_type = "Nearest pairs") 


results_mods3_gamma2_cl <- extract_model_results(mod_mods3_gamma2_cl, 
                                                 diversity_index = "Gamma", 
                                                 diversity_type = "Nearest pairs") 



results_mods4_gamma2_cl <- extract_model_results(mod_mods4_gamma2_cl, 
                                                 diversity_index = "Gamma", 
                                                 diversity_type = "Nearest pairs") 

results_gamma2_cl_check <- extract_model_results(mod_gamma2_cl_check, 
                                                 diversity_index = "Gamma", 
                                                 diversity_type = "Nearest pairs") 


mod_results_gamma2_cl <- rbind.data.frame(results_mods1_gamma2_cl, 
                                          results_mods2_gamma2_cl, 
                                          results_mods3_gamma2_cl,
                                          results_mods4_gamma2_cl,
                                          results_gamma2_cl_check)



### Meta-analysis --------------------------------------------------------------
### 3. models using all plot pairs with rarefied estimates (q=0, q=2) ###

## Effects on Alpha-Diversity (q = 0) ------------------------------------------

### Data Preparation -----------------------------------------------------------
# Filter alpha-diversity data

frag_cont_rare_div2 %>% 
  filter(diversity_index == "alpha" & q_order == "q = 0") -> rare0_alpha2

transformed_data_rare0_alpha2 <- rare0_alpha2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011)

rare0_alpha2_logR <- escalc(measure="ROM", 
                            m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                            sd1i = continuous_sd, sd2i = fragmented_sd,
                            n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                            data = transformed_data_rare0_alpha2)


# Add habitat amount information
rare0_alpha2_logR <- left_join(rare0_alpha2_logR, landscape_info, by = "refshort")

### Model Fitting --------------------------------------------------------------
# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare0_alpha2 <- rma.uni(yi, vi,
                                    data = rare0_alpha2_logR,
                                    method = "REML") # with random effect 

summary(mod_overall_rare0_alpha2) # reported results 

# Tables for plotting later

overall_rare0_alpha2_res <- overall_rma_results(mod_overall_rare0_alpha2)
overall_rare0_alpha2_est <- overall_rare0_alpha2_res %>% 
  dplyr::mutate(diversity = "Alpha", 
                type = "All pairs (q = 0)",
                .before = Estimate)
overall_rare0_alpha2_dat <- mod_overall_rare0_alpha2$data %>% 
  dplyr::mutate(diversity = "Alpha",
                type = "All pairs (q = 0)",
                .before = yi)


### Effects of Moderators -------------------------------------------------------
## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare0_alpha2 <- rma.uni(yi, vi,
                                  data = rare0_alpha2_logR,
                                  mods = ~ ha_diff,
                                  method = "REML") # with random effect 

summary(mod_mods1_rare0_alpha2) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare0_alpha2 <- rma.uni(yi, vi,
                                  data = rare0_alpha2_logR,
                                  mods = ~ ha_class -1,
                                  method = "REML") # with random effect 

summary(mod_mods2_rare0_alpha2) # reported results 

## Comparing South America with other continents 

rare0_alpha2_logR$continent_bin <- ifelse(rare0_alpha2_logR$continent == "South America", "South America", "Others")


mod_mods3_rare0_alpha2 <- rma.uni(yi, vi,
                                  data = rare0_alpha2_logR,
                                  mods = ~ continent_bin - 1,
                                  method = "REML") # with random effect 

summary(mod_mods3_rare0_alpha2) # reported results 

## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare0_alpha2 <- rma.uni(yi, vi,
                                  data = rare0_alpha2_logR,
                                  mods = ~ time_since_fragmentation - 1,
                                  method = "REML") # with random effect 

summary(mod_mods4_rare0_alpha2) # reported results 

### Sensitivity Analysis -------------------------------------------------------
# influential study diagnostics 

mod_infl_rare0_alpha2 <- influence.rma.uni(mod_overall_rare0_alpha2)
plot(mod_infl_rare0_alpha2) # study 10 is an outlier

# double check whether the removel of study 10 affects the results

mod_overall_rare0_alpha2_rem10 <- rma.uni(yi, vi,
                                    data = rare0_alpha2_logR[-10,],
                                    method = "REML") # with random effect 

mod_overall_rare0_alpha2_rem10 # same result. The estimate changed from 0.11 to 0.09 (still significant)

# calculate the inverse of effective sample size
rare0_alpha2_logR$inv_n <- with(rare0_alpha2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare0_alpha2_logR$sqrt_inv_n <- with(rare0_alpha2_logR, sqrt(inv_n))

# Time-lag bias test

rare0_alpha2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare0_alpha2_logR$refshort))
rare0_alpha2_logR$year_c <- as.vector(scale(rare0_alpha2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare0_alpha2_check <- rma(yi, vi,
                              mods = ~1 + sqrt_inv_n + year_c,
                              data = rare0_alpha2_logR,
                              method = "REML")
summary(mod_rare0_alpha2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare0_alpha2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Alpha", 
         type = "All pairs (q = 0)",
         refshort = rare0_alpha2_logR$refshort,
         continent_bin = rare0_alpha2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare0_alpha2
leave1out_mod_rare0_alpha2

### Output Preparation ---------------------------------------------------------
## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type

results_mods1_rare0_alpha2 <- extract_model_results(mod_mods1_rare0_alpha2, 
                                                    diversity_index = "Alpha (q=0)", 
                                                    diversity_type = "All pairs") 


results_mods2_rare0_alpha2 <-extract_model_results(mod_mods2_rare0_alpha2, 
                                                   diversity_index = "Alpha (q=0)", 
                                                   diversity_type = "All pairs") 


results_mods3_rare0_alpha2 <- extract_model_results(mod_mods3_rare0_alpha2, 
                                                    diversity_index = "Alpha (q=0)", 
                                                    diversity_type = "All pairs") 



results_mods4_rare0_alpha2 <- extract_model_results(mod_mods4_rare0_alpha2, 
                                                    diversity_index = "Alpha (q=0)", 
                                                    diversity_type = "All pairs") 

results_rare0_alpha2_check <- extract_model_results(mod_rare0_alpha2_check, 
                                                    diversity_index = "Alpha (q=0)", 
                                                    diversity_type = "All pairs") 


mod_results_rare0_alpha2 <- rbind.data.frame(results_mods1_rare0_alpha2, 
                                             results_mods2_rare0_alpha2, 
                                             results_mods3_rare0_alpha2,
                                             results_mods4_rare0_alpha2,
                                             results_rare0_alpha2_check)



## Effects on Alpha-Diversity (q = 2) ------------------------------------------

### Data Preparation ------------------------------------------------------------
# Filter alpha-diversity data

frag_cont_rare_div2 %>% 
  filter(diversity_index == "alpha" & q_order == "q = 2") -> rare2_alpha2

transformed_data_rare2_alpha2 <- rare2_alpha2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log response ratio (LRR) (Hedges et al., 1999; Lajeunesse, 2011)
rare2_alpha2_logR <- escalc(measure="ROM", 
                            m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                            sd1i = continuous_sd, sd2i = fragmented_sd,
                            n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                            data = transformed_data_rare2_alpha2)

# combine matrices to add habitat amount

rare2_alpha2_logR <- left_join(rare2_alpha2_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare2_alpha2 <- rma.uni(yi, vi,
                                    data = rare2_alpha2_logR,
                                    method = "REML") # with random effect 

summary(mod_overall_rare2_alpha2) # reported results 



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare2_alpha2, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)



# tables for plotting later

overall_rare2_alpha2_res <- mod_results(mod_overall_rare2_alpha2, mod = "1", group = "refshort")
overall_rare2_alpha2_est <- overall_rare2_alpha2_res$mod_table %>% 
  dplyr::mutate(diversity = "Alpha", 
                type = "All pairs (q = 2)",
                .before = name)
overall_rare2_alpha2_dat <- overall_rare2_alpha2_res$data %>% 
  dplyr::mutate(diversity = "Alpha",
                type = "All pairs (q = 2)",
                .before = yi)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare2_alpha2 <- rma.uni(yi, vi,
                                  data = rare2_alpha2_logR,
                                  mods = ~ ha_diff,
                                  method = "REML") # with random effect 

summary(mod_mods1_rare2_alpha2) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare2_alpha2 <- rma.uni(yi, vi,
                                  data = rare2_alpha2_logR,
                                  mods = ~ ha_class -1,
                                  method = "REML") # with random effect 

summary(mod_mods2_rare2_alpha2) # reported results 

orchaRd::orchard_plot(mod_mods2_rare2_alpha2, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare2_alpha2_logR$continent_bin <- ifelse(rare2_alpha2_logR$continent == "South America", "South America", "Others")


mod_mods3_rare2_alpha2 <- rma.uni(yi, vi,
                                  data = rare2_alpha2_logR,
                                  mods = ~ continent_bin - 1,
                                  method = "REML") # with random effect 

summary(mod_mods3_rare2_alpha2) # reported results 

orchaRd::orchard_plot(mod_mods3_rare2_alpha2, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare2_alpha2 <- rma.uni(yi, vi,
                                  data = rare2_alpha2_logR,
                                  mods = ~ time_since_fragmentation - 1,
                                  method = "REML") # with random effect 

summary(mod_mods4_rare2_alpha2) # reported results 

orchaRd::orchard_plot(mod_mods4_rare2_alpha2, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare2_alpha2 <- influence.rma.uni(mod_overall_rare2_alpha2)
plot(mod_infl_rare2_alpha2) # study 26 is an outlier

# calculate the inverse of effective sample size
rare2_alpha2_logR$inv_n <- with(rare2_alpha2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare2_alpha2_logR$sqrt_inv_n <- with(rare2_alpha2_logR, sqrt(inv_n))

# Time-lag bias test

rare2_alpha2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare2_alpha2_logR$refshort))
rare2_alpha2_logR$year_c <- as.vector(scale(rare2_alpha2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare2_alpha2_check <- rma(yi, vi,
                              mods = ~1 + sqrt_inv_n + year_c,
                              data = rare2_alpha2_logR,
                              method = "REML")
summary(mod_rare2_alpha2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare2_alpha2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Alpha", 
         type = "All pairs (q = 2)",
         refshort = rare2_alpha2_logR$refshort,
         continent_bin = rare2_alpha2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare2_alpha2
leave1out_mod_rare2_alpha2


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare2_alpha2 <- extract_model_results(mod_mods1_rare2_alpha2, 
                                                    diversity_index = "Alpha (q=2)", 
                                                    diversity_type = "All pairs") 


results_mods2_rare2_alpha2 <-extract_model_results(mod_mods2_rare2_alpha2, 
                                                   diversity_index = "Alpha (q=2)", 
                                                   diversity_type = "All pairs") 


results_mods3_rare2_alpha2 <- extract_model_results(mod_mods3_rare2_alpha2, 
                                                    diversity_index = "Alpha (q=2)", 
                                                    diversity_type = "All pairs") 



results_mods4_rare2_alpha2 <- extract_model_results(mod_mods4_rare2_alpha2, 
                                                    diversity_index = "Alpha (q=2)", 
                                                    diversity_type = "All pairs") 

results_rare2_alpha2_check <- extract_model_results(mod_rare2_alpha2_check, 
                                                    diversity_index = "Alpha (q=2)", 
                                                    diversity_type = "All pairs") 


mod_results_rare2_alpha2 <- rbind.data.frame(results_mods1_rare2_alpha2, 
                                             results_mods2_rare2_alpha2, 
                                             results_mods3_rare2_alpha2,
                                             results_mods4_rare2_alpha2,
                                             results_rare2_alpha2_check)


### Rarefied BETA diversity (all fragment / plot pairs) - q = 0


# required data frames

frag_cont_rare_div2 %>% 
  filter(diversity_index == "beta" & buffer == 2000 & q_order == "q = 0") -> rare0_beta2

transformed_data_rare0_beta2 <- rare0_beta2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare0_beta2_logR <- escalc(measure="ROM", 
                           m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                           sd1i = continuous_sd, sd2i = fragmented_sd,
                           n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                           data = transformed_data_rare0_beta2)


# combine matrices to add habitat amount

rare0_beta2_logR <- left_join(rare0_beta2_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare0_beta2 <- rma.uni(yi, vi,
                                   data = rare0_beta2_logR,
                                   method = "REML") # with random effect 

summary(mod_overall_rare0_beta2) # reported results 

# tables for plotting later

overall_rare0_beta2_res <- mod_results(mod_overall_rare0_beta2, mod = "1", group = "refshort")
overall_rare0_beta2_est <- overall_rare0_beta2_res$mod_table %>% 
  dplyr::mutate(diversity = "Beta", 
                type = "All pairs (q = 0)",
                .before = name)
overall_rare0_beta2_dat <- overall_rare0_beta2_res$data %>% 
  dplyr::mutate(diversity = "Beta",
                type = "All pairs (q = 0)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare0_beta2, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare0_beta2 <- rma.uni(yi, vi,
                                 data = rare0_beta2_logR,
                                 mods = ~ ha_diff,
                                 method = "REML") # with random effect 

summary(mod_mods1_rare0_beta2) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare0_beta2 <- rma.uni(yi, vi,
                                 data = rare0_beta2_logR,
                                 mods = ~ ha_class -1,
                                 method = "REML") # with random effect 

summary(mod_mods2_rare0_beta2) # reported results 

orchaRd::orchard_plot(mod_mods2_rare0_beta2, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare0_beta2_logR$continent_bin <- ifelse(rare0_beta2_logR$continent == "South America", "South America", "Others")


mod_mods3_rare0_beta2 <- rma.uni(yi, vi,
                                 data = rare0_beta2_logR,
                                 mods = ~ continent_bin - 1,
                                 method = "REML") # with random effect 

summary(mod_mods3_rare0_beta2) # reported results 

orchaRd::orchard_plot(mod_mods3_rare0_beta2, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare0_beta2 <- rma.uni(yi, vi,
                                 data = rare0_beta2_logR,
                                 mods = ~ time_since_fragmentation - 1,
                                 method = "REML") # with random effect 

summary(mod_mods4_rare0_beta2) # reported results 

orchaRd::orchard_plot(mod_mods4_rare0_beta2, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare0_beta2 <- influence.rma.uni(mod_overall_rare0_beta2)
plot(mod_infl_rare0_beta2) # study 26 is an outlier

# calculate the inverse of effective sample size
rare0_beta2_logR$inv_n <- with(rare0_beta2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare0_beta2_logR$sqrt_inv_n <- with(rare0_beta2_logR, sqrt(inv_n))

# Time-lag bias test

rare0_beta2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare0_beta2_logR$refshort))
rare0_beta2_logR$year_c <- as.vector(scale(rare0_beta2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare0_beta2_check <- rma(yi, vi,
                             mods = ~1 + sqrt_inv_n + year_c,
                             data = rare0_beta2_logR,
                             method = "REML")
summary(mod_rare0_beta2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare0_beta2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Beta", 
         type = "All pairs (q = 0)",
         refshort = rare0_beta2_logR$refshort,
         continent_bin = rare0_beta2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare0_beta2
leave1out_mod_rare0_beta2


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare0_beta2 <- extract_model_results(mod_mods1_rare0_beta2, 
                                                   diversity_index = "Beta (q=0)", 
                                                   diversity_type = "All pairs") 


results_mods2_rare0_beta2 <-extract_model_results(mod_mods2_rare0_beta2, 
                                                  diversity_index = "Beta (q=0)", 
                                                  diversity_type = "All pairs") 


results_mods3_rare0_beta2 <- extract_model_results(mod_mods3_rare0_beta2, 
                                                   diversity_index = "Beta (q=0)", 
                                                   diversity_type = "All pairs") 



results_mods4_rare0_beta2 <- extract_model_results(mod_mods4_rare0_beta2, 
                                                   diversity_index = "Beta (q=0)", 
                                                   diversity_type = "All pairs") 

results_rare0_beta2_check <- extract_model_results(mod_rare0_beta2_check, 
                                                   diversity_index = "Beta (q=0)", 
                                                   diversity_type = "All pairs") 


mod_results_rare0_beta2 <- rbind.data.frame(results_mods1_rare0_beta2, 
                                            results_mods2_rare0_beta2, 
                                            results_mods3_rare0_beta2,
                                            results_mods4_rare0_beta2,
                                            results_rare0_beta2_check)




### Rarefied BETA diversity (all fragment / plot pairs) - q = 2


# required data frames

frag_cont_rare_div2 %>% 
  filter(diversity_index == "beta" & buffer == 2000 & q_order == "q = 2") -> rare2_beta2

transformed_data_rare2_beta2 <- rare2_beta2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare2_beta2_logR <- escalc(measure="ROM", 
                           m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                           sd1i = continuous_sd, sd2i = fragmented_sd,
                           n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                           data = transformed_data_rare2_beta2)


# combine matrices to add habitat amount

rare2_beta2_logR <- left_join(rare2_beta2_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare2_beta2 <- rma.uni(yi, vi,
                                   data = rare2_beta2_logR,
                                   method = "REML") # with random effect 

summary(mod_overall_rare2_beta2) # reported results 

# tables for plotting later

overall_rare2_beta2_res <- mod_results(mod_overall_rare2_beta2, mod = "1", group = "refshort")
overall_rare2_beta2_est <- overall_rare2_beta2_res$mod_table %>% 
  dplyr::mutate(diversity = "Beta", 
                type = "All pairs (q = 2)",
                .before = name)
overall_rare2_beta2_dat <- overall_rare2_beta2_res$data %>% 
  dplyr::mutate(diversity = "Beta",
                type = "All pairs (q = 2)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare2_beta2, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare2_beta2 <- rma.uni(yi, vi,
                                 data = rare2_beta2_logR,
                                 mods = ~ ha_diff,
                                 method = "REML") # with random effect 

summary(mod_mods1_rare2_beta2) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare2_beta2 <- rma.uni(yi, vi,
                                 data = rare2_beta2_logR,
                                 mods = ~ ha_class -1,
                                 method = "REML") # with random effect 

summary(mod_mods2_rare2_beta2) # reported results 

orchaRd::orchard_plot(mod_mods2_rare2_beta2, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare2_beta2_logR$continent_bin <- ifelse(rare2_beta2_logR$continent == "South America", "South America", "Others")


mod_mods3_rare2_beta2 <- rma.uni(yi, vi,
                                 data = rare2_beta2_logR,
                                 mods = ~ continent_bin - 1,
                                 method = "REML") # with random effect 

summary(mod_mods3_rare2_beta2) # reported results 

orchaRd::orchard_plot(mod_mods3_rare2_beta2, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare2_beta2 <- rma.uni(yi, vi,
                                 data = rare2_beta2_logR,
                                 mods = ~ time_since_fragmentation - 1,
                                 method = "REML") # with random effect 

summary(mod_mods4_rare2_beta2) # reported results 

orchaRd::orchard_plot(mod_mods4_rare2_beta2, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare2_beta2 <- influence.rma.uni(mod_overall_rare2_beta2)
plot(mod_infl_rare2_beta2) # study 26 is an outlier

# calculate the inverse of effective sample size
rare2_beta2_logR$inv_n <- with(rare2_beta2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare2_beta2_logR$sqrt_inv_n <- with(rare2_beta2_logR, sqrt(inv_n))

# Time-lag bias test

rare2_beta2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare2_beta2_logR$refshort))
rare2_beta2_logR$year_c <- as.vector(scale(rare2_beta2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare2_beta2_check <- rma(yi, vi,
                             mods = ~1 + sqrt_inv_n + year_c,
                             data = rare2_beta2_logR,
                             method = "REML")
summary(mod_rare2_beta2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare2_beta2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Beta", 
         type = "All pairs (q = 2)",
         refshort = rare2_beta2_logR$refshort,
         continent_bin = rare2_beta2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare2_beta2
leave1out_mod_rare2_beta2

## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare2_beta2 <- extract_model_results(mod_mods1_rare2_beta2, 
                                                   diversity_index = "Beta (q=2)", 
                                                   diversity_type = "All pairs") 


results_mods2_rare2_beta2 <-extract_model_results(mod_mods2_rare2_beta2, 
                                                  diversity_index = "Beta (q=2)", 
                                                  diversity_type = "All pairs") 


results_mods3_rare2_beta2 <- extract_model_results(mod_mods3_rare2_beta2, 
                                                   diversity_index = "Beta (q=2)", 
                                                   diversity_type = "All pairs") 



results_mods4_rare2_beta2 <- extract_model_results(mod_mods4_rare2_beta2, 
                                                   diversity_index = "Beta (q=2)", 
                                                   diversity_type = "All pairs") 

results_rare2_beta2_check <- extract_model_results(mod_rare2_beta2_check, 
                                                   diversity_index = "Beta (q=2)", 
                                                   diversity_type = "All pairs") 


mod_results_rare2_beta2 <- rbind.data.frame(results_mods1_rare2_beta2, 
                                            results_mods2_rare2_beta2, 
                                            results_mods3_rare2_beta2,
                                            results_mods4_rare2_beta2,
                                            results_rare2_beta2_check)

### Rarified GAMMA diversity (all fragment / plot pairs) - q = 0


# required data frames

frag_cont_rare_div2 %>% 
  filter(diversity_index == "gamma" & buffer == 2000 & q_order == "q = 0") -> rare0_gamma2

transformed_data_rare0_gamma2 <- rare0_gamma2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare0_gamma2_logR <- escalc(measure="ROM", 
                            m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                            sd1i = continuous_sd, sd2i = fragmented_sd,
                            n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                            data = transformed_data_rare0_gamma2)


# combine matrices to add habitat amount

rare0_gamma2_logR <- left_join(rare0_gamma2_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare0_gamma2 <- rma.uni(yi, vi,
                                    data = rare0_gamma2_logR,
                                    method = "REML") # with random effect 

summary(mod_overall_rare0_gamma2) # reported results 

# tables for plotting later

overall_rare0_gamma2_res <- mod_results(mod_overall_rare0_gamma2, mod = "1", group = "refshort")
overall_rare0_gamma2_est <- overall_rare0_gamma2_res$mod_table %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "All pairs (q = 0)",
                .before = name)
overall_rare0_gamma2_dat <- overall_rare0_gamma2_res$data %>% 
  dplyr::mutate(diversity = "Gamma",
                type = "All pairs (q = 0)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare0_gamma2, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare0_gamma2 <- rma.uni(yi, vi,
                                  data = rare0_gamma2_logR,
                                  mods = ~ ha_diff,
                                  method = "REML") # with random effect 

summary(mod_mods1_rare0_gamma2) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare0_gamma2 <- rma.uni(yi, vi,
                                  data = rare0_gamma2_logR,
                                  mods = ~ ha_class -1,
                                  method = "REML") # with random effect 

summary(mod_mods2_rare0_gamma2) # reported results 

orchaRd::orchard_plot(mod_mods2_rare0_gamma2, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare0_gamma2_logR$continent_bin <- ifelse(rare0_gamma2_logR$continent == "South America", "South America", "Others")


mod_mods3_rare0_gamma2 <- rma.uni(yi, vi,
                                  data = rare0_gamma2_logR,
                                  mods = ~ continent_bin - 1,
                                  method = "REML") # with random effect 

summary(mod_mods3_rare0_gamma2) # reported results 

orchaRd::orchard_plot(mod_mods3_rare0_gamma2, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare0_gamma2 <- rma.uni(yi, vi,
                                  data = rare0_gamma2_logR,
                                  mods = ~ time_since_fragmentation - 1,
                                  method = "REML") # with random effect 

summary(mod_mods4_rare0_gamma2) # reported results 

orchaRd::orchard_plot(mod_mods4_rare0_gamma2, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare0_gamma2 <- influence.rma.uni(mod_overall_rare0_gamma2)
plot(mod_infl_rare0_gamma2) # study 26 is an outlier

# calculate the inverse of effective sample size
rare0_gamma2_logR$inv_n <- with(rare0_gamma2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare0_gamma2_logR$sqrt_inv_n <- with(rare0_gamma2_logR, sqrt(inv_n))

# Time-lag bias test

rare0_gamma2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare0_gamma2_logR$refshort))
rare0_gamma2_logR$year_c <- as.vector(scale(rare0_gamma2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare0_gamma2_check <- rma(yi, vi,
                              mods = ~1 + sqrt_inv_n + year_c,
                              data = rare0_gamma2_logR,
                              method = "REML")
summary(mod_rare0_gamma2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare0_gamma2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Gamma", 
         type = "All pairs (q = 0)",
         refshort = rare0_gamma2_logR$refshort,
         continent_bin = rare0_gamma2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare0_gamma2
leave1out_mod_rare0_gamma2


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare0_gamma2 <- extract_model_results(mod_mods1_rare0_gamma2, 
                                                    diversity_index = "Gamma (q=0)", 
                                                    diversity_type = "All pairs") 


results_mods2_rare0_gamma2 <-extract_model_results(mod_mods2_rare0_gamma2, 
                                                   diversity_index = "Gamma (q=0)", 
                                                   diversity_type = "All pairs") 


results_mods3_rare0_gamma2 <- extract_model_results(mod_mods3_rare0_gamma2, 
                                                    diversity_index = "Gamma (q=0)", 
                                                    diversity_type = "All pairs") 



results_mods4_rare0_gamma2 <- extract_model_results(mod_mods4_rare0_gamma2, 
                                                    diversity_index = "Gamma (q=0)", 
                                                    diversity_type = "All pairs") 

results_rare0_gamma2_check <- extract_model_results(mod_rare0_gamma2_check, 
                                                    diversity_index = "Gamma (q=0)", 
                                                    diversity_type = "All pairs") 


mod_results_rare0_gamma2 <- rbind.data.frame(results_mods1_rare0_gamma2, 
                                             results_mods2_rare0_gamma2, 
                                             results_mods3_rare0_gamma2,
                                             results_mods4_rare0_gamma2,
                                             results_rare0_gamma2_check)

### Rarefied GAMMA diversity (all fragment / plot pairs) - q = 2


# required data frames

frag_cont_rare_div2 %>% 
  filter(diversity_index == "gamma" & buffer == 2000 & q_order == "q = 2") -> rare2_gamma2

transformed_data_rare2_gamma2 <- rare2_gamma2 %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare2_gamma2_logR <- escalc(measure="ROM", 
                            m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                            sd1i = continuous_sd, sd2i = fragmented_sd,
                            n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                            data = transformed_data_rare2_gamma2)


# combine matrices to add habitat amount

rare2_gamma2_logR <- left_join(rare2_gamma2_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare2_gamma2 <- rma.uni(yi, vi,
                                    data = rare2_gamma2_logR,
                                    method = "REML") # with random effect 

summary(mod_overall_rare2_gamma2) # reported results 

# tables for plotting later

overall_rare2_gamma2_res <- mod_results(mod_overall_rare2_gamma2, mod = "1", group = "refshort")
overall_rare2_gamma2_est <- overall_rare2_gamma2_res$mod_table %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "All pairs (q = 2)",
                .before = name)
overall_rare2_gamma2_dat <- overall_rare2_gamma2_res$data %>% 
  dplyr::mutate(diversity = "Gamma",
                type = "All pairs (q = 2)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare2_gamma2, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)



### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare2_gamma2 <- rma.uni(yi, vi,
                                  data = rare2_gamma2_logR,
                                  mods = ~ ha_diff,
                                  method = "REML") # with random effect 

summary(mod_mods1_rare2_gamma2) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare2_gamma2 <- rma.uni(yi, vi,
                                  data = rare2_gamma2_logR,
                                  mods = ~ ha_class -1,
                                  method = "REML") # with random effect 

summary(mod_mods2_rare2_gamma2) # reported results 

orchaRd::orchard_plot(mod_mods2_rare2_gamma2, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare2_gamma2_logR$continent_bin <- ifelse(rare2_gamma2_logR$continent == "South America", "South America", "Others")


mod_mods3_rare2_gamma2 <- rma.uni(yi, vi,
                                  data = rare2_gamma2_logR,
                                  mods = ~ continent_bin - 1,
                                  method = "REML") # with random effect 

summary(mod_mods3_rare2_gamma2) # reported results 

orchaRd::orchard_plot(mod_mods3_rare2_gamma2, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare2_gamma2 <- rma.uni(yi, vi,
                                  data = rare2_gamma2_logR,
                                  mods = ~ time_since_fragmentation - 1,
                                  method = "REML") # with random effect 

summary(mod_mods4_rare2_gamma2) # reported results 

orchaRd::orchard_plot(mod_mods4_rare2_gamma2, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare2_gamma2 <- influence.rma.uni(mod_overall_rare2_gamma2)
plot(mod_infl_rare2_gamma2) # study 26 is an outlier

# calculate the inverse of effective sample size
rare2_gamma2_logR$inv_n <- with(rare2_gamma2_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare2_gamma2_logR$sqrt_inv_n <- with(rare2_gamma2_logR, sqrt(inv_n))

# Time-lag bias test

rare2_gamma2_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare2_gamma2_logR$refshort))
rare2_gamma2_logR$year_c <- as.vector(scale(rare2_gamma2_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare2_gamma2_check <- rma(yi, vi,
                              mods = ~1 + sqrt_inv_n + year_c,
                              data = rare2_gamma2_logR,
                              method = "REML")
summary(mod_rare2_gamma2_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare2_gamma2) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Gamma", 
         type = "All pairs (q = 2)",
         refshort = rare2_gamma2_logR$refshort,
         continent_bin = rare2_gamma2_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare2_gamma2
leave1out_mod_rare2_gamma2

## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare2_gamma2 <- extract_model_results(mod_mods1_rare2_gamma2, 
                                                    diversity_index = "Gamma (q=2)", 
                                                    diversity_type = "All pairs") 


results_mods2_rare2_gamma2 <-extract_model_results(mod_mods2_rare2_gamma2, 
                                                   diversity_index = "Gamma (q=2)", 
                                                   diversity_type = "All pairs") 


results_mods3_rare2_gamma2 <- extract_model_results(mod_mods3_rare2_gamma2, 
                                                    diversity_index = "Gamma (q=2)", 
                                                    diversity_type = "All pairs") 



results_mods4_rare2_gamma2 <- extract_model_results(mod_mods4_rare2_gamma2, 
                                                    diversity_index = "Gamma (q=2)", 
                                                    diversity_type = "All pairs") 

results_rare2_gamma2_check <- extract_model_results(mod_rare2_gamma2_check, 
                                                    diversity_index = "Gamma (q=2)", 
                                                    diversity_type = "All pairs") 


mod_results_rare2_gamma2 <- rbind.data.frame(results_mods1_rare2_gamma2, 
                                             results_mods2_rare2_gamma2, 
                                             results_mods3_rare2_gamma2,
                                             results_mods4_rare2_gamma2,
                                             results_rare2_gamma2_check)



### Rarefied ALPHA diversity (nearest pairs) - q = 0


# required data frames

frag_cont_rare_div2_cl %>% 
  filter(diversity_index == "alpha" & buffer == 2000 & q_order == "q = 0") -> rare0_alpha2_cl

transformed_data_rare0_alpha2_cl <- rare0_alpha2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare0_alpha2_cl_logR <- escalc(measure="ROM", 
                               m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                               sd1i = continuous_sd, sd2i = fragmented_sd,
                               n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                               data = transformed_data_rare0_alpha2_cl)


# combine matrices to add habitat amount

rare0_alpha2_cl_logR <- left_join(rare0_alpha2_cl_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare0_alpha2_cl <- rma.uni(yi, vi,
                                       data = rare0_alpha2_cl_logR,
                                       method = "REML") # with random effect 

summary(mod_overall_rare0_alpha2_cl) # reported results 

# tables for plotting later

overall_rare0_alpha2_cl_res <- mod_results(mod_overall_rare0_alpha2_cl, mod = "1", group = "refshort")
overall_rare0_alpha2_cl_est <- overall_rare0_alpha2_cl_res$mod_table %>% 
  dplyr::mutate(diversity = "Alpha", 
                type = "Nearest pairs (q = 0)",
                .before = name)
overall_rare0_alpha2_cl_dat <- overall_rare0_alpha2_cl_res$data %>% 
  dplyr::mutate(diversity = "Alpha",
                type = "Nearest pairs (q = 0)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare0_alpha2_cl, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare0_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare0_alpha2_cl_logR,
                                     mods = ~ ha_diff,
                                     method = "REML") # with random effect 

summary(mod_mods1_rare0_alpha2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare0_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare0_alpha2_cl_logR,
                                     mods = ~ ha_class -1,
                                     method = "REML") # with random effect 

summary(mod_mods2_rare0_alpha2_cl) # reported results 

orchaRd::orchard_plot(mod_mods2_rare0_alpha2_cl, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare0_alpha2_cl_logR$continent_bin <- ifelse(rare0_alpha2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_rare0_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare0_alpha2_cl_logR,
                                     mods = ~ continent_bin - 1,
                                     method = "REML") # with random effect 

summary(mod_mods3_rare0_alpha2_cl) # reported results 

orchaRd::orchard_plot(mod_mods3_rare0_alpha2_cl, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare0_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare0_alpha2_cl_logR,
                                     mods = ~ time_since_fragmentation - 1,
                                     method = "REML") # with random effect 

summary(mod_mods4_rare0_alpha2_cl) # reported results 

orchaRd::orchard_plot(mod_mods4_rare0_alpha2_cl, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare0_alpha2_cl <- influence.rma.uni(mod_overall_rare0_alpha2_cl)
plot(mod_infl_rare0_alpha2_cl) # study 26 is an outlier

# calculate the inverse of effective sample size
rare0_alpha2_cl_logR$inv_n <- with(rare0_alpha2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare0_alpha2_cl_logR$sqrt_inv_n <- with(rare0_alpha2_cl_logR, sqrt(inv_n))

# Time-lag bias test

rare0_alpha2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare0_alpha2_cl_logR$refshort))
rare0_alpha2_cl_logR$year_c <- as.vector(scale(rare0_alpha2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare0_alpha2_cl_check <- rma(yi, vi,
                                 mods = ~1 + sqrt_inv_n + year_c,
                                 data = rare0_alpha2_cl_logR,
                                 method = "REML")
summary(mod_rare0_alpha2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare0_alpha2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Alpha", 
         type = "Nearest pairs (q = 0)",
         refshort = rare0_alpha2_cl_logR$refshort,
         continent_bin = rare0_alpha2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare0_alpha2_cl
leave1out_mod_rare0_alpha2_cl

## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare0_alpha2_cl <- extract_model_results(mod_mods1_rare0_alpha2_cl, 
                                                       diversity_index = "Alpha (q=0)", 
                                                       diversity_type = "Nearest pairs") 


results_mods2_rare0_alpha2_cl <-extract_model_results(mod_mods2_rare0_alpha2_cl, 
                                                      diversity_index = "Alpha (q=0)", 
                                                      diversity_type = "Nearest pairs") 


results_mods3_rare0_alpha2_cl <- extract_model_results(mod_mods3_rare0_alpha2_cl, 
                                                       diversity_index = "Alpha (q=0)", 
                                                       diversity_type = "Nearest pairs") 



results_mods4_rare0_alpha2_cl <- extract_model_results(mod_mods4_rare0_alpha2_cl, 
                                                       diversity_index = "Alpha (q=0)", 
                                                       diversity_type = "Nearest pairs") 

results_rare0_alpha2_cl_check <- extract_model_results(mod_rare0_alpha2_cl_check, 
                                                       diversity_index = "Alpha (q=0)", 
                                                       diversity_type = "Nearest pairs") 


mod_results_rare0_alpha2_cl <- rbind.data.frame(results_mods1_rare0_alpha2_cl, 
                                                results_mods2_rare0_alpha2_cl, 
                                                results_mods3_rare0_alpha2_cl,
                                                results_mods4_rare0_alpha2_cl,
                                                results_rare0_alpha2_cl_check)



### Rarefied ALPHA diversity (nearest pairs) - q = 2


# required data frames

frag_cont_rare_div2_cl %>% 
  filter(diversity_index == "alpha" & buffer == 2000 & q_order == "q = 2") -> rare2_alpha2_cl

transformed_data_rare2_alpha2_cl <- rare2_alpha2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare2_alpha2_cl_logR <- escalc(measure="ROM", 
                               m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                               sd1i = continuous_sd, sd2i = fragmented_sd,
                               n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                               data = transformed_data_rare2_alpha2_cl)


# combine matrices to add habitat amount

rare2_alpha2_cl_logR <- left_join(rare2_alpha2_cl_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare2_alpha2_cl <- rma.uni(yi, vi,
                                       data = rare2_alpha2_cl_logR,
                                       method = "REML") # with random effect 

summary(mod_overall_rare2_alpha2_cl) # reported results 

# tables for plotting later

overall_rare2_alpha2_cl_res <- mod_results(mod_overall_rare2_alpha2_cl, mod = "1", group = "refshort")
overall_rare2_alpha2_cl_est <- overall_rare2_alpha2_cl_res$mod_table %>% 
  dplyr::mutate(diversity = "Alpha", 
                type = "Nearest pairs (q = 2)",
                .before = name)
overall_rare2_alpha2_cl_dat <- overall_rare2_alpha2_cl_res$data %>% 
  dplyr::mutate(diversity = "Alpha",
                type = "Nearest pairs (q = 2)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare2_alpha2_cl, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare2_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare2_alpha2_cl_logR,
                                     mods = ~ ha_diff,
                                     method = "REML") # with random effect 

summary(mod_mods1_rare2_alpha2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare2_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare2_alpha2_cl_logR,
                                     mods = ~ ha_class -1,
                                     method = "REML") # with random effect 

summary(mod_mods2_rare2_alpha2_cl) # reported results 

orchaRd::orchard_plot(mod_mods2_rare2_alpha2_cl, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare2_alpha2_cl_logR$continent_bin <- ifelse(rare2_alpha2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_rare2_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare2_alpha2_cl_logR,
                                     mods = ~ continent_bin - 1,
                                     method = "REML") # with random effect 

summary(mod_mods3_rare2_alpha2_cl) # reported results 

orchaRd::orchard_plot(mod_mods3_rare2_alpha2_cl, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare2_alpha2_cl <- rma.uni(yi, vi,
                                     data = rare2_alpha2_cl_logR,
                                     mods = ~ time_since_fragmentation - 1,
                                     method = "REML") # with random effect 

summary(mod_mods4_rare2_alpha2_cl) # reported results 

orchaRd::orchard_plot(mod_mods4_rare2_alpha2_cl, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare2_alpha2_cl <- influence.rma.uni(mod_overall_rare2_alpha2_cl)
plot(mod_infl_rare2_alpha2_cl) # study 26 is an outlier

# calculate the inverse of effective sample size
rare2_alpha2_cl_logR$inv_n <- with(rare2_alpha2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare2_alpha2_cl_logR$sqrt_inv_n <- with(rare2_alpha2_cl_logR, sqrt(inv_n))

# Time-lag bias test

rare2_alpha2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare2_alpha2_cl_logR$refshort))
rare2_alpha2_cl_logR$year_c <- as.vector(scale(rare2_alpha2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare2_alpha2_cl_check <- rma(yi, vi,
                                 mods = ~1 + sqrt_inv_n + year_c,
                                 data = rare2_alpha2_cl_logR,
                                 method = "REML")
summary(mod_rare2_alpha2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)


# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare2_alpha2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Alpha", 
         type = "Nearest pairs (q = 2)",
         refshort = rare2_alpha2_cl_logR$refshort,
         continent_bin = rare2_alpha2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare2_alpha2_cl
leave1out_mod_rare2_alpha2_cl


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare2_alpha2_cl <- extract_model_results(mod_mods1_rare2_alpha2_cl, 
                                                       diversity_index = "Alpha (q=2)", 
                                                       diversity_type = "Nearest pairs") 


results_mods2_rare2_alpha2_cl <-extract_model_results(mod_mods2_rare2_alpha2_cl, 
                                                      diversity_index = "Alpha (q=2)", 
                                                      diversity_type = "Nearest pairs") 


results_mods3_rare2_alpha2_cl <- extract_model_results(mod_mods3_rare2_alpha2_cl, 
                                                       diversity_index = "Alpha (q=2)", 
                                                       diversity_type = "Nearest pairs") 



results_mods4_rare2_alpha2_cl <- extract_model_results(mod_mods4_rare2_alpha2_cl, 
                                                       diversity_index = "Alpha (q=2)", 
                                                       diversity_type = "Nearest pairs") 

results_rare2_alpha2_cl_check <- extract_model_results(mod_rare2_alpha2_cl_check, 
                                                       diversity_index = "Alpha (q=2)", 
                                                       diversity_type = "Nearest pairs") 


mod_results_rare2_alpha2_cl <- rbind.data.frame(results_mods1_rare2_alpha2_cl, 
                                                results_mods2_rare2_alpha2_cl, 
                                                results_mods3_rare2_alpha2_cl,
                                                results_mods4_rare2_alpha2_cl,
                                                results_rare2_alpha2_cl_check)


### Rarefied BETA diversity (nearest pairs) - q = 0


# required data frames

frag_cont_rare_div2_cl %>% 
  filter(diversity_index == "beta" & buffer == 2000 & q_order == "q = 0") -> rare0_beta2_cl

transformed_data_rare0_beta2_cl <- rare0_beta2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare0_beta2_cl_logR <- escalc(measure="ROM", 
                              m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                              sd1i = continuous_sd, sd2i = fragmented_sd,
                              n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                              data = transformed_data_rare0_beta2_cl)


# combine matrices to add habitat amount

rare0_beta2_cl_logR <- left_join(rare0_beta2_cl_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare0_beta2_cl <- rma.uni(yi, vi,
                                      data = rare0_beta2_cl_logR,
                                      method = "REML") # with random effect 

summary(mod_overall_rare0_beta2_cl) # reported results 

# tables for plotting later

overall_rare0_beta2_cl_res <- mod_results(mod_overall_rare0_beta2_cl, mod = "1", group = "refshort")
overall_rare0_beta2_cl_est <- overall_rare0_beta2_cl_res$mod_table %>% 
  dplyr::mutate(diversity = "Beta", 
                type = "Nearest pairs (q = 0)",
                .before = name)
overall_rare0_beta2_cl_dat <- overall_rare0_beta2_cl_res$data %>% 
  dplyr::mutate(diversity = "Beta",
                type = "Nearest pairs (q = 0)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare0_beta2_cl, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare0_beta2_cl <- rma.uni(yi, vi,
                                    data = rare0_beta2_cl_logR,
                                    mods = ~ ha_diff,
                                    method = "REML") # with random effect 

summary(mod_mods1_rare0_beta2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare0_beta2_cl <- rma.uni(yi, vi,
                                    data = rare0_beta2_cl_logR,
                                    mods = ~ ha_class -1,
                                    method = "REML") # with random effect 

summary(mod_mods2_rare0_beta2_cl) # reported results 

orchaRd::orchard_plot(mod_mods2_rare0_beta2_cl, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare0_beta2_cl_logR$continent_bin <- ifelse(rare0_beta2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_rare0_beta2_cl <- rma.uni(yi, vi,
                                    data = rare0_beta2_cl_logR,
                                    mods = ~ continent_bin - 1,
                                    method = "REML") # with random effect 

summary(mod_mods3_rare0_beta2_cl) # reported results 

orchaRd::orchard_plot(mod_mods3_rare0_beta2_cl, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare0_beta2_cl <- rma.uni(yi, vi,
                                    data = rare0_beta2_cl_logR,
                                    mods = ~ time_since_fragmentation - 1,
                                    method = "REML") # with random effect 

summary(mod_mods4_rare0_beta2_cl) # reported results 

orchaRd::orchard_plot(mod_mods4_rare0_beta2_cl, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare0_beta2_cl <- influence.rma.uni(mod_overall_rare0_beta2_cl)
plot(mod_infl_rare0_beta2_cl) # study 26 is an outlier

# calculate the inverse of effective sample size
rare0_beta2_cl_logR$inv_n <- with(rare0_beta2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare0_beta2_cl_logR$sqrt_inv_n <- with(rare0_beta2_cl_logR, sqrt(inv_n))

# Time-lag bias test

rare0_beta2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare0_beta2_cl_logR$refshort))
rare0_beta2_cl_logR$year_c <- as.vector(scale(rare0_beta2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare0_beta2_cl_check <- rma(yi, vi,
                                mods = ~1 + sqrt_inv_n + year_c,
                                data = rare0_beta2_cl_logR,
                                method = "REML")
summary(mod_rare0_beta2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare0_beta2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Beta", 
         type = "Nearest pairs (q = 0)",
         refshort = rare0_beta2_cl_logR$refshort,
         continent_bin = rare0_beta2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare0_beta2_cl
leave1out_mod_rare0_beta2_cl


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare0_beta2_cl <- extract_model_results(mod_mods1_rare0_beta2_cl, 
                                                      diversity_index = "Beta (q=0)", 
                                                      diversity_type = "Nearest pairs") 


results_mods2_rare0_beta2_cl <-extract_model_results(mod_mods2_rare0_beta2_cl, 
                                                     diversity_index = "Beta (q=0)", 
                                                     diversity_type = "Nearest pairs") 


results_mods3_rare0_beta2_cl <- extract_model_results(mod_mods3_rare0_beta2_cl, 
                                                      diversity_index = "Beta (q=0)", 
                                                      diversity_type = "Nearest pairs") 



results_mods4_rare0_beta2_cl <- extract_model_results(mod_mods4_rare0_beta2_cl, 
                                                      diversity_index = "Beta (q=0)", 
                                                      diversity_type = "Nearest pairs") 

results_rare0_beta2_cl_check <- extract_model_results(mod_rare0_beta2_cl_check, 
                                                      diversity_index = "Beta (q=0)", 
                                                      diversity_type = "Nearest pairs") 


mod_results_rare0_beta2_cl <- rbind.data.frame(results_mods1_rare0_beta2_cl, 
                                               results_mods2_rare0_beta2_cl, 
                                               results_mods3_rare0_beta2_cl,
                                               results_mods4_rare0_beta2_cl,
                                               results_rare0_beta2_cl_check)


### Rarefied BETA diversity (nearest pairs) - q = 2


# required data frames

frag_cont_rare_div2_cl %>% 
  filter(diversity_index == "beta" & buffer == 2000 & q_order == "q = 2") -> rare2_beta2_cl

transformed_data_rare2_beta2_cl <- rare2_beta2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare2_beta2_cl_logR <- escalc(measure="ROM", 
                              m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                              sd1i = continuous_sd, sd2i = fragmented_sd,
                              n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                              data = transformed_data_rare2_beta2_cl)


# combine matrices to add habitat amount

rare2_beta2_cl_logR <- left_join(rare2_beta2_cl_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare2_beta2_cl <- rma.uni(yi, vi,
                                      data = rare2_beta2_cl_logR,
                                      method = "REML") # with random effect 

summary(mod_overall_rare2_beta2_cl) # reported results 

# tables for plotting later

overall_rare2_beta2_cl_res <- mod_results(mod_overall_rare2_beta2_cl, mod = "1", group = "refshort")
overall_rare2_beta2_cl_est <- overall_rare2_beta2_cl_res$mod_table %>% 
  dplyr::mutate(diversity = "Beta", 
                type = "Nearest pairs (q = 2)",
                .before = name)
overall_rare2_beta2_cl_dat <- overall_rare2_beta2_cl_res$data %>% 
  dplyr::mutate(diversity = "Beta",
                type = "Nearest pairs (q = 2)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare2_beta2_cl, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)



### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare2_beta2_cl <- rma.uni(yi, vi,
                                    data = rare2_beta2_cl_logR,
                                    mods = ~ ha_diff,
                                    method = "REML") # with random effect 

summary(mod_mods1_rare2_beta2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare2_beta2_cl <- rma.uni(yi, vi,
                                    data = rare2_beta2_cl_logR,
                                    mods = ~ ha_class -1,
                                    method = "REML") # with random effect 

summary(mod_mods2_rare2_beta2_cl) # reported results 

orchaRd::orchard_plot(mod_mods2_rare2_beta2_cl, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare2_beta2_cl_logR$continent_bin <- ifelse(rare2_beta2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_rare2_beta2_cl <- rma.uni(yi, vi,
                                    data = rare2_beta2_cl_logR,
                                    mods = ~ continent_bin - 1,
                                    method = "REML") # with random effect 

summary(mod_mods3_rare2_beta2_cl) # reported results 

orchaRd::orchard_plot(mod_mods3_rare2_beta2_cl, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare2_beta2_cl <- rma.uni(yi, vi,
                                    data = rare2_beta2_cl_logR,
                                    mods = ~ time_since_fragmentation - 1,
                                    method = "REML") # with random effect 

summary(mod_mods4_rare2_beta2_cl) # reported results 

orchaRd::orchard_plot(mod_mods4_rare2_beta2_cl, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare2_beta2_cl <- influence.rma.uni(mod_overall_rare2_beta2_cl)
plot(mod_infl_rare2_beta2_cl) # study 26 is an outlier

# calculate the inverse of effective sample size
rare2_beta2_cl_logR$inv_n <- with(rare2_beta2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare2_beta2_cl_logR$sqrt_inv_n <- with(rare2_beta2_cl_logR, sqrt(inv_n))

# Time-lag bias test

rare2_beta2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare2_beta2_cl_logR$refshort))
rare2_beta2_cl_logR$year_c <- as.vector(scale(rare2_beta2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare2_beta2_cl_check <- rma(yi, vi,
                                mods = ~1 + sqrt_inv_n + year_c,
                                data = rare2_beta2_cl_logR,
                                method = "REML")
summary(mod_rare2_beta2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare2_beta2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Beta", 
         type = "Nearest pairs (q = 2)",
         refshort = rare2_beta2_cl_logR$refshort,
         continent_bin = rare2_beta2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare2_beta2_cl
leave1out_mod_rare2_beta2_cl


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare2_beta2_cl <- extract_model_results(mod_mods1_rare2_beta2_cl, 
                                                      diversity_index = "Beta (q=2)", 
                                                      diversity_type = "Nearest pairs") 


results_mods2_rare2_beta2_cl <-extract_model_results(mod_mods2_rare2_beta2_cl, 
                                                     diversity_index = "Beta (q=2)", 
                                                     diversity_type = "Nearest pairs") 


results_mods3_rare2_beta2_cl <- extract_model_results(mod_mods3_rare2_beta2_cl, 
                                                      diversity_index = "Beta (q=2)", 
                                                      diversity_type = "Nearest pairs") 



results_mods4_rare2_beta2_cl <- extract_model_results(mod_mods4_rare2_beta2_cl, 
                                                      diversity_index = "Beta (q=2)", 
                                                      diversity_type = "Nearest pairs") 

results_rare2_beta2_cl_check <- extract_model_results(mod_rare2_beta2_cl_check, 
                                                      diversity_index = "Beta (q=2)", 
                                                      diversity_type = "Nearest pairs") 


mod_results_rare2_beta2_cl <- rbind.data.frame(results_mods1_rare2_beta2_cl, 
                                               results_mods2_rare2_beta2_cl, 
                                               results_mods3_rare2_beta2_cl,
                                               results_mods4_rare2_beta2_cl,
                                               results_rare2_beta2_cl_check)



### Rarefied GAMMA diversity (nearest pairs) - q = 0


# required data frames

frag_cont_rare_div2_cl %>% 
  filter(diversity_index == "gamma" & buffer == 2000 & q_order == "q = 0") -> rare0_gamma2_cl

transformed_data_rare0_gamma2_cl <- rare0_gamma2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare0_gamma2_cl_logR <- escalc(measure="ROM", 
                               m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                               sd1i = continuous_sd, sd2i = fragmented_sd,
                               n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                               data = transformed_data_rare0_gamma2_cl)


# combine matrices to add habitat amount

rare0_gamma2_cl_logR <- left_join(rare0_gamma2_cl_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare0_gamma2_cl <- rma.uni(yi, vi,
                                       data = rare0_gamma2_cl_logR,
                                       method = "REML") # with random effect 

summary(mod_overall_rare0_gamma2_cl) # reported results 

# tables for plotting later

overall_rare0_gamma2_cl_res <- mod_results(mod_overall_rare0_gamma2_cl, mod = "1", group = "refshort")
overall_rare0_gamma2_cl_est <- overall_rare0_gamma2_cl_res$mod_table %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "Nearest pairs (q = 0)",
                .before = name)
overall_rare0_gamma2_cl_dat <- overall_rare0_gamma2_cl_res$data %>% 
  dplyr::mutate(diversity = "Gamma",
                type = "Nearest pairs (q = 0)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare0_gamma2_cl, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare0_gamma2_cl <- influence.rma.uni(mod_overall_rare0_gamma2_cl)
plot(mod_infl_rare0_gamma2_cl) # study 26 is an outlier

# calculate the inverse of effective sample size
rare0_gamma2_cl_logR$inv_n <- with(rare0_gamma2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare0_gamma2_cl_logR$sqrt_inv_n <- with(rare0_gamma2_cl_logR, sqrt(inv_n))

# Time-lag bias test

rare0_gamma2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare0_gamma2_cl_logR$refshort))
rare0_gamma2_cl_logR$year_c <- as.vector(scale(rare0_gamma2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare0_gamma2_cl_check <- rma(yi, vi,
                                 mods = ~1 + sqrt_inv_n + year_c,
                                 data = rare0_gamma2_cl_logR,
                                 method = "REML")
summary(mod_rare0_gamma2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

### Add moderators 


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare0_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare0_gamma2_cl_logR,
                                     mods = ~ ha_diff,
                                     method = "REML") # with random effect 

summary(mod_mods1_rare0_gamma2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare0_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare0_gamma2_cl_logR,
                                     mods = ~ ha_class -1,
                                     method = "REML") # with random effect 

summary(mod_mods2_rare0_gamma2_cl) # reported results 

orchaRd::orchard_plot(mod_mods2_rare0_gamma2_cl, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare0_gamma2_cl_logR$continent_bin <- ifelse(rare0_gamma2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_rare0_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare0_gamma2_cl_logR,
                                     mods = ~ continent_bin - 1,
                                     method = "REML") # with random effect 

summary(mod_mods3_rare0_gamma2_cl) # reported results 

orchaRd::orchard_plot(mod_mods3_rare0_gamma2_cl, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare0_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare0_gamma2_cl_logR,
                                     mods = ~ time_since_fragmentation - 1,
                                     method = "REML") # with random effect 

summary(mod_mods4_rare0_gamma2_cl) # reported results 

orchaRd::orchard_plot(mod_mods4_rare0_gamma2_cl, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare0_gamma2_cl <- influence.rma.uni(mod_overall_rare0_gamma2_cl)
plot(mod_infl_rare0_gamma2_cl) # study 26 is an outlier

# calculate the inverse of effective sample size
rare0_gamma2_cl_logR$inv_n <- with(rare0_gamma2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare0_gamma2_cl_logR$sqrt_inv_n <- with(rare0_gamma2_cl_logR, sqrt(inv_n))

# Time-lag bias test

rare0_gamma2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare0_gamma2_cl_logR$refshort))
rare0_gamma2_cl_logR$year_c <- as.vector(scale(rare0_gamma2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare0_gamma2_cl_check <- rma(yi, vi,
                                 mods = ~1 + sqrt_inv_n + year_c,
                                 data = rare0_gamma2_cl_logR,
                                 method = "REML")
summary(mod_rare0_gamma2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare0_gamma2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Gamma", 
         type = "Nearest pairs (q = 0)",
         refshort = rare0_gamma2_cl_logR$refshort,
         continent_bin = rare0_gamma2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare0_gamma2_cl
leave1out_mod_rare0_gamma2_cl

## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare0_gamma2_cl <- extract_model_results(mod_mods1_rare0_gamma2_cl, 
                                                       diversity_index = "Gamma (q=0)", 
                                                       diversity_type = "Nearest pairs") 


results_mods2_rare0_gamma2_cl <-extract_model_results(mod_mods2_rare0_gamma2_cl, 
                                                      diversity_index = "Gamma (q=0)", 
                                                      diversity_type = "Nearest pairs") 


results_mods3_rare0_gamma2_cl <- extract_model_results(mod_mods3_rare0_gamma2_cl, 
                                                       diversity_index = "Gamma (q=0)", 
                                                       diversity_type = "Nearest pairs") 



results_mods4_rare0_gamma2_cl <- extract_model_results(mod_mods4_rare0_gamma2_cl, 
                                                       diversity_index = "Gamma (q=0)", 
                                                       diversity_type = "Nearest pairs") 

results_rare0_gamma2_cl_check <- extract_model_results(mod_rare0_gamma2_cl_check, 
                                                       diversity_index = "Gamma (q=0)", 
                                                       diversity_type = "Nearest pairs") 


mod_results_rare0_gamma2_cl <- rbind.data.frame(results_mods1_rare0_gamma2_cl, 
                                                results_mods2_rare0_gamma2_cl, 
                                                results_mods3_rare0_gamma2_cl,
                                                results_mods4_rare0_gamma2_cl,
                                                results_rare0_gamma2_cl_check)


### Rarefied GAMMA diversity (nearest pairs) - q = 2


# required data frames

frag_cont_rare_div2_cl %>% 
  filter(diversity_index == "gamma" & buffer == 2000 & q_order == "q = 2") -> rare2_gamma2_cl

transformed_data_rare2_gamma2_cl <- rare2_gamma2_cl %>% 
  dplyr::select(refshort, patch_type, diversity_value, sd, n_pairs) %>% 
  pivot_wider(names_from = patch_type, 
              values_from = c(diversity_value, sd, n_pairs),
              names_glue = "{patch_type}_{.value}")


# Calculate log transformed ratio of means (Hedges et al., 1999; Lajeunesse, 2011)

rare2_gamma2_cl_logR <- escalc(measure="ROM", 
                               m1i = continuous_diversity_value, m2i = fragmented_diversity_value,
                               sd1i = continuous_sd, sd2i = fragmented_sd,
                               n1i = continuous_n_pairs, n2i = fragmented_n_pairs,
                               data = transformed_data_rare2_gamma2_cl)


# combine matrices to add habitat amount

rare2_gamma2_cl_logR <- left_join(rare2_gamma2_cl_logR, landscape_info, by = "refshort")


# Meta-Analysis via Linear (Mixed-Effects) Models

mod_overall_rare2_gamma2_cl <- rma.uni(yi, vi,
                                       data = rare2_gamma2_cl_logR,
                                       method = "REML") # with random effect 

summary(mod_overall_rare2_gamma2_cl) # reported results 

# tables for plotting later

overall_rare2_gamma2_cl_res <- mod_results(mod_overall_rare2_gamma2_cl, mod = "1", group = "refshort")
overall_rare2_gamma2_cl_est <- overall_rare2_gamma2_cl_res$mod_table %>% 
  dplyr::mutate(diversity = "Gamma", 
                type = "Nearest pairs (q = 2)",
                .before = name)
overall_rare2_gamma2_cl_dat <- overall_rare2_gamma2_cl_res$data %>% 
  dplyr::mutate(diversity = "Gamma",
                type = "Nearest pairs (q = 2)",
                .before = yi)



## simple plot of the overall effects:
# POSITIVE = continuous landscapes have larger diversity than fragmented landscapes 
# NEGATIVE = fragmented landscapes have larger diversity than continuous landscapes 

orchaRd::orchard_plot(mod_overall_rare2_gamma2_cl, mod = "1", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Add moderators 


## Habitat amount (difference between fragmented and continuous landscapes)

mod_mods1_rare2_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare2_gamma2_cl_logR,
                                     mods = ~ ha_diff,
                                     method = "REML") # with random effect 

summary(mod_mods1_rare2_gamma2_cl) # reported results 

## Habitat amount (quantiles)

mod_mods2_rare2_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare2_gamma2_cl_logR,
                                     mods = ~ ha_class -1,
                                     method = "REML") # with random effect 

summary(mod_mods2_rare2_gamma2_cl) # reported results 

orchaRd::orchard_plot(mod_mods2_rare2_gamma2_cl, mod = "ha_class", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing South America with other continents 

rare2_gamma2_cl_logR$continent_bin <- ifelse(rare2_gamma2_cl_logR$continent == "South America", "South America", "Others")


mod_mods3_rare2_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare2_gamma2_cl_logR,
                                     mods = ~ continent_bin - 1,
                                     method = "REML") # with random effect 

summary(mod_mods3_rare2_gamma2_cl) # reported results 

orchaRd::orchard_plot(mod_mods3_rare2_gamma2_cl, mod = "continent_bin", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


## Comparing the effect sizes based on time since fragmentation (intermediate vs. long)  

mod_mods4_rare2_gamma2_cl <- rma.uni(yi, vi,
                                     data = rare2_gamma2_cl_logR,
                                     mods = ~ time_since_fragmentation - 1,
                                     method = "REML") # with random effect 

summary(mod_mods4_rare2_gamma2_cl) # reported results 

orchaRd::orchard_plot(mod_mods4_rare2_gamma2_cl, mod = "time_since_fragmentation", group = "refshort", xlab = "Log response ratio (LRR)",
                      transfm = "none", twig.size = 0.5, trunk.size = 1)


### Sensitivity Analysis 

# influential study diagnostics 

mod_infl_rare2_gamma2_cl <- influence.rma.uni(mod_overall_rare2_gamma2_cl)
plot(mod_infl_rare2_gamma2_cl) # study 26 is an outlier

# calculate the inverse of effective sample size
rare2_gamma2_cl_logR$inv_n <- with(rare2_gamma2_cl_logR, (continuous_n_pairs + fragmented_n_pairs)/(continuous_n_pairs * fragmented_n_pairs))
rare2_gamma2_cl_logR$sqrt_inv_n <- with(rare2_gamma2_cl_logR, sqrt(inv_n))

# Time-lag bias test

rare2_gamma2_cl_logR$year <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rare2_gamma2_cl_logR$refshort))
rare2_gamma2_cl_logR$year_c <- as.vector(scale(rare2_gamma2_cl_logR$year, scale = F))

## All-in publication bias test (multi-moderator: sample size and year of publication) - based on Nakagawa et al. 2021

mod_rare2_gamma2_cl_check <- rma(yi, vi,
                                 mods = ~1 + sqrt_inv_n + year_c,
                                 data = rare2_gamma2_cl_logR,
                                 method = "REML")
summary(mod_rare2_gamma2_cl_check) # there is no evidence of small-study effect and time lag bias (multi-moderator bias test)

# Leave-one-out
# note, I added the column refshort, but the row value represents the estimate after "excluding" that refshort

leave1out(mod_overall_rare2_gamma2_cl) %>% 
  as.data.frame() %>% 
  mutate(diversity = "Gamma", 
         type = "Nearest pairs (q = 2)",
         refshort = rare2_gamma2_cl_logR$refshort,
         continent_bin = rare2_gamma2_cl_logR$continent_bin, 
         .before = estimate) -> leave1out_mod_rare2_gamma2_cl
leave1out_mod_rare2_gamma2_cl


## Combined moderators and sensitivity analysis to a final data.frame per diversity index/type


results_mods1_rare2_gamma2_cl <- extract_model_results(mod_mods1_rare2_gamma2_cl, 
                                                       diversity_index = "Gamma (q=2)", 
                                                       diversity_type = "Nearest pairs") 


results_mods2_rare2_gamma2_cl <-extract_model_results(mod_mods2_rare2_gamma2_cl, 
                                                      diversity_index = "Gamma (q=2)", 
                                                      diversity_type = "Nearest pairs") 


results_mods3_rare2_gamma2_cl <- extract_model_results(mod_mods3_rare2_gamma2_cl, 
                                                       diversity_index = "Gamma (q=2)", 
                                                       diversity_type = "Nearest pairs") 



results_mods4_rare2_gamma2_cl <- extract_model_results(mod_mods4_rare2_gamma2_cl, 
                                                       diversity_index = "Gamma (q=2)", 
                                                       diversity_type = "Nearest pairs") 

results_rare2_gamma2_cl_check <- extract_model_results(mod_rare2_gamma2_cl_check, 
                                                       diversity_index = "Gamma (q=2)", 
                                                       diversity_type = "Nearest pairs") 


mod_results_rare2_gamma2_cl <- rbind.data.frame(results_mods1_rare2_gamma2_cl, 
                                                results_mods2_rare2_gamma2_cl, 
                                                results_mods3_rare2_gamma2_cl,
                                                results_mods4_rare2_gamma2_cl,
                                                results_rare2_gamma2_cl_check)


### save results in data.frame to export and prepare figures 

overall_combined_est <- rbind.data.frame(
  overall_alpha2_est,
  overall_beta2_est,
  overall_gamma2_est,
  overall_alpha2_cl_est,
  overall_beta2_cl_est,
  overall_gamma2_cl_est,
  overall_rare0_alpha2_est,
  overall_rare0_beta2_est,
  overall_rare0_gamma2_est,
  overall_rare0_alpha2_cl_est,
  overall_rare0_beta2_cl_est,
  overall_rare0_gamma2_cl_est,
  overall_rare2_alpha2_est,
  overall_rare2_beta2_est,
  overall_rare2_gamma2_est,
  overall_rare2_alpha2_cl_est,
  overall_rare2_beta2_cl_est,
  overall_rare2_gamma2_cl_est)

overall_combined_dat <- rbind.data.frame(
  overall_alpha2_dat,
  overall_beta2_dat,
  overall_gamma2_dat,
  overall_alpha2_cl_dat,
  overall_beta2_cl_dat,
  overall_gamma2_cl_dat,
  overall_rare0_alpha2_dat,
  overall_rare0_beta2_dat,
  overall_rare0_gamma2_dat,
  overall_rare0_alpha2_cl_dat,
  overall_rare0_beta2_cl_dat,
  overall_rare0_gamma2_cl_dat,
  overall_rare2_alpha2_dat,
  overall_rare2_beta2_dat,
  overall_rare2_gamma2_dat,
  overall_rare2_alpha2_cl_dat,
  overall_rare2_beta2_cl_dat,
  overall_rare2_gamma2_cl_dat)

overall_combined_dat$scale <- 1 / sqrt(overall_combined_dat$vi)


mod_results <- rbind.data.frame(
  mod_results_alpha2,
  mod_results_beta2,
  mod_results_gamma2,
  mod_results_alpha2_cl,
  mod_results_beta2_cl,
  mod_results_gamma2_cl,
  mod_results_rare0_alpha2,
  mod_results_rare0_beta2,
  mod_results_rare0_gamma2,
  mod_results_rare0_alpha2_cl,
  mod_results_rare0_beta2_cl,
  mod_results_rare0_gamma2_cl,
  mod_results_rare2_alpha2,
  mod_results_rare2_beta2,
  mod_results_rare2_gamma2,
  mod_results_rare2_alpha2_cl,
  mod_results_rare2_beta2_cl,
  mod_results_rare2_gamma2_cl)


# write.csv(mod_results, "05_model_results/meta_an_mod_results.csv")

### Figures

### Figure 2


## set desired order of diversity analysis

overall_combined_est$type_ord <- factor(overall_combined_est$type,
                                        levels = c("Nearest pairs (q = 2)", "All pairs (q = 2)",
                                                   "Nearest pairs (q = 0)", "All pairs (q = 0)",
                                                   "Nearest pairs", "All pairs"))
overall_combined_dat$type_ord <- factor(overall_combined_dat$type,
                                        levels = c("Nearest pairs (q = 2)", "All pairs (q = 2)",
                                                   "Nearest pairs (q = 0)", "All pairs (q = 0)",
                                                   "Nearest pairs", "All pairs"))


## plot 

### order by diversity

overall_combined_est$type_div <- factor(overall_combined_est$type,
                                        levels = c("Nearest pairs (q = 2)", "Nearest pairs (q = 0)",  "Nearest pairs",
                                                   "All pairs (q = 2)", "All pairs (q = 0)", "All pairs"))
overall_combined_dat$type_div <- factor(overall_combined_dat$type,
                                        levels = c("Nearest pairs (q = 2)", "Nearest pairs (q = 0)",  "Nearest pairs",
                                                   "All pairs (q = 2)", "All pairs (q = 0)", "All pairs"))


## plot 

overall_combined_est %>%
  ggplot(aes(x = estimate, y = type_div, group = diversity, color = diversity)) +
  geom_rect(aes(ymin = 0.5, ymax = 1.5, xmin = -Inf, xmax = Inf), color = "white", fill = "lightgrey", alpha = 0.03) + 
  geom_rect(aes(ymin = 2.5, ymax = 3.5, xmin = -Inf, xmax = Inf), color = "white", fill = "lightgrey", alpha = 0.03)+
  geom_rect(aes(ymin = 4.5, ymax = 5.5, xmin = -Inf, xmax = Inf), color = "white", fill = "lightgrey", alpha = 0.03)+
  ggbeeswarm::geom_quasirandom(data = overall_combined_dat, 
                               aes(x = yi, y = type_div, group = diversity, color = diversity, size = scale),
                               alpha = 0.3, dodge.width=.9) +
  geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL), position=position_dodge(width=.9), height = 0, linewidth = 1, color = "black") +
  geom_vline(xintercept = 0, linetype = 3, colour = "black") +
  geom_point(aes(fill = diversity), size = 3, shape = 21, color = "black", position=position_dodge(width=.9),  show.legend = F) +
  theme_bw() +
  scale_color_manual(values = c("#ddaa33", "#bb5566", "#004488"))+
  scale_fill_manual(values = c("#ddaa33", "#bb5566", "#004488"))+
  xlab("Log response ratio (LRR)") + 
  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position.inside = c(1, 0),
        legend.justification = c(1, 0),
        legend.title = element_text(size = 9),
        legend.direction="horizontal",
        legend.background = element_blank(),
        axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Figure2


ggsave("03_simple_unedited_figures/Figure2.pdf",
       Figure2,
       dpi = 300,
       height = 14.85,
       width = 21,
       units = "cm")



ggsave("03_simple_unedited_figures/Figure2.png",
       Figure2,
       dpi = 300,
       height = 14.85,
       width = 21,
       units = "cm")

ggsave("03_simple_unedited_figures/Figure2.svg",
       Figure2,
       dpi = 300,
       height = 14.85,
       width = 21,
       units = "cm")



### Back transform LRR to estimate % differences between the mean diversity in continuous and fragmented landscapes

overall_combined_est %>% 
  dplyr::mutate(proportional_change = (exp(estimate) - 1) * 100) -> overall_combined_est_p


overall_combined_est_p %>% 
  dplyr::select(diversity, type, proportional_change)

overall_combined_est_p %>% 
  dplyr::select(diversity, type, proportional_change) %>% 
  dplyr::filter(diversity == "Alpha")



overall_combined_est_p %>% 
  dplyr::select(diversity, proportional_change) %>% 
  group_by(diversity) %>% 
  dplyr::summarise(mean(proportional_change))



### Figure "leave one out"

leave1oust_combined$metric

leave1oust_combined <- rbind.data.frame(
  leave1out_mod_alpha2,
  leave1out_mod_beta2,
  leave1out_mod_gamma2,
  leave1out_mod_alpha2_cl,
  leave1out_mod_beta2_cl,
  leave1out_mod_gamma2_cl,
  leave1out_mod_rare0_alpha2,
  leave1out_mod_rare0_beta2,
  leave1out_mod_rare0_gamma2,
  leave1out_mod_rare0_alpha2_cl,
  leave1out_mod_rare0_beta2_cl,
  leave1out_mod_rare0_gamma2_cl,
  leave1out_mod_rare2_alpha2,
  leave1out_mod_rare2_beta2,
  leave1out_mod_rare2_gamma2,
  leave1out_mod_rare2_alpha2_cl,
  leave1out_mod_rare2_beta2_cl,
  leave1out_mod_rare2_gamma2_cl)

leave1oust_combined_summ <- summarySE(leave1oust_combined, 
                                      measurevar = "estimate",
                                      groupvars = c("continent_bin", "type", "diversity"))


leave1oust_combined_summ %>% 
  ggplot(aes(x = continent_bin, y = estimate)) +
  geom_linerange(aes(ymax = estimate + sd, ymin = estimate - sd, color = continent_bin), size = 4, alpha = 0.3)+
  geom_pointrange(aes(ymax = estimate + ci, ymin = estimate - ci, color = continent_bin), size = 0.4) + 
  facet_grid(type~diversity, scales = "free")+
  geom_hline(data = overall_combined_est, mapping = aes(yintercept = estimate), linetype=3) +
  scale_color_manual(values = c("#1b7837", "#762a83")) + 
  # scale_fill_manual(values = c("#1b7837", "#762a83")) + 
  coord_flip() + 
  ylab("Estimate value after 'leave-one-out' ") +
  guides(fill = "none", colour = "none", size = "none") +
  theme_bw()+
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.title = element_text(size = 9),
        legend.direction="horizontal",
        legend.background = element_blank(),
        axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold", angle=0),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) -> figure_leave1out

ggsave("03_simple_unedited_figures/figure_leave1out.pdf",
       figure_leave1out,
       dpi = 300,
       height = 15,
       width = 30,
       units = "cm")


ggsave("03_simple_unedited_figures/figure_leave1out.png",
       figure_leave1out,
       dpi = 300,
       height = 15,
       width = 30,
       units = "cm")


ggsave("03_simple_unedited_figures/figure_leave1out.svg",
       figure_leave1out,
       dpi = 300,
       height = 15,
       width = 30,
       units = "cm")

