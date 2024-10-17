# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: Estimate pairwise diversity using all plot pairs
# Author: Thiago Gonçalves Souza
# Date: [October 10, 2024]
# Notes: R version 4.4.1 [2024-06-14 ucrt]
# ==============================================================================


# Using renv for package management
renv::restore() # it requires installing the package renv

# Load required packages ----------------------------------------------------

library(vegan)
library(tidyverse)
library(ggplot2)


# Data --------------------------------------------------------------------

landfrag_abundance <- read.csv("data/landfrag_37sset_abundance.csv")
landfrag_landscape <- read.csv("data/landfrag_37sset_landscapes.csv")

# Calculate median habitat amount per site / patch type / buffer size 

hamount_np_med <- landfrag_landscape %>% 
  dplyr::filter(buffer %in% c(1000, 1200, 1400, 1600, 1800, 2000)) %>%
  dplyr::mutate(buffer = as.character(buffer)) %>% 
  dplyr::select(refshort, patch_type, buffer, pland, np) %>%  # pland = habitat amount (it varies between zero and 100)
  dplyr::group_by(refshort, patch_type, buffer) %>% 
  dplyr::summarise(habitat_amount = median(pland, na.rm = TRUE),
                   number_patches = median(np, na.rm = TRUE),
                   .groups = "keep")

## NOTE: always double check whether the analysis are selecting the appropriate buffer size


hamount_np_med <- hamount_np_med %>%
  ungroup() %>% 
  dplyr::filter(buffer == 2000) %>% 
  dplyr::select(-buffer)



# Calculating pairwise diversity ("diversity of 2") per site ---------------------------------------------------------------------
# this includes all pairs within a give study / landscape type

study_list <- unique(landfrag_abundance$refshort)
diversity_of_2 <- data.frame()


for (i in study_list){
  
  cat("Processing:", i, "\n")
  
  # Subsetting main matrix by "ith" study
  
  abund_mat = landfrag_abundance %>% filter(refshort == i)
  abund_vec = abund_mat %>% dplyr::select(fragment_id, plot_id, patch_type, scientific_name, abundance)
  abund_vec2 = abund_vec %>% 
    dplyr::group_by(fragment_id, plot_id, patch_type, scientific_name) %>% 
    dplyr::summarise(abundance = sum(abundance)) %>% 
    ungroup()
  study_i_comp = abund_vec2 %>% 
    dplyr::select(fragment_id, plot_id, patch_type, scientific_name, abundance) %>% 
    pivot_wider(names_from = scientific_name,
                values_from = abundance,
                values_fill = 0,
                values_fn = sum) 
  
  # Creating matrices for continuous and fragmented patches
  
  ## fragments 
  
  frag_comm = study_i_comp %>% filter(patch_type == "fragmented") %>%
    dplyr::mutate(frag_plot = interaction(fragment_id, plot_id), .before = patch_type) %>% 
    dplyr::select(-patch_type, -fragment_id, -plot_id) %>% 
    tibble::column_to_rownames("frag_plot")
  
  # Remove sites/plots with only one species
  
  frag_remove_1sp = specnumber(frag_comm) == 1
  
  frag_comm <- frag_comm[!frag_remove_1sp,]
  frag_comm <- frag_comm[,colSums(frag_comm)>0]
  
  ## continuous
  cont_comm = study_i_comp %>% filter(patch_type == "continuous") %>%
    dplyr::mutate(frag_plot = interaction(fragment_id, plot_id), .before = patch_type) %>% 
    dplyr::select(-patch_type, -fragment_id, -plot_id) %>% 
    tibble::column_to_rownames("frag_plot")
  
  # Remove sites/plots with only one species
  
  cont_remove_1sp = specnumber(cont_comm) == 1
  
  cont_comm <- cont_comm[!cont_remove_1sp,]
  cont_comm <- cont_comm[,colSums(cont_comm)>0]
  
  # Select every combination (without repetition) - fragmented landscape
  
  pair_frag_patches = as.data.frame(t(combn(rownames(frag_comm), 2)))
  colnames(pair_frag_patches) = c("pair1", "pair2")
  
  frag_alpha2_all = NULL
  frag_beta2_all = NULL  
  frag_gamma2_all = NULL
  
  for(j in 1:nrow(pair_frag_patches)){
    print(j)
    frag_par_n <- as.character(pair_frag_patches[j,])
    frag_comm_of_2 = frag_comm[rownames(frag_comm)%in%frag_par_n,] 
    frag_comm_of_2 = frag_comm_of_2[,colSums(frag_comm_of_2)>0]
    frag_gamma2 = specnumber(colSums(frag_comm_of_2))
    frag_alpha2 = mean(specnumber(frag_comm_of_2))
    frag_beta = frag_gamma2/frag_alpha2
    
    frag_diversity_of_2_temp = data.frame(
      refshort = i,
      patch_type = rep("fragmented", 3),
      diversity_index = c("alpha", "gamma", "beta"),
      diversity_value = c(frag_alpha2, frag_gamma2, frag_beta))
    
    # save values of all pairs to calculate SD 
    frag_alpha2_all[j] = frag_alpha2
    frag_beta2_all[j] = frag_beta
    frag_gamma2_all[j] = frag_gamma2
    
    
  }
  
  # calculate standard deviation of diversity metrics in fragmented landscapes
  frag_diversity_of_2_temp$sd = c(
    sd(frag_alpha2_all, na.rm = TRUE), 
    sd(frag_gamma2_all, na.rm = TRUE), 
    sd(frag_beta2_all, na.rm = TRUE))
  
  frag_diversity_of_2_temp$n_pairs = c(length(frag_alpha2_all), length(frag_gamma2_all), length(frag_beta2_all))
  
  
  # Select every combination (without repetition) - continuous landscape
  pair_cont_patches = as.data.frame(t(combn(rownames(cont_comm), 2)))
  colnames(pair_cont_patches) = c("pair1", "pair2")
  
  # create object to save individual values for all pairs 
  cont_alpha2_all = NULL
  cont_beta2_all = NULL  
  cont_gamma2_all = NULL
  
  for(k in 1:nrow(pair_cont_patches)){
    print(k)
    cont_par_n <- as.character(pair_cont_patches[k,])
    cont_comm_of_2 = cont_comm[rownames(cont_comm)%in%cont_par_n,] 
    cont_comm_of_2 = cont_comm_of_2[,colSums(cont_comm_of_2)>0]
    cont_gamma2 = specnumber(colSums(cont_comm_of_2))
    cont_alpha2 = mean(specnumber(cont_comm_of_2))
    cont_beta = cont_gamma2/cont_alpha2
    cont_diversity_of_2_temp = data.frame(
      refshort = i,
      patch_type = rep("continuous", 3),
      diversity_index = c("alpha", "gamma", "beta"),
      diversity_value = c(cont_alpha2, cont_gamma2, cont_beta))
    
    # save values of all pairs to calculate SD 
    cont_alpha2_all[k] = cont_alpha2
    cont_beta2_all[k] = cont_beta
    cont_gamma2_all[k] = cont_gamma2
    
  }
  
  cont_diversity_of_2_temp$sd = c(
    sd(cont_alpha2_all, na.rm = TRUE), 
    sd(cont_gamma2_all, na.rm = TRUE), 
    sd(cont_beta2_all, na.rm = TRUE))
  
  cont_diversity_of_2_temp$n_pairs = c(length(cont_alpha2_all), length(cont_gamma2_all), length(cont_beta2_all))
  
  ## Save all values in a data.frame 
  
  diversity_of_2 <- rbind.data.frame(diversity_of_2, frag_diversity_of_2_temp, cont_diversity_of_2_temp)
}


## Combine the results of the pairwise diversity with the data.frame with habitat amount

diversity_of_2 <- left_join(diversity_of_2, hamount_np_med, by = c("refshort", "patch_type"), relationship = "many-to-many")

diversity_of_2$sd <- ifelse(is.na(diversity_of_2$sd), 0, diversity_of_2$sd) # replace NA by 0 (for studies with only one pair of plots within continuous landscapes)



## export results to a csv file 
## note, the data is still available in the folder processed_data

# write.csv(diversity_of_2, "processed_data/diversity_of_2.csv")



