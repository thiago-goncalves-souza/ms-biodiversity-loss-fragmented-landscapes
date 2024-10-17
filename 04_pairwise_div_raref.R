# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: Estimate pairwise diversity using q = 0 and q = 2 (all pairs and nearest pairs)
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
library(fields)
library(iNEXT)
library(tidyr)
library(tibble)

### Calculating pairwise diversity ("diversity of 2") per site ---------------------------------------------------------------------
### using rarefaction with q = 0 and q = 2 -----------------------
### Diversity of 2 using all plot pairs 

study_list <- unique(landfrag_abundance$refshort)
frag_cont_raref_2 = data.frame()

for (i in study_list){
  print(i)
  
  # Subsetting main matrix by "ith" study
  
  abund_mat = landfrag_abundance %>% filter(refshort == i)
  abund_vec = abund_mat %>% dplyr::select(fragment_id, plot_id, patch_type, scientific_name, abundance)
  abund_vec2 = abund_vec %>% 
    group_by(fragment_id, plot_id, patch_type, scientific_name) %>% 
    dplyr::summarise(abundance = sum(abundance), .groups = "keep") %>% 
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
  
  frag_remove_1sp = specnumber(frag_comm) == 1
  frag_comm <- frag_comm[!frag_remove_1sp,]  # remove sites with 1 species
  
  frag_remove_2ind = rowSums(frag_comm) > 2 
  frag_comm <- frag_comm[frag_remove_2ind,] # remove sites with less than 3 individuals
  
  frag_comm <- frag_comm[,colSums(frag_comm)>0]
  
  ## continuous
  cont_comm = study_i_comp %>% filter(patch_type == "continuous") %>%
    dplyr::mutate(frag_plot = interaction(fragment_id, plot_id), .before = patch_type) %>% 
    dplyr::select(-patch_type, -fragment_id, -plot_id) %>% 
    tibble::column_to_rownames("frag_plot")
  
  cont_remove_1sp = specnumber(cont_comm) == 1
  
  cont_comm <- cont_comm[!cont_remove_1sp,]  # remove sites with 1 species
  
  cont_remove_2ind = rowSums(cont_comm) > 2
  cont_comm <- cont_comm[cont_remove_2ind,] # remove sites with less than 3 individuals
  
  cont_comm <- cont_comm[,colSums(cont_comm)>0]
  
  
  # Data frame with all patches
  
  frag_pair_patches = as.data.frame(t(combn(rownames(frag_comm), 2)))
  colnames(frag_pair_patches) = c("pair1", "pair2")
  
  cont_pair_patches = as.data.frame(t(combn(rownames(cont_comm), 2)))
  colnames(cont_pair_patches) = c("pair1", "pair2")
  
  gamma_of_2_frag_targets = NULL
  alpha_of_2_frag_targets_p1 = NULL
  alpha_of_2_frag_targets_p2 = NULL
  
  # k = 1
  for(k in 1:nrow(frag_pair_patches)){
    paste(i, "-", "frag", print(k))
    
    pair_df = frag_pair_patches[k,]
    pair = as.character(pair_df[1,])
    frag_pair_comm_all = frag_comm[rownames(frag_comm)%in%pair,]
    frag_pair_comm = frag_pair_comm_all[,colSums(frag_pair_comm_all)>0, drop=FALSE]
    gamma_of_2_frag_targets[k] = sum(t(colSums(frag_pair_comm)))
    alpha_of_2_frag_targets_p1[k] = rowSums(frag_pair_comm)[[1]]
    alpha_of_2_frag_targets_p2[k] = rowSums(frag_pair_comm)[[2]]
  }
  
  gamma_of_2_cont_targets = NULL
  alpha_of_2_cont_targets_p1 = NULL
  alpha_of_2_cont_targets_p2 = NULL
  
  # j = 1
  for(j in 1:nrow(cont_pair_patches)){
    print(paste(i, "-", "cont", print(j)))
    
    pair_df = cont_pair_patches[j,]
    pair = as.character(pair_df[1,])
    cont_pair_comm_all = cont_comm[rownames(cont_comm)%in%pair,]
    cont_pair_comm = cont_pair_comm_all[,colSums(cont_pair_comm_all)>0, drop=FALSE]
    gamma_of_2_cont_targets[j] = sum(t(colSums(cont_pair_comm)))
    alpha_of_2_cont_targets_p1[j] = rowSums(cont_pair_comm)[[1]]
    alpha_of_2_cont_targets_p2[j] = rowSums(cont_pair_comm)[[2]]
  }
  
  gamma2_target = min(gamma_of_2_cont_targets, gamma_of_2_frag_targets)
  alpha2_target = min(c(alpha_of_2_cont_targets_p1, alpha_of_2_cont_targets_p2, alpha_of_2_frag_targets_p1, alpha_of_2_frag_targets_p2))
  
  # w = 1
  
  gamma_2_frag_q0 = NULL
  gamma_2_frag_q2 = NULL 
  alpha_2_frag_q0_p1 = NULL
  alpha_2_frag_q0_p2 = NULL
  alpha_2_frag_q2_p1 = NULL
  alpha_2_frag_q2_p2 = NULL
  beta_2_frag_q0 = NULL
  beta_2_frag_q2 = NULL 
  
  for(w in 1:nrow(frag_pair_patches)){
    print(paste(i, "-", "frag_pair", print(w)))
    
    pair_df = frag_pair_patches[w,]
    pair = as.character(pair_df[1,])
    frag_pair_comm_all = frag_comm[rownames(frag_comm)%in%pair,]
    frag_pair_comm = frag_pair_comm_all[,colSums(frag_pair_comm_all)>0]
    
    alpha_2_frag_q0_p1[w] = (estimateD(t(frag_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_frag_q0_p2[w] = (estimateD(t(frag_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    alpha_2_frag_q2_p1[w] = (estimateD(t(frag_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_frag_q2_p2[w] = (estimateD(t(frag_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    gamma_2_frag_q0[w] = estimateD(colSums(frag_pair_comm), q=0, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    gamma_2_frag_q2[w] = estimateD(colSums(frag_pair_comm), q=2, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    
    beta_2_frag_q0[w] = gamma_2_frag_q0[w] / mean(alpha_2_frag_q0_p1[w], alpha_2_frag_q0_p2[w])
    beta_2_frag_q2[w] = gamma_2_frag_q2[w] / mean(alpha_2_frag_q2_p1[w], alpha_2_frag_q2_p2[w])
    
  }
  
  gamma_2_cont_q0 = NULL
  gamma_2_cont_q2 = NULL 
  alpha_2_cont_q0_p1 = NULL
  alpha_2_cont_q0_p2 = NULL
  alpha_2_cont_q2_p1 = NULL
  alpha_2_cont_q2_p2 = NULL
  beta_2_cont_q0 = NULL
  beta_2_cont_q2 = NULL 
  
  # z = 1
  for(z in 1:nrow(cont_pair_patches)){
    print(paste(i, "-", "cont_pair", print(z)))
    pair_df = cont_pair_patches[z,]
    pair = as.character(pair_df[1,])
    cont_pair_comm_all = cont_comm[rownames(cont_comm)%in%pair,]
    cont_pair_comm = cont_pair_comm_all[,colSums(cont_pair_comm_all)>0]
    
    
    alpha_2_cont_q0_p1[z] = (estimateD(t(cont_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_cont_q0_p2[z] = (estimateD(t(cont_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    alpha_2_cont_q2_p1[z] = (estimateD(t(cont_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_cont_q2_p2[z] = (estimateD(t(cont_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    gamma_2_cont_q0[z] = estimateD(colSums(cont_pair_comm), q=0, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    gamma_2_cont_q2[z] = estimateD(colSums(cont_pair_comm), q=2, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    
    beta_2_cont_q0[z] = gamma_2_cont_q0[z] / mean(alpha_2_cont_q0_p1[z], alpha_2_cont_q0_p2[z])
    beta_2_cont_q2[z] = gamma_2_cont_q2[z] / mean(alpha_2_cont_q2_p1[z], alpha_2_cont_q2_p2[z])
    
    
  }
  
  
  # Calculate mean alpha and gamma of 2 
  
  alpha_2_frag_q0_mean = mean(c(alpha_2_frag_q0_p1, alpha_2_frag_q0_p2))
  alpha_2_frag_q2_mean = mean(c(alpha_2_frag_q2_p1, alpha_2_frag_q2_p2))
  
  gamma_2_frag_q0_mean = mean(gamma_2_frag_q0)
  gamma_2_frag_q2_mean = mean(gamma_2_frag_q2)
  
  alpha_2_cont_q0_mean = mean(c(alpha_2_cont_q0_p1, alpha_2_cont_q0_p2))
  alpha_2_cont_q2_mean = mean(c(alpha_2_cont_q2_p1, alpha_2_cont_q2_p2))
  
  gamma_2_cont_q0_mean = mean(gamma_2_cont_q0)
  gamma_2_cont_q2_mean = mean(gamma_2_cont_q2) 
  
  # Calculate mean beta of 2 
  
  beta_2_frag_q0_mean = gamma_2_frag_q0_mean / alpha_2_frag_q0_mean 
  beta_2_cont_q0_mean = gamma_2_cont_q0_mean / alpha_2_cont_q0_mean
  
  beta_2_frag_q2_mean = gamma_2_frag_q2_mean / alpha_2_frag_q2_mean 
  beta_2_cont_q2_mean = gamma_2_cont_q2_mean / alpha_2_cont_q2_mean 
  
  # Calculate sd alpha, beta and gamma of 2 
  
  alpha_2_frag_q0_sd = sd(c(alpha_2_frag_q0_p1, alpha_2_frag_q0_p2))
  alpha_2_frag_q2_sd = sd(c(alpha_2_frag_q2_p1, alpha_2_frag_q2_p2))
  
  gamma_2_frag_q0_sd = sd(gamma_2_frag_q0)
  gamma_2_frag_q2_sd = sd(gamma_2_frag_q2)
  
  
  beta_2_frag_q0_sd = sd(beta_2_frag_q0)
  beta_2_frag_q2_sd = sd(beta_2_frag_q2)
  
  alpha_2_cont_q0_sd = sd(c(alpha_2_cont_q0_p1, alpha_2_cont_q0_p2))
  alpha_2_cont_q2_sd = sd(c(alpha_2_cont_q2_p1, alpha_2_cont_q2_p2))
  
  gamma_2_cont_q0_sd = sd(gamma_2_cont_q0)
  gamma_2_cont_q2_sd = sd(gamma_2_cont_q2) 
  
  beta_2_cont_q0_sd = sd(beta_2_cont_q0)
  beta_2_cont_q2_sd = sd(beta_2_cont_q2)
  
  # Calculate sd alpha, beta and gamma of 2 
  
  frag_npairs = length(gamma_2_frag_q0) # alpha, gamma and beta = same number of pairs 
  cont_npairs = length(gamma_2_cont_q0) # alpha, gamma and beta = same number of pairs 
  
  # Combine all metrics in a temporary file and rbind the results from each loop
  
  frag_cont_raref_2_temp = data.frame(refshort = i,
                                      patch_type = rep(c("continuous", "fragmented"), each = 6),
                                      diversity_index = rep(rep(c("alpha", "gamma", "beta"), each = 2), 2),
                                      q_order = rep(c("q = 0", "q = 2"), 6),
                                      diversity_value = c(alpha_2_cont_q0_mean, alpha_2_cont_q2_mean,
                                                          gamma_2_cont_q0_mean, gamma_2_cont_q2_mean,
                                                          beta_2_cont_q0_mean, beta_2_cont_q2_mean,
                                                          alpha_2_frag_q0_mean, alpha_2_frag_q2_mean,
                                                          gamma_2_frag_q0_mean, gamma_2_frag_q2_mean,
                                                          beta_2_frag_q0_mean, beta_2_frag_q2_mean),
                                      sd = c(alpha_2_cont_q0_sd, alpha_2_cont_q2_sd,
                                             gamma_2_cont_q0_sd, gamma_2_cont_q2_sd,
                                             beta_2_cont_q0_sd, beta_2_cont_q2_sd,
                                             alpha_2_frag_q0_sd, alpha_2_frag_q2_sd,
                                             gamma_2_frag_q0_sd, gamma_2_frag_q2_sd,
                                             beta_2_frag_q0_sd, beta_2_frag_q2_sd),
                                      n_pairs = c(rep(cont_npairs, 6), rep(frag_npairs, 6)))
  
  frag_cont_raref_2 = rbind.data.frame(frag_cont_raref_2, frag_cont_raref_2_temp)
  
}

## join with median habitat amount data 

frag_cont_raref_2 <- left_join(frag_cont_raref_2, hamount_np_med, by = c("refshort", "patch_type"), relationship = "many-to-many")
frag_cont_raref_2$sd <- ifelse(is.na(frag_cont_raref_2$sd), 0, frag_cont_raref_2$sd) # studies with only 1 pair have a NA sd, so I'm replacing it by 0


## export results to a csv file 
## note, the data is still available in the folder processed_data

# write.csv(frag_cont_raref_2, "processed_data/frag_cont_raref_2.csv")


### rarefied with q = 0 and q = 1 
### Diversity of 2 -- closest pairs 

study_list <- unique(landfrag_abundance$refshort)

frag_cont_raref_close_2 = data.frame()


for (i in study_list){
  print(i)
  
  # Subsetting main matrix by "ith" study
  
  abund_mat = landfrag_abundance %>% filter(refshort == i)
  abund_vec = abund_mat %>% dplyr::select(fragment_id, plot_id, patch_type, scientific_name, abundance)
  abund_vec2 = abund_vec %>% 
    group_by(fragment_id, plot_id, patch_type, scientific_name) %>% 
    dplyr::summarise(abundance = sum(abundance), .groups = "keep") %>% 
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
  
  frag_remove_1sp = specnumber(frag_comm) == 1
  frag_comm <- frag_comm[!frag_remove_1sp,]  # remove sites with 1 species
  
  frag_remove_2ind = rowSums(frag_comm) > 2 
  frag_comm <- frag_comm[frag_remove_2ind,] # remove sites with less than 3 individuals
  
  frag_comm <- frag_comm[,colSums(frag_comm)>0]
  
  frag_coords = landfrag_coords %>% filter(refshort == i & patch_type == "fragmented") %>%
    dplyr::mutate(frag_plot = interaction(fragment_id, plot_id), .before = patch_type) %>%
    dplyr::select(frag_plot, latitude, longitude) %>% 
    unique() %>% 
    dplyr::filter(frag_plot %in% rownames(frag_comm)) %>% 
    tibble::column_to_rownames("frag_plot")
  
  ## continuous
  cont_comm = study_i_comp %>% filter(patch_type == "continuous") %>%
    dplyr::mutate(frag_plot = interaction(fragment_id, plot_id), .before = patch_type) %>% 
    dplyr::select(-patch_type, -fragment_id, -plot_id) %>% 
    tibble::column_to_rownames("frag_plot")
  
  cont_remove_1sp = specnumber(cont_comm) == 1
  
  cont_comm <- cont_comm[!cont_remove_1sp,]  # remove sites with 1 species
  
  cont_remove_2ind = rowSums(cont_comm) > 2
  cont_comm <- cont_comm[cont_remove_2ind,] # remove sites with less than 3 individuals
  
  cont_comm <- cont_comm[,colSums(cont_comm)>0]
  
  cont_coords = landfrag_coords %>% filter(refshort == i & patch_type == "continuous") %>%
    dplyr::mutate(frag_plot = interaction(fragment_id, plot_id), .before = patch_type) %>%
    dplyr::select(frag_plot, latitude, longitude) %>% 
    unique() %>% 
    dplyr::filter(frag_plot %in% rownames(cont_comm)) %>% 
    tibble::column_to_rownames("frag_plot")
  
  
  # Distance between every two points (fragments and continuous)
  
  fdist_pts = rdist.earth(frag_coords[,c("longitude", "latitude")], miles = FALSE)
  diag(fdist_pts) <- NA
  fdist_pts_df = as.data.frame(fdist_pts)
  colnames(fdist_pts_df) = rownames(frag_coords)
  rownames(fdist_pts_df) = rownames(frag_coords)
  fclosest = fdist_pts_df %>% 
    dplyr::mutate(ID=rownames(.)) %>% 
    tidyr::gather('closest','dist',-ID) %>% 
    dplyr::filter(!is.na(dist) & dist > 0) %>% 
    dplyr::group_by(ID) %>% 
    dplyr::arrange(dist) %>% 
    dplyr::slice(1)
  fclosest_df = fclosest[,-3]
  fclosest_df = fclosest_df[!duplicated(apply(fclosest_df, 1, function(x) paste(sort(x), collapse=""))),]
  
  cdist_pts = rdist.earth(cont_coords[,c("longitude", "latitude")], miles = FALSE)
  diag(cdist_pts) <- NA
  cdist_pts_df = as.data.frame(cdist_pts)
  colnames(cdist_pts_df) = rownames(cont_coords)
  rownames(cdist_pts_df) = rownames(cont_coords)
  cclosest = cdist_pts_df %>% 
    dplyr::mutate(ID=rownames(.)) %>% 
    tidyr::gather('closest','dist',-ID) %>% 
    dplyr::filter(!is.na(dist) & dist > 0) %>% 
    dplyr::group_by(ID) %>% 
    dplyr::arrange(dist) %>% 
    dplyr::slice(1)
  cclosest_df = cclosest[,-3]
  cclosest_df = cclosest_df[!duplicated(apply(cclosest_df, 1, function(x) paste(sort(x), collapse=""))),]
  
  # Alpha and Gamma Diversity of 2 
  
  frag_pair_patches = fclosest_df # closest fragments 
  colnames(frag_pair_patches) = c("pair1", "pair2")
  
  cont_pair_patches = cclosest_df # closest continuous
  colnames(cont_pair_patches) = c("pair1", "pair2")
  
  gamma_of_2_frag_targets = NULL
  alpha_of_2_frag_targets_p1 = NULL
  alpha_of_2_frag_targets_p2 = NULL
  
  # k = 1
  for(k in 1:nrow(frag_pair_patches)){
    paste(i, "-", "frag", print(k))
    
    pair_df = frag_pair_patches[k,]
    pair = as.character(pair_df[1,])
    frag_pair_comm_all = frag_comm[rownames(frag_comm)%in%pair,]
    frag_pair_comm = frag_pair_comm_all[,colSums(frag_pair_comm_all)>0, drop=FALSE]
    gamma_of_2_frag_targets[k] = sum(t(colSums(frag_pair_comm)))
    alpha_of_2_frag_targets_p1[k] = rowSums(frag_pair_comm)[[1]]
    alpha_of_2_frag_targets_p2[k] = rowSums(frag_pair_comm)[[2]]
  }
  
  gamma_of_2_cont_targets = NULL
  alpha_of_2_cont_targets_p1 = NULL
  alpha_of_2_cont_targets_p2 = NULL
  
  # j = 1
  for(j in 1:nrow(cont_pair_patches)){
    print(paste(i, "-", "cont", print(j)))
    
    pair_df = cont_pair_patches[j,]
    pair = as.character(pair_df[1,])
    cont_pair_comm_all = cont_comm[rownames(cont_comm)%in%pair,]
    cont_pair_comm = cont_pair_comm_all[,colSums(cont_pair_comm_all)>0, drop=FALSE]
    gamma_of_2_cont_targets[j] = sum(t(colSums(cont_pair_comm)))
    alpha_of_2_cont_targets_p1[j] = rowSums(cont_pair_comm)[[1]]
    alpha_of_2_cont_targets_p2[j] = rowSums(cont_pair_comm)[[2]]
  }
  
  gamma2_target = min(gamma_of_2_cont_targets, gamma_of_2_frag_targets)
  alpha2_target = min(c(alpha_of_2_cont_targets_p1, alpha_of_2_cont_targets_p2, alpha_of_2_frag_targets_p1, alpha_of_2_frag_targets_p2))
  
  # w = 1
  
  gamma_2_frag_q0 = NULL
  gamma_2_frag_q2 = NULL 
  alpha_2_frag_q0_p1 = NULL
  alpha_2_frag_q0_p2 = NULL
  alpha_2_frag_q2_p1 = NULL
  alpha_2_frag_q2_p2 = NULL
  beta_2_frag_q0 = NULL
  beta_2_frag_q2 = NULL 
  
  for(w in 1:nrow(frag_pair_patches)){
    print(paste(i, "-", "frag_pair", print(w)))
    
    pair_df = frag_pair_patches[w,]
    pair = as.character(pair_df[1,])
    frag_pair_comm_all = frag_comm[rownames(frag_comm)%in%pair,]
    frag_pair_comm = frag_pair_comm_all[,colSums(frag_pair_comm_all)>0]
    
    alpha_2_frag_q0_p1[w] = (estimateD(t(frag_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_frag_q0_p2[w] = (estimateD(t(frag_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    alpha_2_frag_q2_p1[w] = (estimateD(t(frag_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_frag_q2_p2[w] = (estimateD(t(frag_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    gamma_2_frag_q0[w] = estimateD(colSums(frag_pair_comm), q=0, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    gamma_2_frag_q2[w] = estimateD(colSums(frag_pair_comm), q=2, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    
    beta_2_frag_q0[w] = gamma_2_frag_q0[w] / mean(alpha_2_frag_q0_p1[w], alpha_2_frag_q0_p2[w])
    beta_2_frag_q2[w] = gamma_2_frag_q2[w] / mean(alpha_2_frag_q2_p1[w], alpha_2_frag_q2_p2[w])
    
  }
  
  gamma_2_cont_q0 = NULL
  gamma_2_cont_q2 = NULL 
  alpha_2_cont_q0_p1 = NULL
  alpha_2_cont_q0_p2 = NULL
  alpha_2_cont_q2_p1 = NULL
  alpha_2_cont_q2_p2 = NULL
  beta_2_cont_q0 = NULL
  beta_2_cont_q2 = NULL 
  
  # z = 1
  for(z in 1:nrow(cont_pair_patches)){
    print(paste(i, "-", "cont_pair", print(z)))
    pair_df = cont_pair_patches[z,]
    pair = as.character(pair_df[1,])
    cont_pair_comm_all = cont_comm[rownames(cont_comm)%in%pair,]
    cont_pair_comm = cont_pair_comm_all[,colSums(cont_pair_comm_all)>0]
    
    
    alpha_2_cont_q0_p1[z] = (estimateD(t(cont_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_cont_q0_p2[z] = (estimateD(t(cont_pair_comm), q=0, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    alpha_2_cont_q2_p1[z] = (estimateD(t(cont_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[1]
    alpha_2_cont_q2_p2[z] = (estimateD(t(cont_pair_comm), q=2, datatype="abundance", base="size", level=alpha2_target) %>% pull(qD))[2]
    
    gamma_2_cont_q0[z] = estimateD(colSums(cont_pair_comm), q=0, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    gamma_2_cont_q2[z] = estimateD(colSums(cont_pair_comm), q=2, datatype="abundance", base="size", level=gamma2_target) %>% pull(qD)
    
    beta_2_cont_q0[z] = gamma_2_cont_q0[z] / mean(alpha_2_cont_q0_p1[z], alpha_2_cont_q0_p2[z])
    beta_2_cont_q2[z] = gamma_2_cont_q2[z] / mean(alpha_2_cont_q2_p1[z], alpha_2_cont_q2_p2[z])
    
    
  }
  
  
  # Calculate mean alpha and gamma of 2 
  
  alpha_2_frag_q0_mean = mean(c(alpha_2_frag_q0_p1, alpha_2_frag_q0_p2))
  alpha_2_frag_q2_mean = mean(c(alpha_2_frag_q2_p1, alpha_2_frag_q2_p2))
  
  gamma_2_frag_q0_mean = mean(gamma_2_frag_q0)
  gamma_2_frag_q2_mean = mean(gamma_2_frag_q2)
  
  alpha_2_cont_q0_mean = mean(c(alpha_2_cont_q0_p1, alpha_2_cont_q0_p2))
  alpha_2_cont_q2_mean = mean(c(alpha_2_cont_q2_p1, alpha_2_cont_q2_p2))
  
  gamma_2_cont_q0_mean = mean(gamma_2_cont_q0)
  gamma_2_cont_q2_mean = mean(gamma_2_cont_q2) 
  
  # Calculate mean beta of 2 
  
  beta_2_frag_q0_mean = gamma_2_frag_q0_mean / alpha_2_frag_q0_mean 
  beta_2_cont_q0_mean = gamma_2_cont_q0_mean / alpha_2_cont_q0_mean
  
  beta_2_frag_q2_mean = gamma_2_frag_q2_mean / alpha_2_frag_q2_mean 
  beta_2_cont_q2_mean = gamma_2_cont_q2_mean / alpha_2_cont_q2_mean 
  
  # Calculate sd alpha, beta and gamma of 2 
  
  alpha_2_frag_q0_sd = sd(c(alpha_2_frag_q0_p1, alpha_2_frag_q0_p2))
  alpha_2_frag_q2_sd = sd(c(alpha_2_frag_q2_p1, alpha_2_frag_q2_p2))
  
  gamma_2_frag_q0_sd = sd(gamma_2_frag_q0)
  gamma_2_frag_q2_sd = sd(gamma_2_frag_q2)
  
  
  beta_2_frag_q0_sd = sd(beta_2_frag_q0)
  beta_2_frag_q2_sd = sd(beta_2_frag_q2)
  
  alpha_2_cont_q0_sd = sd(c(alpha_2_cont_q0_p1, alpha_2_cont_q0_p2))
  alpha_2_cont_q2_sd = sd(c(alpha_2_cont_q2_p1, alpha_2_cont_q2_p2))
  
  gamma_2_cont_q0_sd = sd(gamma_2_cont_q0)
  gamma_2_cont_q2_sd = sd(gamma_2_cont_q2) 
  
  beta_2_cont_q0_sd = sd(beta_2_cont_q0)
  beta_2_cont_q2_sd = sd(beta_2_cont_q2)
  
  # Calculate sd alpha, beta and gamma of 2 
  
  frag_npairs = length(gamma_2_frag_q0) # alpha, gamma and beta = same number of pairs 
  cont_npairs = length(gamma_2_cont_q0) # alpha, gamma and beta = same number of pairs 
  
  # Combine all metrics in a temporary file and rbind the results from each loop
  
  frag_cont_raref_close_2_temp = data.frame(refshort = i,
                                            patch_type = rep(c("continuous", "fragmented"), each = 6),
                                            diversity_index = rep(rep(c("alpha", "gamma", "beta"), each = 2), 2),
                                            q_order = rep(c("q = 0", "q = 2"), 6),
                                            diversity_value = c(alpha_2_cont_q0_mean, alpha_2_cont_q2_mean,
                                                                gamma_2_cont_q0_mean, gamma_2_cont_q2_mean,
                                                                beta_2_cont_q0_mean, beta_2_cont_q2_mean,
                                                                alpha_2_frag_q0_mean, alpha_2_frag_q2_mean,
                                                                gamma_2_frag_q0_mean, gamma_2_frag_q2_mean,
                                                                beta_2_frag_q0_mean, beta_2_frag_q2_mean),
                                            sd = c(alpha_2_cont_q0_sd, alpha_2_cont_q2_sd,
                                                   gamma_2_cont_q0_sd, gamma_2_cont_q2_sd,
                                                   beta_2_cont_q0_sd, beta_2_cont_q2_sd,
                                                   alpha_2_frag_q0_sd, alpha_2_frag_q2_sd,
                                                   gamma_2_frag_q0_sd, gamma_2_frag_q2_sd,
                                                   beta_2_frag_q0_sd, beta_2_frag_q2_sd),
                                            n_pairs = c(rep(cont_npairs, 6), rep(frag_npairs, 6)))
  
  frag_cont_raref_close_2 = rbind.data.frame(frag_cont_raref_close_2, frag_cont_raref_close_2_temp)
  
  
}

## join with median habitat amount data 

frag_cont_raref_close_2 <- left_join(frag_cont_raref_close_2, hamount_np_med, by = c("refshort", "patch_type"), relationship = "many-to-many")
frag_cont_raref_close_2$sd <- ifelse(is.na(frag_cont_raref_close_2$sd), 0, frag_cont_raref_close_2$sd) # studies with only 1 pair have a NA sd, so I'm replacing it by 0



## export results to a csv file 
## note, the data is still available in the folder processed_data

# write.csv(frag_cont_raref_close_2, "processed_data/frag_cont_raref_close_2.csv")


