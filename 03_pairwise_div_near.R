# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: Estimate pairwise diversity using the nearest plot pairs
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

# Data --------------------------------------------------------------------

landfrag_abundance <- read.csv("data/landfrag_37sset_abundance.csv")
landfrag_landscape <- read.csv("data/landfrag_37sset_landscapes.csv")

# create an object with coordinates to calculate distance between plots or fragments  

landfrag_coords <- landfrag_landscape %>% 
  dplyr::select(refshort, fragment_id, plot_id, patch_type, latitude, longitude)


# Calculating pairwise diversity ("diversity of 2") per site ---------------------------------------------------------------------
# this includes only the nearest pairs within a give study / landscape type


study_list <- unique(landfrag_abundance$refshort)
diversity_of_2_close <- data.frame()

for (i in study_list){
  print(i)
  
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
  
  frag_remove_1sp = specnumber(frag_comm) == 1
  
  frag_comm <- frag_comm[!frag_remove_1sp,]
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
  
  cont_comm <- cont_comm[!cont_remove_1sp,]
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
    dplyr::filter(!is.na(dist)) %>% 
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
    dplyr::filter(!is.na(dist)) %>% 
    dplyr::group_by(ID) %>% 
    dplyr::arrange(dist) %>% 
    dplyr::slice(1)
  cclosest_df = cclosest[,-3]
  cclosest_df = cclosest_df[!duplicated(apply(cclosest_df, 1, function(x) paste(sort(x), collapse=""))),]
  
  # Create the community matrix for the closest pairs 
  
  ## FRAGMENTED  
  
  # create object to save individual values for all pairs 
  frag_alpha2_cl = NULL
  frag_beta2_cl = NULL  
  frag_gamma2_cl = NULL
  
  for(j in 1:nrow(fclosest_df)){
    print(j)
    fpar_n <- as.character(fclosest_df[j,])
    fcomm_of_2_close = frag_comm[rownames(frag_comm) %in% fpar_n,] 
    frag_gamma2_close = specnumber(colSums(fcomm_of_2_close))
    frag_alpha2_close = mean(specnumber(fcomm_of_2_close))
    frag_beta2_close = frag_gamma2_close/frag_alpha2_close
    frag_diversity_of_2_close_temp = data.frame(
      refshort = i,
      patch_type = rep("fragmented", 3),
      diversity_index = c("alpha", "gamma", "beta"),
      diversity_value = c(frag_alpha2_close, frag_gamma2_close, frag_beta2_close)
    )
    
    # save values of all pairs to calculate SD
    frag_alpha2_cl[j] = frag_alpha2_close
    frag_beta2_cl[j] = frag_beta2_close
    frag_gamma2_cl[j] = frag_gamma2_close
    
  }
  
  
  # calculate standard deviation of diversity metrics in fragmented landscapes
  frag_diversity_of_2_close_temp$sd = c(
    sd(frag_alpha2_cl, na.rm = TRUE), 
    sd(frag_gamma2_cl, na.rm = TRUE), 
    sd(frag_beta2_cl, na.rm = TRUE))
  
  frag_diversity_of_2_close_temp$n_pairs = c(length(frag_alpha2_cl), length(frag_gamma2_cl), length(frag_beta2_cl))
  
  ## CONTINUOUS 
  
  # create object to save individual values for all pairs 
  
  cont_alpha2_cl = NULL
  cont_beta2_cl = NULL  
  cont_gamma2_cl = NULL
  
  
  for(k in 1:nrow(cclosest_df)){
    print(k)
    cpar_n <- as.character(cclosest_df[k,])
    ccomm_of_2_close = cont_comm[rownames(cont_comm) %in% cpar_n,] 
    cont_gamma2_close = specnumber(colSums(ccomm_of_2_close))
    cont_alpha2_close = mean(specnumber(ccomm_of_2_close))
    cont_beta2_close = cont_gamma2_close/cont_alpha2_close
    
    cont_diversity_of_2_close_temp = data.frame(
      refshort = i,
      patch_type = rep("continuous", 3),
      diversity_index = c("alpha", "gamma", "beta"),
      diversity_value = c(cont_alpha2_close, cont_gamma2_close, cont_beta2_close)
    )
    
    # save values of all pairs to calculate SD
    cont_alpha2_cl[k] = cont_alpha2_close
    cont_beta2_cl[k] = cont_beta2_close
    cont_gamma2_cl[k] = cont_gamma2_close
    
  }
  
  cont_diversity_of_2_close_temp$sd = c(
    sd(cont_alpha2_cl, na.rm = TRUE), 
    sd(cont_gamma2_cl, na.rm = TRUE), 
    sd(cont_beta2_cl, na.rm = TRUE))
  
  cont_diversity_of_2_close_temp$n_pairs = c(length(cont_alpha2_cl), length(cont_gamma2_cl), length(cont_beta2_cl))
  
  diversity_of_2_close <- rbind.data.frame(diversity_of_2_close, frag_diversity_of_2_close_temp, cont_diversity_of_2_close_temp)    
  
}

## join with median habitat amount data 

diversity_of_2_close <- left_join(diversity_of_2_close, hamount_np_med, by = c("refshort", "patch_type"), relationship = "many-to-many")
diversity_of_2_close$sd <- ifelse(is.na(diversity_of_2_close$sd), 0, diversity_of_2_close$sd) # studies with only 1 pair have a NA sd, so I'm replacing it by 0


## export results to a csv file 
## note, the data is still available in the folder processed_data

# write.csv(diversity_of_2_close, "processed_data/diversity_of_2_close.csv")


