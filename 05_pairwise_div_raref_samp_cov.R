# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: Estimate alpha diversity based on study design
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
library(mobr)
library(iNEXT)
library(purrr)
library(sf)
library(fields)
library(usedist)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
source("Engel_novel_code.R") # from Engel et al. 2021, 10.1002/ecs2.3745

# Data -------------------------------------------------------------------------

landfrag_abundance <- read.csv("data/landfrag_37sset_abundance.csv")
diversity_raref <- read.csv("processed_data/frag_cont_raref_2.csv")

### For loops --- estimated alpha, beta and gamma = WITH SAMPLE COVERAGE -------
## q = 0, q = 2 ---- running with all pairs only 

gamma_of_2_frag_values = data.frame()
gamma_of_2_cont_values = data.frame()

frag_cont_comp = data.frame()

frag_studies <- unique(landfrag_abundance$refshort)

for (i in frag_studies){
  print(i)
  
  # Subsetting main matrix by "ith" study
  
  abund_mat = landfrag_abundance %>% filter(refshort == i)
  abund_vec = abund_mat %>% dplyr::select(fragment_plot_comb, patch_type, scientific_name, abundance)
  abund_vec2 = abund_vec %>% 
    group_by(fragment_plot_comb, patch_type, scientific_name) %>% 
    dplyr::summarise(abundance = sum(abundance)) %>% 
    ungroup()
  study_i_comp = abund_vec2 %>% 
    dplyr::select(fragment_plot_comb, patch_type, scientific_name, abundance) %>% 
    pivot_wider(names_from = scientific_name,
                values_from = abundance,
                values_fill = 0,
                values_fn = sum) 
  
  # Creating matrices for continuous and fragmented patches
  
  frag_comm = study_i_comp %>% filter(patch_type == "fragmented") %>%
    dplyr::select(-patch_type) %>% 
    tibble::column_to_rownames("fragment_plot_comb")
  
  frag_remove_1sp = specnumber(frag_comm) == 1
  
  frag_comm <- frag_comm[!frag_remove_1sp,]
  
  cont_comm = study_i_comp %>% filter(patch_type == "continuous") %>%
    dplyr::select(-patch_type) %>% 
    tibble::column_to_rownames("fragment_plot_comb")
  
  cont_remove_1sp = specnumber(cont_comm) == 1
  
  cont_comm <- cont_comm[!cont_remove_1sp,]
  cont_comm <- cont_comm[,colSums(cont_comm)>0]
  
  # Defining target coverage by comparing fragments and continuous patches
  
  target_coverage = min(C_target(frag_comm), C_target(cont_comm))
  target_coverage_g = min(C_target(t(colSums(frag_comm))), C_target(t(colSums(cont_comm))))
  
  # Changing from data.frame to list as required by iNEXT functions 
  
  new_frag_list <- c(rownames(frag_comm), rownames(cont_comm)) # just in case there are fragments with only one species
  
  frag_sites = study_i_comp %>% dplyr::filter(fragment_plot_comb %in% new_frag_list) %>%  dplyr::select(-patch_type) %>% column_to_rownames("fragment_plot_comb") 
  frag_sites2  = as.list(as.data.frame(t(frag_sites)))
  frag_ish = lapply(frag_sites2, function(x) {sort(x [x>0], decreasing = TRUE)})
  
  # Alpha diversity by patch (based on target_coverage) -> q = 0
  inext_rare <- estimateD(frag_ish, q=0, datatype="abundance", base="coverage", level=target_coverage)
  alpha_rare <- inext_rare %>% dplyr::select(Assemblage, qD) %>% dplyr::rename(fragment_plot_comb = Assemblage)
  alpha_rare_id = study_i_comp %>% dplyr::filter(fragment_plot_comb %in% new_frag_list) %>% dplyr::select(fragment_plot_comb, patch_type) %>% left_join(alpha_rare, by = "fragment_plot_comb")
  aver_alpha = alpha_rare_id %>% dplyr::select(-fragment_plot_comb) %>% dplyr::group_by(patch_type) %>% dplyr::summarise(alpha0 = mean(qD))
  alpha_C_cont_q0 = aver_alpha %>% dplyr::filter(patch_type=="continuous") %>% pull(alpha0)
  alpha_C_frag_q0 = aver_alpha %>% dplyr::filter(patch_type=="fragmented") %>% pull(alpha0)
  
  # Alpha diversity by patch (based on target_coverage) -> q = 2
  inext_rare2 <- estimateD(frag_ish, q=2, datatype="abundance", base="coverage", level=target_coverage)
  alpha_rare2 <- inext_rare2 %>% dplyr::select(Assemblage, qD) %>% dplyr::rename(fragment_plot_comb = Assemblage)
  alpha_rare_id2 = study_i_comp %>% dplyr::filter(fragment_plot_comb %in% new_frag_list) %>% dplyr::select(fragment_plot_comb, patch_type) %>% left_join(alpha_rare2, by = "fragment_plot_comb")
  aver_alpha2 = alpha_rare_id2 %>% dplyr::select(-fragment_plot_comb) %>% dplyr::group_by(patch_type) %>% dplyr::summarise(alpha2 = mean(qD))
  alpha_C_cont_q2 = aver_alpha2 %>% dplyr::filter(patch_type=="continuous") %>% pull(alpha2)
  alpha_C_frag_q2 = aver_alpha2 %>% dplyr::filter(patch_type=="fragmented") %>% pull(alpha2)
  
  # Gamma diversity of 2 
  
  frag_pair_patches = as.data.frame(t(combn(rownames(frag_comm), 2)))
  colnames(frag_pair_patches) = c("pair1", "pair2")
  
  cont_pair_patches = as.data.frame(t(combn(rownames(cont_comm), 2)))
  colnames(cont_pair_patches) = c("pair1", "pair2")
  
  gamma_of_2_frag_targets = NULL
  
  # k = 1
  for(k in 1:nrow(frag_pair_patches)){
    paste(i, "-", "frag", print(k))
    
    pair_df = frag_pair_patches[k,]
    pair = as.character(pair_df[1,])
    frag_pair_comm_all = frag_comm[rownames(frag_comm)%in%pair,]
    frag_pair_comm = frag_pair_comm_all[,colSums(frag_pair_comm_all)>0, drop=FALSE]
    gamma_of_2_frag_targets[k] = C_target(t(colSums(frag_pair_comm)))
  }
  
  gamma_of_2_cont_targets = NULL
  # j = 1
  for(j in 1:nrow(cont_pair_patches)){
    print(paste(i, "-", "cont", print(j)))
    
    pair_df = cont_pair_patches[j,]
    pair = as.character(pair_df[1,])
    cont_pair_comm_all = cont_comm[rownames(cont_comm)%in%pair,]
    cont_pair_comm = cont_pair_comm_all[,colSums(cont_pair_comm_all)>0, drop=FALSE]
    gamma_of_2_cont_targets[j] = C_target(t(colSums(cont_pair_comm)))
  }
  
  gamma2_target = min(gamma_of_2_cont_targets, gamma_of_2_frag_targets)
  
  # w = 1
  for(w in 1:nrow(frag_pair_patches)){
    print(paste(i, "-", "frag_pair", print(w)))
    
    pair_df = frag_pair_patches[w,]
    pair = as.character(pair_df[1,])
    frag_pair_comm_all = frag_comm[rownames(frag_comm)%in%pair,]
    frag_pair_comm = frag_pair_comm_all[,colSums(frag_pair_comm_all)>0]
    PIE_2_frag = calc_PIE(colSums(frag_pair_comm))
    gamma_2_frag_q0 = estimateD(colSums(frag_pair_comm), q=0, datatype="abundance", base="coverage", level=gamma2_target) %>% pull(qD)
    gamma_2_frag_q2 = estimateD(colSums(frag_pair_comm), q=2, datatype="abundance", base="coverage", level=gamma2_target) %>% pull(qD)
    gamma_frag_df_temp = data.frame(refshort = i,
                                    patch_type = "fragmented",
                                    pair1 = pair[1],
                                    pair2 = pair[2],
                                    gamma_2_q0 = gamma_2_frag_q0,
                                    gamma_2_q2 = gamma_2_frag_q2,
                                    PIE_2 = PIE_2_frag
    ) 
    gamma_of_2_frag_values = rbind.data.frame(gamma_of_2_frag_values, gamma_frag_df_temp)
  }
  
  # z = 1
  for(z in 1:nrow(cont_pair_patches)){
    print(paste(i, "-", "cont_pair", print(z)))
    
    pair_df = cont_pair_patches[z,]
    pair = as.character(pair_df[1,])
    cont_pair_comm_all = cont_comm[rownames(cont_comm)%in%pair,]
    cont_pair_comm = cont_pair_comm_all[,colSums(cont_pair_comm_all)>0]
    PIE_2_cont = calc_PIE(colSums(cont_pair_comm)) 
    gamma_2_cont_q0 = estimateD(colSums(cont_pair_comm), q=0, datatype="abundance", base="coverage", level=gamma2_target) %>% pull(qD)
    gamma_2_cont_q2 = estimateD(colSums(cont_pair_comm), q=2, datatype="abundance", base="coverage", level=gamma2_target) %>% pull(qD)
    gamma_cont_df_temp = data.frame(refshort = i,
                                    patch_type = "continuous",
                                    pair1 = pair[1],
                                    pair2 = pair[2],
                                    gamma_2_q0 = gamma_2_cont_q0,
                                    gamma_2_q2 = gamma_2_cont_q2,
                                    PIE_2 = PIE_2_cont 
    ) 
    gamma_of_2_cont_values = rbind.data.frame(gamma_of_2_cont_values, gamma_cont_df_temp)
  }
  
  
  # Gamma diversity by landscape type (based on target_coverage) -> q = 0 
  
  all_frags = frag_comm %>% 
    dplyr::summarise(across(everything(), sum)) %>% 
    t.data.frame()
  all_conts = cont_comm %>% 
    dplyr::summarise(across(everything(), sum))%>% 
    t.data.frame()
  
  # BetaC (based on target_coverage)
  
  frag_data = study_i_comp %>% 
    dplyr::filter(patch_type == "fragmented") %>%
    dplyr::filter(fragment_plot_comb %in% new_frag_list) %>% 
    dplyr::select(-fragment_plot_comb, -patch_type)%>% 
    dplyr::select(where(~ is.numeric(.x) && sum(.x) != 0))
  
  cont_data = study_i_comp %>% filter(patch_type == "continuous") %>%
    dplyr::filter(fragment_plot_comb %in% new_frag_list) %>% 
    dplyr::select(-fragment_plot_comb, -patch_type) %>% 
    dplyr::select(where(~ is.numeric(.x) && sum(.x) != 0))
  
  beta_C_frag = as.numeric(beta_C(frag_data, target_coverage))
  beta_C_cont = as.numeric(beta_C(cont_data, target_coverage))
  
  # Combine all metrics in a temporary file and rbind the results from each loop
  
  frag_cont_comp_temp = data.frame(refshort = i,
                                   patch_type = rep(c("continuous", "fragmented"), each = 3),
                                   diversity_index = rep(c("alphaC", "alphaC2", "beta_C"), 2),
                                   diversity_value = c(alpha_C_cont_q0,
                                                       alpha_C_cont_q2,
                                                       beta_C_cont,
                                                       alpha_C_frag_q0,
                                                       alpha_C_frag_q2,
                                                       beta_C_frag))
  
  frag_cont_comp = rbind.data.frame(frag_cont_comp, frag_cont_comp_temp)
  
}



## save values from fragmented landscapes 

gamma_of_2_frag_values %>% 
  dplyr::select(refshort, patch_type, gamma_2_q0, gamma_2_q2) %>% 
  pivot_longer(!c(refshort, patch_type), names_to = "diversity_index", values_to = "diversity_value") %>% 
  group_by(refshort, patch_type, diversity_index) %>% 
  dplyr::summarise(diversity_value = mean(diversity_value)) -> gamma_of_C2_frag_summ


## save values from continuous landscapes 

gamma_of_2_cont_values %>% 
  dplyr::select(refshort, patch_type, gamma_2_q0, gamma_2_q2) %>% 
  pivot_longer(!c(refshort, patch_type), names_to = "diversity_index", values_to = "diversity_value") %>% 
  group_by(refshort, patch_type, diversity_index) %>% 
  dplyr::summarise(diversity_value = mean(diversity_value)) -> gamma_of_C2_cont_summ

## combine them in the same data.frame
gamma_C <- rbind.data.frame(gamma_of_C2_frag_summ, gamma_of_C2_cont_summ) 


## Combine gamma (q = 0, q = 2) with coverage-based alpha (q = 0, q = 2) and beta 

div_rareC_ABG <- rbind.data.frame(frag_cont_comp, gamma_C) %>% 
  ungroup() %>% 
  dplyr::arrange(refshort)

## export results to a csv file 
## note, the data is still available in the folder processed_data

# write.csv(div_rareC_ABG, "processed_data/div_rareC.csv")



### rename to match other results from individual-based rarefaction to compare alpha and gamma diversity

div_rareC_ABG %>% 
  dplyr::filter(diversity_index != "beta_C") %>% 
  mutate(
    diversity_index_reformatted = case_when(
      grepl("alphaC", diversity_index) ~ "alpha",
      grepl("gamma_2", diversity_index) ~ "gamma"
    ),
    q_order = case_when(
      grepl("C2|q2", diversity_index) ~ "q = 2",
      grepl("C$|q0", diversity_index) ~ "q = 0"
    )
  ) %>%
  dplyr::select(-diversity_index) %>%
  rename(diversity_index = diversity_index_reformatted) %>% 
  dplyr::arrange(refshort, patch_type, diversity_index) -> div_rareC_ABG_reform


### correlation between individual- and coverage-based rarefaction (q = 0)

### alpha diversity 

div_rareC_ABG_reform %>% 
  filter(q_order == "q = 0" & diversity_index == "alpha") -> alpha_q0_coverage

diversity_raref %>% 
  filter(q_order == "q = 0" & diversity_index == "alpha") -> alpha_q0_individual

cor(alpha_q0_coverage$diversity_value, alpha_q0_individual$diversity_value)
plot(alpha_q0_coverage$diversity_value, alpha_q0_individual$diversity_value)


### gamma diversity 
div_rareC_ABG_reform %>% 
  filter(q_order == "q = 0" & diversity_index == "gamma") -> gamma_q0_coverage

diversity_raref %>% 
  filter(q_order == "q = 0" & diversity_index == "gamma") -> gamma_q0_individual

cor(gamma_q0_coverage$diversity_value, gamma_q0_individual$diversity_value)
plot(gamma_q0_coverage$diversity_value, gamma_q0_individual$diversity_value)



### correlation between individual- and coverage-based rarefaction (q = 2)

### alpha diversity 

div_rareC_ABG_reform %>% 
  filter(q_order == "q = 2" & diversity_index == "alpha") -> alpha_q2_coverage

diversity_raref %>% 
  filter(q_order == "q = 2" & diversity_index == "alpha") -> alpha_q2_individual

cor(alpha_q2_coverage$diversity_value, alpha_q2_individual$diversity_value)
plot(alpha_q2_coverage$diversity_value, alpha_q2_individual$diversity_value)


### gamma diversity 
div_rareC_ABG_reform %>% 
  filter(q_order == "q = 2" & diversity_index == "gamma") -> gamma_q2_coverage

diversity_raref %>% 
  filter(q_order == "q = 2" & diversity_index == "gamma") -> gamma_q2_individual

cor(gamma_q2_coverage$diversity_value, gamma_q2_individual$diversity_value)
plot(gamma_q2_coverage$diversity_value, gamma_q2_individual$diversity_value)

