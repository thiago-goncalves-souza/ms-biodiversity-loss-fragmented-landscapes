
# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: Estimate alpha diversity based on study design
# Author: Thiago Gonçalves Souza
# Date: [October 10, 2024]
# Notes: R version 4.4.1 [2024-06-14 ucrt]
# ==============================================================================


# Using renv for package management
renv::restore() # it requires installing the package renv

### Load required packages

library(tidyverse)
library(ggplot2)
library(vegan)
library(iNEXT)

### Load required data

# Species abundances in plots / fragments

landfrag_abundance <- read.csv("data/landfrag_37sset_abundance.csv")

# General information of each study, including landscape variables 
# note: this is a long format including different buffer sizes. For the paper, we used buffer == 2000

landfrag_landscape <- read.csv("data/landfrag_37sset_landscapes.csv") %>% 
  dplyr::filter(buffer == 2000)


landfrag_short <- landfrag_landscape %>% 
  dplyr::select(refshort, fragment_id, patch_type) %>% 
  unique()


### 1. Study Design: steps i and ii from the Extended Data Figure 16b ------------------------------------------------------------------


## Calculate Nstd to estimate a standardized number of species (ALPHAstd) per sampling unit
## this includes two sampling designs: standardized sample = standardized_fragment and standardized subsample = standardized_subsamples


std_frag <- landfrag_abundance %>% 
  dplyr::filter(sampling_design != "pooled") %>% 
  dplyr::select(-fragment_plot_comb, -patch_type, -sampling_effort, -sampling_design)

std_frag %>% 
  group_by(refshort, fragment_id, plot_id) %>% 
  mutate(richness = length(unique(scientific_name))) %>% 
  dplyr::select(refshort, fragment_id, plot_id, abundance, richness) %>% 
  group_by(refshort, fragment_id, plot_id) %>% 
  dplyr::summarise(total_abundance = sum(abundance),
                   total_richness = mean(richness),
                   .groups = "keep") %>% # this "mean" of richness is just a trick (it will not changing nothing)
  ungroup() %>% 
  group_by(refshort, fragment_id) %>% 
  dplyr::summarise(Nstd = round(mean(total_abundance)),
                   alpha_Sstd = round(mean(total_richness)), 
                   .groups = "keep") -> std_frag2


### Estimate Gamma Diversity
## note: estimating GAMMA diversity is problematic because of the reasons expressed in the paper and supplementary materials 

landfrag_abundance %>% 
  dplyr::filter(sampling_design != "pooled") %>% 
  dplyr::select(refshort, patch_type, scientific_name, abundance) %>% 
  group_by(refshort, patch_type) %>% 
  dplyr::summarise(gamma = length(unique(scientific_name)), 
                   .groups = "keep") -> std_sampdesign_gamma


### Estimate Beta Diversity
## note: estimating BETA diversity is problematic because of the reasons expressed in the paper and supplementary materials 

landfrag_short %>% 
  ungroup() %>% 
  inner_join(std_frag2, by = c("refshort", "fragment_id")) %>% 
  dplyr::select(-fragment_id, -Nstd) %>% 
  group_by(refshort, patch_type) %>% 
  dplyr::summarise(alpha_Sstd = mean(alpha_Sstd), .groups = "keep") %>% 
  left_join(std_sampdesign_gamma, by = c("refshort", "patch_type")) %>% 
  mutate(beta = gamma / alpha_Sstd) -> std_sampdesign_abg # alpha, beta and gamma



### 1. Study Design: step (iii) from the Extended Data Figure 16b ------------------------------------------------------------------

## Calculate Nstd and Sstd by for the pooled design using rarefaction


# species list of pooled design studies

pooled_frag_sp_list <- landfrag_abundance %>% 
  ungroup() %>% 
  dplyr::filter(sampling_design == "pooled") %>% 
  dplyr::select(refshort, fragment_id, scientific_name, abundance) %>% 
  group_by(refshort, fragment_id, scientific_name) %>% 
  summarise(abundance = sum(abundance))


# Number of individuals per fragment 

pooled_frag <- landfrag_abundance %>% 
  ungroup() %>% 
  dplyr::filter(sampling_design == "pooled") %>% 
  dplyr::select(-fragment_plot_comb, -patch_type, -sampling_design, -scientific_name) %>% 
  dplyr::mutate(sampling_effort = as.numeric(as.character(sampling_effort))) %>% 
  group_by(refshort, fragment_id) %>% 
  dplyr::summarise(total_abundance = sum(abundance),
                   sampling_effort = mean(sampling_effort),
                   .groups = "keep") %>% 
  ungroup() %>% 
  group_by(refshort) %>% 
  dplyr::mutate(rel_sampl_effort = sampling_effort / min(sampling_effort),
                Nstd = total_abundance / rel_sampl_effort) %>% 
  dplyr::select(refshort, fragment_id, Nstd)


#  studied patches 

pooled_frag_patches <-  landfrag_abundance %>% 
  ungroup() %>% 
  dplyr::filter(sampling_design == "pooled") %>% 
  dplyr::select(refshort, fragment_id, patch_type) %>% 
  unique()


# loops to estimated alpha, beta and gamma diversity with rarefaction analysis 
# note: estimating BETA and GAMMA diversity is problematic because of the reasons expressed in the paper and supplementary materials 

std_pooled <- data.frame()
gamma_pooled <- data.frame()

pooled_studies <- unique(pooled_frag_sp_list$refshort)

for(i in pooled_studies){
  print(i)
  ifrag = pooled_frag %>% filter(refshort == i) %>% 
    mutate(Nstd = round(Nstd)) # round abundance in case it is not integer
  ifrag2 = ifrag %>% 
    filter(Nstd > 2)
  icomm = pooled_frag_sp_list %>% filter(refshort == i) %>% 
    ungroup() %>% 
    dplyr::select(-refshort) %>% 
    pivot_wider(names_from = scientific_name,
                values_from = abundance, 
                values_fill = 0,
                values_fn = sum) %>% 
    column_to_rownames("fragment_id")
  min_sample = min(ifrag2$Nstd)
  
  icomm = icomm[rownames(icomm) %in% ifrag2$fragment_id,] 
  remove_1sp = specnumber(icomm) > 1
  icomm = icomm[remove_1sp,] 
  
  icomm_list = as.list(as.data.frame(t(icomm)))
  icomm_list2 = lapply(icomm_list, function(x) {sort(x [x>0], decreasing = TRUE)})
  
  # Alpha diversity by patch q = 0
  inext_rare <- estimateD(icomm_list2, q=0, datatype="abundance", base="size", level=min_sample)
  df = inext_rare %>% dplyr::select(Assemblage, qD) %>% 
    dplyr::rename(fragment_id = Assemblage,
                  Sstd = qD)
  df_temp = left_join(ifrag2, df, by = "fragment_id")
  
  std_pooled = rbind.data.frame(std_pooled, df_temp)
  
  ### gamma diversity
  
  pooled_frag_patches
  ipatch = pooled_frag_sp_list %>% 
    left_join(pooled_frag_patches, by = c("refshort", "fragment_id")) %>% 
    filter(refshort == i) %>% 
    ungroup() %>% 
    dplyr::select(patch_type, scientific_name, abundance) %>% 
    group_by(patch_type, scientific_name) %>% 
    dplyr::summarise(abundance= sum(abundance, na.rm = TRUE),
                     .groups = "keep") %>% 
    pivot_wider(names_from = scientific_name,
                values_from = abundance, 
                values_fill = 0,
                values_fn = sum) %>% 
    column_to_rownames("patch_type")
  min_sample_patches = min(rowSums(ipatch))
  
  ipatch_list = as.list(as.data.frame(t(ipatch)))
  ipatch_list2 = lapply(ipatch_list, function(x) {sort(x [x>0], decreasing = TRUE)})
  
  # Gamma diversity by patch type - q = 0
  inext_patcht_rare <- estimateD(ipatch_list2, q=0, datatype="abundance", base="size", level=min_sample_patches)
  df_patcht = inext_patcht_rare %>% dplyr::select(Assemblage, qD) %>% 
    dplyr::rename(patch_type = Assemblage,
                  gamma = qD)
  gamma_pooled_temp = data.frame(
    refshort = i,
    df_patcht
  )
  gamma_pooled = rbind.data.frame(gamma_pooled, gamma_pooled_temp)
  
}

std_pooled %>% 
  left_join(pooled_frag_patches, by = c("refshort", "fragment_id")) %>% 
  ungroup() %>% 
  dplyr::select(refshort, patch_type, Sstd) %>% 
  group_by(refshort, patch_type) %>% 
  dplyr::summarise(alpha_Sstd = mean(Sstd, na.rm = TRUE), .groups = "keep") %>% 
  left_join(gamma_pooled, by = c("refshort", "patch_type")) %>% 
  mutate(beta = gamma / alpha_Sstd) -> pooled_abg

### Last, combine the results using the three sampling designs
### std_sampdesign_abg stands for Alpha, Beta, and Gamma using the "old" way, ignoring sampling effort and distance decay
### note: at least alpha diversity is standardized by sampling effort and study design (see methods)
  
combined_abg_old <- rbind.data.frame(pooled_abg, std_sampdesign_abg) %>% 
  rename(alpha = alpha_Sstd) %>% 
  pivot_longer(cols = alpha:beta,
               names_to = "diversity_index",
               values_to = "observed_value")


### Figures -------------------------------------------------------------------------------------------------------------

## Extended Data Fig. 16. 

## first, calculate mean and SE for each landscape type (patch_type variable)
combined_abg_summary <- combined_abg_old %>% 
  dplyr::select(-refshort) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean_value = mean(observed_value, na.rm = TRUE),
                   se_value = sd(observed_value, na.rm = TRUE) / sqrt(length(unique(combined_abg_old$refshort)))) 


combined_abg_old %>% 
  ggplot() + 
  theme_bw()+
  geom_line(aes(x = patch_type, y = observed_value, group = refshort), colour = "darkgrey") + 
  geom_point(aes(x = patch_type, y = observed_value), size = 1.4, shape = 21, fill = "grey") + 
  facet_wrap(~ diversity_index, scales = "free_y", 
             labeller = label_parsed) +
  geom_pointrange(data = combined_abg_summary, 
                  aes(x = patch_type,
                      y = mean_value,
                      ymax = mean_value+se_value,
                      ymin = mean_value-se_value,
                      fill = patch_type),
                  size = 0.7,
                  lwd=1,
                  shape = 21,
                  colour = "black",
                  position = position_dodge2(width = 3))+
  scale_fill_manual(values = c("#2166ac", "#b2182b"))+
  labs(title = "Standardized alpha, but unstandardized beta and gamma diversity",x = "Landscape type", y = "Diversity value")+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(face = "bold", size = rel(1.5)),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black")) ->  Extended_Data_Fig_16

Extended_Data_Fig_16


ggsave("figures/Extended_Data_Fig_16.pdf", 
       Extended_Data_Fig_16,
       width = 22,
       height = 26,
       units = "cm",
       dpi = 300)






