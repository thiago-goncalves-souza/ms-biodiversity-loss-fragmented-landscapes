# ==============================================================================
# Data Analysis for Gonçalves-Souza et al. Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes
# Analysis: Figures 2 and 3
# Author: Thiago Gonçalves Souza
# Date: [October 10, 2024]
# Notes: R version 4.4.1 [2024-06-14 ucrt]
# ==============================================================================


# Using renv for package management
renv::restore() # it requires installing the package renv

# Load required packages ----------------------------------------------------

library(cowplot)
library(tidyverse)
library(ggplot2)
library(svglite)

# Data --------------------------------------------------------------------

# Diversity of 2 (pair = all pairs)
frag_cont_div2 <- read.csv("processed_data/diversity_of_2.csv", row.names = "X") 

# Diversity of 2 (pair = nearest fragments)

frag_cont_div2_cl <- read.csv("processed_data/diversity_of_2_close.csv", row.names = "X") 

# Diversity of 2 (pair = all pairs) - rarefied estimates (q = 0 and q = 2)

frag_cont_rar_div2 <- read.csv("processed_data/frag_cont_raref_2.csv", row.names = "X") 

# Diversity of 2 (pair = closest fragments) - rarefied estimates  (q = 0 and q = 2)

frag_cont_rar_div2_cl <- read.csv("processed_data/frag_cont_raref_close_2.csv", row.names = "X") 


# Figure 2 - models 1 (all pairs) and 2 (closest pairs)

# add alpha, beta and gamma symbols

greek_letters <- list(
  alpha = bquote(alpha),
  beta = bquote(beta),
  gamma = bquote(gamma)
)

levels(frag_cont_div2$diversity_index) <- greek_letters[levels(frag_cont_div2$diversity_index)] 


### Figure 2, Panel A

# calculate averages and SE/SD
frag_cont_div2_summary <- frag_cont_div2 %>% 
  dplyr::select(-refshort) %>% 
  dplyr::group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean_value = mean(diversity_value, na.rm = TRUE),
                   se_value = sd(diversity_value, na.rm = TRUE) / sqrt(length(unique(frag_cont_div2$refshort))),
                   sd_value = sd(diversity_value, na.rm = TRUE),
                   .groups = "keep") 

frag_cont_div2 %>% 
  dplyr::filter(diversity_index %in% c("alpha", "beta", "gamma")) %>% 
  ggplot() + 
  theme_bw()+
  geom_line(aes(x = patch_type, y = diversity_value, group = refshort), colour = "darkgrey") + 
  geom_point(aes(x = patch_type, y = diversity_value), size = 1.4, shape = 21, fill = "grey") + 
  facet_wrap(~ diversity_index, scales = "free_y", 
             labeller = label_parsed) +
  geom_pointrange(data = frag_cont_div2_summary, 
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
  labs(x = "Landscape type", y = "Diversity value")+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_text(face = "bold", size = rel(1.5)),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black")) ->  Fig2_panel_A

Fig2_panel_A


### Figure 2, Panel B

frag_cont_div2_cl_summary <- frag_cont_div2_cl %>% 
  dplyr::select(-refshort) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean_value = mean(diversity_value, na.rm = TRUE),
                   se_value = sd(diversity_value, na.rm = TRUE) / sqrt(length(unique(frag_cont_div2_cl$refshort)))) 

frag_cont_div2_cl %>% 
  dplyr::filter(diversity_index %in% c("alpha", "beta", "gamma")) %>% 
  ggplot() + 
  theme_bw()+
  geom_line(aes(x = patch_type, y = diversity_value, group = refshort), colour = "darkgrey") + 
  geom_point(aes(x = patch_type, y = diversity_value), size = 1.4, shape = 21, fill = "grey") + 
  facet_wrap(~ diversity_index, scales = "free_y", labeller = label_parsed) +
  geom_pointrange(data = frag_cont_div2_cl_summary, 
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
  labs(x = "Landscape type", y = "Diversity value")+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_blank()) ->  Fig2_panel_B


## combined panels 

combined_plot = plot_grid(Fig2_panel_A, Fig2_panel_B, 
                          ncol = 1, 
                          labels = c('a', 'b'),
                          align = "vh")


combined_plot

## use ggsave if you want to export it 

# suggested code
ggsave("figures/Figure2AB.pdf", 
       combined_plot,
       width = 22,
       height = 26,
       units = "cm",
       dpi = 300)


## Figure 3, Panel A (all pairs, q = 0)

frag_cont_rar_div2_q0_summary <- frag_cont_rar_div2 %>% 
  dplyr::filter(q_order == "q = 0") %>% 
  dplyr::select(-refshort, -q_order) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean_value = mean(diversity_value, na.rm = TRUE),
                   se_value = sd(diversity_value, na.rm = TRUE) / sqrt(length(unique(frag_cont_rar_div2$refshort))),
                   .groups = "keep") 


frag_cont_rar_div2 %>%
  filter(q_order == "q = 0") %>% 
  ggplot() + 
  theme_bw()+
  geom_line(aes(x = patch_type, y = diversity_value, group = refshort), colour = "darkgrey") + 
  geom_point(aes(x = patch_type, y = diversity_value), size = 1.4, shape = 21, fill = "grey") + 
  facet_wrap(~ diversity_index, scales = "free_y", labeller = label_parsed) +
  geom_pointrange(data = frag_cont_rar_div2_q0_summary, 
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
  labs(x = "Landscape type", y = "Diversity value (rarefied)")+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.text =  element_text(size = 16)) -> Fig_div2_raref_q0

Fig_div2_raref_q0

## Figure 3, Panel B (all pairs, q = 2)

frag_cont_rar_div2_q2_summary <- frag_cont_rar_div2 %>% 
  filter(q_order == "q = 2") %>% 
  dplyr::select(-refshort, -q_order) %>% 
  group_by(patch_type, diversity_index) %>% 
  dplyr::summarise(mean_value = mean(diversity_value, na.rm = TRUE),
                   se_value = sd(diversity_value, na.rm = TRUE) / sqrt(length(unique(frag_cont_rar_div2$refshort))),
                   .groups = "keep") 


frag_cont_rar_div2 %>%
  filter(q_order == "q = 2") %>% 
  ggplot() + 
  theme_bw()+
  geom_line(aes(x = patch_type, y = diversity_value, group = refshort), colour = "darkgrey") + 
  geom_point(aes(x = patch_type, y = diversity_value), size = 1.4, shape = 21, fill = "grey") + 
  facet_wrap(~ diversity_index, scales = "free_y", labeller = label_parsed) +
  geom_pointrange(data = frag_cont_rar_div2_q2_summary, 
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
  labs(x = "Landscape type", y = "Diversity value (rarefied)")+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_blank()) -> Fig_div2_raref_q2

Fig_div2_raref_q2

## combined panels 

combined_plot_raref = plot_grid(Fig_div2_raref_q0, Fig_div2_raref_q2, 
                                ncol = 1, 
                                labels = c('A) q = 0', 'B) q = 2'),
                                align = "vh")


## use ggsave if you want to export it 
# suggested code


ggsave("figures/Figure3AB.pdf", 
       combined_plot_raref,
       width = 22,
       height = 26,
       units = "cm",
       dpi = 300)

