#### Libraries ####
library(multcomp)
library(tidyverse)
library(microbiome)
library(phyloseq)
library(biomeUtils)
library(glmmTMB)
library(emmeans)
library(patchwork)

source('~/R/R Functions/diagPlots.R')

#### Functions ####
process_outemmeans <- function(emmeans_in){
  broom::tidy(emmeans_in, conf.int = TRUE) %>%
    separate(time_treat, into = c('time', 'anti', 'health')) %>%
    mutate(time = factor(time, levels = c('before', 'after')),
           anti = factor(anti, levels = c('N', 'A')),
           health = factor(health, levels = c('H', 'D')))
}


#### Data ####
tank_data <- read_rds("../intermediate_files/full_tank_microbiome.rds") %>%
  rarefy_even_depth(verbose = TRUE, rngseed = 5222024) %>%
  identity()

summarize_phyloseq(tank_data)

metadata <- sample_data(tank_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac) %>%
  mutate(treatment = str_c(anti, health, sep = '_') %>% factor,
         time_treat = str_c(time, anti, health, sep = '_') %>% factor)


#### Calculate alpha metrics ####
alpha_metrics <- alpha(tank_data, zeroes = TRUE, 
      index = c('Observed', 'shannon', 'camargo', 'gini')) %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(calculatePD(tank_data, justDF = TRUE) %>%
              as_tibble(rownames = 'sample_id') %>%
              select(sample_id, PD),
            by = 'sample_id') %>%
  rename(richness = observed,
         diversity = diversity_shannon,
         evenness = evenness_camargo,
         dominance = dominance_gini,
         phylogenetic = PD) %>%
  left_join(metadata, by = 'sample_id')

#
#### Basic Plot ####
alpha_metrics %>% 
  pivot_longer(cols = c(richness:phylogenetic),
               names_to = 'metric') %>%
  
  ggplot(aes(x = time, y = value, fill = health, shape = anti)) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5)) +
  
  facet_wrap(~metric, scales = 'free_y') +
  
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 1, fill = 'black'))) +
  labs(x = NULL, 
       y = 'Chao1',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment')


#### Richness ####
ggplot(alpha_metrics, aes(x = richness)) +
  geom_histogram(bins = 30)

richness_model_pois <- glmmTMB(richness ~ time_treat + 
                            (1 | tank) + 
                            (1 | geno / fragment), 
                          family = 'poisson',
                          data = alpha_metrics)


car::Anova(richness_model_pois)
emmeans(richness_model_pois, ~time_treat, type = 'response') %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))

ref_grid(richness_model_pois) %>%
  add_grouping('anti', 'time_treat', c('A', NA, 'N', 'A', 'N')) %>%
  emmeans(~anti, type = 'response') %>%
  contrast('revpairwise')

ref_grid(richness_model_pois) %>%
  add_grouping('disease', 'time_treat', c(NA, 'D', 'H', NA, 'H')) %>%
  emmeans(~disease, type = 'response') %>%
  contrast('revpairwise')

richness_plot <- emmeans(richness_model_pois, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  rename(response = rate) %>%
  ggplot(aes(x = time, y = response, 
             ymin = response - std.error, 
             ymax = response + std.error,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  # facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
  labs(x = NULL, 
       y = 'Richness',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())

#### Evenness ####
ggplot(alpha_metrics, aes(x = evenness)) +
  geom_histogram(bins = 30)

evenness_model <- glmmTMB(evenness ~ time_treat + 
          (1 | tank) + 
          (1 | geno / fragment), 
        family = 'beta_family',
        data = alpha_metrics)


car::Anova(evenness_model)
emmeans(evenness_model, ~time_treat) %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))



evenness_plot <- emmeans(evenness_model, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  ggplot(aes(x = time, y = response, 
             ymin = response - std.error, 
             ymax = response + std.error,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  # facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
  labs(x = NULL, 
       y = 'Evenness',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())

#### Diversity ####
ggplot(alpha_metrics, aes(x = diversity)) +
  geom_histogram(bins = 30)

diversity_model <- glmmTMB(diversity ~ time_treat + 
          (1 | tank) + 
          (1 | geno / fragment), 
        family = Gamma(link = 'log'),
        data = alpha_metrics)

car::Anova(diversity_model)
emmeans(diversity_model, ~time_treat) %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))

ref_grid(diversity_model) %>%
  add_grouping('anti', 'time_treat', c('A', NA, 'N', 'A', 'N')) %>%
  emmeans(~anti, type = 'response') %>%
  contrast('revpairwise')


diversity_plot <- emmeans(diversity_model, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  ggplot(aes(x = time, y = response, 
             ymin = response - std.error, 
             ymax = response + std.error,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  # facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
  labs(x = NULL, 
       y = 'Diversity',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())


#### Dominance ####
ggplot(alpha_metrics, aes(x = dominance)) +
  geom_histogram(bins = 30)

dominance_model <- glmmTMB(dominance ~ time_treat + 
          (1 | tank) + 
          (1 | geno), 
        family = 'beta_family',
        data = alpha_metrics)

diag.plots(dominance_model, col.nos = c('tank', 'anti', 'geno', 'health', 'time'), data = alpha_metrics)

# dominance_model2 <- glmmTMB(dominance ~ time_treat + 
#                              (1 | tank) + 
#                              (1 | geno/fragment), 
#                            dispformula = ~time,
#                            family = 'beta_family',
#                            data = alpha_metrics)
# 
# MuMIn::AICc(dominance_model, dominance_model2)
# AIC(dominance_model, dominance_model2)
# BIC(dominance_model, dominance_model2)


car::Anova(dominance_model)
emmeans(dominance_model, ~time_treat) %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))

ref_grid(dominance_model) %>%
  add_grouping('anti', 'time_treat', c('A', NA, 'N', 'A', 'N')) %>%
  emmeans(~anti, type = 'response') %>%
  contrast('pairwise')


dominance_plot <- emmeans(dominance_model, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  ggplot(aes(x = time, y = response, 
             ymin = response - std.error, 
             ymax = response + std.error,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  # facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
  labs(x = NULL, 
       y = 'Dominance',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())
dominance_plot


#### Phylogenetic Diveristy ####
ggplot(alpha_metrics, aes(x = phylogenetic)) +
  geom_histogram(bins = 30)

phylo_model_norm <- glmmTMB(phylogenetic ~ time_treat + 
          (1 | tank) + 
          (1 | geno / fragment), 
        family = 'gaussian',
        data = alpha_metrics)

phylo_model_gamma <- glmmTMB(phylogenetic ~ time_treat + 
                              (1 | tank) + 
                              (1 | geno / fragment), 
                            family = Gamma(link = 'log'),
                            data = alpha_metrics)

MuMIn::AICc(phylo_model_norm, phylo_model_gamma)
AIC(phylo_model_norm, phylo_model_gamma)
BIC(phylo_model_norm, phylo_model_gamma)

diag.plots(phylo_model_norm, col.nos = c('tank', 'anti', 'geno', 'health', 'time'), data = alpha_metrics)


car::Anova(phylo_model_norm)
emmeans(phylo_model_norm, ~time_treat) %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))

phylo_plot <- emmeans(phylo_model_norm, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  rename(response = estimate) %>%
  ggplot(aes(x = time, y = response, 
             ymin = response - std.error, 
             ymax = response + std.error,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = 'black'))) +
  # facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
  labs(x = NULL, 
       y = 'Phylogenetic Diversity',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())

#### Join Together ####
list(phylo_plot, evenness_plot,
     richness_plot, diversity_plot, 
     dominance_plot, guide_area()) %>%
  wrap_plots(ncol = 2) +
  plot_layout(guides = 'collect', axes = 'collect',
              axis_titles = 'collect') +
  plot_annotation(tag_levels = 'A')
ggsave('../Results/Fig3_alpha_diversity.png', height = 8, width = 6)
