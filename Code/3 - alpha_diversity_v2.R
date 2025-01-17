#### Libraries ####
library(multcomp)
library(tidyverse)
library(microbiome)
library(phyloseq)
library(biomeUtils)
library(lmerTest)
library(emmeans)
library(patchwork)

source('~/R/R Functions/diagPlots.R')

#
#### Functions ####
process_outemmeans <- function(emmeans_in){
  broom::tidy(emmeans_in, conf.int = TRUE) %>%
    separate(time_treat, into = c('time', 'anti', 'health')) %>%
    mutate(time = factor(time, levels = c('before', 'after')),
           anti = factor(anti, levels = c('N', 'A')),
           health = factor(health, levels = c('H', 'D')))
}

rarify_sample <- function(otu_mat, sample_depth){
  out <- vegan::rrarefy(otu_mat, sample = sample_depth)
  out[,colSums(out) > 0]
}

rarify_alpha <- function(phylo, sample_depth, ...){
  rare_data <- rarefy_even_depth(phylo, sample.size = sample_depth,
                    trimOTUs = TRUE, verbose = FALSE)
  
  normal_alpha <- alpha(rare_data, ...) %>%
    as_tibble(rownames = 'sample_id')
    
  faiths_div <- calculatePD(rare_data, justDF = TRUE) %>%
    as_tibble(rownames = 'sample_id') %>%
    select(sample_id, PD)
  
  full_join(normal_alpha, faiths_div, 
            by = 'sample_id')
}

#
#### Data ####
raw_tank_data <- read_rds("../intermediate_files/full_tank_microbiome.rds")

summarize_phyloseq(raw_tank_data)

minimum_reads <- 1056

metadata <- sample_data(raw_tank_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac) %>%
  mutate(treatment = str_c(anti, health, sep = '_') %>% factor,
         time_treat = str_c(time, anti, health, sep = '_') %>% factor)

#### Rarefaction Curves & Good's Coverage ####
rare_data <- otu_table(raw_tank_data) %>%
  as.data.frame %>%
  vegan::rarecurve(step = 10, tidy = TRUE) %>%
  rename(sample_id = Site,
         n_reads = Sample,
         n_asv = Species) %>%
  
  left_join(metadata, 
            by = 'sample_id') 

# minimum_reads <- 3000

rare_data %>%
  mutate(treatment = case_when(anti == 'N' ~ health,
                               TRUE ~ str_c('H_A'))) %>%
  
  ggplot(aes(x = n_reads, y = n_asv, group = sample_id, 
             colour = treatment, linetype = time, fill = anti)) +
  
  geom_line() +
  geom_vline(xintercept = minimum_reads) +
  
  scale_x_continuous(labels = scales::comma_format()) +
  scale_colour_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                    type = "continuous"),
                                           "#80D1E9"),
                                         c('H', 'D', 'H_A')),
                      breaks = c('D', 'H_A', 'H'), 
                      labels = c('H' = 'Healthy', 'H_A' = 'Antibiotic', 'D' = 'Diseased'),
                      drop = FALSE) +
  scale_linetype_manual(values = c('before' = 'dashed', 'after' = 'solid'),
                      breaks = c('before', 'after'), 
                      labels = c('before' = 'Before', 'after' = 'After'),
                      drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(line.width = 5))) +
  labs(x = 'Number of Reads',
       y = 'ASV Richness',
       colour = 'Health\nState \n&\n Antibiotic\nTreatment',
       linetype = 'Dose\nTiming') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())
ggsave('../Results/rarefaction_curves.png', height = 7, width = 7)
ggsave('../Results/Figure S2.tiff', height = 7, width = 7)

goods_coverage <- otu_table(raw_tank_data) %>%
  as.data.frame %>%
  as_tibble(rownames = 'sample_id') %>%
  pivot_longer(cols = where(is.numeric),
               names_to = 'asv_id', 
               values_to = 'n_reads') %>%
  summarise(n_seqs = sum(n_reads),
            goods_coverage = 100 * (1 - sum(n_reads == 1) / sum(n_reads)),
            .by = sample_id) %>%
  left_join(metadata, 
            by = 'sample_id')
# 20/99
filter(goods_coverage, n_seqs > minimum_reads)
# minimum_reads <- 5118
goods_coverage %>%
  arrange(goods_coverage) %>%
  mutate(treatment = case_when(anti == 'N' ~ health,
                               TRUE ~ str_c('H_A'))) %>% 
  ggplot(aes(x = n_seqs, y = goods_coverage, 
             colour = health, shape = anti, 
             fill = interaction(time == 'before', health))) +
  geom_point(size = 2, stroke = 1) +
  # geom_vline(xintercept = minimum_reads) +
  
  scale_x_continuous(labels = scales::comma_format(), limits = c(0, NA)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  
  scale_colour_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                    type = "continuous"),
                                           "#80D1E9"),
                                         c('H', 'D', 'H_A')),
                      breaks = c('D', 'H'), 
                      labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                      drop = FALSE) +
  scale_fill_manual(values = set_names(c("#80D1E9", wesanderson::wes_palette("Zissou1", 2, 
                                                                             type = "continuous"),
                                         'white', 'white'),
                                       c('FALSE.H_A', 'FALSE.H', 'FALSE.D', 'TRUE.H_A', 'TRUE.H')),
                    breaks = c('TRUE.H', 'FALSE.H'),
                    labels = c('TRUE.H' = 'Before', 'FALSE.H' = "After")) +
  
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 'circle'), stroke = 0.1),
         fill = guide_legend(override.aes = list(size = 4, shape = c('Before' = 'circle open', 
                                                                     'After' = 'circle'), 
                                                 linewidth = 2, colour = 'black',
                                                 stroke = 0.1)),
         shape = guide_legend(override.aes = list(size = 4, fill = c("#3A9AB2", "#80D1E9"), 
                                                  stroke = 0.1))) +
  labs(x = 'Number of Reads', 
       y = "Good's Coverage",
       colour = 'Health\nState',
       fill = 'Dose\nTiming',
       shape = 'Antibiotic\nTreatment') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())
ggsave('../Results/Figure S1.tiff', height = 7, width = 7)

#### Calculate alpha metrics ####
if(file.exists(str_c('../intermediate_files/alpha_metrics_', minimum_reads, '.csv'))){
  alpha_metrics <- read_csv(str_c('../intermediate_files/alpha_metrics_', minimum_reads, '.csv'),
                            show_col_types = FALSE)
} else {
  library(furrr)
  plan(multisession)
  alpha_metrics <- future_map_dfr(1:1000, 
                           ~rarify_alpha(raw_tank_data, minimum_reads,
                                         zeroes = TRUE, 
                                         index = c('Observed', 'shannon', 
                                                   'camargo', 'gini')),
                           seed = TRUE) %>%
    summarise(across(where(is.numeric), mean),
              .by = sample_id) %>%
    
    rename(richness = observed,
           diversity = diversity_shannon,
           evenness = evenness_camargo,
           dominance = dominance_gini,
           phylogenetic = PD) %>%
    left_join(metadata, by = 'sample_id')
  write_csv(alpha_metrics, str_c('../intermediate_files/alpha_metrics_', minimum_reads, '.csv'))
}


#
#### Basic Plot ####
alpha_metrics %>% 
  pivot_longer(cols = c(richness:dominance),
               names_to = 'metric') %>%
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health)) %>%
  
  ggplot(aes(x = time, y = value, fill = health, shape = anti)) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5)) +
  
  facet_wrap(~metric, scales = 'free_y') +
  
  scale_fill_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"), "#80D1E9"),
                                       c('H', 'D', 'H_A')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 1, fill = c("#3A9AB2", "#80D1E9")))) +
  labs(x = NULL, 
       y = 'Chao1',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment')


#### Richness ####
ggplot(alpha_metrics, aes(x = richness)) +
  geom_histogram(bins = 30)

richness_model <- lmer(richness ~ time_treat + 
                            (1 | tank) + 
                            (1 | geno / fragment), 
                          data = alpha_metrics)

diag.plots(richness_model, data = alpha_metrics, col.nos = c('time_treat'))

anova(richness_model, ddf = 'S')
emmeans(richness_model, ~time_treat, type = 'response') %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))

ref_grid(richness_model) %>%
  add_grouping('anti', 'time_treat', c('A', NA, 'N', 'A', 'N')) %>%
  emmeans(~anti, type = 'response') %>%
  contrast('revpairwise')

ref_grid(richness_model) %>%
  add_grouping('disease', 'time_treat', c(NA, 'D', 'H', NA, 'H')) %>%
  emmeans(~disease, type = 'response') %>%
  contrast('revpairwise')

richness_plot <- emmeans(richness_model, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) %>%
  # rename(response = rate) %>%
  ggplot(aes(x = time, y = estimate, 
             # ymin = response - std.error, 
             # ymax = response + std.error,
             ymin = conf.low, ymax = conf.high,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                         "#80D1E9"),
                                       c('H', 'D', 'H_A')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = c("#3A9AB2", "#80D1E9")))) +
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

evenness_model <- lmer(evenness ~ time_treat + 
          (1 | tank) + 
          (1 | geno / fragment), 
        data = alpha_metrics)

diag.plots(evenness_model, data = alpha_metrics, col.nos = c('time_treat'))

anova(evenness_model, ddf = 'S')
emmeans(evenness_model, ~time_treat) %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))



evenness_plot <- emmeans(evenness_model, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) %>%
  ggplot(aes(x = time, y = estimate, 
             # ymin = response - std.error, 
             # ymax = response + std.error,
             ymin = conf.low, ymax = conf.high,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"),
                                         "#80D1E9"),
                                       c('H', 'D', 'H_A')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = c("#3A9AB2", "#80D1E9")))) +
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

diversity_model <- lmer(diversity ~ time_treat + 
          (1 | tank) + 
          (1 | geno / fragment), 
        data = alpha_metrics)

diag.plots(diversity_model, data = alpha_metrics, col.nos = c('time_treat'))
anova(diversity_model, ddf = "S")
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
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) %>%
  ggplot(aes(x = time, y = estimate, 
             # ymin = response - std.error, 
             # ymax = response + std.error,
             ymin = conf.low, ymax = conf.high,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"),
                                         "#80D1E9"),
                                       c('H', 'D', 'H_A')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = c("#3A9AB2", "#80D1E9")))) +
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

dominance_model <- lmer(dominance ~ time_treat + 
          (1 | tank) + 
          (1 | geno/fragment), 
        data = alpha_metrics)

diag.plots(dominance_model, data = alpha_metrics, col.nos = c('time_treat'))

anova(dominance_model, ddf = "S")
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
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) %>%
  ggplot(aes(x = time, y = estimate, 
             # ymin = response - std.error, 
             # ymax = response + std.error,
             ymin = conf.low, ymax = conf.high,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"),
                                         "#80D1E9"),
                                       c('H', 'D', 'H_A')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = c("#3A9AB2", "#80D1E9")))) +
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

phylo_model <- lmer(phylogenetic ~ time_treat + 
          (1 | tank) + 
          (1 | geno / fragment), 
        data = alpha_metrics)

diag.plots(phylo_model, col.nos = c('time_treat'), data = alpha_metrics)


anova(phylo_model, ddf = 'S')
emmeans(phylo_model, ~time_treat) %>%
  contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                         'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                         'time' = c(1/2, 0, 1/2, -1/2, -1/2)))

phylo_plot <- emmeans(phylo_model, ~time_treat, type = 'response') %>%
  process_outemmeans %>%
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) %>%
  # rename(response = estimate) %>%
  ggplot(aes(x = time, y = estimate, 
             # ymin = response - std.error, 
             # ymax = response + std.error,
             ymin = conf.low, ymax = conf.high,
             fill = health, shape = anti)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), 
             size = 2) +
  
  scale_fill_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"),
                                         "#80D1E9"),
                                       c('H', 'D', 'H_A')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_x_discrete(labels = ~str_to_sentence(.) %>% str_c(., '\nDisease Dose')) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         shape = guide_legend(override.aes = list(size = 4, fill = c("#3A9AB2", "#80D1E9")))) +
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
  plot_annotation(tag_levels = 'A') &
  theme(legend.box = "horizontal")
ggsave('../Results/Fig4_alpha_diversity.png', height = 8, width = 6)
ggsave('../Results/Fig4_r5.tiff', height = 8, width = 6, dpi = 'print')

#### Significance Table ####
# model <- phylo_model
model_2_table <- function(model){
  main_effect <- anova(model, ddf = 'S') %>%
    as_tibble() %>%
    select(-`Sum Sq`, -`Mean Sq`) %>%
    rename(n_df = NumDF,
           d_df = DenDF,
           f_value = `F value`,
           p_value = `Pr(>F)`)
  
  individual_effects <- emmeans(model, ~time_treat) %>%
    contrast(method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                           'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                           'time' = c(1/2, 0, 1/2, -1/2, -1/2))) %>%
    as_tibble %>%
    select(contrast, estimate, p.value) %>%
    mutate(direction = if_else(estimate < 0, '-', '+'),
           .keep = 'unused') %>%
    pivot_wider(names_from = contrast,
                values_from = c(direction, p.value), 
                names_vary = 'slowest')
  
  bind_cols(main_effect, individual_effects)
}

tribble(
  ~'metric', ~'model',
  'richness', richness_model,
  'evenness', evenness_model,
  'diversity', diversity_model,
  'dominance', dominance_model,
  'phylo', phylo_model
) %>%
  rowwise(metric) %>%
  reframe(model_2_table(model)) %>%
  write_csv('../Results/alpha_diversity_metrics.csv')
  
