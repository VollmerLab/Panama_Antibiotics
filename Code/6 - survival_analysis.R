
#### Libraries ####
library(tidyverse)
library(survival)
library(survminer)
library(ggsurvfit)

#### Functions ####
survDiff_pairwise <- function(keep_treats, RHO = 0){
  tmp_out <- survdiff(Surv(time, infected) ~ treatment, 
                      data = filter(survival_data, 
                                    treatment %in% keep_treats),
                      rho = RHO)
  
  tibble(df = length(keep_treats) - 1, chisq = tmp_out$chisq, p.value = tmp_out$pvalue)
}

#### Data ####
survival_data <- read_csv('../Data/su17_tank_surv_data.csv', show_col_types = FALSE) %>%
  pivot_longer(cols = c(where(is.numeric), -Block),
               names_to = 'time',
               values_to = 'infected', 
               names_transform = as.numeric,
               values_transform = as.integer, 
               values_drop_na = TRUE) %>% 
  mutate(time = time / 24,
         antibiotics = if_else(str_detect(antibiotics, '^no'), 'Untreated', 'Antibiotic') %>%
           factor(levels = c('Untreated', 'Antibiotic')),
         exposure = str_to_sentence(exposure) %>% factor(levels = c('Healthy', 'Disease')),
         treatment = str_c(antibiotics, exposure, sep = '_') %>% fct_relevel('Untreated_Healthy', after = 0),
         Block = str_c(exposure, Block, sep = '_') %>% factor(),
         genotype = factor(genotype),
         id = str_c(Block, genotype, antibiotics, sep = '_'),
         dummy = 1) %>%
  
  #Add a single infection to allow model fitting - check
  # mutate(infected = case_when(genotype == 'Black' & Block == 'Healthy_2' & time == max(time) &
  #                             antibiotics == 'Antibiotic' & exposure == 'Healthy' ~ 1,
  #                             TRUE ~ infected)) %>%
  identity()


survival_data %>%
  filter(antibiotics == 'Antibiotic' & exposure == 'Healthy', infected == 1)

survival_data %>%
  count(antibiotics, exposure, time, infected) %>%
  pivot_wider(names_from = infected, values_from = n,
              values_fill = 0L) %>%
  filter(time > 5)

#### Plot IDs ####
survival_data %>%
  ggplot(aes(x = time, y = infected)) +
  geom_point() +
  stat_summary_bin(bins = 30) +
  facet_grid(antibiotics ~ exposure)

#### Binomial Model ####
filter(survival_data, time == max(time)) %>%
  group_by(antibiotics, exposure) %>%
  summarise(n_infected = sum(infected),
            total = n(),
            .groups = 'drop') %>%
  glm(cbind(n_infected, total - n_infected) ~ antibiotics * exposure, data = .,
      family = 'binomial') %>%
  car::Anova()

#### K-M Model ####
km_model <- survfit(Surv(time, infected) ~ antibiotics + exposure, 
                    data = survival_data)

summary(km_model, times = 7)
1-0.6284
survdiff(Surv(time, infected) ~ antibiotics + exposure, 
         data = survival_data, rho = 0)


survdiff(Surv(time, infected) ~ antibiotics, 
         data = survival_data, rho = 0)

survdiff(Surv(time, infected) ~ exposure, 
         data = survival_data, rho = 0)

expand_grid(treatment1 = unique(survival_data$treatment),
            treatment2 = unique(survival_data$treatment)) %>%
  filter(as.character(treatment1) < as.character(treatment2)) %>%
  rowwise(treatment1, treatment2) %>%
  reframe(survDiff_pairwise(c(treatment1, treatment2), RHO = 0))  %>%
  mutate(padj = p.adjust(p.value, 'holm')) %>%
  
  filter(treatment1 == 'Untreated_Disease' | treatment2 == 'Untreated_Disease')

sig_groups <- tribble(
  ~antibiotics, ~exposure, ~.group,
  "Untreated",   'Healthy', 'B',
  "Untreated",   'Disease', 'C',
  "Antibiotic",  'Healthy', 'A',
  "Antibiotic",  'Disease', 'B'
) %>%
  full_join(broom::tidy(km_model) %>%
              filter(time == max(time)) %>%
              select(strata, estimate) %>% 
              separate(strata, into = c('antibiotics', 'exposure'), sep = ', ') %>%
              mutate(antibiotics = str_remove(antibiotics, 'antibiotics='),
                     exposure = str_remove(exposure, 'exposure=')),
            by = c('antibiotics', 'exposure'))


broom::tidy(km_model) %>% 
  separate(strata, into = c('antibiotics', 'exposure'), sep = ', ') %>%
  mutate(antibiotics = str_remove(antibiotics, 'antibiotics='),
         exposure = str_remove(exposure, 'exposure=')) %>% 
  mutate(exposure = case_when(exposure == 'Healthy' & antibiotics == 'Antibiotic' ~ 'H_A',
                              exposure == 'Disease' & antibiotics == 'Antibiotic' ~ 'D_A',
                              TRUE ~ exposure)) %>%
  ggplot(aes(x = time, y = estimate, colour = exposure,
             ymin = estimate - std.error, ymax = estimate + std.error,
             fill = exposure, linetype = antibiotics)) +
  pammtools::geom_stepribbon(alpha = 0.25, colour = NA, show.legend = FALSE,
                             direction = 'mid') +
  pammtools::geom_stephazard(show.legend = TRUE, direction = 'mid',
                             linewidth = 0.75) +
  
  geom_text(data = sig_groups, aes(x = Inf, y = estimate, label = .group),
            inherit.aes = FALSE, hjust = 1.2) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_linetype_manual(values = c('Antibiotic' = 'dotdash', 'Untreated' = 'solid')) +
  scale_colour_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"), 
                                           "#80D1E9", "#FFA19F"),
                                         c('Healthy', 'Disease', 'H_A', 'D_A')),
                      breaks = c('Disease', 'Healthy'), 
                      labels = c('Healthy' = 'Control', 'Disease' = 'Diseased'),
                      drop = FALSE) +
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"),
                                         c('Healthy', 'Disease')),
                      breaks = c('Disease', 'Healthy'), 
                      labels = c('Healthy' = 'Control', 'Disease' = 'Diseased'),
                      drop = FALSE) +
  # guides(colour = guide_legend(override.aes = list(pch = 1))) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5)),
         linetype = guide_legend(override.aes = list(linewidth = 1.5))) +
  labs(x = 'Time Post-Exposure (D)',
       y = 'Fragment Survival (%)',
       colour = 'Exposure',
       fill = 'Exposure',
       linetype = 'Antibiotic\nTreatment') +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank(),
        axis.text = element_text(colour = 'black', size = 12),
        axis.title = element_text(colour = 'black', size = 14),
        legend.text = element_text(colour = 'black', size = 12),
        legend.title = element_text(colour = 'black', size = 14))
ggsave('../Results/Fig1_survival.png', height = 7, width = 10)

broom::tidy(km_model) %>%
  filter(time == max(time))
