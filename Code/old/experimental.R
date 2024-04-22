#### Libraries ####
library(tidyverse)
library(glmmTMB)


#### Data ####
tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
                      show_col_types = FALSE) %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time2 = if_else(time == 0, 'before', 'after'),
         time2 = factor(time2, levels = c('before', 'after')))

#### Mega Model ####
subset_asvs <- nest(tank_data, data = -asv_id) %>%
  sample_n(100) %>%
  unnest(data)

nt <- min(parallel::detectCores(), 5)
complete_model <- glmmTMB(log2_cpm_norm ~ time * anti * exposure + 
                            (1 + time * exposure * anti | asv_id) + 
                            (1 | tank) + (1 | geno / fragment),
                          dispformula = ~asv_id,
                          data = tank_data,
                          control = glmmTMBControl(parallel = nt),
                          verbose = TRUE)
summary(complete_model)

library(emmeans)
emmeans(complete_model, ~time * anti * exposure,
        at = list(time = c(0, 2, 8))) %>%
  broom::tidy() %>%
  filter(!(time == 0 & exposure == 'N')) %>%
  ggplot(aes(x = time, y = estimate, 
             ymin = estimate - std.error,
             ymax = estimate + std.error, 
             colour = exposure,
             shape = anti)) +
  geom_errorbar(position = position_dodge(0.5), width = 0.1) +
  geom_point(size = 2, position = position_dodge(0.5))


coef(complete_model)$cond$asv_id %>%
  as_tibble(rownames = 'asv_id') %>%
  
  left_join(select(tank_data, asv_id, family, genus, species) %>%
              distinct,
            by = 'asv_id') %>%
  
  pivot_longer(cols = where(is.numeric),
               names_to = 'term',
               values_to = 'coef') %>%
  filter(asv_id == 'ASV97')


filter(term == 'time:antiA') %>%
  arrange(-abs(coef))

tank_data %>%
  filter(asv_id == 'ASV97') %>%
  ggplot(aes(x = time, y = log2_cpm_norm, colour = exposure, shape = anti)) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5))


coef(complete_model)$cond$asv_id %>%
  scale %>%
  dist %>%
  hclust(method = 'average') %>%
  plot()


the_pca <- coef(complete_model)$cond$asv_id %>%
  prcomp()

the_pca$x %>%
  as_tibble(rownames = 'asv_id') %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_segment(data = the_pca$rotation,
               aes(xend = 0, yend = 0, x = PC1 * 10, y = PC2 * 10), colour = 'red') +
  geom_text(data = as_tibble(the_pca$rotation, rownames = 'term'),
            aes(label = term, x = PC1 * 10, y = PC2 * 10), colour = 'red') 


filter(tank_data, exposure == 'D') %>%
  group_by(asv_id) %>%
  filter(var(log2_cpm_norm) < 0.1)

complete_model_disease <- glmmTMB(log2_cpm_norm ~ time * (anti + health) + 
                                    (1 + time * (anti + health) | asv_id) + 
                                    (1 | tank) + (1 | geno / fragment),
                                  dispformula = ~asv_id,
                                  data = filter(tank_data, exposure == 'D'),
                                  control = glmmTMBControl(parallel = nt),
                                  verbose = TRUE)
summary(complete_model_disease)



the_pca <- coef(complete_model_disease)$cond$asv_id %>%
  prcomp()

the_pca$x %>%
  as_tibble(rownames = 'asv_id') %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_segment(data = the_pca$rotation,
               aes(xend = 0, yend = 0, x = PC1 * 10, y = PC2 * 10), colour = 'red') +
  geom_text(data = as_tibble(the_pca$rotation, rownames = 'term'),
            aes(label = term, x = PC1 * 10, y = PC2 * 10), colour = 'red') 



coef(complete_model_disease)$cond$asv_id %>%
  as_tibble(rownames = 'asv_id') %>%
  
  left_join(select(tank_data, asv_id, family, genus, species) %>%
              distinct,
            by = 'asv_id') %>%
  
  pivot_longer(cols = where(is.numeric),
               names_to = 'term',
               values_to = 'coef') %>%
  filter(term == 'healthH') %>%
  arrange(-abs(coef)) %>%
  filter(coef < 0)

coef(complete_model_disease)$cond$asv_id %>%
  as_tibble(rownames = 'asv_id') %>%
  
  left_join(select(tank_data, asv_id, family, genus, species) %>%
              distinct,
            by = 'asv_id') %>%
  
  pivot_longer(cols = where(is.numeric),
               names_to = 'term',
               values_to = 'coef') %>%
  filter(asv_id %in% c('ASV25', 'ASV8'))


tank_data %>%
  filter(asv_id == 'ASV25',
         exposure == 'D') %>%
  mutate(time2 = if_else(time == 0, 'before', 'after'),
         time2 = factor(time2, levels = c('before', 'after'))) %>%
  ggplot(aes(x = time2, y = log2_cpm_norm, colour = health, shape = anti)) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5))

filter(tank_data, asv_id == 'ASV1') %>%
  count(anti, health, time, exposure) %>%
  pivot_wider(names_from = health,
              values_from = n,
              values_fill = 0)

library(lmerTest)
library(emmeans)
tmp <- tank_data %>%
  filter(exposure == 'D') %>%
  nest_by(asv_id) %>%
  filter(asv_id %in% c('ASV25', 'ASV8')) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ time2 * anti + health + 
                             (1 | tank) + (1 | geno / fragment),
                           data = data)),
         antibiotic_effect = list(emmeans(model, ~anti) %>%
                                    contrast('revpairwise')),
         outcome_effect = list(emmeans(model, ~health, at = list(time = 'after')) %>%
                                 contrast('revpairwise')),
         
         plot_data = list(emmeans(model, ~time2 * anti + health) %>% broom::tidy(conf.int = TRUE)))

tmp$antibiotic_effect
tmp$outcome_effect

tmp %>%
  select(asv_id, plot_data) %>%
  unnest(plot_data) %>%
  mutate(time2 = factor(time2, levels = c('before', 'after'))) %>%
  filter(!(time2 == 'before' & health == 'D'), !(health == 'D' & anti == 'A')) %>%
  ggplot(aes(x = time2, y = estimate, 
             ymin = conf.low, ymax = conf.high,
             colour = health, shape = anti)) +
  geom_errorbar(position = position_dodge(0.25), width = 0.1) +
  geom_point(position = position_dodge(0.25)) +
  facet_wrap(~asv_id)


tmp <- tank_data %>%
  filter(exposure == 'D') %>%
  nest_by(asv_id) %>%
  filter(asv_id %in% c('ASV25', 'ASV8')) %>%
  mutate(model = list(lmer(log2_cpm_norm ~ time2 + anti + health +
                             (1 | tank) + (1 | geno / fragment),
                           data = data)),
         antibiotic_effect = list(emmeans(model, ~anti, lmer.df = 'satterthwaite') %>%
                                    contrast('revpairwise')),
         outcome_effect = list(emmeans(model, ~health, lmer.df = 'satterthwaite') %>%
                                 contrast('revpairwise')),
         
         plot_data = list(emmeans(model, ~time2 + anti + health) %>% broom::tidy(conf.int = TRUE)))

anova(tmp$model[[2]], ddf = 'Satterthwaite')
tmp$antibiotic_effect[[2]]
tmp$outcome_effect[[2]]

tmp %>%
  select(asv_id, plot_data) %>%
  unnest(plot_data) %>%
  mutate(time2 = factor(time2, levels = c('before', 'after'))) %>%
  filter(!(time2 == 'before' & health == 'D'), !(health == 'D' & anti == 'A')) %>%
  ggplot(aes(x = time2, y = estimate, 
             ymin = conf.low, ymax = conf.high,
             colour = health, shape = anti)) +
  geom_errorbar(position = position_dodge(0.25), width = 0.1) +
  geom_point(position = position_dodge(0.25)) +
  facet_wrap(~asv_id)




#### ASV Model ####
filter(tank_data, asv_id == 'ASV1') %>%
  count(time, anti, exposure, health) %>%
  pivot_wider(names_from = health, values_from = n,
              values_fill = 0)
group_by(anti, exposure)

asv_models <- tank_data %>%
  nest_by(asv_id) %>%
  filter(asv_id %in% c('ASV25', 'ASV8')) %>%
  
  # ungroup %>%
  # sample_n(10) %>%
  # rowwise %>%
  
  mutate(model = list(glmmTMB(log2_cpm_norm ~ time * anti * exposure + 
                                (1 | tank) + (1 | geno / fragment),
                              data = data)),
         
         summary(model)$coef$cond %>%
           as_tibble(rownames = 'term') %>%
           filter(term != '(Intercept)') %>%
           mutate(term = str_remove_all(term, 'N') %>% str_replace_all(':', 'X')) %>%
           select(term, `z value`, `Pr(>|z|)`) %>%
           rename(z = `z value`,
                  p = `Pr(>|z|)`) %>%
           pivot_wider(names_from = term, 
                       values_from = c(z, p))) %>%
  ungroup


#### Upset ####
asv_models %>%
  select(starts_with('p'))


library(emmeans)
asv_models %>%
  rowwise(asv_id) %>%
  reframe(emmeans(model, ~time * exposure * anti, data = data,
                  at = list(time = c(0, 2, 8))) %>%
            broom::tidy()) %>%
  filter(!(exposure == 'N' & time == 0)) %>%
  ggplot(aes(x = time, y = estimate, 
             ymin = estimate - std.error,
             ymax = estimate + std.error,
             colour = exposure,
             shape = anti)) +
  geom_errorbar(position = position_dodge(0.5), width = 0.1) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~asv_id)




asv_models %>%
  rowwise(asv_id) %>%
  reframe(emmeans(model, ~anti | exposure * time, data = data,
                  at = list(time = c(0, 2, 8))) %>%
            contrast('revpairwise') %>%
            broom::tidy()) %>%
  filter(!(time == 0 & exposure == 'N')) %>%
  ggplot(aes(x = time, y = estimate, 
             ymin = estimate - std.error,
             ymax = estimate + std.error,
             colour = exposure)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_errorbar(position = position_dodge(0.5), width = 0.1) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~asv_id)



asv_models %>%
  rowwise(asv_id) %>%
  reframe(emmeans(model, ~exposure | anti * time, data = data,
                  at = list(time = c(2, 8))) %>%
            contrast('pairwise') %>%
            broom::tidy()) %>%
  ggplot(aes(x = time, y = estimate, 
             ymin = estimate - std.error,
             ymax = estimate + std.error,
             colour = anti)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_errorbar(position = position_dodge(0.5), width = 0.1) +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~asv_id)


#### Just disease exposed ####
asv_models <- tank_data %>%
  filter(asv_id %in% c('ASV25', 'ASV8')) %>%
  mutate(time_anti_expose_health = str_c(time, anti, exposure, health, sep = '_')) %>% 
  nest_by(asv_id) %>%
  
  # ungroup %>%
  # sample_n(10) %>%
  # rowwise %>%
  
  mutate(model = list(glmmTMB(log2_cpm_norm ~ time_anti_expose_health + 
                                (1 | tank) + (1 | geno / fragment),
                              data = data)),
         
         summary(model)$coef$cond %>%
           as_tibble(rownames = 'term') %>%
           filter(term != '(Intercept)') %>%
           mutate(term = str_remove_all(term, 'N') %>% str_replace_all(':', 'X')) %>%
           select(term, `z value`, `Pr(>|z|)`) %>%
           rename(z = `z value`,
                  p = `Pr(>|z|)`) %>%
           pivot_wider(names_from = term, 
                       values_from = c(z, p))) %>%
  ungroup


asv_models %>%
  rowwise(asv_id) %>%
  reframe(emmeans(model, ~time_anti_expose_health, data = data) %>%
            broom::tidy()) %>%
  separate(time_anti_expose_health, into = c('time', 'anti', 'exposure', 'health')) %>%
  mutate(time = as.integer(time)) %>%
  ggplot(aes(x = time, y = estimate, 
             ymin = estimate - std.error,
             ymax = estimate + std.error,
             colour = health,
             fill = exposure,
             shape = anti)) +
  geom_errorbar(position = position_dodge(0.5), width = 0.1) +
  geom_point(position = position_dodge(0.5), size = 2) +
  scale_shape_manual(values = c('A' = 'triangle filled', 'N' = 'triangle down filled')) +
  scale_colour_manual(values = c('H' = 'blue', 'D' = 'red')) +
  scale_fill_manual(values = c('N' = 'blue', 'D' = 'red')) +
  facet_wrap(~asv_id)
