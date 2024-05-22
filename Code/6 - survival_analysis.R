
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
  scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                  type = "continuous"),
                                         c('Healthy', 'Disease')),
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



# ph_model <- coxph(Surv(time, infected) ~ antibiotics * exposure, 
#       data = survival_data)
# car::Anova(ph_model)
# emmeans(ph_model, ~antibiotics * exposure) %>%
#   contrast('pairwise')

broom::tidy(km_model) %>% 
  mutate(strata = str_remove(strata, fixed('strata(treatment)='))) %>% 
  ggplot(aes(time, estimate, color = strata, shape = strata)) + 
  geom_point(size = 1) + 
  geom_step(linetype='dashed') +
  scale_color_brewer(palette = 'Set1') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x='Day', 
       y='Survival Probability', 
       color='Treatment', 
       shape='Treatment') + 
  theme_classic() 


#### Survivorship Model ####
full_model <- gam(time ~  antibiotics * exposure +
                    s(genotype, bs = 're', by = dummy), 
                  weights = infected,
                  family = cox.ph(),
                  data = survival_data,
                  method = 'REML')
summary(full_model)
anova(full_model)
emmeans(full_model, ~exposure * antibiotics) %>%
  contrast('pairwise')

emmeans(full_model, ~exposure * antibiotics) %>%
  cld(Letters = LETTERS) %>%
  broom::tidy() %>%
  mutate(.group = str_trim(.group)) %>%
  ggplot(aes(x = estimate, y = antibiotics,
             xmin = estimate - std.error, 
             xmax = estimate + std.error,
             colour = exposure)) +
  geom_errorbar(width = 0.1, position = position_dodge(0.25),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.25)) +
  geom_text(aes(x = Inf, label = .group), hjust = 1,
            position = position_dodge(0.25), 
            show.legend = FALSE)
  



survival_data %>%
  select(antibiotics, exposure, treatment) %>%
  distinct %>%
  mutate(dummy = 0,
         Block = 1,
         genotype = 1) %>%
  expand_grid(time = seq(0, 7.5, length.out = 1000)) %>%
  bind_cols(., predict(full_model, newdata = ., se.fit = TRUE, type = 'response')) %>%
  
  ggplot(aes(x = time, y = fit, colour = exposure,
             ymin = fit - se.fit, ymax = fit + se.fit,
             fill = exposure, linetype = antibiotics)) +
  pammtools::geom_stepribbon(alpha = 0.5, colour = NA, show.legend = FALSE) +
  pammtools::geom_stephazard(show.legend = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = 'Time Post-Exposure (D)',
       y = 'Fragment Survival (%)') +
  theme_classic() 


#### Experimental ####

full_model <- gam(time ~  s(antibiotics, bs = 're') +
                    s(exposure, bs = 're') +
                    s(treatment, bs = 're') +
                    s(genotype, bs = 're', by = dummy), 
                  weights = infected,
                  family = cox.ph(),
                  data = survival_data,
                  method = 'REML')

no_interaction <- gam(time ~  s(antibiotics, bs = 're') +
                        s(exposure, bs = 're') +
                        s(genotype, bs = 're', by = dummy), 
                      weights = infected,
                      family = cox.ph(),
                      data = survival_data,
                      method = 'REML')


no_exposure <- gam(time ~  s(antibiotics, bs = 're') +
                     s(genotype, bs = 're', by = dummy), 
                   weights = infected,
                   family = cox.ph(),
                   data = survival_data,
                   method = 'REML')

no_antibiotic <- gam(time ~  s(exposure, bs = 're') +
                       s(genotype, bs = 're', by = dummy), 
                     weights = infected,
                     family = cox.ph(),
                     data = survival_data,
                     method = 'REML')


anova(no_interaction, full_model, test = 'Chisq') #non-significant interaction
anova(no_exposure, no_interaction, test = 'Chisq') #significant exposure
anova(no_antibiotic, no_interaction, test = 'Chisq') #significant antibiotic


survival_data %>%
  select(antibiotics, exposure, treatment) %>%
  distinct %>%
  mutate(dummy = 0,
         Block = 1,
         genotype = 1) %>%
  expand_grid(time = seq(0, 7.5, length.out = 1000)) %>%
  bind_cols(., predict(full_model, newdata = ., se.fit = TRUE, type = 'response')) %>%
  
  ggplot(aes(x = time, y = fit, colour = exposure,
             ymin = fit - se.fit, ymax = fit + se.fit,
             fill = exposure, linetype = antibiotics)) +
  pammtools::geom_stepribbon(alpha = 0.5, colour = NA, show.legend = FALSE) +
  pammtools::geom_stephazard(show.legend = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = 'Time Post-Exposure (D)',
       y = 'Fragment Survival (%)') +
  theme_classic() 

survival_data %>%
  select(antibiotics, exposure, treatment) %>%
  distinct %>%
  mutate(dummy = 0,
         Block = 1,
         genotype = 1) %>%
  bind_cols(., predict(full_model, newdata = ., se.fit = TRUE)) %>%
  ggplot(aes(x = antibiotics, colour = exposure, 
             y = fit, ymin = fit - se.fit, ymax = fit + se.fit)) +
  geom_errorbar(width = 0.1, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5))


#### Simpler ####
library(flexsurv)
library(emmeans)

lung.model <- coxph(Surv(time, infected) ~ antibiotics * exposure, survival_data)
lung.model <- flexsurvreg(Surv(time, infected) ~ antibiotics * exposure, 
                          dist = 'gengamma',
                          data = survival_data)

emmeans(lung.model, ~antibiotics)

#### Bayes ####
library(rstanarm)
library(tidybayes)
tst <- stan_surv(Surv(time, infected) ~ antibiotics * exposure + (1 | genotype), 
                 data = survival_data,
                 basehaz = 'weibull',
                 cores = 4, 
                 chains = 4)

plot(tst, plotfun = 'basehaz')
ps_check(tst)
plot(tst)

as_tibble(tst) %>%
  transmute(Untreated_Healthy = `(Intercept)`,
            Antibiotic_Healthy = `(Intercept)` + antibioticsAntibiotic,
            Untreated_Disease = `(Intercept)` + exposureDisease,
            Antibiotic_Disease = `(Intercept)` + exposureDisease + `antibioticsAntibiotic:exposureDisease`) %>%
  pivot_longer(cols = everything(),
               names_to = c('antibiotics', 'exposure'),
               names_pattern = '(.*)_(.*)',
               values_to = 'hazard') %>%
  group_by(antibiotics, exposure) %>%
  point_interval() %>%
  mutate(antibiotics = fct_relevel(antibiotics, 'Untreated'),
         exposure = fct_relevel(exposure, 'Healthy')) %>%
  ggplot(aes(x = antibiotics, colour = exposure, 
             y = hazard, ymin = .lower, ymax = .upper)) +
  geom_errorbar(width = 0.1, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5))


tst

bind_cols(survival_model_minimal$family$data[c('tr', 'h', 'q', 'km')]) %>%
  mutate(se = sqrt(q),
         h = -h,
         km = -km) %>%
  
  ggplot(aes(x = tr, y = exp(h), ymin = exp(h - 1.96*se), ymax = exp(h + 1.96*se))) 




pred_data <- survival_data %>%
  select(antibiotics, exposure, treatment, genotype) %>%
  distinct %>%
  mutate(id = row_number())
  
pred_haz <- posterior_survfit(tst,
                              newdata = pred_data,
                              type = 'surv',
                              times = 0,
                              standardise = FALSE) %>%
  as_tibble

tail(pred_haz)

full_join(pred_data, pred_haz, by = 'id') %>%
  group_by(antibiotics, exposure, time) %>%
  summarise(fit = mean(median),
            se.fit = sd(median)) %>%
  
  ggplot(aes(x = time, y = fit, colour = exposure,
            ymin = fit - se.fit, ymax = fit + se.fit,
            fill = exposure, linetype = antibiotics)) +
  pammtools::geom_stepribbon(alpha = 0.5, colour = NA, show.legend = FALSE) +
  pammtools::geom_stephazard(show.legend = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = 'Time Post-Exposure (D)',
       y = 'Fragment Survival (%)') +
  theme_classic() 


# visualise survival curve for each plastic
# plot predictions for each plastic averaged over the random effect of each replicate for each treatment
# hoping this is equivalent to looking at reform=NA in brms and linear mixed effects models

library(emmeans)

emmeans(tst, ~antibiotics * exposure)

pred_data <- survival_data %>%
  select(antibiotics, exposure, treatment) %>%
  distinct

pred_out <- posterior_survfit(tst,
                              newdata = pred_data,
                              times = 0,
                              standardise = TRUE,
                              extrapolate = TRUE,
                              dynamic = TRUE)

d_preds <- select(d, plastic, sample) %>%
  distinct() %>%
  mutate(id = 1:n(),
         plastic2 = plastic) %>%
  nest_legacy(-c(plastic2)) %>%
  mutate(., preds = map(
    data,
    ~ posterior_survfit(
      mod,
      newdata = .x,
      times = 0,
      standardise = TRUE,
      extrapolate = TRUE,
      dynamic = TRUE
    )
  ))

d_preds <- unnest(d_preds, preds) %>%
  select(-data)

ggplot(d_preds,
       aes(
         time,
         col = plastic2,
         fill = plastic2,
         group = plastic2
       )) +
  geom_line(aes(y = median)) +
  geom_ribbon(
    aes(time, ymin = ci_lb, ymax = ci_ub),
    col = NA,
    alpha = 0.2,
    show.legend = FALSE
  ) +
  ylim(c(0, 1)) +
  theme_bw()


#### Fit as time-to-event regression ####
tst_data <- survival_data %>%
  group_by(time, antibiotics, exposure, treatment) %>%
  summarise(n = n(),
            s = sum(infected),
            .groups = 'drop') %>%
  nest(data = -time) %>%
  mutate(tm = time - (time - lag(time, 1, default = 0)) / 2) %>%
  unnest(data)

f <- function(x, d) pmin((x), 0)

# https://dpananos.github.io/posts/2024-01-20-logistic-survival/
fit <- glm(cbind(s, n - s) ~ treatment * (tm + I(f(tm - 2)^2) + I(f(tm - 2)^3)),
           data = tst_data, 
           family = binomial())

fit <- gam(cbind(s, n - s) ~ treatment + s(tm, by = treatment),
           data = tst_data, 
           family = 'binomial')

preds <- tst_data %>% 
  bind_cols(predict(fit, newdata = ., se.fit = T)) %>% 
  mutate(p = plogis(fit),
         p.low = plogis(fit - se.fit),
         p.high = plogis(fit + se.fit)) %>% 
  group_by(treatment) %>% 
  arrange(treatment, time) %>% 
  mutate(S = cumprod(1 - p),
         S.low = cumprod(1 - p.low),
         S.high = cumprod(1 - p.high))


preds %>% 
  ggplot(aes(tm, p, color = treatment)) +
  geom_line() + 
  scale_color_brewer(palette = 'Set1') +
  scale_y_continuous(labels = scales::percent) + 
  labs(x='Day', y='Hazard', 
       color='Treatment', shape='Treatment', 
       title='Hazard (Risk Estimate From Logistic Regression)')



broom::tidy(km_model) %>% 
  mutate(strata = str_remove(strata, fixed('strata(treatment)='))) %>% 
  ggplot(aes(time, estimate, color = strata, shape = strata)) + 
  geom_point(size = 1) + 
  geom_step(linetype='dashed') +
  geom_line(data = preds, aes(x = time, y = S, colour = treatment),
            inherit.aes = FALSE) +
  scale_color_brewer(palette = 'Set1') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x='Day', 
       y='Survival Probability', 
       color='Treatment', 
       shape='Treatment') + 
  theme_classic() 

# https://dpananos.github.io/posts/2023-12-02-gen-gamma/#fnref1
library(cmdstanr)
stan_data <- list(
  n = nrow(survival_data),
  n_trt = length(unique(survival_data$treatment)),
  time = survival_data$time,
  trt = survival_data$treatment,
  censored = 1 - survival_data$infected,
  
  nt = 1000,
  pred_time = seq(0, 7.5, length.out=1000),
  
  mu_df = 10,
  mu_loc = 0,
  mu_scale = 1,
  
  sigma_df = 30,
  sigma_loc = 1,
  sigma_scale = 1,
  
  k_df = 30,
  k_loc = 1,
  k_scale = 1
)

model <- cmdstan_model('genGamma.stan')
fit <- model$sample(stan_data, refresh = 100, parallel_chains = 4)

library(tidybayes)
fit %>% 
  spread_draws(survival_curve[i, trt]) %>% 
  mutate(time = stan_data$pred_time[i]) %>% 
  group_by(time, trt) %>% 
  mean_qi %>% 
  mutate(trt = as.character(trt)) %>% 
  ggplot(aes(time, survival_curve, ymin = survival_curve.lower, ymax = survival_curve.upper, color=trt, fill=trt)) + 
  geom_line() +
  # geom_step(
  #   data=bkm,
  #   aes(time, 1-estimate),
  #   inherit.aes = F
  # ) + 
  geom_ribbon(alpha=0.5, size=0) + 
  scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1')+
  labs(
    x='Time',
    y = expression(1-S(t))
  ) + 
  scale_y_continuous(labels=scales::percent)
