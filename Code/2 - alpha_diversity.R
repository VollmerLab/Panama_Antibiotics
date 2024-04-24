#### Libraries ####
library(tidyverse)
library(magrittr)
library(phyloseq)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(ComplexUpset)
library(broom)

alpha <- 0.05

#### Functions ####
t.test_jds <- function(x, y = NULL, 
                       alternative = c("two.sided", "less", "greater"),
                       mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
                       ...){
  alternative <- match.arg(alternative)
  
  if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if(!missing(conf.level) &&
     (length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if( !is.null(y) ) {
    dname <- paste(deparse(substitute(x)),"and",
                   deparse(substitute(y)))
    if(paired)
      xok <- yok <- complete.cases(x,y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    dname <- deparse(substitute(x))
    if (paired) stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if (paired) {
    x <- x-y
    y <- NULL
  }
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  if(is.null(y)) {
    if(nx < 2) stop("not enough 'x' observations")
    df <- nx-1
    stderr <- sqrt(vx/nx)
    if(stderr < 10 *.Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    tstat <- (mx-mu)/stderr
    method <- if(paired) "Paired t-test" else "One Sample t-test"
    estimate <-
      setNames(mx, if(paired)"mean of the differences" else "mean of x")
  } else {
    ny <- length(y)
    if(nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if(var.equal && nx+ny < 3) stop("not enough observations")
    my <- mean(y)
    vy <- var(y)
    method <- paste(if(!var.equal)"Welch", "Two Sample t-test")
    estimate <- c(mx,my)
    names(estimate) <- c("mean of x","mean of y")
    if(var.equal) {
      df <- nx+ny-2
      v <- 0
      if(nx > 1) v <- v + (nx-1)*vx
      if(ny > 1) v <- v + (ny-1)*vy
      v <- v/df
      stderr <- sqrt(v*(1/nx+1/ny))
    } else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))
    }
    if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df) )
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
  attr(cint,"conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, stderr_pair = c(stderrx, stderry),
               method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}

safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)

#### Data ####
tank_data <- read_rds("../intermediate_files/full_tank_microbiome.rds")

#### Calculate alpha metrics ####
alpha_diversity <- tank_data %>%
  microbiome::alpha(zeroes = FALSE) %>%
  as_tibble(rownames = 'sample_id') %>%
  pivot_longer(cols = -sample_id,
               names_to = 'alpha_metric') %>%
  left_join(sample_data(tank_data) %>%
              as_tibble(rownames = 'sample_id'),
            by = 'sample_id') %>%
  mutate(time = time_fac, .keep = 'unused')

sample_data(tank_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  filter(time == 0) %>%
  count(geno, anti)

#### Graphs ####
alpha_diversity %>%
  # filter(alpha_metric == 'observed') %>%
  ggplot(aes(x = time, y = value, colour = interaction(exposure, anti), 
             group = interaction(time, exposure, anti))) +
  geom_boxplot(position = position_dodge(0.5)) +
  facet_wrap(~alpha_metric, scales = 'free_y') +
  labs(colour = 'Treatment')

alpha_diversity %>%
  # filter(alpha_metric == 'observed') %>%
  ggplot(aes(x = time, y = value, colour = interaction(exposure, anti), 
             group = interaction(time, exposure, anti))) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5)) +
  facet_wrap(~alpha_metric, scales = 'free_y') +
  labs(colour = 'Treatment')

#### Model with/without antibiotics before disease dose ####
predose_effects <- alpha_diversity %>%
  filter(time == 'before',
         exposure == 'D') %>%
  nest_by(alpha_metric) %>%
  mutate(spread = abs(min(data$value) - max(data$value))) %>%
  reframe(ttest = list(with(data, t.test_jds(value[anti == 'N'], value[anti == 'A']))),
          
          broom::tidy(ttest)) %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr'),
         q.value = safe_qvalue(p.value))

filter(predose_effects, fdr < alpha)

#### Model outcome/antibiotics after disease dose ####
postDose_effects <- alpha_diversity %>%
  filter(time == 'after', exposure == 'D') %>%
  mutate(treatment = str_c(anti, health, sep = '_')) %>%
  nest_by(alpha_metric) %>%
  reframe(treatment_model = list(aov(value ~ treatment, data = data)),
          tidy(treatment_model) %>%
            mutate(dDF = df[term == 'Residuals']) %>%
            filter(term != 'Residuals') %>%
            select(-term),
          
          posthoc = list(emmeans(treatment_model, ~treatment)),
          
          planned_contrast = list(contrast(posthoc,
                                           method = list(disease.v.avg = c(-1/2, 1, -1/2),
                                                         disease.v.h = c(0, 1, -1),
                                                         disease.v.anti = c(-1, 1, 0)), 
                                           side = '>')),
          
          tidy(planned_contrast) %>%
            select(contrast, p.value) %>%
            pivot_wider(names_from = contrast, values_from = p.value,
                        names_prefix = 'p.value_')) %>%
  
  mutate(fdr = p.adjust(p.value, method = 'fdr'),
         q.value = safe_qvalue(p.value)) %>%
  mutate(across(starts_with('p.value_disease'), 
                ~p.adjust(., method = 'fdr'),
                .names = 'fdr_{.col}')) %>%
  rename_with(~str_replace_all(., '_p.value_', '_'))


filter(postDose_effects, fdr < alpha)



#### Model only disease exposed through time ####
alpha_models <- alpha_diversity %>%
  filter(exposure == 'D') %>%
  nest_by(alpha_metric) %>%
  
  mutate(model = list(lmer(value ~ time + anti + health + 
                             (1 | tank) + (1 | geno / fragment),
                           data = data)),
         
         anova(model, ddf = 'Satterthwaite') %>%
           as_tibble(rownames = 'term') %>%
           rename(nDF = NumDF,
                  dDF = DenDF,
                  fvalue = `F value`,
                  pvalue = `Pr(>F)`) %>%
           select(-contains(' ')) %>%
           pivot_wider(names_from = term,
                       values_from = -term),
         
         
         effect_time = list(emmeans(model, ~time, 
                                    lmer.df = 'satterthwaite') %>%
                              contrast('revpairwise')),
         effect_anti = list(emmeans(model, ~anti, 
                                          lmer.df = 'satterthwaite') %>%
                                    contrast('revpairwise')),
         effect_health = list(emmeans(model, ~health, 
                                       lmer.df = 'satterthwaite') %>%
                                 contrast('revpairwise')),
         
         across(starts_with('effect'), 
                ~broom::tidy(.) %>%
                  pull(estimate),
                .names = "coef_{.col}"),
         
         plot_data = list(emmeans(model, ~time + anti + health) %>% 
                            broom::tidy(conf.int = TRUE))) %>%
  ungroup %>%
  rename_with(~str_replace_all(., c('coef_effect' = 'coef')))

#### Upset Plot ####
alpha_models %>%
  select(alpha_metric, starts_with('pvalue'), starts_with('coef')) %>%
  pivot_longer(cols = c(starts_with('pvalue'), starts_with('coef')),
               names_to = c('.value', 'metric'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(name = case_when(metric == 'time' & pvalue > 0.05 ~ 'no_tme',
                            metric == 'time' & pvalue < 0.05 & coef < 0 ~ 'before',
                            metric == 'time' & pvalue < 0.05 & coef > 0 ~ 'after',
                            
                            metric == 'anti' & pvalue > 0.05 ~ 'no_anti',
                            metric == 'anti' & pvalue < 0.05 & coef < 0 ~ 'anti',
                            metric == 'anti' & pvalue < 0.05 & coef > 0 ~ 'untreated',
                            
                            metric == 'health' & pvalue > 0.05 ~ 'no_health',
                            metric == 'health' & pvalue < 0.05 & coef < 0 ~ 'disease',
                            metric == 'health' & pvalue < 0.05 & coef > 0 ~ 'healthy'),
         .keep = 'unused') %>%
  mutate(value = TRUE) %>%
  pivot_wider(values_fill = FALSE) %>%
  select(-starts_with('no')) %>%
  mutate(before = FALSE) %>%
  select(alpha_metric, before, after, 
         disease, healthy,
         untreated, anti) %>%
  
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        intersections = list('before',
                           'after',
                           'disease',
                           'healthy',
                           'untreated',
                           'anti',
                           'Outside of known sets'),
        sort_intersections = FALSE,
        sort_sets = FALSE)
ggsave('../Results/alpha_diversity_metic_associations.png')

#### Individual Metric Plots ####
model <- alpha_models$model[[1]]
pvalues <- c(FALSE, TRUE, FALSE)
pvalue_names <- c('time', 'anti', 'health')

plot_metric <- function(model, pvalues, pvalue_names){
  sig_vars <- pvalue_names[pvalues]
  
  plot_data <- emmeans(model, as.formula(str_c('~', sig_vars))) %>%
    broom::tidy(conf.int = TRUE)
  
  
}

library(patchwork)
tmp <- alpha_models %>%
  select(alpha_metric, starts_with('pvalue')) %>%
  pivot_longer(cols = starts_with('pvalue'),
               names_to = 'term',
               values_to = 'pvalue',
               names_prefix = 'pvalue_') %>%
  filter(pvalue < 0.05) %>%
  group_by(alpha_metric) %>%
  summarise(sig_vars = str_c(term, collapse = ', ')) %>%
  left_join(select(alpha_models, alpha_metric, plot_data),
            by = 'alpha_metric') %>%
  rowwise %>%
  mutate(plot = list(plot_data %>%
                       mutate(time = factor(time, levels = c('before', 'after'))) %>%
                       ggplot(aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high, 
                                  fill = health, shape = anti)) +
                       geom_errorbar(width = 0.1, position = position_dodge(0.5)) +
                       geom_point(size = 3, position = position_dodge(0.5)) +
                       scale_shape_manual(values = c('N' = 'square filled',
                                                     'A' = 'diamond filled')) +
                       guides(shape = guide_legend(override.aes = list(fill = 'black')),
                              fill = guide_legend(override.aes = list(shape = 'circle filled'))) +
                       labs(y = alpha_metric,
                            title = sig_vars))) %>%
  group_by(sig_vars) %>%
  summarise(plot = list(wrap_plots(plot, nrow = 1)))

wrap_plots(tmp$plot, guides = 'collect', byrow = TRUE, ncol = 1)
ggsave('../Results/alpha_diversity_values.png', height = 15, width = 15)

#### Model just pre-treatment - i.e. antibiotic effect ####
alpha_diversity %>%
  filter(exposure == 'pre') %>%
  group_by(alpha_metric) %>%
  summarise(broom::tidy(t.test(value ~ anti)),
            .groups = 'drop') %>%
  filter(p.value < 0.05)

library(glmmTMB)

tst <- alpha_diversity %>%
  filter(exposure == 'pre') %>%
  nest_by(alpha_metric) %>%
  mutate(model = list(glmmTMB(value ~ anti + (1 | tank) + (1 | geno),
                                data = data)))

tst$model[[5]] %>% summary %>% coef

antibiotic_effects <- alpha_diversity %>%
  filter(exposure == 'pre') %>%
  nest_by(alpha_metric) %>%
  mutate(model = list(lmer(value ~ anti + (1 | tank) + (1 | geno),
                           data = data)),
         anova(model, ddf = 'Kenward-Roger') %>%
           as_tibble,
         
         VarCorr(model) %>%
           as_tibble %>%
           mutate(var = sdcor^2) %>%
           mutate(pct_var = var / sum(var)) %>%
           filter(grp != 'Residual') %>%
           select(grp, pct_var) %>%
           pivot_wider(names_from = grp,
                       values_from = pct_var,
                       names_prefix = 'pctVar_')) %>%
  
  janitor::clean_names() %>%
  ungroup %>%
  mutate(fdr = p.adjust(pr_f, 'fdr'))

tank_significanct_metrics <- filter(antibiotic_effects, pr_f < 0.05) %>%
  pull(alpha_metric)

filter(antibiotic_effects, alpha_metric %in% tank_significanct_metrics) %>%
  rowwise(alpha_metric) %>%
  reframe(emmeans(model, ~anti) %>% 
            broom::tidy()) %>%
  ggplot(aes(x = anti, y = estimate, ymin = estimate - std.error, 
             ymax = estimate + std.error)) +
  geom_pointrange() 

#### Model metrics by exposure, time, and antibiotic treatment ####
exposure_models <- alpha_diversity %>%
  filter(exposure != 'pre') %>%
  nest_by(alpha_metric) %>%
  mutate(model = list(lmer(value ~ time * anti * exposure + 
                             (1 | tank) + (1 | geno / fragment),
                           data = data)),
         anova(model, ddf = 'Kenward-Roger') %>%
           as_tibble(rownames = 'term') %>%
           mutate(term = str_remove_all(term, 'as.factor\\(|\\)') %>%
                    str_replace_all(':', 'X')) %>%
           janitor::clean_names() %>%
           pivot_wider(names_from = 'term',
                       values_from = -term,
                       names_glue = "{term}.{.value}"),
         
         VarCorr(model) %>%
           as_tibble %>%
           mutate(var = sdcor^2) %>%
           mutate(pct_var = var / sum(var)) %>%
           filter(grp != 'Residual') %>%
           select(grp, pct_var) %>%
           pivot_wider(names_from = grp,
                       values_from = pct_var,
                       names_prefix = 'pctVar_')) %>%
  
  ungroup %>%
  mutate(across(contains('pr_f'), ~p.adjust(., 'fdr'),
                .names = '{.col}_fdr')) %>%
  rename_with(.cols = ends_with('fdr'),
              ~str_remove(., 'pr_f_'))

significant_exposure_metrics <- exposure_models %>%
  # select(ends_with('pr_f'))
  filter(if_any(c(contains('pr_f'), -contains('time')), ~. < 0.05)) %>%
  pull(alpha_metric)

#### Subset plots to only significantly differing effects ####

alpha_diversity %>%
  filter(alpha_metric %in% 
           union(tank_significanct_metrics, 
                 significant_exposure_metrics)) %>%
  # filter(alpha_metric == 'observed') %>%
  ggplot(aes(x = time, y = value, colour = interaction(exposure, anti), 
             group = interaction(time, exposure, anti))) +
  stat_summary(fun.data = mean_se, position = position_dodge(0.5)) +
  facet_wrap(~alpha_metric, scales = 'free_y') +
  labs(colour = 'Treatment')
ggsave('../Results/alpha_diversity..png', height = 7, width = 7)
