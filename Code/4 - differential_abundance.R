#### Libraries ####
library(multcomp)
library(tidyverse)
library(multidplyr)
library(broom)
library(lmerTest)
library(qvalue)
library(emmeans)
library(ComplexUpset)
library(cowplot)

alpha <- 0.05
fdr_or_qvalue <- 'fdr'
refit_models <- TRUE

#### Functions ####
safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)

#### Data ####
tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
                      show_col_types = FALSE) %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac)



#### Just time 0 between antibiotics ####
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

predose_effects <- tank_data %>%
  filter(time == 'before',
         exposure == 'D') %>%
  nest_by(across(domain:species), asv_id) %>%
  mutate(spread = abs(min(data$log2_cpm_norm) - max(data$log2_cpm_norm))) %>%
  filter(spread > 1) %>%
  reframe(ttest = list(with(data, t.test_jds(log2_cpm_norm[anti == 'N'], log2_cpm_norm[anti == 'A']))),
          
          broom::tidy(ttest)) %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr'),
         q.value = safe_qvalue(p.value))

predose_effects %>%
  filter(q.value < 0.05) %>%
  select(domain:species, asv_id, estimate) %>%
  distinct 

#### Just After Dose comparison ####
# postDose_effects <- tank_data %>%
#   filter(time == 'after') %>%
#   nest_by(across(domain:species), asv_id) %>%
#   mutate(spread = abs(min(data$log2_cpm_norm) - max(data$log2_cpm_norm))) %>%
#   filter(spread > 1) %>%
#   reframe(model = list(aov(log2_cpm_norm ~ anti * exposure, data = data)),
#           tidy(model) %>%
#             select(term, p.value) %>%
#             filter(term != 'Residuals') %>%
#             pivot_wider(names_from = term, values_from = 'p.value',
#                         names_prefix = 'p.value_') %>%
#             rename_with(~str_replace_all(., ':', 'X'))) %>%
#   mutate(across(starts_with('p.value'), ~p.adjust(., method = 'fdr'), .names = 'fdr_{.col}'),
#          across(starts_with('p.value'), ~safe_qvalue(.), .names = 'q.value_{.col}')) %>%
#   rename_with(~str_replace_all(., '_p.value_', '_'))

make_plot_data <- function(preT.test, anova_model, posthoc){
  emmeans(anova_model, ~exposure * anti) %>%
    cld(Letters = LETTERS) %>%
    broom::tidy() %>%
    mutate(time = 'after',
           .group = str_trim(.group)) %>%
    bind_rows(tibble(anti = c('N', 'A'),
                     estimate = preT.test$estimate,
                     std.error = preT.test$stderr_pair,
                     .group = if_else(preT.test$p.value < alpha, '*', ''),
                     time = 'before',
                     exposure = 'pre')) %>%
    mutate(time = factor(time, levels = c('before', 'after')))
}


#### Just disease dose after ####
postDose_effects <- tank_data %>%
  filter(time == 'after', exposure == 'D') %>%
  mutate(treatment = str_c(anti, health, sep = '_')) %>%
  nest_by(across(domain:species), asv_id) %>%
  reframe(treatment_model = list(aov(log2_cpm_norm ~ treatment, data = data)),
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


#### Combine ####
tst <- full_join(select(predose_effects,
                 asv_id, ttest, p.value, fdr) %>%
            rename(p.value_preDose = p.value,
                   fdr_preDose = fdr),
          rename(postDose_effects,
                 p.value_treatment = p.value,
                 fdr_treatment = fdr),
          by = 'asv_id') 

tst %>%
  select(asv_id, starts_with(fdr_or_qvalue)) %>%
  mutate(across(starts_with(fdr_or_qvalue), ~. < alpha)) %>%
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames)
ggsave('../Results/asv_upset.png', 
       height = 12, width = 10, bg = 'white')

# preT.test <- tst$ttest[[1]]; posthoc <- tst$posthoc[[1]]; contrast <- tst$planned_contrast[[1]]
make_plot_data <- function(preT.test, posthoc, contrast){
  posthoc %>%
    cld(Letters = LETTERS) %>%
    broom::tidy() %>%
    mutate(time = 'after',
           .group = str_trim(.group)) %>%
    separate(treatment, into = c('anti', 'health')) %>%
    mutate(simple_sig = if_else(all(tidy(contrast)$p.value < alpha), '*', '')) %>%
    bind_rows(tibble(anti = c('N', 'A'),
                     estimate = preT.test$estimate,
                     std.error = preT.test$stderr_pair,
                     simple_sig = if_else(preT.test$p.value < alpha, '*', ''),
                     time = 'before', health = 'H',
                     exposure = 'pre')) %>%
    mutate(time = factor(time, levels = c('before', 'after')))
}

tst %>%
  filter(asv_id %in% c('ASV25', 'ASV8')) %>%
  rowwise(asv_id) %>%
  reframe(make_plot_data(ttest, posthoc, planned_contrast)) %>%

  ggplot(aes(x = time, y = estimate, 
             ymin = estimate - std.error,
             ymax = estimate + std.error,
             shape = anti,
             colour = health)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), 
             size = 4) +
  
  geom_text(data = . %>% filter((time == 'after' & health == 'D') | 
                                  (time == 'before' & anti == 'A')),
            aes(y = Inf, label = simple_sig),
            colour = 'black', vjust = 1, size = 8) +
  
  facet_wrap(~asv_id) +
  labs(x = NULL, 
       y = 'log2CPM')



tst %>% 
  # filter(asv_id == 'ASV8') %>%
  filter(fdr_treatment < alpha) %>%
  filter(fdr_disease.v.avg < alpha, 
         p.value_disease.v.h < alpha, 
         p.value_disease.v.anti < alpha) %>% 
  mutate(name = case_when(!is.na(species) & !is.na(family) ~ 
                            str_c(family, '\n', 
                                  species, ' (', asv_id, ')'),
                          !is.na(genus) & !is.na(family) ~ 
                            str_c(family, '\n', 
                                  genus, ' (', asv_id, ')'),
                          !is.na(family) ~ 
                            str_c(family, '\n', 
                                  '(', asv_id, ')'),
                          !is.na(order) ~ 
                            str_c(order, '\n', 
                                  '(', asv_id, ')'),
                          !is.na(phylum) ~ 
                            str_c(phylum, '\n', 
                                  '(', asv_id, ')'),
                          TRUE ~ asv_id)) %>%
  # filter(!is.na(estimate)) %>%
  rowwise(name) %>%
  reframe(make_plot_data(ttest, posthoc, planned_contrast)) %>%
  
  ggplot(aes(x = time, y = estimate, 
             ymin = estimate - std.error,
             ymax = estimate + std.error,
             shape = anti,
             colour = health)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), 
             size = 4) +
  
  geom_text(data = . %>% filter((time == 'after' & health == 'D') | 
                                  (time == 'before' & anti == 'A')),
            aes(y = Inf, label = simple_sig),
            colour = 'black', vjust = 1, size = 8) +
  
  facet_wrap(~name, scales = 'free_y') +
  labs(x = NULL, 
       y = 'log2CPM')
ggsave('../Results/asvs_changing_postExposure.png', height = 15, width = 15)


full_join(predose_effects,
          postDose_effects,
          by = join_by(domain, phylum, class, order, 
                       family, genus, species, asv_id)) %>%
  select(-ttest, -estimate1, -estimate2, -statistic, -parameter:-alternative,
         -model) %>%
  rename_with(~str_c('preDose_', .), .cols = all_of(c('estimate', 'p.value', 'fdr', 'q.value'))) %>%
  write_csv('../Results/individual_asv_results.csv')

#### Model ASVs ####
if(file.exists('../intermediate_files/asv_models.rds.xz') & !refit_models){
  asv_models <- read_rds('../intermediate_files/asv_models.rds.xz')
} else {
  
  cluster <- new_cluster(parallel::detectCores() - 1)
  cluster_library(cluster, c('dplyr', 'lmerTest', 'emmeans', 
                             'stringr', 'tidyr'))
  
  asv_models <- tank_data %>%
    filter(exposure == 'D') %>%
    nest_by(across(domain:species), asv_id) %>%
    # filter(asv_id %in% c('ASV25', 'ASV8')) %>%
    
    partition(cluster) %>%
    mutate(model = list(lmer(log2_cpm_norm ~ time * anti + health + 
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
    collect %>%
    ungroup %>%
    mutate(across(starts_with('pvalue'), 
                  ~p.adjust(., 'fdr'),
                  .names = 'fdr_{.col}')) %>%
    rename_with(~str_replace(., 'fdr_pvalue', 'fdr')) %>%
    mutate(across(starts_with('pvalue'), safe_qvalue,
                  .names = 'qvalue_{.col}')) %>%
    rename_with(~str_replace_all(., c('coef_effect' = 'coef',
                                      'qvalue_pvalue' = 'qvalue')))
  
  cluster <- NULL
  write_rds(asv_models, '../intermediate_files/asv_models.rds.xz', compress = 'xz')
}

#### Write Results ####
asv_models %>%
  select(-where(is.list)) %>%
  mutate(association_time = if_else(coef_time < 0, 'before', 'after'),
         association_anti = if_else(coef_anti < 0, 'untreated', 'antibiotic'),
         association_health = if_else(coef_time < 0, 'healthy', 'disease')) %>%
  write_csv('../Results/individual_asv_results.csv')


#### Upset Plot ####
colour_options <- read_rds('../../Panama_Tank_Field/intermediate_files/asv_colors.rds')
microbe_colors <- set_names(colour_options$color_palette$hex,
                            colour_options$color_palette$group)
levels(colour_options$asv_clumping$Top_order)


base_plot_data <- asv_models %>%
  select(domain:species, asv_id, starts_with(fdr_or_qvalue), starts_with('coef')) %>%
  pivot_longer(cols = c(starts_with(fdr_or_qvalue), starts_with('coef')),
               names_to = c('.value', 'metric'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(name = case_when(metric == 'time' & !!sym(fdr_or_qvalue) > 0.05 ~ 'no_tme',
                          metric == 'time' & !!sym(fdr_or_qvalue) < 0.05 & coef < 0 ~ 'before',
                          metric == 'time' & !!sym(fdr_or_qvalue) < 0.05 & coef > 0 ~ 'after',
                          
                          metric == 'anti' & !!sym(fdr_or_qvalue) > 0.05 ~ 'no_anti',
                          metric == 'anti' & !!sym(fdr_or_qvalue) < 0.05 & coef > 0 ~ 'anti',
                          metric == 'anti' & !!sym(fdr_or_qvalue) < 0.05 & coef < 0 ~ 'untreated',
                          
                          metric == 'health' & !!sym(fdr_or_qvalue) > 0.05 ~ 'no_health',
                          metric == 'health' & !!sym(fdr_or_qvalue) < 0.05 & coef > 0 ~ 'disease',
                          metric == 'health' & !!sym(fdr_or_qvalue) < 0.05 & coef < 0 ~ 'healthy'),
         .keep = 'unused') %>%
  mutate(value = TRUE) %>%
  pivot_wider(values_fill = FALSE) %>%
  select(-starts_with('no')) %>%
  mutate(anti = FALSE) %>%
  select(domain:species, asv_id, 
         before, after, 
         disease, healthy,
         untreated, anti) %>%
  
  left_join(mutate(colour_options$asv_clumping,
                   genus = str_remove_all(genus, '\\<i\\>|\\</i\\>')),
            by = c('order', 'genus')) %>%
  mutate(group = case_when(!is.na(group) ~ group,
                           is.na(group) & order == 'Oceanospirillales' ~ 'Oceanospirillales-Other',
                           is.na(group) & order == 'Vibrionales' ~ 'Vibrionales-Other',
                           is.na(group) & order == 'Alteromonadales' ~ 'Alteromonadales-Other',
                           is.na(group) & order == 'Verrucomicrobiales' ~ 'Verrucomicrobiales-Other',
                           is.na(group) & order == 'Campylobacterales' ~ 'Campylobacterales-Other',
                           TRUE ~ 'Other-Other'),
         group = factor(group, levels = levels(colour_options$asv_clumping$group))) %>%
  # filter(!is.na(group)) %>%
  select(-Top_order, -Top_genus) %>%
  rename(the_colour = group) 


base_plot <- base_plot_data %>%
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames,
        base_annotations = list(
          'Intersection size' = intersection_size(
            counts = TRUE,
          ) + 
            scale_y_continuous(labels=scales::comma_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            labs(tag = 'B') +
            theme_classic() +
            theme(legend.position = 'none',
                  panel.background = element_rect(colour = 'black'),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text = element_text(colour = 'black', size = 10),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                          face = 'bold', hjust = 0))
        ),
        annotations = list(
          'Order' = ggplot(mapping = aes(fill = the_colour)) +
            geom_bar(fill = 'white', colour = 'black', stat = 'count', 
                     position = 'fill', show.legend = FALSE) +
            geom_bar(stat = 'count', position = 'fill', show.legend = FALSE) +
            scale_y_continuous(labels=scales::percent_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            scale_fill_manual(values = microbe_colors) +
            labs(y = 'Relative Abundance',
                 tag = 'A') +
            theme_classic() +
            theme(legend.position = 'none',
                  panel.background = element_rect(colour = 'black'),
                  axis.ticks.y = element_line(colour = 'black'),
                  axis.minor.ticks.y.left = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text = element_text(colour = 'black', size = 10),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                          face = 'bold', hjust = 0))
        ),
        
        sort_sets = FALSE, 
        set_sizes=(
          upset_set_size(geom = geom_bar(fill = c('gray65', 'gray65', '#F21A00', '#3B9AB2', 'purple')),
                         position = 'right') + #'gray65'
            geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = 'count',
                      colour = 'white', fontface = 'bold', size = 4) +
            # + annotate(geom='text', label='@', x='Drama', y=850, color='white', size=3) +
            # expand_limits(y = 200) +
            labs(tag = 'D') +
            theme_classic() +
            theme(panel.background = element_rect(colour = 'black'),
                  axis.ticks.x = element_line(colour = 'black'),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text.x = element_text(angle = 0,  colour = 'black', size = 10),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.tag = element_text(colour = 'black', size = 18,
                                          face = 'bold', hjust = 0))
        ),
        
        matrix = (
          intersection_matrix() +
            labs(tag = 'C') +
            theme_classic() 
        ),

        # sort_intersections=FALSE,
        # intersections=list(c('Healthy', 'Time'),
        #                    'Healthy',
        #                    c('Disease', 'Time'),
        #                    'Disease',
        #                    'Time',
        #                    'Outside of known sets'),
        
        name = NULL,
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(color = c('black', 'black', '#F21A00', '#3B9AB2', 'purple', 'purple'),
                                                                                                    size = 14, face = 'bold'),
                                                                         panel.background = element_rect(colour = 'black'),
                                                                         panel.border = element_rect(colour = 'black', fill = 'transparent'),
                                                                         plot.tag = element_text(colour = 'black', size = 18, 
                                                                                                 face = 'bold', hjust = 0),
                                                                         axis.ticks = element_line(colour = 'black'),
                                                                         axis.line = element_line(colour = "black", 
                                                                                                  linewidth = rel(1))))),
        upset_stripes(colors = NA)
  )

ggdraw(base_plot) +
  draw_plot(colour_options$legend, 
            x = 0.75, y = 0.25, width = 0.25, height = 0.75)
ggsave('../Results/asv_upset.png', 
       height = 12, width = 10, bg = 'white')


#### Who we care about ####
asv_models %>%
  filter(asv_id %in% c('ASV25')) %>%
  pull(effect_health)
