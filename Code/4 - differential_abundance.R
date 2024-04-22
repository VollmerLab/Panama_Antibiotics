#### Libraries ####
library(tidyverse)
library(multidplyr)
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
    mutate(model = list(lmer(log2_cpm_norm ~ time + anti + health + 
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
