
#### Libraries ####
library(multcomp)
library(tidyverse)
library(broom)
library(lmerTest)
library(qvalue)
library(ComplexUpset)
# library(biostat)
library(cowplot)
library(patchwork)
# library(metap)
library(ggtext)

alpha_test <- 0.05
fdr_or_qvalue <- '_fdr.bh'#'_q.value'
y_metric <- 'rclr' #'log2_cpm_norm'

#### Functions ####
safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_) #better method for doing this! - only want to sub 0 when it fails

tidy_model <- function(model, DF_method){
  the_anova <- anova(model, ddf = DF_method)
  
  as_tibble(the_anova) %>%
    rename(ndf = NumDF,
           ddf = DenDF,
           f.value = `F value`,
           p.value = `Pr(>F)`)
}

posthoc_padjust <- function(data, alpha){
  
  #FDR.bh
  fdrBH_step_data <- filter(data, fdr.bh < alpha) %>%
    mutate(across(matches('(disease|time|antibiotic)_p.value'), 
                  list(fdr.bh = ~p.adjust(., 'fdr')))) %>%
    rename_with(~str_replace_all(., '_p.value_', '_')) %>%
    bind_rows(filter(data, fdr.bh >= alpha | is.na(fdr.bh))) 
  
  #FDR.by
  fdrBY_step_data <- filter(fdrBH_step_data, fdr.by < alpha) %>%
    mutate(across(matches('(disease|time|antibiotic)_p.value'), 
                  list(fdr.by = ~p.adjust(., 'BY')))) %>%
    rename_with(~str_replace_all(., '_p.value_', '_')) %>%
    bind_rows(filter(fdrBH_step_data, fdr.by >= alpha | is.na(fdr.by))) 
  
  #Q value
  filter(fdrBY_step_data, q.value < alpha) %>%
    mutate(across(matches('(disease|time|antibiotic)_p.value'), 
                  list(q.value = ~qvalue(., lambda = 0)$qvalues))) %>%
    rename_with(~str_replace_all(., '_p.value_', '_')) %>%
    bind_rows(filter(fdrBY_step_data, q.value >= alpha | is.na(q.value))) %>%
    arrange(str_extract(asv_id, '[0-9]+') %>% as.integer())
  
}

contrast_to_table <- function(the_contrast){
  broom::tidy(the_contrast) %>%
    select(contrast, estimate, std.error, df, statistic, p.value) %>%
    rename(t.value = statistic) %>%
    pivot_wider(names_from = contrast, 
                values_from = -contrast,
                names_glue = '{contrast}_{.value}')
}

#### Data ####
tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
                      show_col_types = FALSE) %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac) %>%
  mutate(treatment = str_c(anti, health, sep = '_') %>% factor,
         time_treat = str_c(time, anti, health, sep = '_') %>% factor)

sebastians_field <- read_csv('../intermediate_files/sebastians_reef_microbiomes.csv', 
                             show_col_types = FALSE)

#### Run analyses ####
if(file.exists('../intermediate_files/all_asv_models.rds.xz')){
  prePost_effects <- read_rds('../intermediate_files/all_asv_models.rds.xz')
} else {
  library(multidplyr)
  cluster <- new_cluster(parallel::detectCores() - 1)
  cluster_library(cluster, c('multcomp', 'emmeans',
                             'tidyverse', 'broom',
                             'lmerTest'))
  cluster_copy(cluster, c('y_metric', 'tidy_model', 'alpha_test'))
  
  prePost_effects <- tank_data %>%
    
    #add a tiny bit of noise to data
    # mutate(!!sym(y_metric) := !!sym(y_metric) + rnorm(nrow(.), 0, 1e-8)) %>%
    
    nest_by(across(domain:species), asv_id) %>%
    partition(cluster) %>%
    
    summarise(model = list(lmer(!!sym(y_metric) ~ time_treat + 
                                  (1 | tank) + 
                                  (1 | geno / fragment),
                                data, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))),
              
              tidy_model(model, 'Satterthwaite'),
              
              #Posthocs
              posthoc = list(emmeans::emmeans(model, 
                                              ~ time_treat, 
                                              data = data, 
                                              lmer.df = 'satterthwaite',
                                              adjust = 'none')),
              contrast = list(emmeans::contrast(posthoc,
                                                method = list('disease' = c(0, 1, -1/2, 0, -1/2),
                                                              'antibiotic' = c(1/2, 0, -1/2, 1/2, -1/2),
                                                              'time' = c(1/2, 0, 1/2, -1/2, -1/2)),
                                                adjust = 'none')),
              
              contrast %>%
                tidy %>%
                select(contrast, estimate, p.value) %>%
                rename(groups = contrast,
                       coef = estimate) %>%
                pivot_wider(names_from = groups, 
                            values_from = c(coef, p.value), 
                            names_glue = "{groups}_{.value}"),
              
              #Plotting
              metrics = list(cld(posthoc,
                                 Letters = LETTERS, 
                                 alpha = alpha_test) %>%
                               tidy(conf.int = TRUE) %>%
                               select(time_treat, estimate, std.error, 
                                      conf.low, conf.high, .group) %>%
                               mutate(.group = str_trim(.group)) %>%
                               rename(mean_est = estimate,
                                      se_est = std.error) %>%
                               separate(time_treat, into = c('time', 'anti', 'health')) %>%
                               mutate(anti = factor(anti, levels = c('N', 'A')),
                                      health = factor(health, levels = c('H', 'D')),
                                      time = factor(time, levels = c('before', 'after'))))) %>%
    collect %>%
    ungroup %>%
    
    mutate(fdr.bh = p.adjust(p.value, 'BH'),
           fdr.by = p.adjust(p.value, 'BY'),
           q.value = safe_qvalue(p.value), 
           .after = 'p.value') %>%
    
    posthoc_padjust(alpha_test) %>%
    mutate(name = case_when(!is.na(species) & !is.na(family) ~ 
                              str_c(family, '<br></br>', 
                                    '<i>', species, '</i> (', asv_id, ')'),
                            !is.na(genus) & !is.na(family) ~ 
                              str_c(family, '<br></br>', 
                                    '<i>', genus, ' sp.</i> (', asv_id, ')'),
                            !is.na(family) ~ 
                              str_c(family, '<br></br>', 
                                    '(', asv_id, ')'),
                            !is.na(order) ~ 
                              str_c(order, '<br></br>', 
                                    '(', asv_id, ')'),
                            !is.na(phylum) ~ 
                              str_c(phylum, '<br></br>', 
                                    '(', asv_id, ')'),
                            TRUE ~ asv_id)) 
  
  cluster <- NULL
  write_rds(prePost_effects, '../intermediate_files/all_asv_models.rds.xz', compress = 'xz')  
}

prePost_effects %>%
  select(-where(is.list), -name) %>%
  write_csv('../Results/individual_asv_results.csv')

#### Numbers for Paper ####
prePost_effects %>%
  filter(asv_id == 'ASV8') %>%
  pull(posthoc) %>%
  pluck(1) 

prePost_effects %>%
  filter(asv_id == 'ASV8') %>%
  pull(contrast) %>%
  pluck(1) %>%
  broom::tidy()
mean(exp(c(1.41, 1.9)))

#### Summarize Results ####
prePost_effects %>%
  select(asv_id, class, order, family, genus, contains(str_remove(fdr_or_qvalue, '^_'))) %>%
  mutate(across(contains(str_remove(fdr_or_qvalue, '^_')), ~replace_na(., 1))) %>%
  rename(!!sym(str_c('global', fdr_or_qvalue)) := !!sym(str_remove(fdr_or_qvalue, '^_'))) %>%
  pivot_longer(cols = contains(str_remove(fdr_or_qvalue, '^_')),
               names_to = 'metric',
               values_to = str_remove(fdr_or_qvalue, '^_'), 
               values_drop_na = TRUE) %>%
  mutate(metric = str_remove(metric, fdr_or_qvalue)) %>%
  group_by(metric) %>%
  summarise(total_asv = n(),
            total_class = n_distinct(class),
            total_family = n_distinct(family),
            sig_class = n_distinct(class[!!sym(str_remove(fdr_or_qvalue, '^_')) < alpha_test]),
            # sig_order = n_distinct(order[!!sym(str_remove(fdr_or_qvalue, '^_')) < alpha_test]),
            sig_family = n_distinct(family[!!sym(str_remove(fdr_or_qvalue, '^_')) < alpha_test]),
            # sig_genus = n_distinct(genus[!!sym(str_remove(fdr_or_qvalue, '^_')) < alpha_test]),
            sig_asv = n_distinct(asv_id[!!sym(str_remove(fdr_or_qvalue, '^_')) < alpha_test]))
  

asv_assocation_directions <- prePost_effects %>%
  select(asv_id, class, order, family, genus, 
         contains(fdr_or_qvalue), contains('coef')) %>%
  pivot_longer(cols = c(contains(fdr_or_qvalue), contains('coef')),
               names_to = c('metric', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  filter(!is.na(!!sym(str_remove(fdr_or_qvalue, '^_')))) %>%
  filter(!!sym(str_remove(fdr_or_qvalue, '^_')) < alpha_test) %>%
  select(order, family, metric, asv_id, coef) %>%
  mutate(association = case_when(coef > 0 & metric == 'disease' ~ 'disease',
                                 coef < 0 & metric == 'disease' ~ 'healthy',
                                 
                                 coef > 0 & metric == 'antibiotic' ~ 'antibiotic',
                                 coef < 0 & metric == 'antibiotic' ~ 'untreated',
                                 
                                 coef > 0 & metric == 'time' ~ 'after',
                                 coef < 0 & metric == 'time' ~ 'before'),
         .keep = 'unused') %>%
  mutate(n = n_distinct(asv_id), .by = c(order, family)) %>%
  pivot_wider(names_from = association, 
              values_from = asv_id,
              values_fn = ~str_c(., collapse = '; ')) %>%
  arrange(-n)

asv_assocation_directions %>%
  filter(!is.na(disease))

asv_assocation_directions %>%
  filter(!is.na(disease) & !is.na(untreated)) %>%
  select(disease, untreated)

prePost_effects %>%
  filter(asv_id %in% c('ASV9', 'ASV93', 'ASV8', 'ASV25'))
  
  

#### Tables ####
prePost_effects %>%
  select(phylum:asv_id, 
         ndf, ddf, f.value, p.value, 
         all_of(str_remove(fdr_or_qvalue, '^_'))) %>%
  write_csv('../Results/TableS1_ASVmainEffects.csv')

prePost_effects %>%
  select(family:asv_id, 
         ndf, ddf, f.value, p.value, 
         all_of(str_remove(fdr_or_qvalue, '^_'))) %>%
  filter(!!sym(str_remove(fdr_or_qvalue, '^_')) < alpha_test) %>%
  arrange(family, str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  mutate(species = case_when(!is.na(species) ~ species,
                             is.na(species) & is.na(genus) ~ '',
                             TRUE ~ str_c(genus, ' sp.')),
         across(c(p.value, fdr.bh), ~scales::pvalue(.)),
         ddf = sprintf('%.1f', ddf),
         f.value = sprintf('%.2f', f.value)) %>%
  select(-genus) %>%
  relocate(asv_id, .before = 'family') %>%
  write_csv('../Results/Table1_ASVmainEffects_Sig.csv')


posthocTables <- prePost_effects %>%
  filter(if_any(all_of(str_remove(fdr_or_qvalue, '^_')), ~. < alpha_test)) %>%
  arrange(str_extract(asv_id, '[0-9]+') %>% as.integer()) %>%
  select(family:asv_id, contrast,
         contains(fdr_or_qvalue)) %>%
  rowwise %>%
  mutate(contrast_to_table(contrast), .keep = 'unused') %>%
  pivot_longer(cols = c(starts_with('disease'), starts_with('antibiotic'), starts_with('time')),
               names_to = c('association', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  relocate(all_of(str_remove(fdr_or_qvalue, '^_')), .after = p.value)
  

posthocTables %>%
  select(-family:-species) %>%
  write_csv('../Results/TableS2_ASVpostHocs.csv')

filter(posthocTables, fdr.bh < alpha_test) %>%
  mutate(species = case_when(!is.na(species) ~ species,
                             is.na(species) & is.na(genus) ~ '',
                             TRUE ~ str_c(genus, ' sp.')),
         across(c(p.value, fdr.bh), ~scales::pvalue(.)),
         df = sprintf('%.1f', df),
         t.value = sprintf('%.2f', t.value),
         estimate = str_c(sprintf('%.1f', estimate), sprintf('%.2f', std.error), sep = ' ± ')) %>%
  select(-genus, -std.error) %>%
  arrange(association) %>%
  relocate(association, .before = family) %>%
  write_csv('../Results/Table2_ASVpostHocsSig.csv')

#### Numbers for Text ####
filter(prePost_effects, asv_id %in% c('ASV8', 'ASV93')) %>%
  select(asv_id, posthoc) %>%
  rowwise(asv_id) %>%
  reframe(emmeans::contrast(posthoc,
                            method = list(anti = c(1/2, 0, 0, 1/2, 0),
                                          untreated = c(0, 0, 1/2, 0, 1/2),
                                          disease = c(0, 1, 0, 0, 0))) %>%
            broom::tidy())

filter(prePost_effects, asv_id %in% c('ASV25', 'ASV9')) %>%
  select(asv_id, posthoc) %>%
  rowwise(asv_id) %>%
  reframe(emmeans::contrast(posthoc,
                            method = list(untreated = c(0, 0, 1/2, 0, 1/2),
                                          disease = c(0, 1, 0, 0, 0))) %>%
            broom::tidy())

#### Just the Pathogens ####
prePost_effects %>%
  filter(asv_id %in% c('ASV25', 'ASV8')) %>% #, 
  # filter(asv_id %in% c('ASV846')) %>%
  select(asv_id, ends_with(fdr_or_qvalue), metrics) 

prePost_effects$model[[which(prePost_effects$asv_id == 'ASV8')]]

#
#### Make Upset Plot ####
colour_options <- read_rds('../intermediate_files/asv_colors.rds')
microbe_colors <- set_names(colour_options$color_palette$hex,
                            colour_options$color_palette$group)
levels(colour_options$asv_clumping$Top_order)


set_size_plot <- prePost_effects %>%
  select(asv_id, ends_with(fdr_or_qvalue),
         ends_with('coef')) %>%
  mutate(across(ends_with(fdr_or_qvalue), ~. < alpha_test),
         across(ends_with('coef'), ~if_else(. < 0, -1L, 1L))) %>%
  pivot_longer(cols = -asv_id,
               names_to = c('association', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(!!sym(str_remove(fdr_or_qvalue, '^_')) := replace_na(!!sym(str_remove(fdr_or_qvalue, '^_')), FALSE)) %>%
  filter(!!sym(str_remove(fdr_or_qvalue, '^_'))) %>%
  mutate(coef = if_else(coef == -1, 'Neg', 'Pos')) %>%
  # group_by(association) %>%
  # summarise(n_sig = n(),
  #           n_pos = sum(coef > 0)) %>%
  mutate(association = fct_infreq(association) %>% fct_rev()) %>%
  
  
  ggplot(aes(y = association)) +
  geom_bar(aes(fill = coef)) +
  geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = 'count',
            colour = 'white', fontface = 'bold', size = 4) +
  scale_fill_manual(values = c('Pos' = 'gray35', 'Neg' = 'gray65')) +
  labs(tag = 'D',
       x = 'Set Size') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.background = element_rect(colour = 'black'),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text.x = element_text(angle = 0,  colour = 'black', size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.tag = element_text(colour = 'black', size = 18,
                                face = 'bold', hjust = 0))



base_plot_data <- prePost_effects %>%
  select(order, genus, asv_id, ends_with(fdr_or_qvalue),
         ends_with('coef')) %>%
  mutate(across(ends_with(fdr_or_qvalue), ~. < alpha_test),
         across(ends_with('coef'), ~if_else(. < 0, -1L, 1L)),
         genus = if_else(is.na(genus), 'Other', genus),
         order = if_else(is.na(order), 'Other', order)) %>%
  
  mutate(Disease = !!sym(str_c('disease', fdr_or_qvalue)) & disease_coef > 0,
         Healthy = !!sym(str_c('disease', fdr_or_qvalue)) & disease_coef < 0,
         
         Untreated = !!sym(str_c('antibiotic', fdr_or_qvalue)) & antibiotic_coef < 0,
         Antibiotic = !!sym(str_c('antibiotic', fdr_or_qvalue)) & antibiotic_coef > 0,
         
         Postdose = !!sym(str_c('time', fdr_or_qvalue)) & time_coef > 0,
         Predose = !!sym(str_c('time', fdr_or_qvalue)) & time_coef < 0,
         .keep = 'unused') %>%
  
  mutate(Antibiotic = Untreated | Antibiotic,
         Disease = Disease | Healthy,
         Time = Predose | Postdose,
         .keep = 'unused') %>%
  
  left_join(mutate(colour_options$asv_clumping,
                   genus = str_remove_all(genus, '\\<i\\>|\\</i\\>')),
            by = c('order', 'genus')) %>%
  mutate(group = case_when(!is.na(group) ~ group,
                           is.na(group) & order == 'Oceanospirillales' ~ 'Oceanospirillales-Other',
                           is.na(group) & order == 'Flavobacteriales' ~ 'Flavobacteriales-Other',
                           is.na(group) & order == 'Alteromonadales' ~ 'Alteromonadales-Other',
                           is.na(group) & order == 'Verrucomicrobiales' ~ 'Verrucomicrobiales-Other',
                           is.na(group) & order == 'Rhodobacterales' ~ 'Rhodobacterales-Other',
                           TRUE ~ 'Other-Other'),
         group = factor(group, levels = levels(colour_options$asv_clumping$group))) %>%
  # filter(!is.na(group)) %>%
  select(-Top_order, -Top_genus) %>%
  rename(the_colour = group) %>%
  mutate(across(where(is.logical), ~replace_na(., FALSE)))

base_plot <- base_plot_data %>%
  
  filter(!if_all(where(is.logical), ~!.)) %>%
  
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
        
        # sort_sets = FALSE, 
        sort_sets = 'descending',
        set_sizes=(
          upset_set_size(geom = geom_bar(fill = 'gray65'),
                         position = 'right') + 
            geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = 'count',
                      colour = 'white', fontface = 'bold', size = 4) +
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
        
        
        name = NULL,
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(color = c('black'),
                                                                                                    size = 14, face = 'bold'),
                                                                         panel.background = element_rect(colour = 'black'),
                                                                         panel.border = element_rect(colour = 'black', 
                                                                                                     fill = 'transparent'),
                                                                         plot.tag = element_text(colour = 'black', size = 18, 
                                                                                                 face = 'bold', hjust = 0),
                                                                         axis.ticks = element_line(colour = 'black'),
                                                                         axis.line = element_line(colour = "black", 
                                                                                                  linewidth = rel(1))))),
        upset_stripes(colors = NA)
  )


base_plot

base_plot[[6]] <- set_size_plot
base_plot

ggdraw(base_plot) +
  draw_plot(colour_options$legend, 
            x = 0.75, y = 0.25, width = 0.25, height = 0.75)
ggsave('../Results/Fig5_asv_upset.png', 
       height = 12, width = 10, bg = 'white')
ggsave('../Results/Fig5.jpg', height = 12, width = 10, dpi = 'print', bg = 'white')

#### Plot by grouping ####
select(prePost_effects, name, ends_with(fdr_or_qvalue), metrics) %>%
  mutate(across(ends_with(fdr_or_qvalue), ~ . < alpha_test)) %>%
  group_by(across(ends_with(fdr_or_qvalue))) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(n)



# data <- tmp$metrics[[4]]; association <- tmp$association[[4]]
plot_asv <- function(data, association, significance = TRUE){
  association_classes <- str_extract_all(association, 
                                     '(Pre|Post)dose|Healthy|Disease|Antibiotic|Untreated') %>%
    unlist
  
  base_plot <- data %>%
    mutate(health = if_else(anti == 'A', 'A', health) %>% fct_relevel('D', after = Inf)) %>%
    ggplot(aes(x = time, y = mean_est, 
               ymin = conf.low,
               ymax = conf.high,
               fill = health,
               shape = anti)) +
    geom_errorbar(width = 0.1,
                  position = position_dodge(0.5),
                  show.legend = FALSE) +
    geom_point(position = position_dodge(0.5), 
               size = 2.5) +
    
    scale_fill_manual(values = set_names(c("#3A9AB2", "#80D1E9", "#F11B00"),
                                         c('H', 'A', 'D')),
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
         y = case_when(y_metric == 'rclr' ~ str_to_upper(y_metric),
                       TRUE ~ y_metric),
         fill = 'Health\nState',
         shape = 'Antibiotic\nTreatment') +
    theme_classic() +
    theme(strip.background = element_blank(),
          panel.background = element_rect(colour = 'black'),
          legend.key = element_blank())
  
  
  if(significance){
    max_y <- max(data$conf.high)
    n_sig <- 0
    
    if(any(str_detect(association_classes, '(Pre|Post)dose'))){
      y_height <- 1.25 * max_y; n_sig <- n_sig + 1
      base_plot <- base_plot +
        ggpubr::geom_bracket(xmin = 0.8, xmax = "after", 
                             y.position = y_height,
                             label = "*", label.size = 6,
                             vjust = 0.5,
                             inherit.aes = FALSE)
    }
    
    if(any(str_detect(association_classes, 'Antibiotic|Untreated'))){
      y_height <- if_else(n_sig > 0, 1.05 * max_y, 1.25 * max_y); n_sig <- n_sig + 1
      
      which_time_anti <- data %>%
        filter(health == 'H') %>%
        group_by(time) %>%
        summarise(t.value = diff(mean_est) / mean(se_est)) %>%
        mutate(p.value = 1 * pt(abs(t.value), Inf, lower.tail = FALSE)) %>%
        filter(p.value < alpha_test) %>%
        pull(time) %>%
        as.character()
      
      if('before' %in% which_time_anti){
        base_plot <- base_plot +
          ggpubr::geom_bracket(xmin = 0.8, xmax = 1.2, 
                               y.position = y_height,
                               label = "#", label.size = 3,
                               vjust = 0.25,
                               inherit.aes = FALSE) 
      }
      
      if('after' %in% which_time_anti){
        base_plot <- base_plot +
          ggpubr::geom_bracket(xmin = 1.8, xmax = 2.05, 
                               y.position = y_height,
                               label = "#", label.size = 3,
                               vjust = 0.25,
                               inherit.aes = FALSE)
      }
    }
    
    if(any(str_detect(association_classes, 'Healthy|Disease'))){
      y_height <- 1.25 * max_y; n_sig <- n_sig + 1
      
      base_plot <- base_plot +
        ggpubr::geom_bracket(xmin = 2.12, xmax = 2.21, 
                             y.position = y_height,
                             label = "+", label.size = 4,
                             vjust = 0.5, #tip.length = 4,
                             size = 0, colour = 'black',
                             inherit.aes = FALSE)
    }
  }
  
  
  base_plot
}

plot_fc <- function(fold_change){
  broom::tidy(fold_change, conf.int = TRUE) %>%
    mutate(contrast = str_to_sentence(contrast) %>%
             factor(levels = c('Time', 'Antibiotic', 'Disease'))) %>%
    ggplot(aes(x = contrast, y = estimate, 
               ymin = conf.low,
               ymax = conf.high)) +
    geom_hline(yintercept = 0) +
    geom_errorbar(width = 0.1) +
    geom_point() +
    geom_text(aes(y = Inf, label = if_else(p.value < 0.05, '*', '')),
              size = 10, vjust = 1) +
    labs(x = NULL, 
         y = "logFC") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text = element_text(colour = 'black'),
          panel.background = element_rect(colour = 'black'),
          legend.key = element_blank())
}


associated_asvs <- select(prePost_effects, domain:asv_id, name, 
              ends_with('coef'), ends_with(fdr_or_qvalue), 
              metrics, contrast) %>%
  filter(if_any(ends_with(fdr_or_qvalue), ~. < alpha_test)) %>%
  mutate(across(ends_with(fdr_or_qvalue), ~. < alpha_test),
         across(ends_with('coef'), ~if_else(. < 0, -1L, 1L))) %>%
  
  mutate(Disease = !!sym(str_c('disease', fdr_or_qvalue)) & disease_coef > 0,
         Healthy = !!sym(str_c('disease', fdr_or_qvalue)) & disease_coef < 0,
         
         Untreated = !!sym(str_c('antibiotic', fdr_or_qvalue)) & antibiotic_coef < 0,
         Antibiotic = !!sym(str_c('antibiotic', fdr_or_qvalue)) & antibiotic_coef > 0,
         
         Postdose = !!sym(str_c('time', fdr_or_qvalue)) & time_coef > 0,
         Predose = !!sym(str_c('time', fdr_or_qvalue)) & time_coef < 0,
         .keep = 'unused') %>%
  
  pivot_longer(cols = where(is.logical),
               names_to = 'association',
               values_to = 'significance') %>%
  filter(significance) %>%
  group_by(across(c(where(is.character), -association))) %>%
  summarise(association = str_c(sort(association), collapse = ' & '),
            n_associations = n(),
            metrics = unique(metrics),
            contrast = unique(contrast),
            .groups = 'drop')

associated_asvs %>%
  filter(str_detect(association, 'Disease')) %>%
  # arrange(n_associations) %>%
  arrange(family, asv_id) %>%
  rowwise %>%
  mutate(plot = list(plot_asv(metrics, association, significance = FALSE) + labs(subtitle = name)),
         plot_fc = list(plot_fc(contrast))) %>%
  select(asv_id, starts_with('plot')) %>%
  pivot_longer(cols = starts_with('plot'),
               names_to = 'type',
               values_to = 'plot') %>%
  pull(plot) %>%
  wrap_plots(ncol = 2, guides = 'collect', widths = c(0.8, 0.2)) +
  plot_annotation(tag_levels = list(c('A', '', 'B', '', 'C', '', 'D', ''))) +
  plot_layout(axes = 'collect', axis_titles = 'collect',
              guides = 'collect') &
  expand_limits(y = c(0, 3.7)) &
  theme(axis.text = element_text(colour = 'black'),
        # plot.subtitle = element_blank(),
        plot.subtitle = element_markdown(colour = 'black', size = 8))


associated_asvs %>%
  filter(str_detect(association, 'Disease')) %>%
  # arrange(n_associations) %>%
  arrange(family, asv_id) %>%
  rowwise %>%
  mutate(plot = list(plot_asv(metrics, association, significance = FALSE) + labs(subtitle = name)),
         plot_fc = list(plot_fc(contrast) +
                          scale_y_continuous(position = "right"))) %>%
  summarise(plot = list((plot + theme(plot.margin = margin(r = 0))) +
                          (plot_fc + theme(plot.margin = margin(l = 0))) +
                          plot_layout(widths = c(0.8, 0.2), tag_level = 'new'))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 1, guides = 'collect') +
  plot_annotation(tag_levels = c('A', '1'), tag_sep = '.') +
  plot_layout(axes = 'collect', axis_titles = 'collect',
              guides = 'collect') &
  # expand_limits(y = c(0, 3.7)) &
  theme(axis.text = element_text(colour = 'black'),
        # plot.subtitle = element_blank(),
        plot.subtitle = element_markdown(colour = 'black', size = 8))
# ggsave('../Results/Fig6_diseaseAssociated_ASVs_v1.png', height = 10, width = 6)

associated_asvs %>%
  filter(str_detect(association, 'Disease')) %>%
  # arrange(n_associations) %>%
  arrange(family, asv_id) %>%
  rowwise(name) %>%
  reframe(broom::tidy(contrast, conf.int = TRUE)) %>%
  mutate(contrast = str_to_sentence(contrast) %>%
           factor(levels = c('Time', 'Antibiotic', 'Disease'))) %>%
  ggplot(aes(x = contrast, y = estimate, 
             colour = name,
             ymin = conf.low,
             ymax = conf.high)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.25),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.25)) +
  geom_text(aes(label = if_else(p.value < 0.05, '*', ''), y = Inf),
            vjust = 1, position = position_dodge(0.25), size = 6,
            show.legend = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(y = 'logFC',
       x = NULL,
       colour = 'ASV') +
  theme_classic() +
  theme(legend.text = element_markdown(colour = 'black', size = 8))
# ggsave('../Results/Fig6_diseaseAssociated_ASVs_v2.png', height = 6, width = 6)

#### Field Models ####
field_associations <- sebastians_field %>%
  filter(asv_id %in% unique(associated_asvs$asv_id)) %>%
  group_by(asv_id) %>%
  filter(!all(rclr == 0)) %>%
  summarise(model = list(lm(rclr ~ health)),
            .groups = 'rowwise') %>%
  mutate(posthoc = list(emmeans::emmeans(model, ~health))) %>%
  reframe(anova(model) %>%
            as_tibble(rownames = 'term') %>%
            mutate(dDF = Df[term == 'Residuals']) %>%
            filter(term == 'health') %>%
            select(Df, dDF, `F value`, `Pr(>F)`),
          
          emmeans::contrast(posthoc, 'pairwise') %>%
            broom::tidy(conf.int = TRUE) %>%
            select(estimate, conf.low, conf.high, p.value) %>%
            mutate(contrast = 'Field')) %>%
  left_join(select(associated_asvs, asv_id, name),
            by = 'asv_id') %>%
  mutate(fdr = p.adjust(p.value, str_remove(fdr_or_qvalue, '^_fdr.') %>% str_to_upper()))

filter(field_associations, fdr < 0.05)

filter(field_associations, fdr >= 0.05)


#### All Contrasts ####
prePost_effects %>%
  filter(if_any(all_of(str_c(c('antibiotic', 'disease'), fdr_or_qvalue)), ~. < alpha_test)) %>%
  mutate(grouping = case_when(!!sym(str_c('antibiotic', fdr_or_qvalue)) < alpha_test & 
                                !!sym(str_c('disease', fdr_or_qvalue)) < alpha_test ~ 'Disease & Antibiotic',
                              !!sym(str_c('disease', fdr_or_qvalue)) < alpha_test ~ 'Disease',
                              !!sym(str_c('antibiotic', fdr_or_qvalue)) < alpha_test ~ 'Antibiotic')) %>%

  filter((grouping == 'Disease' & disease_coef > 0) |
           (grouping == 'Antibiotic' & antibiotic_coef < 0) |
           (grouping == 'Disease & Antibiotic' & antibiotic_coef < 0 & disease_coef > 0)) %>%
  
  rowwise(asv_id, name, grouping) %>%
  mutate(fdr = list(c_across(all_of(str_c(c('disease', 'antibiotic', 'time'), fdr_or_qvalue))))) %>%
  
  reframe(broom::tidy(contrast, conf.int = TRUE) %>%
            mutate(fdr = fdr)) %>%
  select(asv_id, name, grouping, contrast, 
         estimate, conf.low, conf.high, fdr) %>%
  # bind_rows(left_join(field_associations, 
  #                     distinct(., asv_id, grouping),
  #                     by = 'asv_id')) %>%
  mutate(contrast = str_to_title(contrast),
         grouping = factor(grouping, levels = c('Disease', 'Disease & Antibiotic', 'Antibiotic')),
         fill_value = if_else(fdr >= alpha_test, NA_character_, contrast)) %>%
  
  filter(!is.na(grouping)) %>%
  filter(contrast != 'Time') %>%
  
  mutate(fct_sort = case_when(grouping == 'Disease' ~ estimate[contrast == 'Disease'],
                              grouping == 'Antibiotic' ~ -estimate[contrast == 'Antibiotic'],
                              TRUE ~ (estimate[contrast == 'Disease'] + estimate[contrast == 'Antibiotic']) / 2),
       .by = 'asv_id') %>%
  mutate(name = fct_reorder(name, fct_sort)) %>%
  
  ggplot(aes(y = name, x = estimate, 
             colour = contrast,
             fill = fill_value,
             xmin = conf.low,
             xmax = conf.high)) +
  geom_vline(xintercept = 0) +
  geom_errorbar(width = 0.5,
                position = position_dodge(0.75),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.75), 
             shape = 'circle filled', size = 4) +
  scale_colour_manual(values = c('Time' = 'gray50',
                                 'Antibiotic' = "#80D1E9",
                                 'Disease' = "#F11B00",
                                 'Field' = 'forestgreen')) +
  scale_fill_manual(values = c('Time' = 'gray50',
                                 'Antibiotic' = "#80D1E9",
                                 'Disease' = "#F11B00",
                                 'Field' = 'forestgreen'), 
                    na.value = 'white',
                    breaks = c('Disease', 'Antibiotic'),
                    labels = c('Disease' = 'Significant', 'Antibiotic' = "Non-Significant")) +
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 'circle')),
         fill = guide_legend(override.aes = list(size = 4, 
                                                 shape = c('Significant' = 'circle', 
                                                           "Non-Significant" = 'circle open'), 
                                                 fill = 'black'))) +
  facet_grid(grouping ~ ., scales = 'free_y', space = 'free_y',
             labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
  labs(x = 'log<sub>2</sub>(FC)',
       y = NULL,
       colour = 'Contrast',
       fill = 'Significance') +
  theme_classic() +
  theme(axis.text.y = element_markdown(colour = 'black', size = 8),
        panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        legend.key = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_markdown(colour = 'black', size = 12))
ggsave('../Results/Fig6_diseaseAssociated_ASVs_v3.png', height = 6, width = 6)
ggsave('../Results/Fig6_r1.tiff', height = 6, width = 6, dpi = 'print')


#Aureispira question
prePost_effects %>%
  filter(str_detect(genus, 'Aureispira'))  %>%
  select(species, asv_id, fdr.bh)

prePost_effects %>%
  filter(str_detect(genus, 'Aureispira')) %>%
  rowwise(asv_id, name) %>%
  mutate(fdr = list(c_across(all_of(str_c(c('disease', 'antibiotic', 'time'), fdr_or_qvalue))))) %>%
  
  reframe(broom::tidy(contrast, conf.int = TRUE) %>%
            mutate(fdr = fdr)) %>%
  select(asv_id, name, contrast, 
         estimate, conf.low, conf.high, fdr) %>%
  ggplot(aes(y = name, x = estimate, 
             colour = contrast,
             # fill = fill_value,
             xmin = conf.low,
             xmax = conf.high)) +
  geom_vline(xintercept = 0) +
  geom_errorbar(width = 0.5,
                position = position_dodge(0.75),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.75), 
             shape = 'circle filled', size = 4) +
  scale_colour_manual(values = c('time' = 'gray50',
                                 'antibiotic' = "#80D1E9",
                                 'disease' = "#F11B00",
                                 'Field' = 'forestgreen')) +
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 'circle')),
         fill = guide_legend(override.aes = list(size = 4, 
                                                 shape = c('Significant' = 'circle', 
                                                           "Non-Significant" = 'circle open'), 
                                                 fill = 'black'))) 

#### FC vs Control (untreated & healthy) ####
tmp <- prePost_effects %>%
  filter(if_any(all_of(str_c(c('antibiotic', 'disease'), fdr_or_qvalue)), ~. < alpha_test)) %>%
  mutate(grouping = case_when(!!sym(str_c('antibiotic', fdr_or_qvalue)) < alpha_test & 
                                !!sym(str_c('disease', fdr_or_qvalue)) < alpha_test ~ 'Both',
                              !!sym(str_c('disease', fdr_or_qvalue)) < alpha_test ~ 'Disease',
                              !!sym(str_c('antibiotic', fdr_or_qvalue)) < alpha_test ~ 'Antibiotic')) %>%
  rowwise(asv_id, name, grouping) %>%
  mutate(new_contrast = list(emmeans::contrast(posthoc,
                                               method = list(before_anti = c(0, 0, 0, 1, -1),
                                                             after_anti = c(1, 0, -1, 0, 0),
                                                             after_disease = c(0, 1, -1, 0, 0))))) 

tmp %>%
  reframe(broom::tidy(new_contrast, conf.int = TRUE)) %>%
  select(asv_id, name, grouping, contrast, 
         estimate, conf.low, conf.high, p.value) %>%
  bind_rows(left_join(field_associations, 
                      distinct(., asv_id, grouping),
                      by = 'asv_id')) %>%
  mutate(time = str_extract(contrast, 'Field|before|after') %>% str_to_title %>%
           factor(levels = rev(c('Field', 'Before', 'After'))),
         antibiotic = if_else(str_detect(contrast, 'anti'), 'Antibiotic', 'Untreated'),
         disease = if_else(!str_detect(contrast, 'anti'), 'Disease', 'Healthy'),
         control_v = if_else(str_detect(contrast, 'anti'), 'Antibiotic', 'Disease'),
         signigicance_control_v = case_when(p.value > alpha_test ~ NA_character_,
                                            TRUE ~ control_v)) %>%
  # filter(asv_id == 'ASV8') %>%
  
  ggplot(aes(y = name, x = estimate, xmin = conf.low, xmax = conf.high, 
             shape = time, colour = control_v, fill = signigicance_control_v,
             group = time)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_errorbar(width = 0.1, position = position_dodge(0.5),
                show.legend = FALSE) +
  geom_point(position = position_dodge(0.5), size = 2) +
  scale_shape_manual(values = c('Field' = 'square filled', 
                                'Before' = 'triangle down filled', 
                                'After' = 'circle filled'),
                     breaks = c('Field', 'Before', 'After')) +
  scale_fill_manual(na.value = 'white', values = c('Antibiotic' = "#3A9AB2", 'Disease' = "#F11B00")) +
  scale_colour_manual(values = c('Antibiotic' = "#3A9AB2", 'Disease' = "#F11B00")) +
  guides(fill = 'none') +
  facet_grid(grouping ~ ., scales = 'free_y', space = 'free_y') +
  labs(x = 'logFC', y = NULL,
       colour = 'Contrast',
       shape = 'Time') +
  theme_classic() +
  theme(axis.text.y = element_markdown(colour = 'black', size = 8),
        panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        legend.key = element_blank())
ggsave('../Results/Figure6_vsControl.png', height = 10, width = 7)

#### Other Plots ####
associated_asvs %>%
  filter(str_detect(association, 'Healthy')) %>%
  arrange(n_associations) %>%
  rowwise %>%
  mutate(plot = list(plot_asv(metrics, association) + labs(subtitle = name))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 1, guides = 'collect') +
  plot_layout(axes = 'collect', axis_titles = 'collect',
              guides = 'collect') &
  theme(axis.text = element_text(colour = 'black'),
        plot.subtitle = element_markdown(colour = 'black', size = 8))
ggsave('../Results/healthyAssociated_ASVs.png', height = 10, width = 6)

associated_asvs %>%
  filter(str_detect(association, 'Untreated')) %>%
  arrange(n_associations) %>%
  rowwise %>%
  mutate(plot = list(plot_asv(metrics, association) + labs(subtitle = name))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 2) +
  plot_layout(axes = 'collect', axis_titles = 'collect',
              guides = 'collect') &
  theme(axis.text = element_text(colour = 'black'),
        plot.subtitle = element_markdown(colour = 'black', size = 8))
ggsave('../Results/untreatedAssociated_ASVs.png', height = 10, width = 6)

associated_asvs %>%
  filter(str_detect(association, 'Antibiotic')) %>%
  arrange(n_associations) %>%
  rowwise %>%
  mutate(plot = list(plot_asv(metrics, association) + labs(subtitle = name))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 1) +
  plot_layout(axes = 'collect', axis_titles = 'collect',
              guides = 'collect') &
  theme(axis.text = element_text(colour = 'black'),
        plot.subtitle = element_markdown(colour = 'black', size = 8))
ggsave('../Results/antiAssociated_ASVs.png', height = 10, width = 6)

associated_asvs %>%
  filter(str_detect(association, 'Predose')) %>%
  arrange(n_associations) %>%
  rowwise %>%
  mutate(plot = list(plot_asv(metrics, association) + labs(subtitle = name))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 1) +
  plot_layout(axes = 'collect', axis_titles = 'collect',
              guides = 'collect') &
  theme(axis.text = element_text(colour = 'black'),
        plot.subtitle = element_markdown(colour = 'black', size = 8))
ggsave('../Results/preAssociated_ASVs.png', height = 10, width = 6)

associated_asvs %>%
  filter(str_detect(association, 'Postdose')) %>%
  arrange(n_associations) %>%
  rowwise %>%
  mutate(plot = list(plot_asv(metrics, association) + labs(subtitle = name))) %>%
  pull(plot) %>%
  c(list(guide_area())) %>%
  wrap_plots(ncol = 2) +
  plot_layout(axes = 'collect', axis_titles = 'collect',
              guides = 'collect') &
  theme(axis.text = element_text(colour = 'black'),
        plot.subtitle = element_markdown(colour = 'black', size = 8))
ggsave('../Results/postAssociated_ASVs.png', height = 10, width = 6)




