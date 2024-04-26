#### Libraries ####
library(tidyverse)
library(broom)
library(biostat)
library(qvalue)
library(ComplexUpset)
library(cowplot)
library(patchwork)

alpha <- 0.05
fdr_or_qvalue <- 'fdr'
refit_models <- TRUE

#### Functions ####
safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_)

chisq_from_formula = function(formula, data) {
  chisq.test(
    ftable(formula, data)
  )
}

upset_test_jds <- function (data, intersect, ...){
  comparison = do.call(compare_between_intersections_jds, c(list(data, 
                                                                 intersect), list(...)))
  if (nrow(comparison) == 0) {
    stop("No variables to compare")
  }
  significant_hits = comparison[comparison$fdr < 0.05, "variable"]
  if (length(significant_hits)) {
    print(paste(paste(significant_hits, collapse = ", "), 
                "differ significantly between intersections"))
  }
  result = comparison[order(comparison$fdr), ]
  result
}

compare_between_intersections_jds <- function (data, intersect, test = kruskal.test, tests = list(), 
                                               ignore = list(), ignore_mode_columns = TRUE, mode = "exclusive_intersection", 
                                               ...) 
{
  data = upset_data(data, intersect, mode = mode, ...)
  isect = data$with_sizes
  isect = isect[(isect[, paste0("in_", mode)] == 1) & (isect$intersection %in% 
                                                         data$plot_intersections_subset), ]
  modes = c("exclusive_intersection", "inclusive_intersection", 
            "exclusive_union", "inclusive_union")
  if (ignore_mode_columns) {
    ignore = c(ignore, paste0("in_", modes), paste0(modes, 
                                                    "_size"), "exclusive_intersection")
  }
  ignore = c("intersection", ignore)
  candidate_variables = setdiff(setdiff(names(isect), intersect), 
                                ignore)
  results = list()
  for (variable in candidate_variables) {
    if (variable %in% names(tests)) {
      var_test = tests[[variable]]
    }
    else {
      var_test = test
    }
    isect$intersection = as.factor(isect$intersection)
    result = var_test(as.formula(paste(variable, "~ intersection")), 
                      data = isect)
    results[[length(results) + 1]] = list(variable = variable, 
                                          p.value = result$p.value, statistic = result$statistic, 
                                          parameter = result$parameter,
                                          test = result$method)
  }
  if (length(results) == 0) {
    headers = c("variable", "p.value", "statistic", "parameter", "test", 
                "fdr")
    results = data.frame(matrix(ncol = length(headers), 
                                nrow = 0, dimnames = list(NULL, headers)))
  }
  else {
    results = do.call(rbind, lapply(results, data.frame))
    results$fdr = p.adjust(results$p.value, "BH")
    rownames(results) = results$variable
  }
  results
}

#### Data ####
tank_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
                      show_col_types = FALSE) %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac)

#### Run analyses ####
prePost_effects <- tank_data %>%
  filter(exposure == 'D') %>%
  mutate(treatment = str_c(anti, health, sep = '_')) %>%
  
  #add a tiny bit of noise to data
  # mutate(log2_cpm_norm = rnorm(nrow(.), 0, 1e-8)) %>%
  
  nest_by(across(domain:species), asv_id) %>%
  
  #Remove samples which will break tests
  filter(!summarise(data,
                   var = sd(log2_cpm_norm) < 1e-12,
                   .by = c('time', 'treatment')) %>%
           pull(var) %>%
           any) %>%
  
  mutate(predose_data = list(filter(data, time == 'before')),
         postdose_data = list(filter(data, time == 'after'))) %>%
  
  #Test for equal variance in both pre and post data and use appropriate test
  #Test for normality in both pre and post data and use appropriate test - TODO
  
  reframe(var.test(log2_cpm_norm ~ anti, data = predose_data) %>%
            tidy %>%
            select(p.value) %>%
            rename(predoseVAR_p.value = p.value),
          
          car::leveneTest(log2_cpm_norm ~ anti, 
                          data = postdose_data) %>%
            tidy() %>%
            select(p.value) %>%
            rename(postdoseVAR_p.value = p.value),
          
          predose_test = list(with(predose_data, 
                                   t.test(log2_cpm_norm ~ anti, 
                                          var.equal = predoseVAR_p.value > alpha))),
        
          postdose_test = list(oneway.test(log2_cpm_norm ~ treatment, 
                                           data = postdose_data, 
                                           var.equal = postdoseVAR_p.value > alpha)),
          
          #Plotting
          metrics = list(summarise(data, 
                              mean_est = mean(log2_cpm_norm), 
                              se_est = sd(log2_cpm_norm) / sqrt(n()),
                              .by = c('time', 'anti', 'health'))),
          
          tidy(predose_test) %>%
            select(statistic, parameter, p.value) %>%
            rename(predose_t.value = statistic,
                   predose_df = parameter,
                   predose_p.value = p.value),
          
          tidy(postdose_test) %>%
            select(-method) %>%
            rename(postdose_num.df = num.df,
                   postdose_den.df = den.df,
                   postdose_f.value = statistic,
                   postdose_p.value = p.value),
          
          #Posthocs
          posthoc_anova(log2_cpm_norm ~ treatment, 
                        data = postdose_data, 
                        method = if_else(postdoseVAR_p.value > alpha, 
                                         'tukey', 'games-howell'))$output$result %>%
            select(groups, difference, p) %>%
            mutate(groups = str_remove_all(groups, '_') %>% str_replace_all('-','v')) %>%
            rename(coef = difference,
                   p.value = p) %>%
            pivot_wider(names_from = groups, 
                        values_from = c(coef, p.value), 
                        names_glue = "{groups}_{.value}")) %>%
  
  mutate(across(c(ends_with('p.value'), 
                  -matches('[DHNA]v[DHNA]'), 
                  -contains('VAR')),
                list(fdr = ~p.adjust(., 'fdr'), 
                     q.value = ~safe_qvalue(.)))) %>%
  rename_with(~str_replace_all(., '_p.value_', '_')) %>%
  # select(asv_id, postdose_fdr, matches('[DHNA]+v[DHNA]+_p.value')) %>%
  
  group_by(post_sig_fdr = postdose_fdr < alpha) %>%
  mutate(across(matches('[DHNA]+v[DHNA]+_p.value'), 
                list(fdr = ~p.adjust(., 'fdr')))) %>%
  ungroup %>%
  select(-post_sig_fdr) %>%
  rename_with(~str_replace_all(., '_p.value_', '_')) %>%
  
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
                          TRUE ~ asv_id)) 



#### Make Upset Plot ####
colour_options <- read_rds('../intermediate_files/asv_colors.rds')
microbe_colors <- set_names(colour_options$color_palette$hex,
                            colour_options$color_palette$group)
levels(colour_options$asv_clumping$Top_order)


base_plot_data <- prePost_effects %>%
  select(order, genus, asv_id, ends_with(fdr_or_qvalue)) %>%
  mutate(across(ends_with(fdr_or_qvalue), ~. < alpha),
         genus = if_else(is.na(genus), 'Other', genus),
         order = if_else(is.na(order), 'Other', order)) %>%
  # mutate(across(ends_with(fdr_or_qvalue), ~if_else(!postdose_fdr, FALSE, .))) %>%
  
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
  rename(the_colour = group) 


base_plot <- base_plot_data %>%
  select(order:asv_id, the_colour, contains('v'), postdose_fdr, predose_fdr) %>%
  select(-contains('v')) %>%
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
          upset_set_size(geom = geom_bar(fill = c('gray65')),
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

ggdraw(base_plot) +
  draw_plot(colour_options$legend, 
            x = 0.75, y = 0.25, width = 0.25, height = 0.75)
ggsave('../Results/asv_upset.png', 
       height = 12, width = 10, bg = 'white')



prePost_effects %>%
  select(phylum:genus, asv_id, ends_with(fdr_or_qvalue)) %>%
  select(-matches('[NDAH]v[NDAH]')) %>%
  mutate(significnce_group = case_when(predose_fdr < alpha & postdose_fdr < alpha ~ 'both',
                                       predose_fdr < alpha ~ 'pre',
                                       postdose_fdr < alpha ~ 'post',
                                       TRUE ~ 'neither'),
         .keep = 'unused') %>%
  pivot_longer(cols = phylum:genus,
               names_to = 'taxonomic_value',
               values_to = 'name', 
               values_drop_na = TRUE) %>%
  group_by(significnce_group, taxonomic_value, name) %>%
  summarise(n = n_distinct(asv_id), .groups = 'drop') %>%
  pivot_wider(names_from = significnce_group, 
              values_from = n,
              values_fill = 0L) %>%
  
  select(-neither) %>%
  filter(!if_all(where(is.numeric), ~. == 0)) %>%
  
  nest_by(taxonomic_value) %>%
  mutate(column_to_rownames(data, 'name') %>%
           chisq.test(simulate.p.value = TRUE) %>%
           tidy)

#### Just the Pathogens ####
prePost_effects %>%
  filter(asv_id %in% c('ASV25', 'ASV8')) %>%
  # filter(asv_id %in% c('ASV846')) %>%
  select(asv_id, ends_with('fdr'), metrics) %>%
  unnest(metrics) %>%

  ggplot(aes(x = time, y = mean_est, 
             ymin = mean_est - se_est,
             ymax = mean_est + se_est,
             shape = anti,
             colour = health)) +
  geom_errorbar(width = 0.1,
                position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), 
             size = 4) +
  
  geom_text(data = . %>%
              filter((time == 'after' & health == 'D') |
                       (time == 'before' & anti == 'A')) %>%
              mutate(significance = case_when(time == 'after' & postdose_fdr < alpha ~ '*',
                                              time == 'before' & predose_fdr < alpha ~ '*',
                                              TRUE ~ '')),
            aes(y = Inf, label = significance),
            colour = 'black', vjust = 1, size = 8) +
  
  facet_wrap(~asv_id) +
  labs(x = NULL, 
       y = 'log2CPM')

#### Plot by grouping ####
select(prePost_effects, name, ends_with('fdr'), metrics) %>%
  mutate(across(ends_with('fdr'), ~ . < alpha)) %>%
  count(predose_fdr, postdose_fdr, NDvAH_fdr, NHvAH_fdr, NHvND_fdr) %>%
  filter(postdose_fdr) %>%
  arrange(n)

plot_asv <- function(data, main, sub){
  if(is.na(sub)){
    sub <- NULL
  }
  
  data %>%
    ggplot(aes(x = time, y = mean_est, 
               ymin = mean_est - se_est,
               ymax = mean_est + se_est,
               shape = anti,
               fill = health)) +
    geom_errorbar(width = 0.1,
                  position = position_dodge(0.5),
                  show.legend = FALSE) +
    geom_point(position = position_dodge(0.5), 
               size = 2) +
    
    geom_text(data = . %>%
                filter((time == 'after' & health == 'D') |
                         (time == 'before' & anti == 'A')) %>%
                mutate(significance = case_when(time == 'after' & postdose_fdr ~ '*',
                                                time == 'before' & predose_fdr ~ '*',
                                                TRUE ~ '')),
              aes(y = Inf, label = significance),
              colour = 'black', vjust = 1, size = 8) +
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
    facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
    labs(x = NULL, 
         y = 'log2CPM',
         fill = 'Health\nState',
         shape = 'Antibiotic\nTreatment',
         title = main,
         subtitle = sub) +
    theme_classic() +
    theme(strip.background = element_blank(),
          panel.background = element_rect(colour = 'black'),
          legend.key = element_blank())
}

tmp <- select(prePost_effects, name, ends_with('fdr'), metrics) %>%
  filter(if_any(ends_with('fdr'), ~. < alpha)) %>%
  mutate(across(ends_with('fdr'), ~. < alpha)) %>%
  mutate(main_group = case_when(predose_fdr & postdose_fdr ~ 'both',
                               predose_fdr ~ 'pre',
                               postdose_fdr ~ 'post',
                               TRUE ~ 'none'),
         sub_group = case_when(main_group %in% c('both', 'post') & NDvAH_fdr & NHvAH_fdr & NHvND_fdr ~ 
                                 'Disease vs Antibiotic and Healthy vs Disease and Healthy vs Antibiotic',
                               
                               main_group %in% c('both', 'post') & NDvAH_fdr & NHvND_fdr ~ 'Disease vs Antibiotic and Healthy vs Disease',
                               main_group %in% c('both', 'post') & NDvAH_fdr & NHvAH_fdr ~ 'Disease vs Antibiotic and Healthy vs Antibiotic',
                               main_group %in% c('both', 'post') & NHvND_fdr & NHvAH_fdr ~ 'Healthy vs Disease and Healthy vs Antibiotic',
                               
                               main_group %in% c('both', 'post') & NDvAH_fdr ~ 'Disease vs Antibiotic',
                               main_group %in% c('both', 'post') & NHvAH_fdr ~ 'Healthy vs Antibiotic',
                               main_group %in% c('both', 'post') & NHvND_fdr ~ 'Healthy vs Disease',
                               
                               TRUE ~ NA_character_),
         main_group = case_when(main_group == 'both' ~ 'PreDose & PostDose',
                                main_group == 'pre' ~ 'PreDose',
                                main_group == 'post' ~ 'PostDose'),
         main_group = factor(main_group, levels = c('PreDose', 'PostDose', 'PreDose & PostDose')),
         sub_group = fct_relevel(sub_group, 'Disease vs Antibiotic and Healthy vs Antibiotic', after = Inf),
         sub_group = fct_relevel(sub_group, 'Healthy vs Antibiotic', after = 0L)) %>%
  unnest(metrics) %>% 
  nest_by(main_group, sub_group) %>%
  mutate(plot = list(plot_asv(data, main_group, sub_group))) %>%
  ungroup %>%
  arrange(main_group, sub_group)


wrap_plots(tmp$plot, ncol = 1) +
  plot_layout(guides = 'collect')
ggsave('../Results/asvs_changing_postExposure.png', height = 20, width = 17)


prePost_effects %>%
  filter(str_detect(name, 'Keto')) %>%
  select(asv_id, ends_with('fdr'))

prePost_effects %>%
  select(-where(is.list)) %>%
  write_csv('../Results/individual_asv_results.csv')

