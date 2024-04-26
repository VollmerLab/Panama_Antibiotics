#### Libraries ####
library(tidyverse)
library(microbiome)
library(phyloseq)
library(broom)
library(biostat)
library(qvalue)
library(ComplexUpset)
library(cowplot)
library(patchwork)

alpha <- 0.05
fdr_or_qvalue <- 'fdr'

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
tank_data <- read_rds("../intermediate_files/full_tank_microbiome.rds")

#### Calculate alpha metrics ####
alpha_diversity <- microbiome::tax_tibble(tank_data, 'asv') %>%
  colnames %>%
  tibble(taxonomic_level = .) %>%
  filter(!taxonomic_level %in% c('domain', 'asv')) %>%
  rowwise() %>%
  mutate(data = list(aggregate_taxa(tank_data, taxonomic_level))) %>%
  add_row(taxonomic_level = 'asv', data = list(tank_data)) %>%
  rowwise(taxonomic_level) %>%
  
  #Calculate alpha diversity
  reframe(alpha(data, zeroes = FALSE) %>%
            as_tibble(rownames = 'sample_id')) %>%
  rename(richness_observed = observed,
         richness_chao1 = chao1) %>%
  pivot_longer(cols = -c(taxonomic_level, sample_id),
               names_to = 'alpha_metric') %>%
  # filter(alpha_metric %in% c('richness_observed', 'diversity_shannon',
  #                            'evenness_pielou', 'dominance_core_abundance',
  #                            'rarity_rare_abundance')) %>%
  
  left_join(sample_data(tank_data) %>%
              as_tibble(rownames = 'sample_id'),
            by = 'sample_id') %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac)


#### Run analyses ####
prePost_effects <- alpha_diversity %>%
  filter(exposure == 'D') %>%
  mutate(treatment = str_c(anti, health, sep = '_')) %>%
  nest_by(taxonomic_level, alpha_metric) %>%
  filter(!alpha_metric %in% c('diversity_coverage')) %>%
  
  mutate(predose_data = list(filter(data, time == 'before')),
         postdose_data = list(filter(data, time == 'after'))) %>%
  
  #Test for equal variance in both pre and post data and use appropriate test
  #Test for normality in both pre and post data and use appropriate test - TODO
  
  reframe(var.test(value ~ anti, data = predose_data) %>%
            tidy %>%
            select(p.value) %>%
            rename(predoseVAR_p.value = p.value),
          
          car::leveneTest(value ~ anti, 
                          data = postdose_data) %>%
            tidy() %>%
            select(p.value) %>%
            rename(postdoseVAR_p.value = p.value),
          
          predose_test = list(with(predose_data,
                                   t.test(value ~ anti,
                                          var.equal = predoseVAR_p.value > alpha))),
          
          postdose_test = list(oneway.test(value ~ treatment, 
                                           data = postdose_data, 
                                           var.equal = postdoseVAR_p.value > alpha)),

          
          #Plotting
          metrics = list(summarise(data, 
                                   mean_est = mean(value), 
                                   se_est = sd(value) / sqrt(n()),
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
          posthoc_anova(value ~ treatment, 
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
  
  
  group_by(taxonomic_level) %>%
  mutate(across(c(ends_with('p.value'), 
                  -matches('[DHNA]v[DHNA]'), 
                  -contains('VAR')),
                list(fdr = ~p.adjust(., 'fdr'), 
                     q.value = ~safe_qvalue(.)))) %>%
  rename_with(~str_replace_all(., '_p.value_', '_')) %>%
  # select(asv_id, postdose_fdr, matches('[DHNA]+v[DHNA]+_p.value')) %>%
  
  group_by(taxonomic_level, post_sig_fdr = postdose_fdr < alpha) %>%
  mutate(across(matches('[DHNA]+v[DHNA]+_p.value'), 
                list(fdr = ~p.adjust(., 'fdr')))) %>%
  ungroup %>%
  select(-post_sig_fdr) %>%
  rename_with(~str_replace_all(., '_p.value_', '_')) %>%
  
  mutate(alpha_type = str_extract(alpha_metric, '[a-z]+') %>% str_to_sentence(),
         alpha_metric = str_remove(alpha_metric, '^[a-z]+_') %>% 
           str_replace_all('_', ' ') %>%
           str_to_title() %>% str_c(alpha_type, ': ', .),
         alpha_metric = str_replace_all(alpha_metric, c('Dbp' = 'DBP', 
                                                        'Dmn' = 'DMN', 
                                                        'Gini' = 'GINI')),
         taxonomic_level = case_when(taxonomic_level == 'asv' ~ 'ASV',
                                     TRUE ~ str_to_sentence(taxonomic_level)),
         taxonomic_level = factor(taxonomic_level, levels = c('Phylum', 'Class', 
                                                              'Order', 'Family',
                                                              'Genus', 'Species', 'ASV')))

prePost_effects %>%
  filter(predose_fdr < alpha | postdose_fdr < alpha)

#### Upset Plot ####
upset_data <- prePost_effects %>%
  select(taxonomic_level, alpha_type, alpha_metric, ends_with(fdr_or_qvalue)) %>%
  mutate(across(ends_with(fdr_or_qvalue), ~. < alpha)) %>%
  select(taxonomic_level, alpha_type, alpha_metric, contains('v'), postdose_fdr, predose_fdr) %>%
  select(-matches('[NDAH]v[NDAH]'))

metric_colours <- ggnested::nested_palette(upset_data,
                                           group = "alpha_type",
                                           subgroup = "alpha_metric")


upset_plot <- upset_data %>%
  
  upset(data = ., 
        intersect = select(., where(is.logical)) %>% colnames, 
        sort_sets = FALSE,
        guides = 'over',
        
        annotations = list(
          'taxonomic_level' = ggplot(mapping = aes(fill = taxonomic_level)) +
            geom_bar(fill = 'white', colour = 'black', stat = 'count', 
                     position = 'fill', show.legend = FALSE) +
            geom_bar(stat = 'count', position = 'fill', show.legend = TRUE) +
            scale_y_continuous(labels=scales::percent_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            # scale_fill_manual(values = microbe_colors) +
            labs(y = 'Relative Abundance',
                 tag = 'A',
                 fill = 'Taxonomic\nLevel') +
            theme_classic() +
            theme(legend.position = 'left',
                  panel.background = element_rect(colour = 'black'),
                  axis.ticks.y = element_line(colour = 'black'),
                  axis.minor.ticks.y.left = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title = element_text(colour = 'black', size = 14),
                  axis.text = element_text(colour = 'black', size = 10),
                  plot.tag = element_text(colour = 'black', size = 18, 
                                          face = 'bold', hjust = 0)),
          
          'alpha_metric' = ggplot(mapping = aes(fill = alpha_metric)) +
            geom_bar(fill = 'white', colour = 'black', stat = 'count', 
                     position = 'fill', show.legend = FALSE) +
            geom_bar(stat = 'count', position = 'fill', show.legend = TRUE) +
            scale_y_continuous(labels=scales::percent_format(), 
                               expand = expansion(mult = c(0.01, 0.05))) +
            scale_fill_manual(values = with(metric_colours, set_names(subgroup_colour, alpha_metric))) +
            labs(y = 'Relative Abundance',
                 tag = 'B',
                 fill = 'Alpha\nDiversity') +
            theme_classic() +
            theme(legend.position = 'left',
              panel.background = element_rect(colour = 'black'),
              axis.ticks.y = element_line(colour = 'black'),
              axis.minor.ticks.y.left = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title = element_text(colour = 'black', size = 14),
              axis.text = element_text(colour = 'black', size = 10),
              plot.tag = element_text(colour = 'black', size = 18, 
                                      face = 'bold', hjust = 0))))


ggsave('../Results/alpha_diversity_metic_associations.png',
       plot = upset_plot,
       height = 12, width = 15, bg = 'white')



upset_data %>%
  select(-matches('[NDAH]v[NDAH]'), -alpha_metric) %>%
  mutate(type_by_level = str_c(taxonomic_level, alpha_type, sep = '_')) %>%
  upset_test_jds(intersect = select(., where(is.logical)) %>% colnames, 
             tests = list(taxonomic_level = chisq_from_formula,
                          # alpha_metric = chisq_from_formula,
                          alpha_type = chisq_from_formula,
                          type_by_level = chisq_from_formula))

#### Individual Metrics ####
plot_alpha <- function(data, metric){
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
                mutate(significance = case_when(time == 'after' & postdose_fdr < alpha ~ '*',
                                                time == 'before' & predose_fdr < alpha ~ '*',
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
    facet_wrap( ~ taxonomic_level, scales = 'free_y', nrow = 1) +
    labs(x = NULL, 
         y = metric,
         fill = 'Health\nState',
         shape = 'Antibiotic\nTreatment') +
    theme_classic() +
    theme(strip.background = element_blank(),
          panel.background = element_rect(colour = 'black'),
          legend.key = element_blank())
}

tmp <- prePost_effects %>% 
  select(taxonomic_level, alpha_type, alpha_metric, ends_with('fdr'), metrics) %>%
  unnest(metrics) %>%
  nest(data = -c(alpha_type, alpha_metric)) %>%
  mutate(row_id = row_number(),
         alpha_metric = str_replace(alpha_metric, ': ', '\n')) %>%
  rowwise %>%
  mutate(plot = list(plot_alpha(data, alpha_metric)),
         plot = case_when(row_id %% 2 == 0 ~ list(plot + scale_y_continuous(position = "right")),
                          TRUE ~ list(plot)),
         plot = case_when(row_id == 1 ~ list(plot + theme(axis.text.x = element_blank())),
                          row_id == 21 ~ list(plot + theme(strip.text = element_blank())),
                          TRUE ~ list(plot + 
                                        theme(strip.text = element_blank(),
                                              axis.text.x = element_blank())))) %>%
  ungroup


wrap_plots(tmp$plot) +
  plot_layout(ncol = 1, guides = 'collect')
  
ggsave('../Results/alpha_changing_postExposure.png', height = 15, width = 15)


