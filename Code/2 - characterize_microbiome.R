
library(tidyverse)
library(microshades)
library(phyloseq)
library(speedyseq)
library(microshades)
library(cowplot)
library(ggtext)

#### Function ####
#Update microshade to allow specified list of genera
create_color_dfs_jds <- function(mdf, selected_groups = c("Proteobacteria", "Actinobacteria", 
                                                          "Bacteroidetes", "Firmicutes"), 
                                 desired_genera,
                                 top_n_subgroups = 4, group_level = "Phylum", 
                                 subgroup_level = "Genus", 
                                 cvd = FALSE, top_orientation = FALSE) 
{
  if (top_n_subgroups > 4) {
    stop("'top_n_subgroups' exceeds MAX value 4")
  }
  if (class(mdf) != "data.frame") {
    stop("mdf argument must be a data frame")
  }
  if (!is.null(mdf$group)) {
    stop("'group' column name already exists; consider renaming or removing")
  }
  if (is.null(mdf[[group_level]])) {
    stop("'group_level' does not exist")
  }
  if (is.null(mdf[[subgroup_level]])) {
    stop("'subgroup_level' does not exist")
  }
  
  
  #### Add logic to get desired genera as non-other
  in_top <- filter(desired_genera, !!sym(group_level) %in% selected_groups)
  in_other <- filter(desired_genera, !(!!sym(group_level) %in% selected_groups))
  
  #### Colour the main groups ####
  col_name_group <- paste0("Top_", group_level)
  mdf[[col_name_group]] <- "Other"
  rows_to_change <- mdf[[group_level]] %in% selected_groups
  taxa_names_mdf <- row.names(mdf[rows_to_change, ])
  mdf[taxa_names_mdf, col_name_group] <- as.character(mdf[taxa_names_mdf, 
                                                          group_level])
  mdf[[col_name_group]] <- factor(mdf[[col_name_group]], levels = c("Other", 
                                                                    selected_groups))
  if (sum(selected_groups %in% as.character(unique(mdf[[col_name_group]]))) != 
      length(selected_groups)) {
    stop("some 'selected_groups' do not exist in the dataset. Consider SILVA 138 c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes')")
  }
  mdf_unknown_subgroup <- mdf %>% mutate(`:=`(!!sym(subgroup_level), 
                                              fct_na_value_to_level(!!sym(subgroup_level), "Unknown")))
  col_name_subgroup <- paste0("Top_", subgroup_level)
  subgroup_ranks <- mdf_unknown_subgroup %>% group_by_at(c(paste(subgroup_level), 
                                                           paste(col_name_group))) %>% summarise(rank_abundance = sum(Abundance)) %>% 
    arrange(desc(rank_abundance)) %>% group_by_at(c(paste(col_name_group))) %>% 
    mutate(order = row_number()) %>% ungroup()
  subgroup_ranks[[col_name_subgroup]] <- "Other"
  
  z <- 0; row_to_change_list <- vector('list', length(selected_groups) + 1)
  for(big_group in c(selected_groups, 'Other')){
    z <- z + 1
    if(big_group == 'Other'){
      desired_genera_subset <- in_other
    } else {
      desired_genera_subset <- filter(desired_genera, !!sym(group_level) == big_group)
    }
    
    
    data_subset <- mutate(subgroup_ranks, row = row_number()) %>%
      filter(!!sym(col_name_group) == big_group)
    
    top_4_rows <- filter(data_subset, order <= top_n_subgroups) %>%
      pull(row)
    
    if(nrow(desired_genera_subset) > 0){
      desired_rows <- filter(data_subset, genus %in% desired_genera_subset$genus) %>%
        pull(row)
      
      n_desire <- length(desired_rows)
      if(n_desire < top_n_subgroups){
        bonus_rows <- sort(top_4_rows[!top_4_rows %in% desired_rows])[1:(top_n_subgroups - n_desire)]
        keep_rows <- sort(c(desired_rows, bonus_rows))
      } else {
        keep_rows <- sort(c(desired_rows))
      }
      
      
    } else {
      keep_rows <- top_4_rows
    }
    
    row_to_change_list[[z]] <- keep_rows
  }
  
  
  # rows_to_change <- subgroup_ranks$order <= top_n_subgroups
  rows_to_change <- unlist(row_to_change_list)
  
  subgroup_ranks[rows_to_change, col_name_subgroup] <- as.vector(subgroup_ranks[rows_to_change, 
                                                                                subgroup_level])
  group_info <- subgroup_ranks %>% mutate(group = paste(!!sym(col_name_group), 
                                                        !!sym(col_name_subgroup), sep = "-"))
  
  group_info <- arrange(group_info, !!sym(col_name_group), order) %>%
    mutate(is_other = !!sym(col_name_subgroup) == 'Other') %>%
    group_by(!!sym(col_name_group), is_other) %>%
    mutate(order = row_number()) %>%
    ungroup 
  
  
  group_info$order[group_info[[col_name_subgroup]] == "Other"] <- top_n_subgroups + 1
  
  group_info_to_merge <- group_info[, c(col_name_group, subgroup_level, 
                                        col_name_subgroup, "group")]
  mdf_group <- mdf_unknown_subgroup %>% left_join(group_info_to_merge, 
                                                  by = c(col_name_group, subgroup_level))
  prep_cdf <- group_info %>% 
    select(all_of(c("group", "order", col_name_group, col_name_subgroup))) %>% 
    arrange(!!sym(col_name_group), order)
  
  
  num_group_colors <- length(selected_groups) + 1
  hex_df <- default_hex(num_group_colors, cvd)
  cdf <- prep_cdf %>% group_by_at(c(paste(col_name_group))) %>% 
    tidyr::nest() %>% arrange(!!sym(col_name_group))
  if ("Other" %in% mdf[[col_name_group]]) {
    start <- 1
  }
  else {
    start <- 2
    num_group_colors <- num_group_colors - 1
  }
  for (i in 1:num_group_colors) {
    cdf$data[[i]]$hex <- hex_df[1:length(cdf$data[[i]]$group), 
                                start]
    start = start + 1
  }
  cdf <- cdf %>% ungroup() %>% arrange(desc(row_number())) %>% 
    tidyr::unnest(data) %>% select(!!sym(col_name_group), 
                                   !!sym(col_name_subgroup), group, hex) %>% mutate_all(as.character)
  cdf <- cdf %>% filter(!is.na(hex))
  if (top_orientation) {
    level_assign = unique(cdf$group)
  }
  else {
    level_assign = unique(rev(cdf$group))
  }
  mdf_group$group <- factor(mdf_group$group, levels = level_assign)
  list(mdf = mdf_group, cdf = cdf)
}


custom_legend_jds <- function(mdf, cdf, group_level = "Phylum", subgroup_level = "Genus", 
                              x = "Sample", y = "Abundance", legend_key_size = 0.4, legend_text_size = 10, 
                              legend_orientation = "vertical", compression_vect = NULL) 
{
  if (is.null(mdf[[group_level]])) {
    stop("mdf 'group_level' does not exist")
  }
  if (is.null(mdf[[subgroup_level]])) {
    stop("mdf 'subgroup_level' does not exist")
  }
  if (is.null(cdf$hex)) {
    stop("cdf 'hex' does not exist")
  }
  if (!(legend_orientation %in% c("vertical", "horizontal"))) {
    stop("legend orientation must be \"vertical\" or \"horizontal\"")
  }
  col_name_group <- paste0("Top_", group_level)
  col_name_subgroup <- paste0("Top_", subgroup_level)
  group_level_names <- unique(cdf[[col_name_group]])
  for (i in 1:length(group_level_names)) {
    if (i == 1) {
      complete_legend <- individual_legend_jds(mdf, cdf, group_level_names[i], 
                                               col_name_group, col_name_subgroup, legend_key_size = legend_key_size, 
                                               legend_text_size = legend_text_size)
    }
    else {
      new_legend <- individual_legend_jds(mdf, cdf, group_level_names[i], 
                                          col_name_group, col_name_subgroup, legend_key_size = legend_key_size, 
                                          legend_text_size = legend_text_size)
      complete_size <- i - 1
      new_size <- 1
      if (legend_orientation == "vertical") {
        complete_legend <- plot_grid(complete_legend, 
                                     new_legend, ncol = 1, 
                                     rel_heights = c(complete_size, new_size))
      }
      else if (legend_orientation == "horizontal") {
        complete_legend <- plot_grid(complete_legend, 
                                     new_legend, ncol = 2, 
                                     rel_widths = c(complete_size, new_size))
      }
    }
  }
  
  if(!all(is.null(compression_vect))){
    complete_legend <- plot_grid(complete_legend, NULL, ncol = 1, rel_heights = compression_vect)
  }
  complete_legend
  
}

individual_legend_jds <- function (mdf, cdf, group_name, col_name_group = "Top_Phylum", 
                                   col_name_subgroup = "Top_Genus", x = "Sample", y = "Abundance", 
                                   legend_key_size = 0.4, legend_text_size = 10) 
{
  select_mdf <- mdf %>% filter(!!sym(col_name_group) == group_name)
  select_cdf <- cdf %>% filter(!!sym(col_name_group) == group_name)
  select_plot <- ggplot(select_mdf, aes(x = .data[[x]], y = .data[[y]], 
                                        fill = .data[[col_name_subgroup]], text = .data[[col_name_subgroup]])) + 
    geom_col(position = "fill") + 
    scale_fill_manual(name = group_name, 
                      values = select_cdf$hex, 
                      breaks = select_cdf[[col_name_subgroup]]) + 
    theme(legend.justification = "left",
          legend.title = element_text(face = "bold"),
          legend.key.size = unit(legend_key_size, "lines"), 
          legend.key = element_rect(colour = 'black'),
          text = element_text(size = legend_text_size),
          legend.text = element_markdown())
  legend <- get_legend(select_plot)
}

plot_microshades_jds <- function(mdf_group, cdf, group_label = "Phylum Genus", x = "Sample", y = "Abundance"){
  if (class(mdf_group) != "data.frame") {
    stop("mdf_group argument must be a data frame")
  }
  if ((sum(is.na(cdf$hex)) > 0) || (sum(is.na(cdf$group)) > 
                                    0)) {
    stop("cdf does not contain complete color information - missing hex or group info")
  }
  plot <- mdf_group %>% 
    ggplot(aes(x = !!sym(x), y = !!sym(y), fill = group)) +
    # geom_col(colour = 'black', fill = 'white') +
    stat_summary(fun = sum, geom = "col",
                 colour = "black", fill = "white") +
    geom_col() +
    scale_fill_manual(name = group_label, values = cdf$hex, breaks = cdf$group) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot
}


# data <- microbiome_data
make_microbe_plot <- function(data, highest_taxon = 'order', lower_taxon = 'genus',
                              top_choices = NULL, sample_font_size, sample_vjust, ...){
  
  grouping_levels <- c('health', 'time', 'anti')
  
  sample_sizes <- sample_data(data) %>%
    as_tibble %>%
    mutate(time = time_fac) %>%
    summarise(n = n(),
              n_frag = n_distinct(fragment),
              .by = all_of(grouping_levels)) %>%
    # count(across(all_of(grouping_levels))) %>%
    mutate(health = str_replace_all(health, c('H' = 'Healthy', 'D' = 'Diseased')),
           time = str_replace_all(time, c('before' = 'Before', 'after' = 'After')),
           anti = str_replace_all(anti, c('A' = 'Antibiotics', 'N' = 'Untreated'))) %>%
    mutate(Sample = str_c(!!!syms(grouping_levels[c(3, 1)]), sep = '_'),
           Sample = fct_relevel(Sample, 'Untreated_Diseased', after = Inf),
           time = factor(time, levels = c('Before', 'After'))) %>%
    select(Sample, time, n, n_frag)
  
  update_data <- data %>%
    tax_glom(lower_taxon) %>%
    psmelt() %>% 
    filter(Abundance > 0) %>%
    as_tibble() %>% 
    
    filter(time_fac == 'before' | (time_fac == 'after' & exposure == 'D')) %>%
    
    mutate(health = str_replace_all(health, c('H' = 'Healthy', 'D' = 'Diseased')),
           time = str_replace_all(time_fac, c('before' = 'Before', 'after' = 'After')),
           anti = str_replace_all(anti, c('A' = 'Antibiotics', 'N' = 'Untreated'))) %>%
    group_by(health, time, anti) %>%
    mutate(total = sum(Abundance)) %>%
    ungroup() %>%
    group_by(!!!syms(grouping_levels), !!sym(lower_taxon)) %>%
    reframe(!!!syms(grouping_levels), 
            !!sym(highest_taxon), !!sym(lower_taxon), 
            total, rel_abun = sum(Abundance)/total) %>%
    distinct() %>%
    rename(Abundance = rel_abun) %>%
    mutate(Sample = str_c(!!!syms(grouping_levels[c(3, 1)]), sep = '_'),
           Sample = fct_relevel(Sample, 'Untreated_Diseased', after = Inf),
           Sample = factor(Sample, levels = c('Untreated_Healthy',
                                              'Antibiotics_Healthy',
                                              'Untreated_Diseased')),
           time = factor(time, levels = c('Before', 'After'))) %>%
    as.data.frame() %>%
    mutate(genus = str_c('<i>', genus, '</i>'))
  
  
  desired_genera <- c('Thalassotalea','Cysteiniphilum',
                      'Vibrio')
  
  
  if(is.null(top_choices)){
    top_choices <- group_by(update_data, !!sym(highest_taxon)) %>%
      summarise(count = sum(Abundance),
                n_low = n_distinct(!!sym(lower_taxon))) %>%
      arrange(-count) %>%
      filter(n_low > 3) %>%
      mutate(row_id = row_number()) %>%
      slice(1:5) %>%
      
      pull(!!sym(highest_taxon))
  }
  
  
  
  if(!is.null(desired_genera)){
    desired_taxa <- data %>%
      tax_glom(lower_taxon) %>%
      tax_table %>%
      as.data.frame() %>%
      as_tibble() %>%
      select(-species) %>%
      distinct %>%
      filter(genus %in% desired_genera) %>%
      select(order, genus) %>%
      
      
      left_join(group_by(update_data, !!sym(lower_taxon)) %>%
                  summarise(n_genus = sum(Abundance)),
                by = lower_taxon) %>%
      left_join(group_by(update_data, !!sym(highest_taxon)) %>%
                  summarise(n_order = sum(Abundance),
                            n_low = n_distinct(!!sym(lower_taxon))),
                by = highest_taxon) %>%
      arrange(-n_order, -n_genus) %>%
      mutate(genus = str_c('<i>', genus, '</i>'))
    
    
    color_objs_ordergenus <- create_color_dfs_jds(update_data, group_level = highest_taxon, 
                                                  desired_genera =  desired_taxa,
                                                  subgroup_level = lower_taxon,
                                                  selected_groups = top_choices, 
                                                  cvd = TRUE)
    
    
  } else {
    color_objs_ordergenus <- create_color_dfs(update_data, group_level = highest_taxon, 
                                              subgroup_level = lower_taxon,
                                              selected_groups = top_choices, 
                                              cvd = TRUE)
  }
  
  
  mdf_ordergenus <- color_objs_ordergenus$mdf
  cdf_ordergenus <- color_objs_ordergenus$cdf
  
  
  if(!is.null(desired_genera)){
    legend_ordergenus <- custom_legend_jds(mdf_ordergenus, cdf_ordergenus,
                                           group_level = highest_taxon,
                                           subgroup_level = lower_taxon, 
                                           ...)
  } else {
    legend_ordergenus <- custom_legend(mdf_ordergenus, cdf_ordergenus,
                                       group_level = highest_taxon,
                                       subgroup_level = lower_taxon, 
                                       ...)
  }
  
  plot_ordergenus_prelim <- plot_microshades_jds(mdf_ordergenus, cdf_ordergenus) + 
    geom_text(data = sample_sizes, 
              aes(x = Sample, y = Inf, label = n), 
              inherit.aes = FALSE, vjust = sample_vjust, size = sample_font_size) +
    scale_x_discrete(labels = ~str_replace_all(., '_', '\n')) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
    facet_grid(~time, scales = "free", space = "free",
               labeller = as_labeller(c('Before' = 'A', 'After' = 'B'))) +
    theme_bw() +
    theme(legend.position = "none", 
          plot.margin = margin(t = 20, r = 20, 
                               b = 6, l = 6),
          strip.placement = 'outside',
          axis.title.x = element_blank(),
          axis.text = element_text(colour = 'black', size = 12),
          strip.background = element_blank(),
          strip.text = element_text(colour = 'black', size = 14),
          axis.title.y = element_text(colour = 'black', size = 16),
          plot.background = element_rect(colour = 'white'))
  

  list(plot = plot_ordergenus_prelim,
       legend = legend_ordergenus,
       color_palette = cdf_ordergenus,
       asv_clumping = as_tibble(mdf_ordergenus) %>%
         select(genus, order, Top_order:group) %>%
         distinct)
  
  
}

# phylo <- microbiome_data
describe_phyloseq <- function(phylo){
  microbiome::tax_tibble(phylo, 'asv_id') %>%
  summarise(across(everything(), 
                   ~n_distinct(., na.rm = TRUE))) %>%
    relocate(asv_id, .after = species)
}

#### Data ####
microbiome_data <- read_rds('../intermediate_files/full_tank_microbiome.rds')

#### Descriptives ####
describe_phyloseq(microbiome_data)
sample_data(microbiome_data)

#### Make Plot ####
microbial_diversity <- microbiome_data %>%
  make_microbe_plot(highest_taxon = 'order', lower_taxon = 'genus',
                    legend_key_size = 1, legend_text_size = 14, 
                    sample_font_size = 5, sample_vjust = 2.5)

microbial_diversity$plot

plot_grid(microbial_diversity$plot + theme(strip.text = element_text(hjust = 0)),
          cowplot::plot_grid(NULL, microbial_diversity$legend, NULL, ncol = 1, rel_heights = c(0.12, 1, 0.12)), 
          rel_widths = c(1, .25))
ggsave('../Results/Fig2_overall_composition.png', height = 10, width = 10, bg = 'white')
ggsave('../Results/Fig2_r4.tiff', height = 10, width = 10, dpi = 'print', bg = 'white')



write_rds(microbial_diversity[-1], '../intermediate_files/asv_colors.rds')

#### Describe post rarity purge ####
analyzed_asvs <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
         show_col_types = FALSE) %>%
  pull(asv_id) %>%
  unique

subset_taxa(microbiome_data, otu_table(microbiome_data) %>% colnames %in% analyzed_asvs) %>%
  describe_phyloseq

nonRare_microbial_diversity <- subset_taxa(microbiome_data, 
                                           otu_table(microbiome_data) %>% 
                                             colnames %in% analyzed_asvs) %>%
  make_microbe_plot(highest_taxon = 'order', lower_taxon = 'genus',
                    legend_key_size = 1, legend_text_size = 14, 
                    top_choices = rev(unique(microbial_diversity$color_palette$Top_order)[-6]))


plot_grid(nonRare_microbial_diversity$plot,
          cowplot::plot_grid(NULL, nonRare_microbial_diversity$legend, NULL, ncol = 1, rel_heights = c(0.12, 1, 0.12)), 
          rel_widths = c(1, .25))
ggsave('../Results/nonRare_composition.png', height = 10, width = 10, bg = 'white')


#### Post Variance Purge ####
analyzed_asvs <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
         show_col_types = FALSE)  %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac) %>%
  filter(exposure == 'D') %>%
  mutate(treatment = str_c(anti, health, sep = '_')) %>% 
  
  
  nest_by(across(domain:species), asv_id) %>%
  
  #Remove samples which will break tests
  filter(!summarise(data,
                    var = sd(rclr) < 1e-12,
                    .by = c('time', 'treatment')) %>%
           pull(var) %>%
           any) %>%
  pull(asv_id) %>%
  unique

subset_taxa(microbiome_data, otu_table(microbiome_data) %>% colnames %in% analyzed_asvs) %>%
  describe_phyloseq



nonRare_variant_microbial_diversity <- subset_taxa(microbiome_data, 
                                           otu_table(microbiome_data) %>% 
                                             colnames %in% analyzed_asvs) %>%
  make_microbe_plot(highest_taxon = 'order', lower_taxon = 'genus',
                    legend_key_size = 1, legend_text_size = 14, 
                    top_choices = rev(unique(microbial_diversity$color_palette$Top_order)[-6]))


plot_grid(nonRare_variant_microbial_diversity$plot,
          cowplot::plot_grid(NULL, nonRare_variant_microbial_diversity$legend, NULL, 
                             ncol = 1, rel_heights = c(0.12, 1, 0.12)), 
          rel_widths = c(1, .25))
ggsave('../Results/nonRare_variant_composition.png', height = 10, width = 10, bg = 'white')
