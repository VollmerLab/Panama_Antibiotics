#### Libraries ####
library(tidyverse)
library(vegan)
library(patchwork)

#### Data ####
tank_data_cpm <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                           show_col_types = FALSE) %>%
  filter(exposure == 'D') %>%
  select(asv_id, sample_id, log2_cpm_norm) %>%
  pivot_wider(names_from = asv_id, values_from = log2_cpm_norm)

tank_metadata <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                          show_col_types = FALSE) %>%
  filter(exposure == 'D') %>%
  select(sample_id, plate:norm.factors) %>%
  distinct %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac)

tank_metadata %>%
  filter(geno == 'BL') %>%
  filter(str_detect(sample_id, 'N_H'))


if(!file.exists('../intermediate_files/taxonomy.csv.gz')){
  taxonomy <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv', 
                       show_col_types = FALSE) %>%
    select(asv_id, domain:species) %>%
    distinct %>%
    mutate(across(everything(), str_replace_na)) %>%
    write_csv('../intermediate_files/taxonomy.csv.gz')
} else {
  taxonomy <- read_csv('../intermediate_files/taxonomy.csv.gz', 
                       show_col_types = FALSE) %>%
    mutate(across(everything(), str_replace_na))
}

#### PERMANOVA ####
if(file.exists('../intermediate_files/adonis.rds.gz')){
  the_adonis <- read_rds('../intermediate_files/adonis.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  adon_y_data <- column_to_rownames(tank_data_cpm, 'sample_id')
  
  the_adonis <- adonis2(vegdist(adon_y_data, binary = FALSE, method = 'bray') ~ 
                          time + health + anti, 
                        permutations = 9999, by = 'terms',
                        data = tank_metadata, parallel = clus)
  
  stopCluster(cl = clus)
  write_rds(the_adonis, '../intermediate_files/adonis.rds.gz')
}

the_adonis

#### VEGAN ####
if(file.exists('../intermediate_files/tank_nmds.rds.gz')){
  the_nmds <- read_rds('../intermediate_files/tank_nmds.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  the_nmds <- column_to_rownames(tank_data_cpm, 'sample_id') %>%
    metaMDS(distance = 'bray', binary = FALSE, trymax = 1000, 
            parallel = clus)
  stopCluster(cl = clus)
  write_rds(the_nmds, '../intermediate_files/field_tank_nmds.rds.gz')
}


health_plot <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(tank_metadata, by = 'sample_id') %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = health)) +
  geom_point(size = 3)

time_plot <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(tank_metadata, by = 'sample_id') %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = time)) +
  geom_point(size = 3)

anti_plot <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(tank_metadata, by = 'sample_id') %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = anti)) +
  geom_point(size = 3)

health_plot + time_plot + anti_plot
ggsave('../Results/nmds_plots.png')

if(file.exists('../intermediate_files/field_tank_asvArrows.rds.gz')){
  asv_fit <- read_rds('../intermediate_files/field_tank_asvArrows.rds.gz')
} else {
  asv_fit <- field_data_cpm %>%
    select(-sample_id) %>%
    envfit(the_nmds, env = .,
           permutations = 9999)
  write_rds(asv_fit, '../intermediate_files/field_tank_asvArrows.rds.gz')
}


#### Plot ####
colony_points_nmds <- scores(the_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(field_metadata, by = 'sample_id')

filter(colony_points_nmds, NMDS1 > 0.8)

health_plot <- colony_points_nmds %>% 
  mutate(health = if_else(health == 'D', 'Diseased', 'Healthy')) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = health, shape = health)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(colour = 'Disease\nState',
       shape = 'Disease\nState') +
  theme_classic() +
  theme(legend.position = 'right',
        panel.background = element_rect(colour = 'black'),
        strip.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10))
# ggsave('../Results/nmds_field_tank.png', height = 7, width = 7)

timepoint_plot <- colony_points_nmds %>% 
  mutate(timepoint = str_replace_all(timepoint, c('W' = 'Jan', 'S' = 'Jul'))) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(data=. %>% select(-timepoint), 
             colour="grey60", size = 0.75) +
  geom_point(aes(colour = site, shape = health), size = 1.5) +
  facet_wrap(~timepoint, labeller = labeller(timepoint = ~str_replace(., '_', ' - '))) +
  guides(shape = 'none') +
  labs(colour = 'Site') +
  theme_classic() +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(colour = 'black', size = 14),
        axis.text = element_text(colour = 'black', size = 10))
