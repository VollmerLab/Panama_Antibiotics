#### Libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(patchwork)

distance_metric <- 'robust.aitchison'
redo_analysis <- FALSE
alpha <- 0.05

#### Data ####
microbiome_data <- read_rds('../intermediate_files/full_tank_microbiome.rds')

tank_full_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
                      show_col_types = FALSE) %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac)

tank_metadata <- tank_full_data %>%
  select(-asv_id, -domain:-species, -log2_cpm_norm:-n_reads) %>%
  distinct %>%
  mutate(treatment = str_c(anti, health, sep = '_'))

tank_data_cpm <- tank_full_data %>%
  select(sample_id, asv_id, log2_cpm_norm) %>%
  pivot_wider(names_from = asv_id, 
              values_from = log2_cpm_norm)

taxonomy <- select(tank_full_data, asv_id, domain:species) %>%
  distinct()

fully_analyzed_asvs <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
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
                    var = sd(log2_cpm_norm) < 1e-12,
                    .by = c('time', 'treatment')) %>%
           pull(var) %>%
           any) %>%
  pull(asv_id) %>%
  unique

#### PERMANOVA ####
if(file.exists('../intermediate_files/adonis_post.rds.gz') & !redo_analysis){
  the_adonis_pre <- read_rds('../intermediate_files/adonis_pre.rds.gz')
  the_adonis_post <- read_rds('../intermediate_files/adonis_post.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  pre_metadata <- filter(tank_metadata, time == 'before')
  pre_distance <- filter(tank_data_cpm, sample_id %in% pre_metadata$sample_id) %>%
    column_to_rownames('sample_id') %>%
    vegdist(., binary = FALSE, method = distance_metric)
  
  # pre_distance <- subset_samples(microbiome_data, otu_table(microbiome_data) %>% 
  #                  rownames %in% pre_metadata$sample_id) %>%
  #   distance(method = distance_metric)
  
  the_adonis_pre <-  adonis2(pre_distance ~ anti, 
            permutations = 9999, by = 'terms',
            data = pre_metadata, parallel = clus)
  write_rds(the_adonis_pre, '../intermediate_files/adonis_pre.rds.gz')
  
  the_beta_pre <- betadisper(pre_distance, pre_metadata$anti,
                             bias.adjust = TRUE) %>%
    anova(permutations = 9999)
  write_rds(the_beta_pre, '../intermediate_files/beta_pre.rds.gz')
  
  post_metadata <- filter(tank_metadata, time == 'after')
  post_distance <- filter(tank_data_cpm, sample_id %in% post_metadata$sample_id) %>%
    column_to_rownames('sample_id') %>%
    vegdist(., binary = FALSE, method = distance_metric)
  # post_distance <- subset_samples(microbiome_data, otu_table(microbiome_data) %>% 
  #                                  rownames %in% post_metadata$sample_id) %>%
  #   distance(method = distance_metric)
  the_adonis_post <- adonis2(post_distance ~ treatment, 
                            permutations = 9999, by = 'terms',
                            data = post_metadata, parallel = clus)
  write_rds(the_adonis_post, '../intermediate_files/adonis_post.rds.gz')
  
  the_beta_post <- betadisper(post_distance, post_metadata$treatment,
                              bias.adjust = TRUE) %>%
    anova(permutations = 9999)
  write_rds(the_beta_post, '../intermediate_files/beta_post.rds.gz')
  
  stopCluster(cl = clus)
  
}

the_adonis_pre
the_beta_pre
the_adonis_post
the_beta_post

#### NMDS ####
if(file.exists('../intermediate_files/post_nmds.rds.gz') & !redo_analysis){
  the_nmds_pre <- read_rds('../intermediate_files/pre_nmds.rds.gz')
  the_nmds_post <- read_rds('../intermediate_files/post_nmds.rds.gz')
} else {
  library(parallel)
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  
  pre_metadata <- filter(tank_metadata, time == 'before')
  the_nmds_pre <- filter(tank_data_cpm, sample_id %in% pre_metadata$sample_id) %>%
    column_to_rownames('sample_id') %>%
    metaMDS(distance = distance_metric, binary = FALSE, trymax = 1000,
            parallel = clus)
  # the_nmds_pre <- subset_samples(microbiome_data, otu_table(microbiome_data) %>% 
  #                  rownames %in% pre_metadata$sample_id) %>%
  #   distance(method = distance_metric) %>%
  #   metaMDS(trymax = 1000, parallel = clus)
  
  
  write_rds(the_nmds_pre, '../intermediate_files/pre_nmds.rds.gz')
  
  
  post_metadata <- filter(tank_metadata, time == 'after')
  the_nmds_post <- filter(tank_data_cpm, sample_id %in% post_metadata$sample_id) %>%
    column_to_rownames('sample_id') %>%
    metaMDS(distance = distance_metric, binary = FALSE, trymax = 1000,
            parallel = clus)
  
  # the_nmds_post <- subset_samples(microbiome_data, otu_table(microbiome_data) %>% 
  #                                   rownames %in% post_metadata$sample_id) %>%
  #   distance(method = distance_metric) %>%
  #   metaMDS(trymax = 1000, parallel = clus)
  
  write_rds(the_nmds_post, '../intermediate_files/post_nmds.rds.gz')
  
  stopCluster(cl = clus)
}

veganCovEllipse<-function(x, se = TRUE, conf = 0.95, npoints = 100){
  #X is a dataframe of 2 coordinates
  
  covariance_mat <- cov.wt(x, wt=rep(1/nrow(x), nrow(x)))
  
  cov <- covariance_mat$cov
  
  if(se){cov <- cov * sum(covariance_mat$wt^2)}
  
  center <- covariance_mat$center
  
  scale <- sqrt(qchisq(conf, 2))
  
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov))) %>% as_tibble()
}

veganGetCentroids <- function(x){
  #X is a dataframe of 2 coordinates
  
  covariance_mat <- cov.wt(x, wt=rep(1/nrow(x), nrow(x)))
  
  cov <- covariance_mat$cov
  enframe(covariance_mat$center) %>%
    pivot_wider()
}


nmds_data <- bind_rows(
  scores(the_nmds_pre)$sites %>%
    as_tibble(rownames = 'sample_id'),
  scores(the_nmds_post)$sites %>%
    as_tibble(rownames = 'sample_id'),
) %>%
  left_join(tank_metadata, by = 'sample_id') %>%
  mutate(time = str_to_sentence(time) %>% str_c('\nDisease Dose') %>% fct_inorder()) 

centroid_data <- nmds_data %>%
  select(time, anti, health, NMDS1, NMDS2) %>%
  nest_by(time, anti, health) %>%
  reframe(veganGetCentroids(data))

ellipse_data <- nmds_data %>%
  select(time, anti, health, NMDS1, NMDS2) %>%
  nest_by(time, anti, health) %>%
  reframe(veganCovEllipse(data, se = FALSE, conf = 0.66))


#### Fit biplot ####
if(file.exists('../intermediate_files/post_asvArrows.rds.gz') & !redo_analysis){
  pre_asv_arrows <- read_rds('../intermediate_files/pre_asvArrows.rds.gz')
  prost_asv_arrows <- read_rds('../intermediate_files/post_asvArrows.rds.gz')
} else {
  
  pre_metadata <- filter(tank_metadata, time == 'before')
  pre_asv_arrows <- filter(tank_data_cpm, sample_id %in% pre_metadata$sample_id) %>%
    select(-sample_id) %>%
    select(all_of(fully_analyzed_asvs)) %>%
    envfit(the_nmds_pre, env = .,
           permutations = 9999)
  write_rds(pre_asv_arrows, '../intermediate_files/pre_asvArrows.rds.gz')
  
  post_metadata <- filter(tank_metadata, time == 'after')
  post_asv_arrows <- filter(tank_data_cpm, sample_id %in% post_metadata$sample_id) %>%
    select(-sample_id) %>%
    select(all_of(fully_analyzed_asvs)) %>%
    envfit(the_nmds_post, env = .,
           permutations = 9999)
  write_rds(post_asv_arrows, '../intermediate_files/post_asvArrows.rds.gz')
}

vectors_to_tibble <- function(vegan_vector){
  as_tibble(vegan_vector$arrows, rownames = 'asv_id') %>%
    mutate(r2 = vegan_vector$r,
           p.value = vegan_vector$pvals)
}

sig_asvs <- bind_rows(before = vectors_to_tibble(pre_asv_arrows$vectors),
          after = vectors_to_tibble(post_asv_arrows$vectors),
          .id = 'time') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(time = str_to_sentence(time) %>% str_c('\nDisease Dose') %>% fct_inorder()) %>%
  mutate(fdr = p.adjust(p.value, 'fdr'), .by = time) %>%
  filter(fdr < alpha) 



bind_rows(before = as_tibble(scores(the_nmds_pre)$species, rownames = 'asv_id'),
         after = as_tibble(scores(the_nmds_post)$species, rownames = 'asv_id'),
         .id = 'time') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  mutate(time = str_to_sentence(time) %>% str_c('\nDisease Dose') %>% fct_inorder()) %>%
  filter(!is.na(family)) %>%
  select(time, family, NMDS1, NMDS2) %>%
  nest_by(time, family) %>%
  filter(nrow(data) > 5) %>%
  reframe(veganGetCentroids(data)) %>%
  
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point() +
  facet_wrap(~time)


bind_rows(before = as_tibble(scores(the_nmds_pre)$species, rownames = 'asv_id'),
          after = as_tibble(scores(the_nmds_post)$species, rownames = 'asv_id'),
          .id = 'time') %>%
  left_join(taxonomy, by = 'asv_id') %>%
  filter(!is.na(family)) %>%
  select(time, family, NMDS1, NMDS2) %>%
  nest_by(time, family) %>%
  filter(nrow(data) > 5) %>%
  reframe(veganGetCentroids(data)) %>%
  pivot_wider(names_from = time, values_from = c(NMDS1, NMDS2)) %>%
  mutate(dist = sqrt((NMDS1_after - NMDS1_before)^2 + (NMDS2_after - NMDS2_before)^2)) %>%
  arrange(-dist)

#### Make Plot ####
nmds_data %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_path(data = ellipse_data, aes(colour = health, linetype = anti),
            linewidth = 1.5, show.legend = FALSE) +
  
  geom_point(aes(fill = health, shape = anti), size = 1.5) +
  
  geom_point(data = centroid_data, aes(fill = health, shape = anti),
             size = 5) +
  
  geom_segment(data = sig_asvs, colour = 'black', xend = 1, yend = 1) +
  
  scale_colour_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  scale_fill_manual(values = set_names(wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                       c('H', 'D')),
                    breaks = c('D', 'H'), 
                    labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                    drop = FALSE) +
  
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  scale_linetype_manual(values = c('N' = 'solid', 'A' = 'dotted'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  
  guides(#colour = guide_legend(override.aes = list(size = 4, shape = 'circle filled')),
         fill = guide_legend(override.aes = list(size = 4, shape = 'circle filled', linewidth = 2)),
         #shape = guide_legend(override.aes = list(size = 4, fill = 'black')),
         linetype = guide_legend(override.aes = list(size = 4, fill = 'black', linewidth = 2))) +
  labs(x = 'NMDS1', 
       y = 'NMDS2',
       colour = 'Health\nState',
       fill = 'Health\nState',
       shape = 'Antibiotic\nTreatment',
       linetype = 'Antibiotic\nTreatment') +
  facet_wrap(~time) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())
ggsave('../Results/asv_nmds.png', height = 5, width = 12)







