#### Libraries ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(patchwork)
library(parallel)

clus <- makeCluster(detectCores() - 1)
clusterEvalQ(clus, library(vegan))

distance_metric <- 'robust.aitchison'
redo_analysis <- TRUE
alpha_test <- 0.05

#### Functions ####
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

safe_qvalue <- possibly(.f = ~qvalue(.)$qvalues, otherwise = NA_real_) #better method for doing this! - only want to sub 0 when it fails

#### Data ####
tank_full_data <- read_csv('../intermediate_files/normalized_tank_asv_counts.csv',
                      show_col_types = FALSE) %>%
  mutate(anti = factor(anti, levels = c('N', 'A')),
         exposure = factor(exposure, levels = c('N', 'D')),
         health = factor(health, levels = c('H', 'D')),
         time_fac = factor(time_fac, levels = c('before', 'after'))) %>%
  select(-time) %>%
  rename(time = time_fac)

tank_metadata <- tank_full_data %>%
  select(-asv_id, -domain:-species, -log2_cpm_norm:-rclr) %>%
  distinct %>%
  mutate(treatment = str_c(anti, health, sep = '_') %>% factor,
         time_treat = str_c(time, anti, health, sep = '_') %>% factor)

tank_data_cpm <- tank_full_data %>%
  select(sample_id, asv_id, n_reads) %>%
  pivot_wider(names_from = asv_id, 
              values_from = n_reads)

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
distiance_mat <- column_to_rownames(tank_data_cpm, 'sample_id') %>%
  vegdist(., binary = FALSE, method = distance_metric)

the_adonis <- adonis2(distiance_mat ~ time_treat, 
                      permutations = 9999, by = 'terms',
                      data = tank_metadata, parallel = clus)

# adonis_posthoc <- pairwise.adonis2(distiance_mat ~ time_treat, 
#                                    permutations = 9999, by = NULL,
#                                    data = tank_metadata, parallel = clus)

# adonis_posthoc <- adonis_posthoc[2:length(adonis_posthoc)] %>%
#   map_dfr(~as_tibble(., rownames = 'term'), .id = 'pair') %>%
#   group_by(pair) %>%
#   mutate(dDF = Df[term == 'Residual']) %>%
#   ungroup() %>%
#   filter(term == 'Model')

the_beta <- betadisper(distiance_mat, tank_metadata$time_treat,
                       bias.adjust = TRUE, sqrt.dist = FALSE,
                       type = 'median') %>%
  anova(permutations = 9999)

# the_beta_tukey <- betadisper(distiance_mat, tank_metadata$time_treat,
#                              bias.adjust = TRUE, type = 'median') %>%
#   TukeyHSD()

stopCluster(cl = clus)

the_adonis
# adonis_posthoc

the_beta
# the_beta_tukey

#### Posthoc replication ####
adonis_posthoc <- function(posthoc, metadata, compositional_data, B = 1000){
  B <- B - 1
  if(posthoc == 'disease'){
    new_meta <- filter(metadata, anti == 'N') %>%
      mutate(indVar = health)
    
  } else if (posthoc == 'antibiotic'){
    new_meta <- filter(metadata, health == 'H') %>%
      mutate(indVar = anti)
    
  } else if (posthoc == 'time'){
    new_meta <- filter(metadata, health == 'H') %>%
      mutate(indVar = time)
    
  } else {
    break('error incorrect posthoc type')
  }
  
  dist_mat <- filter(tank_data_cpm, sample_id %in% new_meta$sample_id) %>%
    column_to_rownames('sample_id') %>%
    vegdist(., binary = FALSE, method = distance_metric)
  
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  posthoc_adon <- adonis2(dist_mat ~ indVar, 
          permutations = B, by = NULL,
          data = new_meta, parallel = clus)
  
  stopCluster(cl = clus)
  
  posthoc_adon %>%
    as_tibble(rownames = 'term') %>%
    mutate(denDF = Df[term == 'Residual'],
           .after = 'Df') %>%
    rename(numDF = Df,
           r2 = R2,
           f.value = 'F',
           p.value = `Pr(>F)`) %>%
    filter(term == 'Model') %>%
    select(-term, -SumOfSqs) %>%
    mutate(contrast = posthoc, .before = numDF)
}

# posthoc <- 'disease'; metadata <- tank_metadata; compositional_data <- tank_data_cpm; B <- 1000
beta_posthoc <- function(posthoc, metadata, compositional_data, B = 1000){
  B <- B - 1
  if(posthoc == 'disease'){
    new_meta <- filter(metadata, anti == 'N') %>%
      mutate(indVar = health)
    
  } else if (posthoc == 'antibiotic'){
    new_meta <- filter(metadata, health == 'H') %>%
      mutate(indVar = anti)
    
  } else if (posthoc == 'time'){
    new_meta <- filter(metadata, health == 'H') %>%
      mutate(indVar = time)
    
  } else {
    break('error incorrect posthoc type')
  }
  
  dist_mat <- filter(tank_data_cpm, sample_id %in% new_meta$sample_id) %>%
    column_to_rownames('sample_id') %>%
    vegdist(., binary = FALSE, method = distance_metric)
  
  clus <- makeCluster(detectCores() - 1)
  clusterEvalQ(clus, library(vegan))
  
  posthoc_adon <- betadisper(dist_mat, pull(new_meta, indVar),
                             bias.adjust = TRUE, type = 'median') %>%
    anova(permutations = B, parallel = clus) 
  stopCluster(cl = clus)
  
  
  posthoc_adon %>%
    as_tibble(rownames = 'term') %>%
    mutate(denDF = Df[term == 'Residuals'],
           .after = 'Df') %>%
    rename(numDF = Df,
           f.value = 'F value',
           p.value = `Pr(>F)`) %>%
    filter(term == 'Groups') %>%
    select(-term, -`Sum Sq`, -`Mean Sq`) %>%
    mutate(contrast = posthoc, .before = numDF)
}
## Time - ignore disease
## Antibiotics - ignore disease
## Disease - ignore antibiotic

planned_posthocs_adon <- map_dfr(c('disease', 'antibiotic', 'time'), 
                                 ~adonis_posthoc(.x, tank_metadata, 
                                                 tank_data_cpm, B = 10000))

planned_posthocs_beta <- map_dfr(c('disease', 'antibiotic', 'time'), 
                                 ~beta_posthoc(.x, tank_metadata, 
                                                 tank_data_cpm, B = 10000))

planned_posthocs_adon
planned_posthocs_beta

#### PCoA ####
the_rda <- dbrda(column_to_rownames(tank_data_cpm, 'sample_id') ~ time_treat, 
                 data = tank_metadata, 
                 distance = distance_metric)

summary(the_rda)
anova(the_rda, by = 'terms')
#adonis2 = anova of RDA


const_var_exp <- summary(the_rda)$concont$importance[2, ]
total_var_exp <- summary(the_rda)$cont$importance[3,]

pcoa_data <- scores(the_rda)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(tank_metadata, by = 'sample_id') %>%
  # mutate(time = str_to_sentence(time) %>% str_c('\nDisease Dose') %>% fct_inorder()) %>%
  rename(Dim1 = dbRDA1, Dim2 = dbRDA2) %>%
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) 

centroid_data <- pcoa_data %>%
  select(time, anti, health, Dim1, Dim2) %>%
  nest_by(time, anti, health) %>%
  reframe(veganGetCentroids(data)) %>%
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) 

ellipse_data <- pcoa_data %>%
  select(time, anti, health, Dim1, Dim2) %>%
  nest_by(time, anti, health) %>%
  reframe(veganCovEllipse(data, se = FALSE, conf = 0.66)) %>%
  mutate(health = case_when(health == 'H' & anti == 'A' ~ 'H_A',
                            TRUE ~ health),
         health = factor(health, levels = c('H', 'H_A', 'D'))) 

dbrda_plot <- pcoa_data %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = health, fill = interaction(time == 'before', health),
             shape = anti)) +
  geom_path(data = ellipse_data, aes(colour = health),
            linewidth = 0.5, show.legend = FALSE, linetype = 'solid') +
  
  geom_point(size = 3, stroke = 1) +
  
  geom_point(data = centroid_data, size = 7, stroke = 2) +
  
  
  scale_colour_manual(values = set_names(c(wesanderson::wes_palette("Zissou1", 2, 
                                                                    type = "continuous"),
                                           "#80D1E9"),
                                         c('H', 'D', 'H_A')),
                      breaks = c('D', 'H'), 
                      labels = c('H' = 'Healthy', 'D' = 'Diseased'),
                      drop = FALSE) +
  scale_fill_manual(values = set_names(c("#80D1E9", wesanderson::wes_palette("Zissou1", 2, 
                                                                type = "continuous"),
                                         'white', 'white'),
                                       c('FALSE.H_A', 'FALSE.H', 'FALSE.D', 'TRUE.H_A', 'TRUE.H')),
                    breaks = c('TRUE.H', 'FALSE.H'),
                    labels = c('TRUE.H' = 'Before', 'FALSE.H' = "After")) +
  
  scale_shape_manual(values = c('N' = 'triangle filled', 'A' = 'triangle down filled'),
                     breaks = c('N', 'A'),
                     labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  
  # scale_linetype_manual(values = c('N' = 'solid', 'A' = 'dotted'),
  #                       breaks = c('N', 'A'),
  #                       labels = c('A'= 'Antibiotic\nTreated', 'N' = 'Untreated')) +
  
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 'circle'), stroke = 0.1),
    fill = guide_legend(override.aes = list(size = 4, shape = c('Before' = 'circle open', 
                                                                'After' = 'circle'), 
                                            linewidth = 2, colour = 'black',
                                            stroke = 0.1)),
    shape = guide_legend(override.aes = list(size = 4, fill = c("#3A9AB2", "#80D1E9"), 
                                             stroke = 0.1))) +
  labs(x = str_c('dbRDA1 (', scales::percent(const_var_exp[1], accuracy = 0.1),')'), 
       y = str_c('dbRDA2 (', scales::percent(const_var_exp[2], accuracy = 0.1),')'),
       colour = 'Health\nState',
       fill = 'Dose\nTiming',
       shape = 'Antibiotic\nTreatment') +
  # facet_wrap(~time) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())
# ggsave('../Results/Fig3_asv_dbrda.png', height = 7, width = 10)
dbrda_plot

data.frame(ord_in$CCA$biplot[, axes])

#### Family/Order Biplot ####
microbiome_data <- read_rds('../intermediate_files/full_tank_microbiome.rds')

microbiome_data <- subset_taxa(microbiome_data, 
                               otu_table(microbiome_data) %>% colnames %in% unique(tank_full_data$asv_id))

postProcess_envfit <- function(the_env, level = 'taxa'){
  as_tibble(the_env$vectors$arrows, rownames = level) %>%
    mutate(r2 = the_env$vectors$r,
           p.value = the_env$vectors$pvals)
}

tax_table(microbiome_data) %>%
  as.data.frame() %>%
  as_tibble(rownames = 'asv_id') %>%
  filter(asv_id %in% str_c('ASV', c(8, 25)))

sample_families <- microbiome_data %>%
  # subset_taxa(., 
  #             otu_table(microbiome_data) %>% 
  #               colnames %in% colnames(tank_data_cpm)[-1]) %>%
  microbiome::aggregate_taxa('family') %>%
  otu_table() %>%
  t %>%
  vegan::decostand(method = 'rclr')


the_env <- envfit(the_rda, sample_families, permutations = 9999) %>%
  postProcess_envfit(level = 'taxa') %>%
  mutate(fdr.bh = p.adjust(p.value, 'BH'),
         fdr.by = p.adjust(p.value, 'BY'),
         q.value = safe_qvalue(p.value), 
         .after = 'p.value') %>%
  rename(Dim1 = dbRDA1, Dim2 = dbRDA2) %>%
  filter(taxa != 'Unknown')

filter(the_env, taxa %in% c('Vibrionaceae', 'Fastidiosibacteraceae'))

filter(the_env, fdr.bh < 0.05, r2 > 0.09)

filter(the_env, fdr.bh < 0.05) %>%
  filter(str_detect(taxa, 'Are'))

filter(the_env, fdr.bh < 0.05) %>%
  filter(Dim1 > 0.9, Dim2 < 0) %>%
  filter(!taxa %in% c('Clade III', 'Verrucomicrobiaceae')) %>%
  arrange(Dim2)

dbrda_plot +
  geom_segment(data = filter(the_env, fdr.bh < 0.05, r2 > 0.09), 
               xend = 0, yend = 0, inherit.aes = FALSE,
               aes(x = 3.5 * Dim1, y = 3.5 * Dim2),
               alpha = 0.25, linetype = 'dashed') +
  geom_text(data = filter(the_env, fdr.bh < 0.05, r2 > 0.09), 
               inherit.aes = FALSE,
               aes(x = 3.5 * Dim1, y = 3.5 * Dim2,
                   label = taxa,
                   hjust = if_else(Dim1 < 0, 1, 0),
                   vjust = case_when(taxa == 'Pasteurellaceae' ~ 0.3,
                                     taxa == 'Maricaulaceae' ~ 0.7,
                                     
                                     # taxa == 'Alteromonadaceae' ~ 0.3,
                                     # taxa == 'Vicingaceae' ~ 0.7,
                                     
                                     taxa == 'Bdellovibrionaceae' ~ -0.4,
                                     taxa == 'Saprospiraceae' ~ 0.9,
                                     
                                     taxa == 'Rubritaleaceae' ~ 0.2,
                                     taxa == 'Lewinellaceae' ~ 0.6,
                                     
                                     taxa == 'Verrucomicrobiaceae' ~ 0.1,
                                     taxa == 'Clade III' ~ 0.9,
                                     
                                     # taxa == 'Cohaesibacteraceae' ~ 0.3,
                                     # taxa == 'Colwelliaceae' ~ 0.7,
                                     
                                     # taxa == 'Stappiaceae' ~ 0.1,
                                     # taxa == 'Pirellulaceae' ~ 1.5,
                                     
                                     # taxa == 'Halobacteriovoraceae' ~ 0.1,
                                     # taxa == 'Roseivirgaceae' ~ 0.8,
                                     
                                     taxa == 'Vibrionaceae' ~ 0.1,
                                     taxa == 'Pseudomonadaceae' ~ 0.9,
                                     
                                     taxa == 'Pseudobdellovibrionaceae' ~ 0.3,
                                     taxa == 'Halieaceae' ~ 0.7,
                                     
                                     taxa == 'AEGEAN-169 marine group' ~ 0.3,
                                     taxa == 'Endozoicomonadaceae' ~ 0.7,
                                     
                                     taxa == 'Moraxellaceae' ~ 0.1,
                                     taxa == 'Staphylococcaceae' ~ 0.1,
                                     taxa == 'Myxococcaceae' ~ 0,
                                     taxa == 'Oceanospirillaceae' ~ 1,
                                     taxa == 'Rickettsiaceae' ~ 0.8,
                                     
                                     taxa == 'Arenicellaceae' ~ -0.8,
                                     taxa == 'Roseobacteraceae' ~ -0.2,
                                     taxa == 'Hyphomonadaceae' ~ 0.5,
                                     taxa == 'Fastidiosibacteraceae' ~ 1.2,
                                     
                                     TRUE ~ 0.5)), 
            check_overlap = FALSE) +
  expand_limits(x = c(-5.5, 6))

ggsave('../Results/Fig4_asv_dbrda.png', height = 9, width = 12) 
ggsave('../Results/Fig4_r1.tiff', height = 9, width = 12, dpi = 'print')


#### NMDS ####
if(file.exists('../intermediate_files/the_nmds.rds')){
  the_nmds <- read_rds('../intermediate_files/the_nmds.rds')
} else {
  the_nmds <- column_to_rownames(tank_data_cpm, 'sample_id') %>%
    metaMDS(distance = distance_metric, 
            binary = FALSE, trymax = 5000,
            parallel = clus)
  write_rds(the_nmds, '../intermediate_files/the_nmds.rds')
}


stressplot(the_nmds_pre)
stressplot(the_nmds)


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

#
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
  
  # geom_segment(data = sig_asvs, colour = 'black', xend = 1, yend = 1) +
  
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





joint_nmds



joint_nmds_data <- scores(joint_nmds)$sites %>%
  as_tibble(rownames = 'sample_id') %>%
  left_join(tank_metadata, by = 'sample_id') %>%
  mutate(time = str_to_sentence(time) %>% str_c('\nDisease Dose') %>% fct_inorder()) 

joint_centroid_data <- joint_nmds_data %>%
  select(time, anti, health, NMDS1, NMDS2) %>%
  nest_by(time, anti, health) %>%
  reframe(veganGetCentroids(data))

joint_ellipse_data <- joint_nmds_data %>%
  select(time, anti, health, NMDS1, NMDS2) %>%
  nest_by(time, anti, health) %>%
  reframe(veganCovEllipse(data, se = FALSE, conf = 0.66))


joint_nmds_data %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_path(data = ellipse_data, aes(colour = health, linetype = anti),
            linewidth = 1.5, show.legend = FALSE) +
  
  geom_point(aes(fill = health, shape = anti), size = 1.5) +
  
  geom_point(data = centroid_data, aes(fill = health, shape = anti),
             size = 5) +
  
  # geom_segment(data = sig_asvs, colour = 'black', xend = 1, yend = 1) +
  
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
  # facet_wrap(~time) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'),
        legend.key = element_blank())
