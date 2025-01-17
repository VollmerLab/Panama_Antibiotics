#### Libraries ####
library(tidyverse)
library(magrittr)
library(phyloseq)
library(microbiome)
library(taxize)
library(metagMisc)
library(edgeR)

remove_samples <- 'none'


#### Functions ####
filter_missingness <- function(data, prop_missing){
  #data = a DGElist object, model_sample = samples to use for the calculation of % missingness
  #prop_missing = maximum percentage of samples which can be 0 and still keep ASV in dataset
  keep <- rowMeans(data$counts == 0) <= prop_missing
  
  message('\n')
  message('ASVs Removed for not being expressed in enough samples: ', scales::comma(sum(!keep)))
  message('ASVs Kept by Filter: ', scales::comma(sum(keep)))
  message('\n')
  data[keep, keep.lib.sizes = FALSE]
}

phyloseq_filter_prevalence_fix <- function (physeq, prev.trh = 0.05, 
                                            abund.trh = NULL, threshold_condition = "OR", 
                                            abund.type = "total") {
  if (prev.trh > 1 | prev.trh < 0) {
    stop("Prevalence threshold should be non-negative value in the range of [0, 1].\n")
  }
  if (!is.null(abund.trh)) {
    if (abund.trh <= 0) {
      stop("Abundance threshold should be non-negative value larger 0.\n")
    }
  }
  prevalenceThreshold <- prev.trh * phyloseq::nsamples(physeq)
  prevdf <- prevalence(physeq, package = 'base')
  if (abund.type == "total") {
    prevdf$AbundFilt <- prevdf$TotalAbundance
  }
  if (abund.type == "mean") {
    prevdf$AbundFilt <- prevdf$MeanAbundance
  }
  if (abund.type == "median") {
    prevdf$AbundFilt <- prevdf$MedianAbundance
  }
  if (is.null(abund.trh)) {
    tt <- prevdf$Prevalence >= prevalenceThreshold
  }
  if (!is.null(abund.trh)) {
    if (threshold_condition == "OR") {
      tt <- (prevdf$Prevalence >= prevalenceThreshold | 
               prevdf$AbundFilt >= abund.trh)
    }
    if (threshold_condition == "AND") {
      tt <- (prevdf$Prevalence >= prevalenceThreshold & 
               prevdf$AbundFilt >= abund.trh)
    }
  }
  keepTaxa <- prevdf$Taxa[tt]
  res <- phyloseq::prune_taxa(taxa = keepTaxa, x = physeq)
  return(res)
}



#### Data ####
tank_microbiome <- read_rds('../intermediate_files/full_tank_microbiome.rds')

cyanos_to_add <- read_rds('../intermediate_files/cyanos_to_add.rds')

# Remove phylogenetic trees from both objects
tank_microbiome <- phyloseq(otu_table(tank_microbiome),
                            sample_data(tank_microbiome),
                            tax_table(tank_microbiome))
cyanos_to_add <- phyloseq(otu_table(cyanos_to_add),
                          sample_data(cyanos_to_add),
                          tax_table(cyanos_to_add))

# Merge the phyloseq objects
tank_microbiome <- merge_phyloseq(tank_microbiome, cyanos_to_add)


#### Normalize ASV counts ####
otu_tmm <- tank_microbiome %>%
  # subset_taxa(rownames(tax_table(microbiome_data)) %in% asvs_to_keep) %>%
  phyloseq_filter_prevalence_fix(prev.trh = 0.1) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  
  #Add any other filtering here
  filter_missingness(0.9) %>%
  
  edgeR::calcNormFactors(method = 'TMMwsp')

rownames(otu_tmm$counts)[rownames(otu_tmm$counts) %in% colnames(otu_table(cyanos_to_add))]

#### Normalize with rCLR - compositional ####
compositional_asvs <- otu_table(tank_microbiome) %>% 
  as.data.frame() %>%
  vegan::decostand(method = 'rclr') %>%
  as_tibble(rownames = 'sample_id') %>%
  pivot_longer(cols = -sample_id,
               names_to = 'asv_id',
               values_to = 'rclr')

#### Output CPM ####
full_data <- full_join(cpm(otu_tmm, log = TRUE, prior.count = 0.5,
                           normalized.lib.sizes = TRUE) %>%
                         as_tibble(rownames = 'asv_id') %>%
                         pivot_longer(cols = -asv_id,
                                      names_to = 'sample_id',
                                      values_to = 'log2_cpm_norm'),
                       
                       cpm(otu_tmm, log = FALSE, 
                           normalized.lib.sizes = TRUE) %>%
                         as_tibble(rownames = 'asv_id') %>%
                         pivot_longer(cols = -asv_id,
                                      names_to = 'sample_id',
                                      values_to = 'cpm_norm'),
                       by = c('asv_id', 'sample_id')) %>%
  
  full_join(cpm(otu_tmm, log = FALSE,
                normalized.lib.sizes = FALSE) %>%
              as_tibble(rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'cpm'),
            by = c('asv_id', 'sample_id')) %>%
  
  full_join(as_tibble(otu_tmm$counts, rownames = 'asv_id') %>%
              pivot_longer(cols = -asv_id,
                           names_to = 'sample_id',
                           values_to = 'n_reads'),
            by = c('asv_id', 'sample_id')) %>%
  
  inner_join(compositional_asvs,
             by = c('asv_id', 'sample_id')) %>%
  
  left_join(sample_data(tank_microbiome) %>%
              as_tibble(rownames = 'sample_id'), 
            by = 'sample_id') %>%
  left_join(as_tibble(otu_tmm$samples, rownames = 'sample_id') %>%
              select(-group),
            by = 'sample_id') %>%
  left_join(tax_table(tank_microbiome) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_id'),
            by = c('asv_id'))
write_csv(full_data, '../intermediate_files/normalized_tank_asv_counts.csv')

