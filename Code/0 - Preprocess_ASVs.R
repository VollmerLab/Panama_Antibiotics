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
top_classification <- function(data, threshold = 80){
  rename_with(data, .cols = domain:species,
              ~str_c(., '_name')) %>%
    pivot_longer(cols = -asv_id,
                 names_to = c('taxon_level', '.value'),
                 names_pattern = '(.*)_(.*)') %>%
    mutate(taxon_level = factor(taxon_level, ordered = TRUE,
                                levels = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'))) %>%
    filter(confidence >= threshold) %>%
    select(-confidence) %>%
    pivot_wider(names_from = taxon_level, 
                values_from = 'name')
}

# taxonomy <- reconciled_higher_taxonomy; lowest_level = 'genus'
synonomize_taxonomy <- function(taxonomy, lowest_level = 'genus'){
  all_levels <- c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'asv_id')
  remove_levels <- all_levels[(which(all_levels == lowest_level) + 1):length(all_levels)]
  
  duplicated_higher_taxonomy <- select(taxonomy, -asv_id, -all_of(remove_levels)) %>%
    distinct %>%
    filter(!is.na(!!sym(lowest_level))) %>%
    group_by(!!sym(lowest_level)) %>%
    filter(n() > 1) %>%
    arrange(!!sym(lowest_level)) %>%
    pull(!!sym(lowest_level)) %>%
    unique
  
  if(lowest_level == 'species'){
    genus_species_mismatch <- select(taxonomy, -asv_id, -all_of(remove_levels)) %>%
      distinct %>%
      filter(!is.na(!!sym(lowest_level))) %>%
      filter(!str_detect(species, genus)) %>%
      pull(species) %>%
      unique %>%
      sort
    
    duplicated_higher_taxonomy <- unique(c(duplicated_higher_taxonomy, genus_species_mismatch))
  }
  
  if(length(duplicated_higher_taxonomy) > 0){
    update_classification <- classification(duplicated_higher_taxonomy, db = 'ncbi')
    
    refresh_taxa <- tibble(!!sym(lowest_level) := names(update_classification)) %>%
      mutate(row = row_number()) %>%
      rowwise(!!sym(lowest_level)) %>%
      reframe(update_classification[[row]]) %>%
      filter(rank %in% all_levels,
             rank != lowest_level) %>%
      select(-id) %>%
      pivot_wider(names_from = rank, 
                  values_from = name) %>%
      mutate(domain = 'Bacteria')
    
    
    
    out <- filter(taxonomy, !(!!sym(lowest_level) %in% refresh_taxa[[lowest_level]])) %>%
      bind_rows(filter(taxonomy, !!sym(lowest_level) %in% refresh_taxa[[lowest_level]]) %>%
                  select(asv_id, !!sym(lowest_level), all_of(remove_levels)) %>%
                  left_join(refresh_taxa, by = lowest_level))
  } else {
    out <- taxonomy
  }
  out
}

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


#### Data ####
microbiome_raw <- read_rds("../Data/field_tank_newPS_deciphersilva.rds") 

microbiome_data <- microbiome_raw %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  subset_taxa(domain == "Bacteria" &
                phylum != "Cyanobacteria" &
                !is.na(phylum) &
                family != "Mitochondria" &
                class != "Chloroplast" & 
                order != "Chloroplast") %>%
  prune_samples(sample_sums(.) > 1000, .) %>%
  subset_samples(site != 'MG') %>% #no disease samples
  prune_samples(!str_detect(rownames(sample_data(.)), 
                            str_c(remove_samples, collapse = '|')), .)

## Reclassify with Decipher
if(!file.exists('../intermediate_files/update_taxonomy.csv')){
  Biostrings::writeXStringSet(file="../intermediate_files/all_asvs.fasta", refseq(microbiome_data))
  #http://www2.decipher.codes/ClassifyOrganisms.html
} else {
  new_taxonomy <- read_csv('../intermediate_files/update_taxonomy.csv', show_col_types = FALSE) %>%
    top_classification() 
}


old_taxonomy <- tax_table(microbiome_data) %>%
  as.data.frame %>%
  as_tibble(rownames = 'asv_id')


updated_taxonomy <- bind_rows(new_taxonomy,
                              filter(old_taxonomy, !asv_id %in% new_taxonomy$asv_id))


reconciled_higher_taxonomy <- synonomize_taxonomy(updated_taxonomy, 'species') %>%
  synonomize_taxonomy('genus') %>%
  synonomize_taxonomy('family') %>%
  synonomize_taxonomy('order') %>%
  synonomize_taxonomy('class') %>%
  synonomize_taxonomy('phylum') 

tax_table(microbiome_data) <- column_to_rownames(reconciled_higher_taxonomy, 'asv_id') %>% as.matrix

metadata <- sample_data(microbiome_data) %>%
  as_tibble(rownames = 'sample_id') %>%
  mutate(health = if_else(health == 'N', 'H', health),
         geno = if_else(geno == 'GE', 'GR', geno),
         sample_id = str_replace(sample_id, 'GE', 'GR'))

sequenced_tank_samples <- filter(metadata, dataset == 'tank') %>%
  mutate(fragment_id = str_remove(sample_id, 'P[0-9]_')) %>%
  pull(fragment_id) %>%
  unique

#### Prep original tank data ####
#Corals put in tank July 11 & dosed with antibiotics
#July 12 dosed with antibiotics again
#July 13 post-antibiotic leisson pic - anti_0
#Homogenates from Sebastian's reef
#disease dosed 12pm July 13
#Post disease pic July 15


original_tank <- read_csv('../Data/su17_tank_surv_data.csv', 
                          show_col_types = FALSE) %>%
  mutate(fragment_id = str_c(str_c('Bin', Block), if_else(antibiotics == 'no_antibiotics', 'N', 'A'),
                             if_else(exposure == 'healthy', 'N', 'D'), 
                             str_sub(genotype, 1, 2) %>% str_to_upper(), sep = '_'),
         .before = everything()) %>%
  filter(fragment_id %in% sequenced_tank_samples) %>%
  filter(!genotype %in% c('Grey', 'Blue')) %>% #Did not sequence these genotypes for 16s
  rename_with(~str_replace_all(., c('21' = 'jul14am',	'30' = 'jul14pm',	
                                    '45' = 'jul15am',	'54' = 'jul15pm',	
                                    '69' = 'jul16am',	'78' = 'jul16pm',	
                                    '93' = 'jul17am',	'102' = 'jul17pm',
                                    '117' = 'jul18am', '126' = 'jul18pm',	
                                    '141' = 'jul19am',	'150' = 'jul19pm',	
                                    '165' = 'jul20am',	'174' = 'jul20pm'))) %>%
  
  pivot_longer(cols = starts_with('jul'),
               names_to = 'date',
               values_to = 'state') %>%
  
  mutate(tank = str_c('Bin', Block),
         health = if_else(state == 0, 'H', 'D'),
         geno = str_sub(genotype, 1, 2) %>% str_to_upper(),
         # date = fct_inorder(date),
         #no time 0_anti in this dataset - must keep the metadata from before - all are healthy
         time = case_when(date %in% c('jul15pm') ~ '2_exp',
                          date %in% c('jul20pm') ~ '8_exp',
                          TRUE ~ 'other'),
         .keep = 'unused') %>%
  filter(time != 'other') %>%
  select(-Condition) %>%
  
  mutate(plate = case_when(time == '2_exp' ~ 'P5',
                           time == '8_exp' ~ 'P6',
                           TRUE ~ NA_character_),
         sample_id = str_c(plate, fragment_id, sep = '_'),
         anti = if_else(antibiotics == 'no_antibiotics', 'N', 'A'),
         exposure = if_else(exposure == 'healthy', 'N', 'D'),
         anti_exposure = str_c(anti, exposure, sep = '_'),
         site = 'tank', season = 'S', dataset = 'tank',
         year = '2017') %>%
  select(sample_id, health, time, tank, anti_exposure, anti,
         exposure, geno, plate, year, season, site, dataset)

#### Swap tank data ####
pre_exposure <- filter(metadata, dataset == 'tank') %>%
  filter(time == '0_anti') %>%
  select(-resist, -anti_health) %>%
  mutate(exposure = 'pre')

updated_tank_data <- filter(metadata, dataset == 'tank', 
                            time != '0_anti') %>%
  select(-health, -resist, -anti_health) %>%
  left_join(select(original_tank, sample_id, health, exposure),
            by = 'sample_id') %>%
  bind_rows(pre_exposure) %>%
  mutate(tank = str_c(tank, anti, exposure, sep = '_'),
         tank = str_remove(tank, '_pre'))

metadata <- select(metadata, -resist, -anti_health) %>%
  filter(dataset != 'tank') %>%
  bind_rows(updated_tank_data)

sample_data(microbiome_data) <- column_to_rownames(metadata, 'sample_id')

#### Subset to just tank data & asvs present in any samples ####
tank_microbiome <- microbiome_data %>%
  subset_samples(dataset == 'tank') %>%
  phyloseq_filter_sample_wise_abund_trim(minabund = 1)

sample_data(tank_microbiome) <- sample_data(tank_microbiome) %>%
  as_tibble(rownames = 'sample_id') %>%
  select(sample_id, health, time) %>%
  separate(sample_id, into = c('plate', 'tank', 'anti', 'exposure', 'geno'),
           sep = '_', remove = FALSE) %>%
  mutate(time = str_extract(time, '[0-9]+') %>% as.integer(),
         time_fac = if_else(time == 0, 'before', 'after'),
         # exposure = if_else(time == 0, 'pre', exposure),
         tank = str_c(tank, anti, exposure, sep = '_'),
         fragment = str_c(geno, tank, sep = '_')) %>%
  column_to_rownames('sample_id')
write_rds(tank_microbiome, '../intermediate_files/full_tank_microbiome.rds')

#### Normalize ASV counts ####
otu_tmm <- tank_microbiome %>%
  # subset_taxa(rownames(tax_table(microbiome_data)) %in% asvs_to_keep) %>%
  phyloseq_filter_prevalence(prev.trh = 0.1) %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  otu_table() %>% 
  t %>% #NOTE: *genus and family do not need the t but ASVs need the t*
  as.data.frame %>%
  as.matrix %>% 
  DGEList(remove.zeros = TRUE) %>%
  
  #Add any other filtering here
  filter_missingness(0.9) %>%
  
  edgeR::calcNormFactors(method = 'TMMwsp')

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
  
  left_join(sample_data(tank_microbiome) %>%
              as_tibble(rownames = 'sample_id'), 
            by = 'sample_id') %>%
  left_join(as_tibble(otu_tmm$samples, rownames = 'sample_id') %>%
              select(-group),
            by = 'sample_id') %>%
  left_join(tax_table(microbiome_data) %>%
              as.data.frame() %>%
              as_tibble(rownames = 'asv_id'),
            by = c('asv_id'))
write_csv(full_data, '../intermediate_files/normalized_tank_asv_counts.csv')
