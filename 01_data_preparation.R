rm(list = ls())
library(tidyverse)
library(phyloseq)
library(skimr)
library(caret)
library(janitor)


# read phyloseq objects for sw (surface water) and ss (surface sediment) -------
phylo_sw <- read_rds('data/ps16S_Lakes380.rds') # surface water data 
phylo_ss <- read_rds('data/ps16S_Lakes380_ss.rds') # surface sediment data 

# check number of distinct lakes-------------
phylo_sw %>% sample_data() %>% as_tibble() %>% distinct(Code)
phylo_ss %>% sample_data() %>% as_tibble() %>% distinct(Code)

# read pathogens list ------
patho <- read_csv('data/putative_patho_Lakes380.csv')

patho %>% 
  separate(
    taxonomy_hit,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";"
  ) %>% 
  distinct(Order, Family, Species) %>% 
  arrange(Order, Family, Species) %>% 
  write_csv('tables/pathogen_list.csv')


# filter pathogens --------
pruned_phylo_sw_raw <- 
  prune_taxa(distinct(patho, query_id) %>% deframe(), phylo_sw) %>% # filter pathogen taxa
  prune_taxa(taxa_sums(.) > 0, .) 

write_rds(pruned_phylo_sw_raw, 'data/pruned_phylo_sw_raw.rds')

pruned_phylo_sw <- 
  pruned_phylo_sw_raw  %>% # remove empty samples
  transform_sample_counts(function(x) {x/sum(x)}) # convert to relative abundance
  
# merge_samples("Code") %>% # sum samples by lake
#   transform_sample_counts(function(x) {x/sum(x)})

pruned_phylo_ss_raw <-
  prune_taxa(distinct(patho, query_id) %>% deframe(), phylo_ss) %>% # filter pathogen taxa
  prune_taxa(taxa_sums(.) > 0, .) 

write_rds(pruned_phylo_ss_raw, 'data/pruned_phylo_ss_raw.rds')


pruned_phylo_ss <- 
  pruned_phylo_ss_raw %>% # remove empty samples
  transform_sample_counts(function(x) {x/sum(x)})  # convert to relative abundance
  # merge_samples("Code") %>% # sum samples by lake
  # transform_sample_counts(function(x) {x/sum(x)})


# check data ---
pruned_phylo_sw 
pruned_phylo_ss 

# otu_all <- merge_phyloseq(pruned_phylo_sw, pruned_phylo_ss)

# melt data and convert to tibble ----
sw_otu_table <-
  psmelt(pruned_phylo_sw) %>% 
  as_tibble() %>% 
  select(Code, OTU, Abundance) %>% 
  pivot_wider(
    names_from = OTU,
    values_from = Abundance,
    values_fn = ~mean(.x, na.rm = T)
  ) %>% 
  drop_na() %>% 
  rename(code = Code)

ss_otu_table <-
  psmelt(pruned_phylo_ss) %>% 
  as_tibble() %>% 
  select(Code, OTU, Abundance) %>% 
  pivot_wider(
    names_from = OTU,
    values_from = Abundance,
    values_fn = ~mean(.x, na.rm = T)
  ) %>% 
  drop_na() %>% 
  rename(code = Code)

# add metadata ---------
t_meta <- read_csv('data/clean_metadata.csv')
sw_pathons_by_lake_metadata <- right_join(t_meta, sw_otu_table, by = 'code')
ss_pathons_by_lake_metadata <- right_join(t_meta, sw_otu_table, by = 'code')


# # relative abundance transformation of otu----------
# rel_otu_sw <- decostand(sw_otu_table[, -1], "total") %>% as_tibble() %>% bind_cols(sw_otu_table[, 1], .)
# rel_otu_ss <- decostand(ss_otu_table[, -1], "total") %>% as_tibble() %>% bind_cols(ss_otu_table[, 1], .)


# get specific metadata---------
t_meta_sw <- t_meta %>% filter(code %in% sw_otu_table$code)
t_meta_ss <- t_meta %>% filter(code %in% ss_otu_table$code)

# write datasets------
# # OTUs
# write_csv(sw_otu_table, 'data/clean_sw_pathogens.csv')
# write_csv(ss_otu_table, 'data/clean_ss_pathogens.csv')

# OTUs rel abund
write_csv(sw_otu_table, 'data/rel_otu_sw.csv')
write_csv(ss_otu_table, 'data/rel_otu_ss.csv')

# meta
write_csv(t_meta_sw, 'data/t_meta_sw.csv')
write_csv(t_meta_ss, 'data/t_meta_ss.csv')

# meta and OTUs
write_csv(sw_pathons_by_lake_metadata, 'data/clean_sw_pathogens_metadata.csv')
write_csv(ss_pathons_by_lake_metadata, 'data/clean_ss_pathogens_metadata.csv')


# # check number of zeros
# sum(melted_sw$Abundance==0)/nrow(melted_sw)*100
# sum(melted_ss$Abundance==0)/nrow(melted_ss)*100
# 
# melted_ss %>% group_by(Code) %>% summarise(N = sum(Abundance)) %>% ggplot() + geom_histogram(aes(N)) + scale_x_log10()
# melted_sw %>% group_by(Code) %>% summarise(N = sum(Abundance)) %>% ggplot() + geom_histogram(aes(N)) + scale_x_log10()
