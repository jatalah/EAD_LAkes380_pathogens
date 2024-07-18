library(tidyverse)
library(caret)
library(janitor)
library(skimr)

source('C:/Users/javiera/OneDrive - Cawthron/Stats/R code/half_detection_limits_JA.R')

# read and prepare metadata----- 
raw_metadata <- read_csv('data/metadata_2021.csv')
skim(raw_metadata)

metadata <- 
  raw_metadata %>% 
  mutate(Forestry =  For.Harv + Exot.For, .keep = 'unused') %>%
  mutate(
    Native  = Alpine.Grass + Indi.Hard + Fernland + Flaxland + Indi.For + Manuka + Matagouri + Sub.Alp.Shrub + Tussock,
    .keep = 'unused'
  ) %>% 
  # select predictor variables based on John's paper (based on correlation and ecological knowledge)
  select(
    Code,
    BulS.TN,
    BulS.Phosphorus,
    BulS.Sulfur,
    High.Prod.Exotic.Grass,
    Low.Prod.Grass,
    Latitude,
    Altitude,
    LakeArea,
    Secchi.Disk,
    WCS.Chla,
    WCS.DOC,
    Max.Depth,
    Forestry,
    Native
  )  %>% 
  # fix half detection limit values 
  mutate(across(c(BulS.Sulfur, BulS.TN, WCS.Chla, WCS.DOC), ~h_detect_lim(.x) %>% as.numeric())) %>% 
  select(-BulS.Sulfur) %>% # r>0.8
  left_join(read_csv('data/road_distance_data.csv')) %>% # add road access data
  clean_names() %>% 
  write_csv('data/clean_metadata.csv')


metadata  <- read_csv('data/clean_metadata.csv')

skim_without_charts(metadata) %>% yank('numeric') %>% as_tibble()

# impute, scale and transform predictors---------
t_para <-
  preProcess(metadata %>% select(-code),
             method = c("bagImpute", "center", "scale", "YeoJohnson"))

t_meta <- 
  predict(t_para, metadata) %>% 
  clean_names() 

# check VIF-------
source('C:/Users/javiera/OneDrive - Cawthron/Stats/R code/vif_simple.R')
vif_func(t_meta)

# check predictors distribution-----
t_meta %>% 
  pivot_longer(-code) %>% 
  ggplot(aes(value)) +
  geom_density() +
  facet_wrap(~name, scale = 'free')

write_csv(t_meta, 'data/clean_metadata.csv')


# summary table of predictor variables------
full_join(
  metadata %>% filter(code %in%  read_csv('data/rel_otu_sw.csv')$code),
  metadata %>% filter(code %in%  read_csv('data/rel_otu_ss.csv')$code),
) %>% skim_without_charts(-code)  %>% 
  yank("numeric") %>% 
  as_tibble() %>% 
  rename(variable = skim_variable) %>% 
  mutate(across(is.numeric, ~round(., 2))) %>% 
  write_csv('tables/predictors_summary.csv')


# separate tables for sediment and water datasets-----------------
metadata %>% filter(code %in%  read_csv('data/rel_otu_sw.csv')$code) %>%
  skim_without_charts(-code)  %>% 
  yank("numeric") %>%
  as_tibble() %>%
  rename(variable = skim_variable) %>%
  mutate(across(is.numeric, ~ round(., 2)))

metadata %>% filter(code %in%  read_csv('data/rel_otu_ss.csv')$code) %>%
  skim_without_charts(-code)  %>% yank("numeric") %>%
  as_tibble() %>%
  rename(variable = skim_variable) %>%
  mutate(across(is.numeric, ~ round(., 2)))


