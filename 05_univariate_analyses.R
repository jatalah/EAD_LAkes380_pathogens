rm(list = ls())
library(tidyverse)
library(caret)
library(vegan)
library(quantreg)
library(ggpubr)
library(broom)
library(lares)
library(performance)
library(MASS)
library(plotmo)
library(gratia)
library(randomForest)
library(DHARMa)

theme_set(theme_minimal())

# read data------------
patho <- read_csv('data/putative_patho_Lakes380.csv')
t_meta_sw <- read_csv('data/t_meta_sw.csv')
rel_otu_sw <- read_csv('data/rel_otu_sw.csv')
t_meta_ss <- read_csv('data/t_meta_ss.csv')
rel_otu_ss <- read_csv('data/rel_otu_ss.csv')

sw <- read_csv('data/clean_sw_pathogens_metadata.csv')
ss <- read_csv('data/clean_ss_pathogens_metadata.csv')

# calculate univariate richness indices from OUTU data.
div_ss <-
  rel_otu_ss %>%
  mutate(S = rowSums(.[, -1] > 0), 
         .keep = 'used') %>%
  bind_cols(t_meta_ss) %>% 
  write_csv('data/div_ss.csv') %>% 
  select(-code)


div_sw <-
  rel_otu_sw %>%
  mutate(S = rowSums(.[,-1] > 0),
         .keep = 'used') %>%
  bind_cols(t_meta_sw) %>% 
  write_csv('data/div_sw.csv') %>% 
  select(-code) 

# summary stats---
div_ss %>% get_summary_stats(S)
div_sw %>% get_summary_stats(S)

# richness vs predictors------
div_ss %>%
  corr_var(S,
           method = "spearman",
           plot = FALSE,
           top = 30)  %>%
  write_csv('tables/cor_richness_ss.csv')

corr_var(div_ss,
         S,
         method = "spearman",
         plot = T,
         top = 30)

# hist(sqrt(div_sw$S))

div_sw %>%
  corr_var(S,
           method = "spearman",
           plot = FALSE) %>% write_csv('tables/cor_richness_sw.csv')

corr_var(div_sw,
         S,
         method = "spearman",
         plot = T,
         top = 30)

# hist(sqrt(div_sw$S))

# Lasso regression ----
lasso_ss <- 
  div_ss %>% 
  mutate(S = log(S)) %>% 
  lasso_vars(S)


print(lasso_ss$coef)
print(lasso_ss$metrics)
print(lasso_ss$plot)


lasso_sw <- lasso_vars(div_sw, S)
print(lasso_sw$coef)
print(lasso_sw$metrics)
print(lasso_sw$plot)

div_sw %>% 
  pivot_longer(cols = names(t_meta_sw[,-1])) %>% 
  ggplot(aes(value, S)) +
  geom_point(alpha = .3) +
  stat_smooth() +
  # scale_y_log10() +
  facet_wrap(~name, scales = 'free')

# library(glmmTMB)
# ss_model <-
#   glmmTMB(
#     S ~ bul_s_tn + bul_s_phosphorus + high_prod_exotic_grass + low_prod_grass +
#       latitude + altitude + lake_area + secchi_disk + wcs_chla + wcs_doc + max_depth +
#       forestry + native + distance_to_road,
#     data = div_ss,
#     family = compois(link = "log")
#   )

ss_model <-
  glm(
    S ~ .,
    data = div_ss,
    family = poisson(link = "log")
  )

summary(ss_model)
glance(ss_model)
check_overdispersion(ss_model)

sw_model <-
  glm(
    S ~ .,
    data = div_sw[-170,],
    family = poisson(link = "log")
  )

summary(sw_model)
glance(sw_model)
check_overdispersion(sw_model)


# Negative binomial -----

ss_nb <- 
glm.nb(S ~ .,
data = div_ss, link = "log"
)

summary(ss_nb)
model_performance(ss_nb)
simulateResiduals(ss_nb) %>% plot()

ss_nb_table <- tidy(ss_nb, conf.int = F, exponentiate = T) 

sw_nb <- 
  glm.nb(S ~ .,
               data = div_sw, link = "log"
  )
summary(sw_nb)

sw_nb_table <- 
tidy(sw_nb, conf.int = F, exponentiate = F) 

glance(sw_nb)
model_performance(sw_nb)
simulateResiduals(sw_nb) %>% plot()


bind_rows(water = sw_nb_table, sediment = ss_nb_table, .id = 'sample') %>% filter(p.value<0.05)

# save tables---
bind_rows(water = sw_nb_table, sediment = ss_nb_table, .id = 'sample') %>%
  mutate(p.value = scales::pvalue(p.value)) %>% 
  mutate(across(where(is.numeric), ~round(.x, 2))) %>% 
  write_csv('tables/richness_glms.csv')

# step_ss_model <- step(ss_model)

# Random forest ---------
rf_ss <- 
train(log10(S) ~ . ,
      data = div_ss[,-c(2:3)],
      method = "rf")

rf_ss

varImp(rf_ss$finalModel) %>% 
  rownames_to_column(var = "predictor") %>% 
  as_tibble() %>% 
  arrange(-Overall)

rf_sw <- 
  train(log10(S) ~ . ,
        data = div_sw[,-c(2:3)],
        method = "rf")


varImp(rf_sw$finalModel) %>% 
  rownames_to_column(var = "predictor") %>% 
  as_tibble() %>% 
  arrange(-Overall)

plotmo(rf_sw, )

# GAMS-----------------
gam_ss <- 
  train(log(S) ~ . ,
        data = div_ss[,-c(2:3)],
        method = "gam")
gam_ss

plotmo(rf_sw$finalModel)

importance(rf_sw$finalModel)

rel_otu_ss %>% 
  mutate(S = rowSums(.[, -1]>0), .keep = 'used') %>% 
  bind_cols(t_meta_ss) %>% 
  pivot_longer(cols = names(t_meta_ss[,-1])) %>% 
  ggplot(aes(value, S)) +
  geom_point(alpha = .3) +
  stat_smooth() +
  scale_y_log10() +
  facet_wrap(~name, scales = 'free')

# logistic models on dominant taxa -----
read_csv('tables/top_20_ss.csv')$ASV[1]
  
rel_otu_sw %>% 
  select(code, read_csv('tables/top_20_sw.csv')$ASV[1]) %>% 
  left_join(t_meta_sw) %>% 
  # left_join(read_csv('data/clean_metadata.csv')) %>% 
  rename(taxa = 2) %>% 
  mutate(taxa = if_else(taxa>0, 1,0)) %>% 
  pivot_longer(cols = -c(1:2)) %>% 
  ggplot(aes(value, taxa)) +
  geom_point(alpha = .3) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  facet_wrap(~name, scales = 'free')


rel_otu_sw %>% 
  select(code, read_csv('tables/top_20_sw.csv')$ASV[2]) %>% 
  left_join(t_meta_sw) %>% 
  rename(taxa = 2) %>% 
  select(-code) %>% 
  mutate(taxa = if_else(taxa>0, 1,0)) %>% 
  glm(taxa~., data = .) %>% 
  summary()


# LASSO models---------
top_20_sw <- read_csv('tables/top_20_sw.csv')$ASV[1:5]

lassos <- 
rel_otu_sw %>% 
  select(code, top_20_sw) %>% 
  left_join(t_meta_sw) %>% 
  pivot_longer(cols = all_of(top_20_sw)) %>% 
  group_by(name) %>% 
  nest() %>% 
  mutate(lasso = map(data, ~lasso_vars(.x, value)),
         coef = map(lasso, ~.x$coef),
         metrics = map(lasso, ~.x$metrics))

lassos %>% 
  select(coef) %>% 
  unnest() %>% 
  rename(ASV = name) %>% 
  left_join(read_csv('tables/top_20_sw.csv'))