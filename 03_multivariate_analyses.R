rm(list = ls())
library(vegan)
library(ggord)
library(ggvegan)
library(broom)
theme_set(theme_minimal(base_size = 9))
load('02_multivariate_analyses.RData')

# read data-----
t_meta_sw <- read_csv('data/t_meta_sw.csv')
rel_otu_sw <- read_csv('data/rel_otu_sw.csv')

t_meta_ss <- read_csv('data/t_meta_ss.csv')
rel_otu_ss <- read_csv('data/rel_otu_ss.csv')

anti_join(t_meta_ss, t_meta_sw)

f_root <- function(x) {sqrt(sqrt(x))}

# Run a multivariate multiple regression based on Bray-Curtis dissimilarities------------
adonis_funct <-
  function(otu, env, type) {
    adonis2(
      f_root(otu[, -1]) ~
        bul_s_tn + bul_s_phosphorus + high_prod_exotic_grass + low_prod_grass + latitude + altitude + lake_area + secchi_disk + wcs_chla + wcs_doc + max_depth + forestry + native + distance_to_road,
      data = env,
      method = 'bray',
      permutations = 999,
      by = type
    )
  }

# DistLMs-----------------
permdist_sw_term <- adonis_funct(otu = rel_otu_sw, env = t_meta_sw, type = 'term')
tidy(permdist_sw_term) %>% write_csv('tables/distlm_sw_term.csv')

permdist_sw <- adonis_funct(otu = rel_otu_sw, env = t_meta_sw, type = 'margin')
tidy(permdist_sw) %>% write_csv('tables/distlm_sw_margin.csv')

permdist_ss <- adonis_funct(otu = rel_otu_ss, env = t_meta_ss, type = 'margin')
tidy(permdist_ss) %>% write_csv('tables/distlm_ss_margin.csv')

permdist_ss_term <- adonis_funct(otu = rel_otu_ss, env = t_meta_ss, type = 'term')
tidy(permdist_ss_term) %>% write_csv('tables/distlm_ss_term.csv')

# redundancy analyses------
rda_func <- function (otu, env, met) {
  dbrda(
    f_root(otu[, -1]) ~ bul_s_tn + bul_s_phosphorus + high_prod_exotic_grass + low_prod_grass + latitude + altitude + lake_area + secchi_disk + wcs_chla + wcs_doc + max_depth + forestry + native + distance_to_road,
    data = env,
    method = met,
    permutations = 999
  )
}

# fit RDAs--------
rda_sw <- rda_func(rel_otu_sw, t_meta_sw, met = 'bray')
rda_sw_jacc <- rda_func(rel_otu_sw, t_meta_sw, met = 'jaccard')
rda_sw_ordisted <- ordistep(rda_sw)

rda_ss <- rda_func(rel_otu_ss, t_meta_ss, met = 'bray')
rda_ss_jacc <- rda_func(rel_otu_ss, t_meta_ss, met = 'jaccard')
rda_ss_ordisted <- ordistep(rda_ss)

# ANOVAs RDA-------
anova_rda_sw_term <- anova(rda_sw, by = 'term', permutations = 999)
anova_rda_sw_term_jacc <- anova(rda_sw_jacc, by = 'term', permutations = 999)

tidy(anova_rda_sw_term) %>% 
  mutate(across(where(is.numeric), ~round(.x, 2)),
         p.value = pvalue(p.value)) %>% 
  write_csv('tables/anova_rda_sw_term.csv', na = "")

tidy(anova_rda_sw_term_jacc) %>% write_csv('tables/anova_rda_sw_term_jacc.csv')

anova_rda_sw_margin <- anova(rda_sw, by = 'margin', permutations = 999)
anova_rda_sw_margin

# anova(rda_sw_ordisted, by = 'term', permutations = 999)

anova_rda_ss_term <- anova(rda_ss, by = 'term', permutations = 999)
anova_rda_ss_term_jacc <- anova(rda_ss_jacc, by = 'term', permutations = 999)

tidy(anova_rda_ss_term) %>% 
  mutate(across(where(is.numeric), ~round(.x, 2)),
         p.value = pvalue(p.value)) %>% 
  write_csv('tables/anova_rda_ss_term.csv', na = "")


tidy(anova_rda_ss_term_jacc) %>% write_csv('tables/anova_rda_ss_term.csv')

anova_rda_ss_margin <- anova(rda_ss, by = 'margin', permutations = 999)
anova_rda_ss_margin

# anova(rda_ss_ordisted, by = 'term', permutations = 999)


# RsquareAdj(rda_sw)$adj.r.squared
RsquareAdj(rda_sw)$r.squared
RsquareAdj(rda_ss)$r.squared

source('C:/Users/javiera/OneDrive - Cawthron/Stats/R code/theme_javier.R')

# biplots figure-----------
lawa_cols <- c("#20a7ad", "#85bb5b", "#ffa827", "#ff8400", "#e85129", "darkred")
tl_sw <- left_join(rel_otu_sw,
                   read_csv('data/sbti_data_clean.csv'),
                   by = c("code" = "lake_id")) %>% select(trophic_level) %>% mutate(
                     trophic_level = fct_relevel(
                       trophic_level,
                       'Microtrophic',
                       "Oligotrophic",
                       "Mesotrophic",
                       "Eutrophic",
                       "Hypertrophic"
                     )
                   )

tl_ss <- left_join(rel_otu_ss,
                   read_csv('data/sbti_data_clean.csv'),
                   by = c("code" = "lake_id")) %>% select(trophic_level) %>% mutate(
                     trophic_level = fct_relevel(
                       trophic_level,
                       'Microtrophic',
                       "Oligotrophic",
                       "Mesotrophic",
                       "Eutrophic",
                       "Hypertrophic"
                     )
                   )

original_names <- c(
  "bul_s_tn",
  "bul_s_phosphorus",
  "high_prod_exotic_grass",
  "low_prod_grass",
  "latitude",
  "altitude",
  "lake_area",
  "secchi_disk",
  "wcs_chla",
  "wcs_doc",
  "max_depth",
  "forestry",
  "native",
  "distance_to_road"
)

# Transformed vector
transformed_names <- c(
  "TN",
  "P",
  "High prod. exotic grass",
  "Low prod. grass",
  "Latitude",
  "Altitude",
  "Area",
  "Secchi disk",
  "Chl-a",
  "DOC",
  "Depth",
  "Forestry",
  "Native",
  "Distance to road"
)

name_changes <- setNames(transformed_names, original_names)

ggpubr::ggarrange(
  ggord(
    rda_sw,
    vec_ext = .7,
    parse = F,
    repel = T,
    max.overlaps = 15,
    alpha = .5,
    coord_fix = F,
    size = 2,
    arrow = .1,
    txt = 3,
    ellipse = F, 
    vec_lab = name_changes,
    grp_title = "Trophic level",
    grp_in = tl_sw$trophic_level
    
  ) +
    theme_javier(base_size = 9) +
    scale_color_manual(values = lawa_cols),
  ggord(
    rda_ss,
    vec_ext = .7,
    parse = F,
    repel = T,
    max.overlaps = 15,
    alpha = .5,
    coord_fix = F, 
    size = 2,
    arrow = .1,
    txt = 3,
    ellipse = F, 
    vec_lab = name_changes,
    grp_title = "Trophic level",
    grp_in = tl_ss$trophic_level
  ) +
    theme_javier(base_size = 9) +
    scale_color_manual(values = lawa_cols),
  labels = 'AUTO', common.legend = T, legend = 'bottom'
) 

ggsave(
  last_plot(),
  filename = 'figures/biplot.svg',
  dpi = 300,
  height = 4,
  width = 8
)

ggsave(
  last_plot(),
  filename = 'figures/biplot.png',
  dpi = 300,
  height = 4,
  width = 8
)

# nMDS -------
mds_func <- function(otu){capscale(sqrt(otu[,-1])~1, distance = 'bray', data = otu)}
mds_ss <- mds_func(otu = rel_otu_ss)

mds_plot_funct <- function(mds, env){
  fortify(mds, display = 'sites') %>% 
    bind_cols(env) %>% 
    ggplot(aes(MDS1, MDS2)) +
    geom_point() +
    scale_color_viridis_c()}

mds_plot_funct(mds_ss, t_meta_ss)


# CAP -----
cap_func <-
  function(data,meta){
  capscale(
    data[, -1] ~ high_prod_exotic_grass,
    # decostand(data[, -1], "pa") ~ high_prod_exotic_grass,
    data = meta,
    method = 'bray',
    permutations = 99
  )}

cap_sw <- cap_func(rel_otu_sw, t_meta_sw)
cap_ss <- cap_func(rel_otu_ss, t_meta_ss)


# CAP plot----
autoplot(cap_sw, layers = 'sites', arrows = TRUE)
autoplot(cap_ss, layers = 'sites', arrows = TRUE)

# ANOVAs
anova(cap_sw, by = "axis")
anova(cap_sw, by = "margin")
anova(cap_ss, by = "axis")
anova(cap_ss, by = "margin")

# Plot CAP1 vs hpg----
summary(cap_ss)$sites %>% 
  as_tibble() %>% 
  bind_cols(t_meta_ss %>% select(high_prod_exotic_grass)) %>%
  ggplot(aes(high_prod_exotic_grass, CAP1)) +
  geom_point(alpha = .5, position = position_jitter(width = .1)) +
  geom_smooth(method = lm)

