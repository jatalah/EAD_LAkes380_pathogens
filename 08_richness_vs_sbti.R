library(tidyverse)
library(ggpubr)
library(broom)
rm(list = ls())
load('07_richness_vs_sbti.RData')
source('theme_javier.R')

sbti_data <- 
read_csv('data/sbti_data_new.csv') %>% 
mutate(trophic_level =
         cut(
           SBTI,
           breaks = c(-Inf, 2, 3, 4, 5, 6, 7),
           labels = c(
             "Microtrophic",
             "Oligotrophic",
             "Mesotrophic",
             "Eutrophic",
             "Supertrophic",
             "Hypertrophic"
           ),
           right = FALSE
         )) %>% 
  write_csv('data/sbti_data_clean.csv')
  
lawa_cols <- c("#20a7ad", "#85bb5b", "#ffa827", "#ff8400", "#e85129", "darkred")

# boxplot---------
sbti_data_richness <- 
bind_rows(Water = read_csv('data/div_sw.csv'),
          Sediment = read_csv('data/div_ss.csv'),
          .id = 'sample') %>%
  left_join(sbti_data, by = c("code" = "lake_id")) %>% 
  as_tibble()

# summary stats --------
sbti_data_richness %>% 
  group_by(trophic_level) %>% 
  get_summary_stats(S, type = "robust") %>% 
  drop_na()


# anova----
sbti_data_richness %>%
  drop_na(trophic_level) %>% 
  filter(sample == "Sediment") %>% 
  car::leveneTest(log(S)~trophic_level, data = .)

# ANOVAS---
aov_water <- aov(S~trophic_level, data = filter(sbti_data_richness, sample == "Water"))
summary(aov_water)
TukeyHSD(aov_water) %>% tidy()


aov_sed <- aov(log10(S)~trophic_level, data = filter(sbti_data_richness, sample == "Sediment"))
summary(aov_sed)
TukeyHSD(aov_sed) %>% tidy()


sbti_data_richness %>%
  drop_na(trophic_level) %>% 
  filter(sample == "Sediment") %>% 
  lm(log10(S)~trophic_level, data = .) %>% 
  anova() %>% 
  tidy()

# Figure 4_ Boxplot_richness ----------
sbti_data_richness %>% 
drop_na(SBTI) %>%
  ggplot(aes(rev(trophic_level), S, fill = rev(trophic_level))) +
  geom_boxplot(alpha = .6,
               width = .5,
               linewidth = .3) +
  coord_flip() +
  scale_fill_manual(values = lawa_cols, guide = 'none') +
  facet_wrap(~sample, scales = "free_x") +
  theme_javier(base_size = 9) +
  labs(y = "Number of bacterial pathogens", x = "Trophic level", fill = NULL) +
  NULL

ggsave(
  last_plot(),
  filename = 'figures/figure_4.svg',
  dpi = 300,
  height = 3,
  width = 6, 
  bg = 'white'
)

# explore oligotrophic sediment sample with high richness values ---
sbti_data_richness %>%
  drop_na(SBTI) %>%
  filter(sample == "Water" &
           trophic_level == "Microtrophic") %>%
  select(S, code) %>%
  arrange(-S) %>%
  left_join(
    .,
    read_csv('data/metadata_2021.csv') %>% select(Lake, Code, Latitude, Longitude),
    by = c("code" = "Code")
  ) %>% 
  write_csv('tables/microtrophic_lakes_water_richness.csv')

