# read pathogens list ------
rm(list = ls())
library(ggpubr)
library(skimr)
library(tidyverse)

patho <- read_csv('data/putative_patho_Lakes380.csv')
rel_otu_sw <- read_csv('data/rel_otu_sw.csv')
rel_otu_ss <- read_csv('data/rel_otu_ss.csv')


# Helper functions to get top 20 species -------------- 
sum_func <- function(d) {
  d %>%
    pivot_longer(names_to = "ASV",
                 values_to = "abundance",
                 cols = -code) %>%
    filter(abundance != 0) %>%
    left_join(patho %>% select(-c(hit_id:consensus_genus)), by = c('ASV' = 'query_id')) %>%
    group_by(code, species) %>%
    summarise(abundance = sum(abundance)) %>%
    group_by(species) %>%
    summarise(
      n_lakes = n_distinct(code),
      # is the number lakes they are present
      mean_rel_abund = sum(abundance) / nrow(d) * 100,
      # average relative abundance as % by lake
      rel_index = n_lakes * mean_rel_abund # rel abundace index
    ) 
}

# data summaries------
sum_func(d = rel_otu_ss) %>%
  skim_without_charts() %>% 
  yank('numeric') %>% 
  as_tibble()

sum_func(d = rel_otu_sw) %>% 
  skim_without_charts() %>% 
  yank('numeric') %>% 
  as_tibble()

# top 20 species---
# sediment
top_20_ss <- 
  sum_func(d = rel_otu_ss)  %>% 
  ungroup() %>% 
  top_n(n = 21, wt = n_lakes) %>% 
  arrange(-n_lakes) %>% 
  filter(species != "uncultured bacterium") %>% 
  write_csv('tables/top_20_ss.csv')

# water column
top_20_sw <- 
sum_func(d = rel_otu_sw)  %>%
  ungroup() %>% 
  top_n(n = 21, wt = n_lakes) %>% 
  arrange(-n_lakes) %>% 
  filter(species != "uncultured bacterium") %>% 
  write_csv('tables/top_20_sw.csv')


# Arrange data by descending order of n_lakes
bar_plot_funct <- function(data) {
  ggplot(data,
         aes(
           x = reorder(species, n_lakes),
           y = n_lakes / 287 * 100,
           fill = mean_rel_abund
         )) +
    geom_col(color = "gray40") +
    coord_flip() +
    theme_minimal(base_size = 10) +
    scale_fill_viridis_c(
      option = 'D',
      alpha = .8,
      trans = 'log10',
      name = "Mean relative abundance (%)"
    ) +
    theme(
      axis.text.y = element_text(face = 'italic'),
      legend.key.size = unit(.4, "cm"),
      plot.margin = unit(c(.5, 0.2, 0, 0.2), "cm") 
    ) +
    labs(y = "Percentage of lakes (%)", x = NULL)
}

figure_3_bar_plots <- 
ggarrange(
  bar_plot_funct(top_20_ss),
  bar_plot_funct(top_20_sw),
  common.legend = T,
  legend = 'bottom',
  labels = c('A. Water', 'B.Sediment'), 
  vjust = 1, 
  font.label = list(face = 'plain', size = 10)
)

figure_3_bar_plots

# write_rds(figure_3_bar_plots, 'figures/figure_3_bar_plots.rds')

# Save figure 3------------
figure_3_bar_plots <- read_rds('figures/rds/figure_3_bar_plots.rds')

ggsave(
  last_plot(),
  filename = 'figures/figure_3.svg',
  width = 7,
  height = 4,
  dpi = 300,
  bg = 'white'
)

# Figure S1 frequency distribution of pathogen taxa by lakes------
bind_rows(
  "Water" = sum_func(d = rel_otu_ss),
  "Sediment" = sum_func(d = rel_otu_sw),
  .id = "sample"
) %>%
  ggplot() +
  geom_histogram(
    aes(n_lakes),
    bins = 15,
    color = 'gray30',
    fill = 'gray90'
  ) +
  labs(x = 'Number of lakes', y = "Number of taxa") +
  facet_wrap( ~ fct_rev(sample))

# Save figure S1---------
ggsave(
  last_plot(),
  filename = 'figures/histograms_txa_per_lake.png',
  width = 6,
  height = 3,
  dpi = 300,
  bg = 'white'
)