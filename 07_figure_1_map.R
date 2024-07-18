library(sf)
library(tidyverse)
library(janitor)
library(rnaturalearth)
library(ggspatial)
library(ggmagnify)
rm(list = ls())

coords <- 
  read_csv('data/metadata_2021.csv') %>% 
  select(Code, Lake, Latitude, Longitude) %>% 
  clean_names()
  # st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

coords_samples <-
  bind_rows(
    'Water' = coords %>% filter(code %in% read_csv('data/t_meta_sw.csv')$code),
    'Sediment' = coords %>% filter(code %in% read_csv('data/t_meta_ss.csv')$code),
    .id = 'type'
  ) 

coords_sample <- 
bind_rows(

# shared
'Both' = semi_join(
  coords %>% filter(code %in% read_csv('data/t_meta_ss.csv')$code),
  coords %>% filter(code %in% read_csv('data/t_meta_sw.csv')$code),
  by = 'code'
),

# sediment
'Sediment' = 
anti_join(
  coords %>% filter(code %in% read_csv('data/t_meta_ss.csv')$code),
  coords %>% filter(code %in% read_csv('data/t_meta_sw.csv')$code),
  by = 'code'
),

# water

'Water' = 
  anti_join(
  coords %>% filter(code %in% read_csv('data/t_meta_sw.csv')$code),
  coords %>% filter(code %in% read_csv('data/t_meta_ss.csv')$code),
  by = 'code'
),
.id = 'type'
 )


write_csv(coords_sample, 'data/coords_samples.csv')

nz_high <- ne_countries(country = 'new zealand', scale = "large", returnclass = "sf")

figure_1_map <- 
ggplot(nz_high) +
  geom_sf(fill = 'gray95') +
  theme_minimal(base_size = 9) +
  coord_sf(xlim = c(166, 178.8),
           ylim = c(-47.35,-34.35)) +
  scale_y_continuous(breaks = seq(-34, -48, by = -4)) + 
  scale_x_continuous(breaks = seq(166, 178, by = 4)) +
  geom_point(
    data = coords_sample,
    aes(longitude, latitude, fill = fct_rev(type)),
    alpha = .8,
    pch = 21,
    color = 'gray60'
  ) +
  theme(legend.position = c(.2, .8)) +
  labs(fill = NULL) +
  scale_fill_viridis_d(option = 'B') +
  annotation_scale(
    location = "br",
    width_hint = 0.2,
    bar_cols = 'white',
    style  = 'ticks'
  ) +
  labs(x = NULL, y = NULL) +
  annotation_north_arrow(
    location = "br",
    pad_y = unit(0.4, "in"),
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  NULL

figure_1_map

write_rds(figure_1_map, 'figures/figure_1_maps.rds')

ggsave(
  figure_1_map,
  filename = 'figures/map_figure.svg',
  dpi = 300,
  height = 4,
  width = 3, 
  bg = 'white'
)