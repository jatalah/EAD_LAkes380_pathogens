rm(list = ls())
library(ggVennDiagram)
library(ggpubr)

rel_otu_sw <- read_csv('data/rel_otu_sw.csv')
rel_otu_ss <- read_csv('data/rel_otu_ss.csv')

taxa_list <- list(ss = names(rel_otu_ss[,-1]), sw = names(rel_otu_sw[,-1]))

ven_taxa <- 
  ggVennDiagram(taxa_list,
              label_alpha = 0.2,
              category.names = c("Sediment", "Water"),
              show_intersect = F, 
              edge_size = .5,
              label_size = 5) + 
  scale_fill_distiller(guide = F, palette = "Reds", direction = 1) +
  scale_color_manual(values = c('gray50', 'grey50'))

ven_taxa

# save plot-----------
ggsave(
  ven_taxa,
  filename = 'figures/venn_diagram,.png',
  width = 4,
  height = 3,
  dpi = 300,
  bg = 'white'
)


t_meta_sw <- read_csv('data/t_meta_sw.csv')
t_meta_ss <- read_csv('data/t_meta_ss.csv')


lakes_list <- list(ss = t_meta_ss$code, sw = t_meta_sw$code)
ven_lakes <- 
  ggVennDiagram(lakes_list,
              label_alpha = 0.2,
              category.names = c("Sediment", "Water"),
              show_intersect = F, 
              edge_size = .5,
              label_size = 5) + 
  scale_fill_distiller(guide = F, palette = "Blues", direction = 1) +
  scale_color_manual(values = c('gray50', 'grey50'))


ggarrange(
  ven_lakes,
  ven_taxa,
  nrow = 1,
  labels = c('A. Number of lakes', "B. Number of pathogens"),
  font.label = list(size = 12, font = 'plain')
)

ggsave(
  last_plot(),
  filename = 'figures/venn_diagram_all.svg',
  width = 8,
  height = 4,
  dpi = 300,
  bg = 'white'
)


t_meta_ss 
t_meta_sw

  
