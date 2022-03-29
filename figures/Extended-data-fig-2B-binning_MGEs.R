# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 2B - Binning MGEs ===================================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(viridis)
library(tidyverse)
library(patchwork)
library(googlesheets4)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('/Users/paolil/polybox/PhD/Exploratorium/sushilab-colors/palettes-paoli.R')

# Load data ==============================================================================

data_tbl = read_sheet('1YqvL7CDaW4fkQl9R9ZAsBKCuFtBCrJ0XV1rjLHA_x4s') %>%
  mutate(`Scaffold Length` = factor(`Scaffold Length`, levels = c("≥10Kb", "≥2Kb", "≥1Kb")),
         `Binning Result` = factor(`Binning Result`, levels = rev(c("MAG", "Bin", "Unbinned"))),
         `Genetic Element` = factor(`Genetic Element`, levels = c("All", "Chromosome", "Plasmid", "Virus", "Unannotated")))

vir_init = viridis(length(unique(data_tbl$`Binning Result`)), begin = .2)
vir_cols = vir_init[1:(length(vir_init) - 1)]
vir_grey = c(vir_cols, DescTools::ColToGrey(vir_init[length(vir_init)])) # Las color as greyscale


p1 = data_tbl  %>%
  ggplot() +
  geom_bar(aes(x = `Scaffold Length`,  y = `Number of scaffolds`), stat = 'identity') +
  facet_wrap(~`Genetic Element`, nrow = 1, scales = "free_y") +
  theme_bw() +
  theme_cell +
  theme(rect = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "mm"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(color = "lightgrey", fill = "lightgrey"),
        strip.text = element_text(size = 6, colour = "white", face = "bold", margin = margin(1, 1, 1, 1, "mm")))

p2 = data_tbl %>%
  ggplot() +
  geom_bar(aes(x = `Scaffold Length`,  y = `Perc. of scaffolds`, fill = `Binning Result`), stat = 'identity') +
  facet_wrap(~`Genetic Element`, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = rev(vir_grey)) +
  theme_bw() +
  theme_cell +
  theme(rect = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "mm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_rect(color = "lightgrey", fill = "lightgrey"),
        strip.text = element_text(size = 6, colour = "white", face = "bold", margin = margin(1, 1, 1, 1, "mm")))


p1 / p2

ggsave(paste0(figures_path_proj, "Figure-SX/Figure-SX-binning_MGEs.pdf"), width = two_col, height = 100, units = col_unit, device = cairo_pdf)
