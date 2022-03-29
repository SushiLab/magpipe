# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 1D - Taxonomic Novelty ==========================================================

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(viridis)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('/Users/paolil/polybox/PhD/Exploratorium/sushilab-colors/palettes-paoli.R')

# Load Summary table ---------------------------------------------------------------------

Summary_table = load_raw_summary() %>%
  left_join(load_general_metadata()) %>%
  mutate(drep_rep_genome = as.logical(drep_rep_genome))

# Format data table ----------------------------------------------------------------------

Summary_table = Summary_table %>%
  mutate(`Novelty Category` = ifelse(!grepl(";s__$", classification), "Known Species", NA)) %>%
  mutate(`Novelty Category` = ifelse(grepl(";s__$", classification), "New Species", `Novelty Category`)) %>%
  mutate(`Novelty Category` = ifelse(grepl(";g__;s__", classification), "New Clade", `Novelty Category`)) %>%
  group_by(drep_cluster) %>%
  mutate(`Cluster Type` = ifelse(all(grepl("_METAG_[A-Z]{8}$", `Bin Id`)), "Int. Only", NA)) %>%
  mutate(`Cluster Type` = ifelse(all(!grepl("_METAG_[A-Z]{8}$", `Bin Id`)), "Ext. Only", `Cluster Type`)) %>%
  mutate(`Cluster Type` = ifelse(is.na(`Cluster Type`), "Mixed", `Cluster Type`)) %>%
  ungroup()

Summary_table %>% pull(closest_placement_ani) %>% as.numeric() %>% summary()
Summary_table %>% filter(`Novelty Category` == "Known Species") %>% pull(closest_placement_ani) %>% as.numeric() %>% summary()
Summary_table %>% filter(`Novelty Category` == "Known Species") %>% pull(fastani_ani) %>% as.numeric() %>% summary()
Summary_table %>% filter(`Novelty Category` == "New Species") %>% pull(closest_placement_ani) %>% as.numeric() %>% summary()
Summary_table %>% filter(`Novelty Category` == "New Species") %>% pull(fastani_ani) %>% as.numeric() %>% summary()
Summary_table %>% filter(`Novelty Category` == "New Clade") %>% pull(closest_placement_ani) %>% as.numeric() %>% summary()
Summary_table %>% filter(`Novelty Category` == "New Clade") %>% pull(fastani_ani) %>% as.numeric() %>% summary()

table(Summary_table$`Novelty Category`)
table(Summary_table$`Cluster Type`)

Cluster_summary = Summary_table %>%
  #filter(grepl("d__Archaea", `classification`)) %>%
  group_by(drep_cluster) %>%
  summarize(`Novelty Category` = `Novelty Category`[drep_rep_genome],
            `Cluster Type` = unique(`Cluster Type`)) %>%
  group_by(`Novelty Category`, `Cluster Type`) %>%
  summarize(n = n(),
            Facet = "Clusters") %>%
  ungroup()
sum(Cluster_summary$n)

Genome_summary = Summary_table %>%
  group_by(`Cluster Type`) %>%
  summarize(`Known Species` = sum(`Novelty Category` == "Known Species"),
            `New Species` = sum(`Novelty Category` == "New Species"),
            `New Clade` = sum(`Novelty Category` == "New Clade"),
            Facet = "Genomes") %>%
  gather(value = n, key = `Novelty Category`, - Facet, - `Cluster Type`)

sum(Genome_summary$n)

# Get the plotting going =================================================================

# Inner layer ----------------------------------------------------------------------------

Figure_2E_1_table = Cluster_summary %>%
  # Add spacers
  rbind(tibble(`Cluster Type` = "Ext. Only", `Novelty Category` = "Spacer", n = round(sum(Cluster_summary %>% filter(`Cluster Type` == "Ext. Only") %>% pull(n))/10), Facet = NA)) %>%
  rbind(tibble(`Cluster Type` = "Mixed", `Novelty Category` = "Spacer", n = round(sum(Cluster_summary %>% filter(`Cluster Type` == "Mixed") %>% pull(n))/10), Facet = NA)) %>%
  rbind(tibble(`Cluster Type` = "Int. Only", `Novelty Category` = "Spacer", n = round(sum(Cluster_summary %>% filter(`Cluster Type` == "Int. Only") %>% pull(n))/10), Facet = NA)) %>%
  group_by(`Cluster Type`) %>%
  mutate(tmp_fraction = n / sum(n)) %>%
  ungroup()

Figure_2E_1_table = Figure_2E_1_table %>%
  mutate(`Cluster Type` = factor(`Cluster Type`, levels = c("Ext. Only", "Mixed", "Int. Only"))) %>%
  mutate(`Novelty Category` = factor(`Novelty Category`, levels = names(novelty_colors_reds))) %>%
  arrange(`Cluster Type`, `Novelty Category`)

# Compute percentages
Figure_2E_1_table$fraction <- Figure_2E_1_table$tmp_fraction / sum(Figure_2E_1_table$tmp_fraction)
# Compute the cumulative percentages (top of each rectangle)
Figure_2E_1_table$ymax <- cumsum(Figure_2E_1_table$fraction)
# Compute the bottom of each rectangle
Figure_2E_1_table$ymin <- c(0, head(Figure_2E_1_table$ymax, n=-1))

Figure_2E_1_table %>%
  ggplot() +
  geom_rect(aes(ymax=ymax, ymin=ymin, xmax=6, xmin=4, fill=`Novelty Category`)) +
  #geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = novelty_colors_reds) +
  coord_polar(theta="y", start = pi/40) +
  xlim(c(0, 10)) +
  theme_void() +
  theme(legend.position = "none")

bar_species = Figure_2E_1_table %>%
  filter(`Novelty Category` != 'Spacer') %>%
  mutate(`Cluster Type` = factor(`Cluster Type`, levels = rev(c("Int. Only", "Mixed", "Ext. Only")))) %>%
  ggplot() +
  geom_col(aes(x=`Cluster Type`, y = n, fill=`Novelty Category`)) +
  scale_fill_manual(values = novelty_colors_reds) +
  scale_y_continuous(labels = c("0", "1k", "2k", "3k", "4k")) +
  ylab("# Species") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        rect = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = unit(6, 'pt')),
        axis.text = element_text(color = 'black', size = unit(6, 'pt')),
        axis.title = element_text(size = unit(6, 'pt')),
        plot.margin = margin(0,0,0,0),
        panel.spacing = margin(0,0,0,0),
        line = element_line(size = unit(0.2, 'mm')))
bar_species

# Outer layer ----------------------------------------------------------------------------

Figure_2E_2_table = Genome_summary %>%
  # Add spacers
  rbind(tibble(`Cluster Type` = "Ext. Only", `Novelty Category` = "Spacer", n = round(sum(Genome_summary %>% filter(`Cluster Type` == "Ext. Only") %>% pull(n))/10), Facet = NA)) %>%
  rbind(tibble(`Cluster Type` = "Mixed", `Novelty Category` = "Spacer", n = round(sum(Genome_summary %>% filter(`Cluster Type` == "Mixed") %>% pull(n))/10), Facet = NA)) %>%
  rbind(tibble(`Cluster Type` = "Int. Only", `Novelty Category` = "Spacer", n = round(sum(Genome_summary %>% filter(`Cluster Type` == "Int. Only") %>% pull(n))/10), Facet = NA)) %>%
  group_by(`Cluster Type`) %>%
  mutate(tmp_fraction = n / sum(n)) %>%
  ungroup()

Figure_2E_2_table = Figure_2E_2_table %>%
  mutate(`Cluster Type` = factor(`Cluster Type`, levels = c("Ext. Only", "Mixed", "Int. Only"))) %>%
  mutate(`Novelty Category` = factor(`Novelty Category`, levels = names(novelty_colors_reds))) %>%
  arrange(`Cluster Type`, `Novelty Category`)

# Compute percentages
Figure_2E_2_table$fraction <- Figure_2E_2_table$tmp_fraction / sum(Figure_2E_2_table$tmp_fraction)
# Compute the cumulative percentages (top of each rectangle)
Figure_2E_2_table$ymax <- cumsum(Figure_2E_2_table$fraction)
# Compute the bottom of each rectangle
Figure_2E_2_table$ymin <- c(0, head(Figure_2E_2_table$ymax, n=-1))

Figure_2E_2_table %>%
  ggplot() +
  geom_rect(aes(ymax=ymax, ymin=ymin, xmax=10, xmin=8, fill=`Novelty Category`)) +
  #geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = novelty_colors_reds) +
  coord_polar(theta="y", start = pi/40) +
  xlim(c(0, 10)) +
  theme_void() +
  theme(legend.position = "none")

bar_genomes = Figure_2E_2_table %>%
  filter(`Novelty Category` != 'Spacer') %>%
  mutate(`Cluster Type` = factor(`Cluster Type`, levels = rev(c("Int. Only", "Mixed", "Ext. Only")))) %>%
  ggplot() +
  geom_col(aes(x=`Cluster Type`, y = n, fill=`Novelty Category`)) +
  scale_fill_manual(values = novelty_colors_reds) +
  scale_y_continuous(labels = c("0", "5k", "10k", "15k")) +
  ylab("# Genomes") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top",
        rect = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = unit(6, 'pt')),
        axis.text = element_text(color = 'black', size = unit(6, 'pt')),
        axis.title = element_text(size = unit(6, 'pt')),
        legend.key.size = unit(2, 'mm'),
        legend.box.margin = margin(0, 0, 0, 0, 'mm'),
        legend.margin = margin(1,1,1,1,'mm'),
        plot.margin = margin(0,0,0,0),
        panel.spacing = margin(0,0,0,0),
        line = element_line(size = unit(0.2, 'mm')))
bar_genomes

# Combined layers ------------------------------------------------------------------------

Figure_2E = ggplot() +
  geom_rect(data = Figure_2E_1_table, aes(ymax=ymax, ymin=ymin, xmax=7, xmin=5, fill=`Novelty Category`)) +
  geom_rect(data = Figure_2E_2_table, aes(ymax=ymax, ymin=ymin, xmax=10, xmin=8, fill=`Novelty Category`)) +
  #geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = novelty_colors_reds) +
  coord_polar(theta="y", start = pi/40) +
  xlim(c(0, 10)) +
  theme_void() +
  #theme_cell +
  theme(legend.title = element_blank(),
        legend.position = c(0.5,0.5),
        legend.key.size = unit(0.2, "mm"),
        legend.text = element_text(size = 6, color = "black"))

ggsave(paste0(figures_path_proj, "Figure-2/Figure-2E.pdf"), Figure_2E, width = 0.3*two_col, height = 50, units = col_unit)

# num clusters

Figure_2E_1_table %>%
  filter(`Novelty Category` != "Spacer") %>%
  group_by(`Cluster Type`) %>%
  summarize(n = sum(n))
Figure_2E_2_table %>%
  filter(`Novelty Category` != "Spacer") %>%
  group_by(`Cluster Type`) %>%
  summarize(n = sum(n))

# Barplot version ------------------------------------------------------------------------

library(patchwork)

Figure_2E_bars =
  bar_genomes + 
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,0,40)) + 
  bar_species + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) + 
  plot_layout(nrow = 1) + 
  plot_annotation(theme = theme(plot.margin = margin(0,0,0,0)))
Figure_2E_bars
ggsave(paste0(figures_path_proj, "Figure-1/Figure-1E-bars.pdf"), Figure_2E_bars, width = 69, height = 25, units = 'mm')



