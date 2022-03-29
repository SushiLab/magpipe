# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 2C - Archaeal Tree ===============================================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(treeio) # YuLab-SMU/treeio
library(ggtree) # YuLab-SMU/treeio
library(tidytree)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Prepare data ---------------------------------------------------------------------------

artree = "data/raw/go_microbiomics/gtdb/go_microbiomics-gtdb-4Chris/gtdbtk.ar122.classify.tree"

# Assumes data prep is loaded

# Archaea ================================================================================


# Inner plot -----------------------------------------------------------------------------

artrnw <- drop.tip(read.newick(artree),
                   genomes_summary %>%
                     filter(!`Representative Genome` |  `GTDB Taxonomy` == "N/A") %>% pull(`Genome Id`)
)
assertthat::assert_that(!any(Species_biosynth %>% pull(`# Genomes`) %>% is.na()))

summary(Species_biosynth$`# Genomes`)
summary(log(Species_biosynth$`# Genomes`))
# After collapsing at 15%, we have a maximum of 1304 Bacteria in a collapsed clade
scaled_colors = 
  var_scale_colors(log(1159), scaling_factor = 10, fun = "viridis", type = "D", begin = 0)# from max(bctree_collapsable_level_nodes$value, na.rm = T), before log
#  var_scale_colors(log(Species_biosynth$`# Genomes`), scaling_factor = 10, fun = "cmocean", type = "ice", begin = .1, end = 0.85, dir = -1)
#  var_scale_colors(log(Species_biosynth$`# Genomes`), scaling_factor = 10, fun = "scico", type = "oslo", begin = 0, end = 0.8)
#  var_scale_colors(log(Species_biosynth$`# Genomes`), scaling_factor = 10, fun = "viridis", type = "D", begin = 0)

artrnw_metadata = 
  as_tibble(artrnw) %>% 
  left_join(Species_biosynth %>% select(label = `Representative Genome`, `# Genomes`) %>% mutate(user_genome = TRUE)) %>%
  mutate(user_genome = ifelse(is.na(user_genome), FALSE, user_genome),
         color = ifelse(is.na(`# Genomes`), "lightgrey", scaled_colors[round(log(`# Genomes`)*10) + 1]))

artree_plot = ggtree(artrnw, layout="fan", open.angle=285, branch.length = "none", size = 0.35, alpha = 1, ) %<+%
  #artree_plot = ggtree(artrnw, layout="fan", open.angle=270, size = 0.35, alpha = 1) %<+%
  artrnw_metadata + # add metadata to the tree data
  aes(color = color) +
  scale_color_identity() +
  coord_polar(theta = 'y', start = 1.5*pi, direction = -1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        line = element_line(size = 0.2))
artree_plot

# Collapsing the clades
# /!\ Uncomment this line and comment the next if you need to update the backbon!!
#artree_collapsable_backbone_nodes = identify_collapsable_backbone(artree_plot, force = T)
artree_collapsable_backbone_nodes = identify_collapsable_backbone(artree_plot)

artree_plot_collapsed_backbone = collapse_backbone(artree_plot, artree_collapsable_backbone_nodes)

artree_collapsable_level_nodes = identify_collapsable_levels(artree_plot_collapsed_backbone, "# Genomes", perc_to_collapse = 0.15) %>%
  mutate(value = round(log(value)*10) + 1)
artree_plot_collapsed = collapse_level(artree_plot_collapsed_backbone, artree_collapsable_level_nodes, scaled_colors)

artree_plot_collapsed

ggsave("data/processed/figures/Figure-2/figure-2-archaeal_backbone.pdf", artree_plot_collapsed, width = 100, height = 100, units = 'mm')

artree_plot_collapsed_no_backbone = remove_backbone_layers(artree_plot_collapsed)
artree_plot_collapsed_no_backbone

ggsave("data/processed/figures/Figure-2/figure-2-archaeal_backbone_coloronly.pdf", artree_plot_collapsed_no_backbone, width = 100, height = 100, units = 'mm')

# Add labels -----------------------------------------------------------------------------

# Identify the taxa to display

artree_plot$data %>% filter(grepl("p__", label))
artree_plot$data %>% filter(grepl("GCA_002789275.1", label)) # -> GCA_002789275.1, Huber
artree_plot$data %>% filter(grepl("GCA_002254415.1", label)) # -> GCA_002254415.1, EX4484-52
artree_plot$data %>% filter(grepl("c__", label))
artree_plot$data%>% filter(label == "1.0:c__Poseidoniia")

MRCA(artree_plot$data, c(2136, 2147, 2154, 2159, 2173, 2187))
artree_plot + geom_cladelabel(node = 1931, label = '1')
offspring(artree_plot$data, 2136)
MRCA(artree_plot$data, c(2136, 2147, 2154, 2159, 2173, 2187, 1935))
artree_plot + geom_cladelabel(node = 212, label = '1')

artree_plot +
  geom_cladelabel(node = 1931, label = '1') +
  geom_cladelabel(node = 1935, label = 'Nano') +
  geom_cladelabel(node = 3381, label = 'Eury') +
  geom_cladelabel(node = 2199, label = 'Halo') +
  geom_cladelabel(node = 2669, label = 'Thermo') +
  geom_cladelabel(node = 3551, label = 'Crena') +
  geom_cladelabel(node = 3541, label = 'Hydro') +
  geom_cladelabel(node = 3549, label = 'Hado') +
  geom_cladelabel(node = 3843, label = 'Asgard')

# leftmost: 2136, 2147, 2154, 2159, 2173, 2187 + (Huber , EX..) --> collapse
# then 1935 -> Nano
# then Asgard --> collapse
# then 3551 -> Crena
# then Hado, Hydro -> collapse 
# then 3381 -> Eury
# then 2199 -> Halo
# Finally 2669 -> Thermoplasmota (-> Highlight Poseido?)
# Needs fixing with collapsed areas not being annotated properly

artree_plot_annot = 
  artree_plot +
  geom_cladelabel(node=1935, label='(c) Nano', fontsize=4) + #Nanoarchaeia 
  geom_cladelabel(node=3551, label='(p) Cren', fontsize=4) + #Crenarchaeota
  geom_cladelabel(node=3381, label='(p) Eury',fontsize=4) + #Euryarchaeota
  geom_cladelabel(node=2199, label='(p) Halo', fontsize=4) + #Halobacterota
  geom_cladelabel(node=2672, label='(c) Poseidoniia', fontsize=4) 
artree_plot_annot

ggsave("data/processed/figures/Figure-2/figure-2-archaeal_backbone_annot.pdf", artree_plot_annot, width = 90, height = 90, units = 'mm')

# Prepare outer barplots -----------------------------------------------------------------

genomes_colors = c("% MAGs" = "#17becf", "% SAGs" = "#dc8b39", "% REFs" = "#7f7f7f")

backbone = artree_plot$data %>% filter(isTip)
circle_ratio = 75/360
total_elt = nrow(backbone) * (1/circle_ratio)

max(artree_plot$data$y, na.rm = T)/250

# BGC Barplot ----------------------------------------------------------------------------

table_figure_bgcs_melt = artree_plot$data %>%
  left_join(Species_biosynth %>% filter(`# Archaea` > 0),
            by = c('label' = 'Representative Genome', '# Genomes' = '# Genomes')) %>%
  filter(isTip) %>%
  mutate(y_r = DescTools::RoundTo(y, 25)) %>%
  group_by(y_r) %>%
  summarize(x = 0,
            y = mean(y),
            `# Biosynthetic Regions` = max_bgcs(`# Biosynthetic Regions`, `# Biosynthetic Products`),
            `RiPPs (Ribosomal Natural Products)` = max_bgcs(`RiPPs (Ribosomal Natural Products)`, `# Biosynthetic Products`),
            `Non-Ribosomal Peptide Synthases` = max_bgcs(`Non-Ribosomal Peptide Synthases`, `# Biosynthetic Products`),
            `Type I Polyketide Synthases` = max_bgcs(`Type I Polyketide Synthases`, `# Biosynthetic Products`),
            `Type II & III Polyketide Synthases` = max_bgcs(`Type II & III Polyketide Synthases`, `# Biosynthetic Products`),
            Terpenes = max_bgcs(Terpenes, `# Biosynthetic Products`),
            Other = max_bgcs(Other, `# Biosynthetic Products`),
            `# Biosynthetic Products` = max_bgcs(`# Biosynthetic Products`, `# Biosynthetic Products`),
            labels = paste(label, collapse = "|")) %>%
  gather(key = Type, value = `Max. Potential`,
         -c(x, y, y_r, labels, `# Biosynthetic Regions`, `# Biosynthetic Products`)) %>%
  mutate(`Max. Potential` = as.numeric(`Max. Potential`)) %>%
  mutate(`Max. Potential` = `Max. Potential` / `# Biosynthetic Products` * `# Biosynthetic Regions`) %>% # Scale by BGCs
  mutate(Type = factor(Type, levels = rev(c("Other", "Terpenes", "Type II & III Polyketide Synthases", "Type I Polyketide Synthases", "Non-Ribosomal Peptide Synthases", "RiPPs (Ribosomal Natural Products)"))),
         `Max. Potential` = ifelse(is.na(`Max. Potential`) & Type == "Other", 0, `Max. Potential`)) 

outer_layer_bgcs = table_figure_bgcs_melt  %>%
  ggplot() +
  geom_col(aes(x = y, fill = Type, y = `Max. Potential`), width = 21) +
  geom_vline(xintercept = 0, size = 0.35*size_converter) +
  geom_vline(xintercept = max(table_figure_bgcs_melt$y, na.rm=T) +  1, size = 0.35*size_converter) +
  scale_fill_manual(values = bgc_colors_serina_old) +
  theme_minimal() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 0, 0, 0, 'mm'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.35*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_polar(start = 1.5*pi, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = c(0, 2, 5, 10, 20, 30), limits = c(-110,35))
outer_layer_bgcs
ggsave("data/processed/figures/Figure-2/figure-2-archaeal_outer_bgcs.pdf", outer_layer_bgcs, width = 140, height = 140, units = 'mm')

# Genome Type Barplot --------------------------------------------------------------------

table_figure_genomes_melt = artree_plot$data %>%
  left_join(Species_biosynth %>% filter(`# Archaea` > 0),
            by = c('label' = 'Representative Genome', '# Genomes' = '# Genomes')) %>%
  filter(isTip) %>%
  mutate(y_r = DescTools::RoundTo(y, 14)) %>%
  group_by(y_r) %>%
  summarize(x = 0,
            y = mean(y),
            `% MAGs` = sum(`# MAGS`, na.rm = T)/sum(`# Genomes`, na.rm = T),
            `% SAGs` = sum(`# SAGS`, na.rm = T)/sum(`# Genomes`, na.rm = T),
            `% REFs` = sum(`# REFG`, na.rm = T)/sum(`# Genomes`, na.rm = T)) %>%
  gather(key = Type, value = `Perc. of Genome Type`,
         -c(x, y, y_r)) %>%
  mutate(Type = factor(Type, levels = rev(c("% MAGs", "% SAGs", "% REFs"))),
         height = 0.9) # 1 makes things buggy..

outer_layer_genomes = table_figure_genomes_melt  %>%
  mutate(`Perc. of Genome Type` = ifelse(`Perc. of Genome Type` == 0, NA, `Perc. of Genome Type`)) %>%
  ggplot() +
  geom_col(aes(x = y, fill = Type, y = height, alpha = `Perc. of Genome Type`), width = 11.3) +
  geom_vline(xintercept = 0, size = 0.35*size_converter) +
  geom_vline(xintercept = max(table_figure_genomes_melt$y) +  1, size = 0.35*size_converter) +
  scale_fill_manual(values = genomes_colors) +
  scale_alpha_continuous(range = c(.1, 1), na.value = 0) +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_polar(start = pi, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = c(0, 2.7), limits = c(-95,20))
outer_layer_genomes

ggsave("data/processed/figures/Figure-S3-GTDB_tree_genome_type/Figure-S3-archaeal_outer_type.pdf", outer_layer_genomes, width = 232.5, height = 232.5, units = 'mm')

