# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 2C - Bacterial Tree =============================================================

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

bctree = "data/raw/go_microbiomics/gtdb/go_microbiomics-gtdb-4Chris/gtdbtk.bac120.classify.tree"

# Assumes data prep is loaded

# Bacteria ===============================================================================


# Inner plot -----------------------------------------------------------------------------

bctrnw <- drop.tip(read.newick(bctree),
                   genomes_summary %>%
                     filter(!`Representative Genome` |  `GTDB Taxonomy` == "N/A") %>% pull(`Genome Id`)
)
assertthat::assert_that(!any(Species_biosynth %>% pull(`# Genomes`) %>% is.na()))

summary(Species_biosynth$`# Genomes`)
summary(log(Species_biosynth$`# Genomes`))

# Redifine the color after the collapse:
max_genome_value = 1159 # from max(bctree_collapsable_level_nodes_raw$value, na.rm = T)
scaled_colors = 
  var_scale_colors(log(max_genome_value), scaling_factor = 10, fun = "viridis", type = "D", begin = 0.1)

bctrnw_metadata = 
  as_tibble(bctrnw) %>% 
  left_join(Species_biosynth %>% select(label = `Representative Genome`, `# Genomes`) %>% mutate(user_genome = TRUE)) %>%
  mutate(user_genome = ifelse(is.na(user_genome), FALSE, user_genome),
         color = ifelse(is.na(`# Genomes`), "lightgrey", scaled_colors[round(log(`# Genomes`)*10) + 1]))

bctree_plot = ggtree(bctrnw, layout="fan", open.angle=90, branch.length = "none", size = 0.35, alpha = 1) %<+%
  #bctree_plot = ggtree(bctrnw, layout="fan", open.angle=180, size = 0.35, alpha = 1) %<+%
  bctrnw_metadata + # add metadata to the tree data
  aes(color = color) +
  scale_color_identity() +
  coord_polar(theta = 'y', start = 0, direction = -1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))
bctree_plot$data

# Collapsing the clades
# /!\ Uncomment this line and comment the next if you need to update the backbone!!
#bctree_collapsable_backbone_nodes = identify_collapsable_backbone(bctree_plot, force = T)
bctree_collapsable_backbone_nodes = identify_collapsable_backbone(bctree_plot)
bctree_plot_collapsed_backbone = collapse_backbone(bctree_plot, bctree_collapsable_backbone_nodes)

bctree_collapsable_level_nodes_raw = identify_collapsable_levels(bctree_plot_collapsed_backbone, "# Genomes", perc_to_collapse = 0.15)
assertthat::are_equal(max(bctree_collapsable_level_nodes_raw$value, na.rm = T), max_genome_value)

bctree_collapsable_level_nodes = bctree_collapsable_level_nodes_raw %>% mutate(value = round(log(value)*10) + 1)
bctree_plot_collapsed = collapse_level(bctree_plot_collapsed_backbone, bctree_collapsable_level_nodes, scaled_colors)

bctree_plot_collapsed$data %>% pull(color) %>% table() %>% sort()

ggsave("data/processed/figures/Figure-2/figure-2-bacterial_backbone.pdf", bctree_plot_collapsed, width = 100, height = 100, units = 'mm')

bctree_plot_collapsed_no_backbone = remove_backbone_layers(bctree_plot_collapsed) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave("data/processed/figures/Figure-2/figure-2-bacterial_backbone_coloronly.pdf", bctree_plot_collapsed_no_backbone, width = 100, height = 100, units = 'mm')

# Add labels -----------------------------------------------------------------------------

# Identify the taxa to display
bctree_plot_annot = bctree_plot
for (n in bctree_plot$data %>% filter(grepl("p__", label)) %>% pull(node)){
  l = bctree_plot_annot$data %>% filter(node == n) %>% pull(label)
  bctree_plot_annot = 
    bctree_plot_annot + geom_cladelabel(node=n, label=l, fontsize=4) 
}
bctree_plot_annot
ggsave("data/processed/recent/figure-2-bacterial_backbone_w_annot.pdf", bctree_plot_annot, width = two_col, height = two_col, units = col_unit)

bctree_plot$data %>% filter(grepl("p__", label)) %>% arrange(y) %>% View()

# Annotations downward
bctree_plot$data %>% filter(grepl("p__Proteobacteria", label))
bctree_plot$data %>% filter(grepl("Acidobacteriota", label))
bctree_plot$data %>% filter(grepl("Myxococcota", label))
bctree_plot$data %>% filter(grepl("Rhodobacteraceae", label))
bctree_plot$data %>% filter(grepl("Gammaproteo", label))
bctree_plot$data %>% offspring(48138) %>% View()
bctree_plot_annot_downward = bctree_plot + 
  geom_cladelabel(node=48138, label="", color = "blue") +
  geom_cladelabel(node=54254, label="", color = "red") +
  geom_cladelabel(node=48140, label="", color = "green") +
  geom_cladelabel(node=57968, label="Pelagibacterales", fontsize=4, hjust = 1) +
  geom_cladelabel(node=55374, label="Rhodobacteraceae", fontsize=4, hjust = 1) +
  geom_cladelabel(node=61608, label="Acido", fontsize=4, hjust = 1, color = "orange")  +
  geom_cladelabel(node=61202, label="Myxo", fontsize=4, hjust = 1, color = "steelblue") 
bctree_plot$data %>% filter(grepl("Bacteroidota", label))
bctree_plot_annot_downward = 
  bctree_plot_annot_downward + 
  geom_cladelabel(node=43514, label="Bacteroidota", color = "steelblue", hjust = 1)
bctree_plot$data %>% filter(grepl("Planctomycetota", label))
bctree_plot_annot_downward = 
  bctree_plot_annot_downward + 
  geom_cladelabel(node=42882, label="Plancto", color = "orange", hjust = 1)
bctree_plot$data %>% filter(grepl("Chloroflexota", label))
bctree_plot_annot_downward = 
  bctree_plot_annot_downward + 
  geom_cladelabel(node=41379, label="Chloroflex", color = "darkgreen", hjust = 1)
ggsave("data/processed/recent/figure-2-bacterial_backbone_w_annot_downward.pdf", bctree_plot_annot_downward, width = two_col, height = two_col, units = col_unit)


# Annotations upward
bctree_plot_annot_upward = bctree_plot
for (n in bctree_plot$data %>% filter(grepl("Firmicutes", label)) %>% pull(node)){
  l = bctree_plot_annot_upward$data %>% filter(node == n) %>% pull(label)
  bctree_plot_annot_upward = 
    bctree_plot_annot_upward + geom_cladelabel(node=n, label="", fontsize=4, hjust = 1) 
}
bctree_plot$data %>% filter(grepl("Cyanobacteria", label))
bctree_plot_annot_upward = bctree_plot_annot_upward + geom_cladelabel(node=35713, label="Cyano", fontsize=4, color = "green", hjust = 1) 
bctree_plot$data %>% filter(grepl("Actinobact", label))
bctree_plot_annot_upward = bctree_plot_annot_upward + geom_cladelabel(node=36842, label="Actino", fontsize=4, color = "steelblue", hjust = 1) 
bctree_plot$data %>% filter(grepl("Patescibacteria", label))
bctree_plot_annot_upward = bctree_plot_annot_upward + geom_cladelabel(node=40285, label="CPR", fontsize=4, color = "blue", hjust = 1) 
ggsave("data/processed/recent/figure-2-bacterial_backbone_w_annot_upward.pdf", bctree_plot_annot_upward, width = two_col, height = two_col, units = col_unit)

# Prepare outer barplots -----------------------------------------------------------------

genomes_colors = c("% MAGs" = "#17becf", "% SAGs" = "#dc8b39", "% REFs" = "#7f7f7f")

backbone = bctree_plot$data %>% filter(isTip)
circle_ratio = 3/4
total_elt = nrow(backbone) * (1/circle_ratio)

# BGC Barplot ----------------------------------------------------------------------------

table_figure_bgcs_melt = bctree_plot$data %>%
  left_join(Species_biosynth %>% filter(`# Bacteria` > 0),
            by = c('label' = 'Representative Genome', '# Genomes' = '# Genomes')) %>%
  filter(isTip) %>%
  mutate(y_r = DescTools::RoundTo(y, 100)) %>%
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
  geom_col(aes(x = y, fill = Type, y = `Max. Potential`), width = 90) +
  geom_vline(xintercept = 0, size = 0.35*size_converter) +
  geom_vline(xintercept = max(table_figure_bgcs_melt$y, na.rm=T) +  1, size = 0.35*size_converter) +
  scale_fill_manual(values = bgc_colors_serina_old) +
  theme_minimal() +
  theme(rect = element_blank(),
        plot.margin = margin(0, 0, 0, 0, 'mm'),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.35*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = c(0, 2, 5, 10, 15, 20, 30), limits = c(-110,35))
outer_layer_bgcs
ggsave("data/processed/figures/Figure-2/figure-2-bacterial_outer_bgcs.pdf", outer_layer_bgcs, width = 140, height = 140, units = 'mm')

# Genome Type Barplot --------------------------------------------------------------------

table_figure_genomes_melt = bctree_plot$data %>%
  left_join(Species_biosynth %>% filter(`# Bacteria` > 0),
            by = c('label' = 'Representative Genome', '# Genomes' = '# Genomes')) %>%
  filter(isTip) %>%
  mutate(y_r = DescTools::RoundTo(y, 100)) %>%
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
  geom_col(aes(x = y, fill = Type, y = height, alpha = `Perc. of Genome Type`), width = 90) +
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
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = c(0, 2.7), limits = c(-95,20))
outer_layer_genomes

ggsave("data/processed/figures/Figure-S3-GTDB_tree_genome_type/Figure-S3-bacterial_outer_type.pdf", outer_layer_genomes, width = 232.5, height = 232.5, units = 'mm')

# Plot legend ----------------------------------------------------------------------------

legend_tbl = tibble(
  color = scaled_colors,
  scaled_value = 1:length(scaled_colors)) %>% 
  filter(scaled_value/10 == round(scaled_value/10) | scaled_value == 1) %>% # round(log(`# Genomes`)*10) + 1
  mutate(genome_value = round(exp(round((scaled_value - 1)/10)))) %>%
  mutate(color = factor(color, levels = rev(color)))

ggplot(legend_tbl) +
  geom_col(aes(x = 1, y = 1, fill  = color)) +
  scale_fill_identity() +
  theme_void()
