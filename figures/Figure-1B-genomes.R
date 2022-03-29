# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 1B - Genome descr. ==============================================================

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('/Users/paolil/polybox/PhD/Exploratorium/sushilab-colors/palettes-paoli.R')

# Load and prepare data ------------------------------------------------------------------

summary_table = load_prettified_summary()

# Utilitary function ---------------------------------------------------------------------

name_formatting <- function(table){
  table = table %>%
    mutate(Dataset = ifelse(Dataset %in% c("TARA", "GORG", "MARD"), c("TARA" = "Delmont 2018 (MAGs)", "GORG" = "GORG (SAGs)", "MARD" = "MarDB")[Dataset], Dataset)) %>%
    mutate(Dataset =
             factor(
               Dataset,
               levels = c("Tara Oceans", "Malaspina", "Biogeotraces", "Bermuda-Atlantic Time-Series", "Hawaiian Ocean Time-Series", "Delmont 2018 (MAGs)", "GORG (SAGs)", "MarDB")
             ))
  if ("Genome Id" %in% names(table)){
    table = mutate(table, genome_type = ifelse(grepl("^MARD_", `Genome Id`), ifelse(grepl("_SAGS_", `Genome Id`), "\n(SAGs)", "\n(Ref. Genomes)"), ""))
  }
  table = table %>%
    mutate(x = paste0(Dataset, genome_type)) %>%
    mutate(x = gsub(" Time-Series", "\nTime-Series", x)) %>%
    mutate(x = gsub("Bermuda-Atlantic", "Bermuda\nAtlantic", x)) %>%
    mutate(x = gsub("Hawaiian Ocean", "Hawaiian\nOcean", x)) %>%
    mutate(x = gsub(" \\(MAGs)", "\n(MAGs)", x)) %>%
    mutate(x = gsub(" \\(SAGs)", "\n(SAGs)", x)) %>%
    mutate(x = factor(x, levels = c("Tara Oceans", "Malaspina", "Biogeotraces", "Hawaiian\nOcean\nTime-Series", "Bermuda\nAtlantic\nTime-Series", "Delmont 2018\n(MAGs)", "GORG\n(SAGs)", "MarDB\n(Ref. Genomes)", "MarDB\n(SAGs)")))
  return(table)
}

# plot quality barplots ------------------------------------------------------------------

Figure_1B_table = 
  table =summary_table %>%
  mutate(Dataset = ifelse(is.na(dataset), gsub("_.*", "", `Genome Id`), dataset)) %>%
  select(`Genome Id`,
         `Mean Completeness`,
         `Mean Contamination`,
         `Q-score`,
         `High Quality`,
         `Good Quality`,
         `Medium Quality`,
         `Low Quality`,
         `Sample`,
         `Dataset`,
         `size fraction`,
         `SpecI Dereplication Cluster`,
         `dRep Dereplication Cluster`,
         `Representative Genome`,
         `GTDB Taxonomy`) %>%
  mutate(Quality = ifelse(`High Quality`, "High Quality", NA)) %>%
  mutate(Quality = ifelse(`Good Quality`, "Good Quality", Quality)) %>%
  mutate(Quality = ifelse(`Medium Quality`, "Medium Quality", Quality)) %>%
  mutate(Quality = ifelse(`Low Quality`, "Low Quality", Quality)) %>%
  mutate(Quality = factor(Quality, levels = rev(c("High Quality", "Good Quality", "Medium Quality", "Low Quality")))) %>%
  name_formatting(.) %>%
  mutate(x = factor(x, levels = levels(x))) %>%
  mutate(facet = ifelse(grepl("_METAG_", `Genome Id`) & !grepl("Delmont", Dataset), "Reconstructed Genomes", "External Genomes")) %>%
  mutate(facet = factor(facet, levels = c("Reconstructed Genomes", "External Genomes")))

names(dataset_colors) = c("Tara Oceans", "Malaspina", "Biogeotraces", "Bermuda-Atlantic Time-Series", "Hawaiian Ocean Time-Series", "GORG (SAGs)", "Delmont 2018 (MAGs)", "MarDB")

Figure_1B = ggplot(Figure_1B_table) +
  geom_bar(aes(x = x, fill = Dataset, alpha = Quality), width=0.8) +
  scale_fill_manual(values = dataset_colors) +
  scale_alpha_manual(values = c(0.4, 0.6, 0.8, 1)) +
  facet_grid(cols = vars(facet), scale = "free_x", space = "free_x") +
  ylab("# Genomes") +
  scale_y_continuous(breaks = c(1000, 5000, 10000, 20000), labels = c("1k", "5k", "10k", "20k"), minor_breaks = c(0, 500, 3000, 7500, 15000), expand = c(0,100), limits = c(0,20000)) +
  theme_bw() +
  theme_cell +
  theme(rect = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = c(.5, .6),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0, 'mm'),
        legend.spacing = unit(0, "mm"),
        legend.box = "horizontal",
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6, margin = margin(0,0,0,0, 'mm')),
        axis.title.y = element_text(size = 6, margin = margin(0,0,0,0, 'mm')),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(1, 'mm'),
        strip.background = element_rect(colour = NA, fill = "black"),
        strip.text = element_text(size = 6, colour = "white", face = "bold", margin = margin(.5,0,.5,0, "mm")))

Figure_1B

ggsave(paste0(figures_path_proj, "Figure-1/Figure-1B.pdf"), Figure_1B, width = 69, height = 33, units = col_unit)



