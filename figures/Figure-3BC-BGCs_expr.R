# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 3BC - Ca. Eudoremicrobium BGC expr. =============================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(patchwork)
library(viridis)
#source('code/figures/figure4D-dataprep.R') # expects it to have been run

source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Data preparation =======================================================================

# functional annotations & metadata ------------------------------------------------------

formatted_metadata = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")
MGs = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-fetchMGs")
BGCs = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-antismash")
metat_featurecounts_summary_processed = read_tsv("data/raw/marine_eremios/marine_eremios-featurecounts-summary_PIAMPJPB_metat-processed.tsv.gz", guess_max = 200000)
samples_clustering = read_tsv("data/raw/marine_eremios/PIAMPJPB-transcriptome_clusters-sample_order_membership.tsv")
serina_bgcs <- googlesheets4::read_sheet("1IqyKwihHOYAG8jY4VX6U8znPViSxNAg9jOSEE57djKM", sheet = "Biosynthetic potential analysis")

# Add the BGC heatmap from Serina---------------------------------------------------------

## Sort and create factor for ggplot2 tick label ordering
dat <- serina_bgcs %>%
  group_by(product_short) %>%
  mutate(core_index = n()) %>%
  ungroup() %>%
  arrange(core_index, desc(order_median_expression_data))

dat$product_short <- factor(dat$product_short, levels = unique(dat$product_short))
dat$genome_short <- factor(dat$genome_short, levels = c("A1", "A2", "B1", "C1", "C2"))

## Set color palette and plot heatmap
pal <- c("#A6D854", "#FFD92F", "#8DA0CB", "#E78AC3", "#FC8D62", "#66C2A5")

ggplot(dat, aes(x = genome_short, 
                y = product_short)) +
  geom_tile(aes(fill = product_class), color = "white", size = 0.5) +
  scale_fill_manual(values = pal) +
  coord_fixed() +
  theme_bw() +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_blank(), 
        legend.position = "none",
        panel.border = element_rect(color = NA, fill= NA),
        axis.text.y = element_text(size = 6, hjust = 1),
        axis.text.x = element_text(size = 6),
        panel.grid.major = element_blank(),
        plot.margin = margin(0,0,0,0),
        axis.ticks.length = unit(0.5, "mm"))

ggsave(filename = "data/processed/figures/Figure-4/Figure-4B.raw.pdf",
       width = 45, height = 80, units = "mm")

# Metat abundances featurecounts ---------------------------------------------------------

metat_featurecounts_summary_processed %>%
  filter(Gene %in% (MGs %>% filter(`mOTUs-10`) %>% pull(bestMGs))) %>%
  ggplot() +
  geom_density(aes(x = value))

metat_featurecounts_summary_spread = metat_featurecounts_summary_processed %>%
  select(Gene, material, value) %>% 
  filter(!is.na(material)) %>%
  spread(Gene, value)

# BGC Analysis ===========================================================================

metat_featurecounts_summary_clustered = metat_featurecounts_summary_processed %>%
  left_join(samples_clustering) %>%
  mutate(material = factor(material, levels = unique(samples_clustering$material))) %>%
  arrange(metat_featurecounts_summary_clustered)

antismash_biosynthetic = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (BGCs %>% filter(biosynthetic_gene | biosynthetic_add_gene) %>% pull(gene))) %>%
  left_join(BGCs, by = c("Gene" = "gene")) %>%
  group_by(biosynthetic_region, material) %>%
  summarize(n_biosynthetic_genes = n(),
            n_inserts = sum(inserts),
            median_inserts = median(inserts),
            value = median(value),
            cluster = unique(cluster)) %>%
  group_by(biosynthetic_region) %>%
  mutate(pvalue = kruskal.test(value ~ cluster)$p.value) %>%
  ungroup()

adjust_pvals = p.adjust(antismash_biosynthetic %>% filter(!duplicated(biosynthetic_region)) %>% pull(pvalue), method = "fdr")
names(adjust_pvals) = antismash_biosynthetic %>% filter(!duplicated(biosynthetic_region)) %>% pull(pvalue) %>% as.character()

antismash_biosynthetic = antismash_biosynthetic %>%
  mutate(pvalue_adj = adjust_pvals[as.character(pvalue)]) %>%
  mutate(significance = ifelse(pvalue_adj <= 0.05, "*", "")) %>%
  mutate(significance = ifelse(pvalue_adj <= 0.01, "**", significance)) %>%
  mutate(significance = ifelse(pvalue_adj <= 0.001, "***", significance)) %>%
  mutate(significance = ifelse(duplicated(biosynthetic_region), NA, significance))


antismash_biosynthetic %>%
  mutate(shape = ifelse(median_inserts > 0, "Evidence for expression", "No evidence")) %>%# View()
  mutate(color = ifelse(abs(value) >= 2, "#f28e2c", "#4e79a7")) %>%
  mutate(biosynthetic_region = factor(biosynthetic_region, levels = unique(dat$region_id))) %>%
  ggplot() +
  geom_boxplot(aes(x = biosynthetic_region, y = value), outlier.shape = NA, size = 0.2) +
  geom_jitter(aes(x = biosynthetic_region, y = value, shape = shape, color = color), width = 0.2, height = 0, alpha = .7, size = 1) +
  geom_text(aes(x = biosynthetic_region, y = max(value) + 1, label = significance), vjust = .75, size = unit(3, 'pt')) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_identity() +
  ylab("Median expression rate per sample (log2)\n(0 Indicates values similar to housekeeping genes)") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = unit(.3, 'mm')),
    legend.title = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    plot.margin = margin(0,0,0,0),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    line = element_line(size = unit(0.3, 'mm'))
  ) +
  coord_flip()

ggsave(filename = "data/processed/figures/Figure-4/Figure-4C.raw.pdf",
       width = 40, height = 79, units = "mm")
