# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 7A - Clustering =====================================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
#source('code/figures/figure4D-dataprep.R')

source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Data preparation =======================================================================

# functional annotations & metadata ------------------------------------------------------

formatted_metadata = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")
metat_featurecounts_summary_processed = read_tsv("data/raw/marine_eremios/marine_eremios-featurecounts-summary_PIAMPJPB_metat-processed.tsv.gz", guess_max = 200000)
MGs = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-fetchMGs")

# Metat abundances featurecounts ---------------------------------------------------------

# Quick diagnostic
metat_featurecounts_summary_processed %>%
  filter(Gene %in% (MGs %>% filter(`mOTUs-10`) %>% pull(bestMGs))) %>%
  ggplot() +
  geom_density(aes(x = value))

metat_featurecounts_summary_spread = metat_featurecounts_summary_processed %>%
  select(Gene, material, value) %>% 
  filter(!is.na(material)) %>%
  spread(Gene, value)

# Gene detection distribution ------------------------------------------------------------

length(unique(metat_featurecounts_summary_processed$Sample))
summary(metat_featurecounts_summary_processed$value)
summary(metat_featurecounts_summary_processed$n_genes_detected)
unique(metat_featurecounts_summary_processed$n_genes_detected) %>% sort

gene_detection = metat_featurecounts_summary_processed %>%
  group_by(Gene) %>%
  summarize(n_sample = sum(inserts > 0)) %>%
  group_by(n_sample) %>%
  summarize(n_genes = n()) 

p1 = gene_detection %>%
  ggplot() +
  geom_col(aes(x = n_sample, y = n_genes)) +
  xlab("Number of samples the genes are detected in") +
  ylab(paste0("Number of genes (out of ", length(unique(metat_featurecounts_summary_processed$Gene)),")")) +
  theme_bw() +
  theme(text = element_text(size = unit(6, 'pt')))

ggsave("data/processed/figures/Figure-S5-GCF_clustering/gene_detection.pdf", p1, width = 70, height = 40, units = "mm")

sample_order = read_tsv("data/raw/marine_eremios/PIAMPJPB-transcriptome_clusters-sample_order_membership.tsv")
p2 = metat_featurecounts_summary_processed %>%
  select(internal_sample_name = Sample, n_genes_detected) %>%
  unique() %>%
  left_join(formatted_metadata %>% select(internal_sample_name, material)) %>%
  left_join(sample_order) %>%
  arrange(cluster) %>%
  ggplot() +
  geom_col(aes(x = material, y = n_genes_detected)) +
  facet_grid(.~cluster, scales = "free_x", space = "free_x") +
  xlab("Metatranscriptomic samples") +
  ylab("Number of genes") +
  theme_bw() +
  theme(text = element_text(size = unit(6, 'pt')),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("data/processed/figures/Figure-S5-GCF_clustering/gene_per_sample.pdf", p2, width = 70, height = 40, units = "mm")


# Genome wide expression patterns - UMAP/HDBSCAN with featurecounts ----------------------

custom.settings = umap::umap.defaults
custom.settings$random_state = 456#26
custom.settings$min_dist = 10**-10 # Minimum output distance: how packed are the points, for viz try not so much, for hdbscan
custom.settings$n_neighbors = 2 # 5 and seed 489 separates 2 clusters
custom.settings$n_components = 2 # How many dimensions in return
custom.settings$knn_repeats = 100 # How many knn starts
custom.settings$n_epochs = 20000 # How many epochs
custom.settings$input = "data"
custom.settings$metric = "euclidean"

metat_featurecounts_summary_umap = umap::umap(metat_featurecounts_summary_spread %>% select(-material),
                                              config = custom.settings)

scores = tibble()
for (k in 2:round(nrow(metat_featurecounts_summary_umap$layout)/2)){
  tmp_hdbscan = dbscan::hdbscan(metat_featurecounts_summary_umap$layout, minPts = k)
  res = tibble(k = k,
               score = sum(tmp_hdbscan$membership_prob),
               n_unassign = sum(tmp_hdbscan$cluster == 0),
               n_clusters = length(unique(tmp_hdbscan$cluster)))
  scores = rbind(scores, res)
  print(paste('For minPts =', res$k, '|', res$n_unassign, 'unassigned,',  'score of', res$score, 'and', res$n_clusters, 'clusters.'))
}

metat_featurecounts_summary_hdbscan = dbscan::hdbscan(metat_featurecounts_summary_umap$layout, minPts = (scores %>% arrange(n_unassign, desc(score)) %>% pull(k))[1])

as_tibble(metat_featurecounts_summary_umap$layout) %>%
  mutate(material = metat_featurecounts_summary_spread$material,
         clusters = as.character(metat_featurecounts_summary_hdbscan$cluster)) %>%
  left_join(formatted_metadata) %>% 
  ggplot() +
  geom_jitter(aes(x = V1, y = V2, color = clusters, fill = clusters, shape = size_fraction), height = 5, width = 5, size = 3, alpha = .8) +
  scale_shape_manual(values = c("0.22-1.6" = 22, "0.22-3" = 21, "5-20" = 23, "0.8->" = 24, "0.8-5" = 25), name = "Size Fraction (um)") +
  xlab("x") +
  ylab("y") +
  theme_bw()

as_tibble(metat_featurecounts_summary_umap$layout) %>%
  mutate(material = metat_featurecounts_summary_spread$material,
         clusters = as.character(metat_featurecounts_summary_hdbscan$cluster)) %>%
  left_join(formatted_metadata) %>% 
  ggplot() +
  geom_jitter(aes(x = V1, y = V2, color = ocean_province, fill = ocean_province, shape = size_fraction), height = 5, width = 5, size = 3, alpha = .8) +
  scale_shape_manual(values = c("0.22-1.6" = 22, "0.22-3" = 21, "5-20" = 23, "0.8->" = 24, "0.8-5" = 25), name = "Size Fraction (um)") +
  scale_color_manual(values = ocean_material2_colors, name = "Ocean Provinces") +
  scale_fill_manual(values = ocean_material2_colors, name = "Ocean Provinces") +
  xlab("x") +
  ylab("y") +
  theme_bw()



vegan::adonis(metat_featurecounts_summary_spread %>% select(-material) ~ 
                as.character(metat_featurecounts_summary_hdbscan$cluster),
              permutations = 10000, by = 'margin',method = "euclidean")

hclust_all_dim = hclust(dist(metat_featurecounts_summary_spread %>% select(-material)))
hclust_all_dim$labels = as.character(metat_featurecounts_summary_spread$material)
plot(hclust_all_dim)

hclust_reduced = hclust(dist(metat_featurecounts_summary_umap$layout))
hclust_reduced$labels = as.character(metat_featurecounts_summary_spread$material)
plot(hclust_reduced)

class(hclust_reduced)
to_tree = ape::as.phylo(hclust_reduced) 
ape::write.tree(phy = to_tree, file="data/raw/marine_eremios/PIAMPJPB-transcriptome_clusters.newick")

tibble(material = metat_featurecounts_summary_spread %>% pull(material),
       cluster = paste("Cluster", metat_featurecounts_summary_hdbscan$cluster)) %>%
  mutate(material = factor(material, levels = hclust_reduced$labels[hclust_reduced$order])) %>%
  arrange(material) %>%
  write_tsv("data/raw/marine_eremios/PIAMPJPB-transcriptome_clusters-sample_order_membership.tsv")
