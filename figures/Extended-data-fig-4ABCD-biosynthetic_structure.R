# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 4A - Structure of the biosynthetic potential ========================

# Libraries ==============================================================================

rm(list = ls())
library(tidyverse)
library(patchwork)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# functions ==============================================================================

source('code/figures/figure2C-functions.R')

# load data ==============================================================================

# get some samples metadata ready
samples_metadata = load_general_metadata() %>%
  mutate(`Internal Sample Name` = ifelse(grepl("TARA_", `Internal Sample Name`), paste0(`Internal Sample Name`, "_METAG"), `Internal Sample Name`))
assertthat::assert_that(length(unique(samples_metadata$Sample)) == 1038)

# load the genome summary (w. motus membership)
genomes_summary = load_prettified_summary() %>%
  left_join(read_tsv("data/raw/go_microbiomics/motus/v1.0/gom_mOTUsv2.mag-memberships.tsv", col_names = c('Genome Id', 'motus_cluster'))) %>%
  mutate(motus_cluster = gsub("ext_mOTU_v26_", "gom_", motus_cluster))

# Load a per bgc summary
bgc_gcf_summary = read_tsv("data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-antismash_bgc_to_gcf-curated-v0.1.tsv", col_names = c('bgc', 'gcf')) %>%
  mutate(gcf = paste0('gcf_', gcf)) %>%
  left_join(read_tsv("data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-antismash_summary-enriched.tsv"), by = c('bgc' = 'region')) %>%
  select(genome, scaffold, bgc, gcf, length, `contig edge`, products, gene = `CDS biosynthetic`) %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, `GTDB Taxonomy`, `dRep Dereplication Cluster`))

# Get the gene information for each bgc
bgc_gcf_gene_summary = bgc_gcf_summary %>%
  separate_rows(gene, sep = ';') %>%
  left_join(read_tsv("data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-biosynth-biosynthgenes-gene_catalog_membership.tsv"))

# Load the GCF class ratio
gcf_classes_distrib = read_tsv("data/raw/go_microbiomics/bgcs/v0.1/gcf_classes_t0.2.tsv") %>%
  mutate(gcf = paste0('gcf_', gcf)) %>%
  gather(key = class, value = freq, -gcf)
length(unique(gcf_classes_distrib$gcf))

# Get the processed abundance data for the bgcs
#gcf_summary_processed = read_tsv("data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-gcf-abundances_processed-archive/go_microbiomics-integrated-cpl50_ctn10-gcf-abundances_processed.tsv.gz")
#gcf_summary_processed = read_tsv("data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-gcf-abundances_processed-archive/go_microbiomics-integrated-cpl50_ctn10-gcf-abundances_processed.tsv.new.gz")
gcf_summary_processed = read_tsv("data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-gcf-abundances_processed-0.1.tsv.gz")
gcf_summary_processed = gcf_summary_processed %>%
  filter(!grepl("_T", sample))
assertthat::assert_that(length(unique(gcf_summary_processed$sample)) == 1038)
assertthat::assert_that(sum(gcf_summary_processed %>% filter(is.na(gcf)) %>% pull(value)) == 0)

# get the motus abundances for genome length estimates
motus_abundances = read_tsv("data/raw/go_microbiomics/motus/v1.0/gom.motus", skip = 2) %>%
  select(all_of(c("#consensus_taxonomy", samples_metadata$`Internal Sample Name`))) %>%
  gather(key = sample, value = motus_count, -`#consensus_taxonomy`) %>%
  mutate(taxonomy = gsub(" \\[.*", "", `#consensus_taxonomy`),
         motus_cluster = gsub(".*\\[|\\]", "", `#consensus_taxonomy`))

# genome sizes
motus_lengths = genomes_summary %>%
  filter(`Mean Completeness` >= 70) %>%
  group_by(motus_cluster) %>%
  summarize(genome_size = mean(100/`Mean Completeness` * `Genome size`))

sample_avg_size = motus_abundances %>%
  left_join(motus_lengths) %>%
  filter(!is.na(genome_size)) %>%
  group_by(sample) %>%
  mutate(motus_count_rel = motus_count / sum(motus_count)) %>%
  ungroup() %>%
  filter(motus_count > 0) %>%
  group_by(sample) %>%
  summarize(avg_size = sum(motus_count_rel * genome_size))

motus_abundances %>%
  group_by(sample) %>%
  summarize(n = sum(motus_count)) %>%
  summary

# Compute distances ======================================================================

gcf_summary_eucl_spread = gcf_summary_processed %>%
  spread(gcf, value)

gcf_summary_eucl_dist = dist(select(gcf_summary_eucl_spread, -sample))

gcf_summary_eucl_dist_matrix = as.matrix(gcf_summary_eucl_dist)
gcf_summary_eucl_dist_matrix_na = gcf_summary_eucl_dist_matrix
gcf_summary_eucl_dist_matrix_na[upper.tri(gcf_summary_eucl_dist_matrix_na, diag = T)] = NA
gcf_summary_eucl_dist_tbl = as_tibble(gcf_summary_eucl_dist_matrix_na)
gcf_summary_eucl_dist_tbl %>%
  gather(key = sample, value = dist) %>%
  filter(!is.na(dist)) %>%
  ggplot() + 
  geom_density(aes(x = dist))

# UMAP dimension reduction & HDBSCAN clustering based on Euclidean dist ==================

eucl.settings = umap::umap.defaults
eucl.settings$input = "dist"
eucl.settings$random_state = NA
eucl.settings$transform_state = NA
eucl.settings$min_dist = 10**-10 # Minimum output distance: how packed are the points, for clustering small values are good
eucl.settings$spread = 1
eucl.settings$n_neighbors = 30 # 
eucl.settings$n_components = 2 # How many dimensions in return
eucl.settings$knn_repeats = 1000 # How many knn starts
eucl.settings$n_epochs = 10000 # How many epochs
#eucl.settings$metric = "cosine"

# Get the umap embedding
rm(.Random.seed, envir=globalenv())
gcf_summary_eucl_umap = umap::umap(gcf_summary_eucl_dist_matrix, config = eucl.settings)
gcf_summary_eucl_umap$config
gcf_summary_eucl_umap_layout = as_tibble(gcf_summary_eucl_umap$layout)

# Optimize clustering
gcf_summary_eucl_hdbscan = hdbscan_auto_optimum(gcf_summary_eucl_umap_layout, hard_min = 30, hard_max = 80)
gcf_summary_eucl_hdbscan = dbscan::hdbscan(gcf_summary_eucl_umap_layout, minPts = 53)

# Look at the results
plot_umap_hdbscan_helper(gcf_summary_eucl_umap_layout, gcf_summary_eucl_hdbscan, gcf_summary_eucl_spread, samples_metadata) %>% 
  ggplot() +
  geom_point(aes(x = V1, y = V2, shape = depth_layer, color = ocean_province, fill = ocean_province), size = 2, alpha = .7) +
  scale_shape_manual(values = c("EPI" = 21, "MES" = 22, "BAT" = 24, "ABY" = 25), name = "Depth Layers") +
  scale_color_manual(values = c(ocean_material2_colors, "Background" = "grey90"), name = "Ocean Provinces") +
  scale_fill_manual(values = c(ocean_material2_colors, "Background" = "grey90"), name = "Ocean Provinces") +
  xlab("x") +
  ylab("y") +
  theme_bw() +
  facet_wrap(~clusters, nrow = 1)

# Prepare/save embedding
umap_eucl_log = paste(c("dist=", gcf_summary_eucl_umap$config$min_dist,
                        "_neighbors=", gcf_summary_eucl_umap$config$n_neighbors,
                        "_seed=", gcf_summary_eucl_umap$config$random_state,
                        "_transform=", gcf_summary_eucl_umap$config$transform_state,
                        "_hdbscan=", length(unique(gcf_summary_eucl_hdbscan$cluster))),
                      collapse = "")
write_tsv(gcf_summary_eucl_umap_layout, paste(c("data/raw/go_microbiomics/bgcs/umap_embedding_eucl", umap_eucl_log, "tsv"), collapse = "."))

# Load pre-computed embedding
gcf_summary_eucl_umap_layout = read_tsv("data/raw/go_microbiomics/bgcs/umap_embedding_eucl.dist=1e-10_neighbors=30_seed=879810110_transform=NA_hdbscan=3.tsv")
gcf_summary_eucl_hdbscan = hdbscan_auto_optimum(gcf_summary_eucl_umap_layout, hard_min = 30, hard_max = 80)

rename_eucl_clusters = c("3" = "1", "2" = "2", "1" = "3")

# Evaluation:
vegan::adonis(gcf_summary_eucl_dist ~ 
                as.character(gcf_summary_eucl_hdbscan$cluster), by = 'margin')

gcf_summary_eucl_rev_eco_summary = gcf_summary_eucl_umap_layout %>%
  mutate(clusters = paste0("Cluster ", rename_eucl_clusters[as.character(gcf_summary_eucl_hdbscan$cluster)])) %>%
  mutate(sample = gcf_summary_eucl_spread$sample) %>%
  left_join(samples_metadata %>% rename(sample = `Internal Sample Name`)) %>%
  mutate(depth_layer = factor(depth_layer, levels = c("EPI", "MES", "BAT", "ABY"))) %>%
  mutate(fraction_facet = ifelse(size_fraction %in% c("<-0.22", "0.1-0.22"), "<0.2", NA)) %>%
  mutate(fraction_facet = ifelse(size_fraction %in% c("0.22-0.45", "0.45-0.8", "0.2-0.8"), "0.2-0.8", fraction_facet)) %>%
  mutate(fraction_facet = ifelse(size_fraction %in% c("0.22-1.6", "0.22-3"), "0.2-3", fraction_facet)) %>%
  mutate(fraction_facet = ifelse(size_fraction == "0.8-20", "0.8-20", fraction_facet)) %>%
  mutate(fraction_facet = ifelse(size_fraction == "0.2<-", ">0.2", fraction_facet)) %>%
  mutate(fraction_facet = factor(fraction_facet, levels = c("<0.2", "0.2-0.8", "0.2-3", "0.8-20", ">0.2")))


# Test for balanced permanova
vegan::betadisper(gcf_summary_eucl_dist, gcf_summary_eucl_rev_eco_summary$clusters)
table(gcf_summary_eucl_rev_eco_summary$clusters)
gcf_summary_eucl_balanced_permanovas_R2 = balanced_anova_R2_estimate(gcf_summary_eucl_rev_eco_summary, gcf_summary_eucl_dist_tbl, gcf_summary_eucl_spread, 75)
gcf_summary_eucl_balanced_permanovas_R2 %>%
  ggplot() +
  geom_density(aes(x = R2)) +
  geom_vline(xintercept = median(gcf_summary_eucl_balanced_permanovas_R2$R2))


# Prepare data for figure ================================================================

# Cluster composition
go_secmet_t1 = gcf_summary_eucl_rev_eco_summary %>%
  group_by(clusters, depth_layer, fraction_facet, ocean_province) %>%
  summarize(n = n())

# GCF abundances, richnesses by class
go_secmet_t2 = gcf_summary_processed %>%
  left_join(gcf_summary_eucl_rev_eco_summary) %>%
  full_join(gcf_classes_distrib %>% filter(gcf %in% unique(gcf_summary_processed$gcf))) %>% 
  mutate(class_value = freq*value) %>%
  group_by(sample, class) %>%
  summarize(abundance = sum(class_value),
            richness = sum(class_value > 0),
            diversity = vegan::diversity(class_value),
            clusters = unique(clusters)) %>%
  ungroup() %>%
  mutate(class = factor(class, levels = c("Non-Ribosomal Peptide Synthetases",
                                          "Type I Polyketide Synthases", 
                                          "Type II/III Polyketide Synthases",
                                          "RiPPs (Ribosomal Natural Products)",
                                          "Terpenes",
                                          "Other")))
names(bgc_colors_serina) = c("Type I Polyketide Synthases", 
                             "RiPPs (Ribosomal Natural Products)",
                             "Terpenes",
                             "Type II/III Polyketide Synthases",
                             "Non-Ribosomal Peptide Synthetases",
                             "Other")

go_secmet_t3 = gcf_summary_eucl_rev_eco_summary %>%
  left_join(sample_avg_size)

# Figure for paper =======================================================================
# A four quadrant 69 by 69 mm

# Panel 1, UMAP --------------------------------------------------------------------------

figure2C1 = gcf_summary_eucl_umap_layout %>% 
  mutate(clusters = paste0("Cluster ", gcf_summary_eucl_hdbscan$cluster)) %>%
  mutate(sample = gcf_summary_eucl_spread$sample) %>%
  left_join(samples_metadata %>% rename(sample = `Internal Sample Name`)) %>%
  ggplot() +  geom_point(aes(x = V1, y = V2, shape = depth_layer, color = ocean_province, fill = ocean_province), size = .5, stroke = .2, alpha = .6) +
  scale_shape_manual(values = c("EPI" = 21, "MES" = 22, "BAT" = 24, "ABY" = 25), name = "Depth Layers") +
  scale_color_manual(values = c(ocean_material2_colors, "Background" = "grey90"), name = "Ocean Provinces") +
  scale_fill_manual(values = c(ocean_material2_colors, "Background" = "grey90"), name = "Ocean Provinces") +
  xlab("x") +
  ylab("y") +
  coord_fixed() +
  theme_bw() +
  theme(rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.ticks.length = unit(0.5, 'mm'),
        axis.ticks = element_line(size = 0.5),
        axis.title = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        plot.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0),
        legend.key.size = unit(1, 'mm'),
        plot.background = element_blank(),
        panel.grid.major = element_line(size = 0.3),
        panel.grid.minor = element_line(size = 0.2),
        legend.position = "none")

figure2C1
ggsave("data/processed/figures/Figure-2/Eucl/Figure-2C1-umap.pdf", figure2C1, width = 70, height = 25, unit = 'mm')

# Panel 2, Context -----------------------------------------------------------------------

fraction_palette = c(">0.2" = "#A8327D",
                     "0.8-20" = "#C43C75",
                     "0.2-3" = "#DD4968",
                     "0.2-0.8" = "#F2605C",
                     "<0.2" = "#FA7F5E")
p1 = go_secmet_t1 %>%
  mutate(clusters = gsub("Cluster ", "", clusters)) %>%
  ggplot() +
  geom_col(aes(x = fraction_facet, y = n, fill = fraction_facet)) +
  scale_fill_manual(values = fraction_palette) +
  scale_y_continuous(breaks = c(0, 200), labels = c(0, 200)) +
  facet_grid(clusters~ "Size fraction", switch = "y") +
  ylab("samples") +
  coord_flip() +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size = unit(1, 'mm')),
        axis.ticks.length = unit(1, 'mm'),
        axis.text.x = element_text(margin = margin(-1,0,0,0)),
        axis.title.x = element_text(margin = margin(-1,0,0,0)),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, .5, .5, .5, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_rect(color = NA, fill = "black"),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
p1

p2 = go_secmet_t1 %>%
  mutate(depth_layer = factor(depth_layer, levels = rev(levels(depth_layer)))) %>%
  group_by(clusters) %>%
  mutate(n = n / sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = depth_layer, y = n, fill = ocean_province)) +
  scale_fill_manual(values = ocean_material2_colors, name = "Ocean Provinces") +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  facet_grid(clusters~ "Biogeography", switch = "y") +
  ylab("%") +
  coord_flip() +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(1, 'mm'),
        axis.ticks = element_line(size = unit(1, 'mm')),
        axis.text.x = element_text(margin = margin(-1,0,0,0)),
        axis.title.x = element_text(margin = margin(-1,0,0,0)),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, 1, .5, 1, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))

p2

p3 = go_secmet_t2%>%
  mutate(class = factor(class, levels = rev(levels(class)))) %>%
  ggplot() +
  geom_boxplot(aes(x = class, y = abundance/1000, fill = class, color = class),
               alpha = .8, size = unit(.25, 'pt'), fatten = 1, outlier.colour = NA) + # FIXME divide by 1000 because I stupidly scaled down the total motus count by 1000
  scale_colour_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  facet_grid(clusters~ "Product class", switch = "y") +
  scale_y_log10(breaks = c(0.000001, 0.0001, 0.01), labels = c("10-6", "10-4", "10-2")) +
  ylab("FPKb/cell") +
  coord_flip() +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = .5),
        axis.ticks.length = unit(.5, 'mm'),
        axis.text.x = element_text(margin = margin(-1,0,0,0)),
        axis.title.x = element_text(margin = margin(-1,0,0,0)),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, 1, .5, 1, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
p3

p4 = go_secmet_t3 %>%
  ggplot() +
  geom_violin(aes(x = "dummy", y = avg_size), fill = "#595959", color = "#595959",
              alpha = .6, size = unit(.25, 'pt')) +
  facet_grid(clusters~"Genome size", switch = "y") +
  ylab("Mbp") +
  scale_y_continuous(limits = c(10**6, 5*10**6),
                     breaks = c(1*10**6, 3*10**6, 5*10**6),
                     labels = c(1, 3, 5)) +
  coord_flip() +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = .5),
        axis.ticks.length = unit(.5, 'mm'),
        axis.text.x = element_text(margin = margin(-1,0,0,0)),
        axis.title.x = element_text(margin = margin(-1,0,0,0)),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, 1, .5, 1, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
p4


go_secmet_t3 %>% 
  left_join(go_secmet_t2 %>% select(-richness, -diversity) %>% spread(class, abundance)) %>%
  select(-sample, -V1, -V2) %>% 
  rename(average_genome_size = avg_size) %>%
  googlesheets4::write_sheet(., ss = "1VIETOjN8M9bjUMftr_2vJ7SnOLPe1yMEqI093FLNAW4", "Structure of the ocean microbiome biosynthetic potential")

go_secmet_t3 = googlesheets4::read_sheet(ss = "1VIETOjN8M9bjUMftr_2vJ7SnOLPe1yMEqI093FLNAW4", "Structure of the ocean microbiome biosynthetic potential") %>% rename(avg_size = average_genome_size)

# Environmental parameters
p_temp = gcf_summary_eucl_rev_eco_summary %>%
  ggplot() +
  geom_boxplot(aes(x = clusters, y = `temperature [°C]`), outlier.size = 1) +
  facet_grid(clusters~ "Temperature", switch = "y", scales = "free") +
  coord_flip() +
  ylab("Temperature (ºC)") +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size = unit(1, 'mm')),
        axis.ticks.length = unit(1, 'mm'),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, .5, .5, .5, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_rect(color = NA, fill = "black"),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
p_depth = gcf_summary_eucl_rev_eco_summary %>%
  ggplot() +
  geom_boxplot(aes(x = clusters, y = as.numeric(depth)), outlier.size = 1) +
  facet_grid(clusters~ "Depth", switch = "y", scales = "free") +
  coord_flip() +
  scale_y_log10() +
  ylab("Depth (m)") +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = .5),
        axis.ticks.length = unit(.5, 'mm'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, 1, .5, 1, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
p_temp | p_depth
ggsave("data/processed/figures/Figure-2/Eucl/Figure-2CX-environment.pdf", width = 113, height = 50, unit = 'mm')

gcf_summary_eucl_rev_eco_summary %>%
  ggplot() +
  geom_violin(aes(x = clusters, y = abs(`latitude`)), outlier.size = 1) +
  facet_grid(clusters~ "Latitude", switch = "y", scales = "free") +
  coord_flip()
gcf_summary_eucl_rev_eco_summary %>%
  ggplot() +
  geom_violin(aes(x = clusters, y = as.numeric(depth)), outlier.size = 1) +
  facet_grid(clusters~ "Depth", switch = "y", scales = "free") +
  coord_flip() +
  #scale_y_log10() +
  ylab("Depth (m)") +
  ylim(0, 250)

go_secmet_t3 %>%
  filter(depth_layer %in% c("EPI", "MES")) %>%
  ggplot() +
  geom_point(aes(x = `temperature [°C]`, y = avg_size), size = .5, alpha = .6) +
  geom_smooth(aes(x = `temperature [°C]`, y = avg_size), method = "lm", size = .5) +
  ylab("Average genome size per sample") +
  scale_y_continuous(breaks = c(2*10**6, 3*10**6, 4*10**6), labels = c("2Mbp", "3Mbp", "4Mbp")) +
  facet_grid(.~fraction_facet) +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = .5),
        axis.ticks.length = unit(.5, 'mm'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, 1, .5, 1, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(1,1,1,1, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
ggsave("data/processed/figures/Figure-2/Eucl/Figure-2CX-genome_size_colder.pdf", width = 113, height = 50, unit = 'mm')


lm_0.2_minus = lm(go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == "<0.2") %>% pull(avg_size) ~
                  go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == "<0.2") %>% pull(`temperature [°C]`))
summary(lm_0.2_minus)
lm_0.2_0.8 = lm(go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == "0.2-0.8") %>% pull(avg_size) ~
                  go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == "0.2-0.8") %>% pull(`temperature [°C]`))
summary(lm_0.2_0.8)
lm_0.2_3 = lm(go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == "0.2-3") %>% pull(avg_size) ~
                  go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == "0.2-3") %>% pull(`temperature [°C]`))
summary(lm_0.2_3)
lm_0.2_plus = lm(go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == ">0.2") %>% pull(avg_size) ~
                  go_secmet_t3 %>% filter(depth_layer %in% c("EPI", "MES") & fraction_facet == ">0.2") %>% pull(`temperature [°C]`))
summary(lm_0.2_plus)

go_secmet_t3 %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = avg_size), size = .3, outlier.size = .5) +
  ylab("Average genome size per sample") +
  scale_y_continuous(breaks = c(2*10**6, 3*10**6, 4*10**6), labels = c("2Mbp", "3Mbp", "4Mbp")) +
  facet_grid(.~fraction_facet, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = .5),
        axis.ticks.length = unit(.5, 'mm'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, 1, .5, 1, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(1,1,1,1, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
ggsave("data/processed/figures/Figure-2/Eucl/Figure-2CX-genome_size_deeper.pdf", width = 113, height = 40, unit = 'mm')

# Statistical significance -- Product classes
go_secmet_t2_spread = go_secmet_t2 %>%
  select(-richness, -diversity) %>%
  unique() %>%
  spread(class, abundance)
pairwise.wilcox.test(go_secmet_t3$`Non-Ribosomal Peptide Synthetases`, go_secmet_t3$clusters, p.adjust.method = 'fdr')
pairwise.wilcox.test(go_secmet_t3$`Type I Polyketide Synthases`, go_secmet_t3$clusters, p.adjust.method = 'fdr')
pairwise.wilcox.test(go_secmet_t3$`Type II/III Polyketide Synthases`, go_secmet_t3$clusters, p.adjust.method = 'fdr')
pairwise.wilcox.test(go_secmet_t3$`RiPPs (Ribosomal Natural Products)`, go_secmet_t3$clusters, p.adjust.method = 'fdr')
pairwise.wilcox.test(go_secmet_t3$Terpenes, go_secmet_t3$clusters, p.adjust.method = 'fdr')
pairwise.wilcox.test(go_secmet_t3$Other, go_secmet_t3$clusters, p.adjust.method = 'fdr')
# Statistical significance -- genome sizes
go_secmet_t3 %>% select(avg_size, clusters) %>% split(.$clusters) %>% map(summary)
pairwise.wilcox.test(go_secmet_t3$avg_size, go_secmet_t3$clusters)
go_secmet_t3 %>% select(avg_size, depth_layer) %>% split(.$depth_layer) %>% map(summary)
pairwise.wilcox.test(go_secmet_t3$avg_size, go_secmet_t3$depth_layer)
kruskal.test(go_secmet_t3$avg_size, go_secmet_t3$depth_layer)
# Statistical significance -- environment
pairwise.wilcox.test(gcf_summary_eucl_rev_eco_summary$`temperature [°C]`, gcf_summary_eucl_rev_eco_summary$clusters)
pairwise.wilcox.test(as.numeric(gcf_summary_eucl_rev_eco_summary$depth), gcf_summary_eucl_rev_eco_summary$clusters)
kruskal.test(as.numeric(go_secmet_t3$depth), go_secmet_t3$clusters)
kruskal.test(as.numeric(go_secmet_t3$`temperature [°C]`), go_secmet_t3$clusters)


figure2C2 = p1 + p2 + p3 + p4 + plot_layout(nrow = 1) + plot_annotation(theme = theme(plot.margin = margin(0,0,0,0)))
figure2C2
ggsave("data/processed/figures/Figure-2/Eucl/Figure-2C2-interpretation.pdf", figure2C2, width = 70, height = 30, unit = 'mm')


# Dominant taxa -----------------------------------------------

motus_abundances_rel = motus_abundances %>%
  group_by(sample) %>%
  mutate(rel_abd = motus_count/sum(motus_count)) %>%
  ungroup() %>%
  filter(rel_abd >= 0.01)

motus_abundances_rel_c1 = motus_abundances_rel %>%
  filter(sample %in% (gcf_summary_eucl_rev_eco_summary %>% filter(clusters == "Cluster 1") %>% pull(sample)))
motus_abundances_rel_c1 %>%
  filter(motus_cluster != -1) %>% 
  group_by(motus_cluster) %>%
  summarize(n = n()/(gcf_summary_eucl_rev_eco_summary %>% filter(clusters == "Cluster 1") %>% nrow),
            avg_rel_abd = mean(rel_abd),
            taxonomy = unique(taxonomy)) %>%
  filter(n > 0.1) %>%
  View("C1")

motus_abundances_rel_c2 = motus_abundances_rel %>%
  filter(sample %in% (gcf_summary_eucl_rev_eco_summary %>% filter(clusters == "Cluster 2") %>% pull(sample)))
motus_abundances_rel_c2 %>%
  filter(motus_cluster != -1) %>% 
  group_by(motus_cluster) %>%
  summarize(n = n()/(gcf_summary_eucl_rev_eco_summary %>% filter(clusters == "Cluster 2") %>% nrow),
            avg_rel_abd = mean(rel_abd),
            taxonomy = unique(taxonomy)) %>%
  filter(n > 0.1) %>%
  View("C2")


# Prepare abundance and prevalence of GCCs -----------------------------------------------

gcc_abd = gcf_summary_processed %>%
  left_join(bgc_clustering %>% select(gcf, gcc) %>% mutate(gcf = paste0("gcf_", gcf)) %>% unique()) %>%
  group_by(gcc, sample) %>%
  summarize(value = sum(value)) %>%
  group_by(gcc) %>%
  summarize(abundance = sum(value),
            prevalence = sum(value > 0)/n())

write_tsv(gcc_abd, "data/raw/go_microbiomics/bgcs/tables-v0.1/gcc_abundance_prevalence.tsv")

# Now, test for novelty within the clusters ----------------------------------------------

# We want to do that fairly, i.e. with the same sampling effort per sample
# Let's approximate the numer of cells with the motus count

motus_abundances_per_sample = motus_abundances %>%
  group_by(sample) %>%
  summarize(motus_sum = sum(motus_count)) %>%
  arrange(motus_sum)

motus_abundances_per_sample %>%
  ggplot() +
  geom_density(aes(x = motus_sum)) +
  geom_vline(xintercept = 2000)

samples_to_exclude = motus_abundances_per_sample %>%
  filter(motus_sum < 2000) %>%
  pull(sample)

# assumes gcf_dist is available
gcf_rarefied = gcf_summary_processed %>%
  filter(!sample %in% samples_to_exclude) %>%
  rowwise() %>%
  mutate(raref = rbinom(1, size = 2000, p = min(value, 1))) %>% # Note that this is stochastic... and a bit slow
  ungroup()

gcf_rarefied_to_plot = gcf_rarefied %>%
  left_join(gcf_summary_eucl_rev_eco_summary) %>%
  left_join(gcf_dist %>% mutate(gcf = paste0('gcf_', gcf))) %>%
  filter(raref > 0) %>%
  right_join(gcf_classes_distrib %>% filter(gcf %in% unique(gcf_rarefied %>% filter(raref > 0) %>% pull(gcf)))) %>%
  group_by(sample, class) %>%
  summarize(richness_all = sum(freq > 0),
            richness_new = sum(freq[median_d >= 0.2] > 0),
            clusters = unique(clusters),
            fraction_facet = unique(fraction_facet),
            depth_layer = unique(depth_layer),
            ocean_province = unique(ocean_province),
            dataset = unique(dataset),
            polar = ifelse(abs(unique(latitude)) >= 60, "polar", "non-polar"))

gcf_rarefied_to_plot %>%
  filter(class %in% c("Terpenes", "Other")) %>%
  ggplot() +
  geom_boxplot(aes(x = polar, y = richness_new, fill = class, color = class), alpha = .6) +
  facet_grid(depth_layer ~ fraction_facet, scales = "free_x", space = "free_x") +
  ylab("Number of new GCFs / 2,000 cells") +
  scale_color_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  #scale_y_log10() +
  theme_bw() +
  theme(axis.title.x = element_blank())

gcf_rarefied_to_plot %>%
  filter(!class %in% c("Terpenes", "Other")) %>%
  ggplot() +
  geom_boxplot(aes(x = polar, y = richness_new, fill = class, color = class), alpha = .6) +
  facet_grid(depth_layer ~ fraction_facet, scales = "free_x", space = "free_x") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  #scale_y_log10() +
  theme_bw() +
  theme(axis.title.x = element_blank())

# Testing for significance =======

#Vir, Epi

gcf_rarefied_latitude_vir = gcf_rarefied_to_plot %>%
  filter(fraction_facet %in% c("<0.2")) %>%
  filter(depth_layer == "EPI")
gcf_rarefied_latitude_vir %>%
  ggplot() +
  geom_boxplot(aes(x = polar, y = richness_new, fill = class, color = class), alpha = .6) +
  facet_grid(depth_layer ~ fraction_facet*class, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  #scale_y_log10() +
  theme_bw() +
  theme(axis.title.x = element_blank())
gcf_rarefied_latitude_vir_terpenes = gcf_rarefied_latitude_vir %>% filter(class == "Terpenes")
wilcox.test(gcf_rarefied_latitude_vir_terpenes$richness_new ~ gcf_rarefied_latitude_vir_terpenes$polar)
gcf_rarefied_latitude_vir_other = gcf_rarefied_latitude_vir %>% filter(class == "Other")
wilcox.test(gcf_rarefied_latitude_vir_other$richness_new ~ gcf_rarefied_latitude_vir_other$polar)
gcf_rarefied_latitude_vir_nrps = gcf_rarefied_latitude_vir %>% filter(class == "Non-Ribosomal Peptide Synthetases")
wilcox.test(gcf_rarefied_latitude_vir_nrps$richness_new ~ gcf_rarefied_latitude_vir_nrps$polar)
gcf_rarefied_latitude_vir_ripps = gcf_rarefied_latitude_vir %>% filter(class == "RiPPs (Ribosomal Natural Products)")
wilcox.test(gcf_rarefied_latitude_vir_ripps$richness_new ~ gcf_rarefied_latitude_vir_ripps$polar)
gcf_rarefied_latitude_vir_t1pks = gcf_rarefied_latitude_vir %>% filter(class == "Type I Polyketide Synthases")
wilcox.test(gcf_rarefied_latitude_vir_t1pks$richness_new ~ gcf_rarefied_latitude_vir_t1pks$polar)
gcf_rarefied_latitude_vir_t23pks = gcf_rarefied_latitude_vir %>% filter(class == "Type II/III Polyketide Synthases")
wilcox.test(gcf_rarefied_latitude_vir_t23pks$richness_new ~ gcf_rarefied_latitude_vir_t23pks$polar)


# All data, testing depth -----------

gcf_rarefied_depth = gcf_rarefied_to_plot %>%
  filter(depth_layer %in% c("EPI", "MES", "BAT")) 
gcf_rarefied_depth %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = class, color = class), alpha = .6) +
  facet_grid(. ~ class, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  #scale_y_log10() +
  theme_bw() +
  theme(axis.title.x = element_blank())

gcf_rarefied_depth_nrps = gcf_rarefied_depth %>% filter(class == "Non-Ribosomal Peptide Synthetases")
pairwise.wilcox.test(gcf_rarefied_depth_nrps$richness_new, gcf_rarefied_depth_nrps$depth_layer)

gcf_rarefied_depth_other = gcf_rarefied_depth %>% filter(class == "Other")
pairwise.wilcox.test(gcf_rarefied_depth_other$richness_new, gcf_rarefied_depth_other$depth_layer)

gcf_rarefied_depth_ripps = gcf_rarefied_depth %>% filter(class == "RiPPs (Ribosomal Natural Products)")
pairwise.wilcox.test(gcf_rarefied_depth_ripps$richness_new, gcf_rarefied_depth_ripps$depth_layer)

gcf_rarefied_depth_terpenes = gcf_rarefied_depth %>% filter(class == "Terpenes")
pairwise.wilcox.test(gcf_rarefied_depth_terpenes$richness_new, gcf_rarefied_depth_terpenes$depth_layer)

gcf_rarefied_depth_t1pks = gcf_rarefied_depth %>% filter(class == "Type I Polyketide Synthases")
pairwise.wilcox.test(gcf_rarefied_depth_t1pks$richness_new, gcf_rarefied_depth_t1pks$depth_layer)

gcf_rarefied_depth_t23pks = gcf_rarefied_depth %>% filter(class == "Type II/III Polyketide Synthases")
pairwise.wilcox.test(gcf_rarefied_depth_t23pks$richness_new, gcf_rarefied_depth_t23pks$depth_layer)

# Small prok, testing depth -----------

gcf_rarefied_depth_prok1 = gcf_rarefied_to_plot %>%
  filter(fraction_facet %in% c("0.2-0.8")) %>%
  filter(depth_layer %in% c("EPI", "MES", "BAT")) 
gcf_rarefied_depth_prok1 %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = class, color = class), alpha = .6) +
  facet_grid(. ~ fraction_facet*class, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  #scale_y_log10() +
  theme_bw() +
  theme(axis.title.x = element_blank())

gcf_rarefied_depth_prok1_nrps = gcf_rarefied_depth_prok1 %>% filter(class == "Non-Ribosomal Peptide Synthetases")
pairwise.wilcox.test(gcf_rarefied_depth_prok1_nrps$richness_new, gcf_rarefied_depth_prok1_nrps$depth_layer)

gcf_rarefied_depth_prok1_other = gcf_rarefied_depth_prok1 %>% filter(class == "Other")
pairwise.wilcox.test(gcf_rarefied_depth_prok1_other$richness_new, gcf_rarefied_depth_prok1_other$depth_layer)

gcf_rarefied_depth_prok1_ripps = gcf_rarefied_depth_prok1 %>% filter(class == "RiPPs (Ribosomal Natural Products)")
pairwise.wilcox.test(gcf_rarefied_depth_prok1_ripps$richness_new, gcf_rarefied_depth_prok1_ripps$depth_layer)

gcf_rarefied_depth_prok1_terpenes = gcf_rarefied_depth_prok1 %>% filter(class == "Terpenes")
pairwise.wilcox.test(gcf_rarefied_depth_prok1_terpenes$richness_new, gcf_rarefied_depth_prok1_terpenes$depth_layer)

gcf_rarefied_depth_prok1_t1pks = gcf_rarefied_depth_prok1 %>% filter(class == "Type I Polyketide Synthases")
pairwise.wilcox.test(gcf_rarefied_depth_prok1_t1pks$richness_new, gcf_rarefied_depth_prok1_t1pks$depth_layer)

gcf_rarefied_depth_prok1_t23pks = gcf_rarefied_depth_prok1 %>% filter(class == "Type II/III Polyketide Synthases")
pairwise.wilcox.test(gcf_rarefied_depth_prok1_t23pks$richness_new, gcf_rarefied_depth_prok1_t23pks$depth_layer)

# Large Prok, testing depth -------------

gcf_rarefied_depth_prok2 = gcf_rarefied_to_plot %>%
  filter(fraction_facet %in% c("0.2-3")) %>%
  filter(depth_layer %in% c("EPI", "MES")) 
gcf_rarefied_depth_prok2 %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = class, color = class), alpha = .6) +
  facet_grid(. ~ fraction_facet*class, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  #scale_y_log10() +
  theme_bw() +
  theme(axis.title.x = element_blank())

gcf_rarefied_depth_prok2_nrps = gcf_rarefied_depth_prok2 %>% filter(class == "Non-Ribosomal Peptide Synthetases")
wilcox.test(gcf_rarefied_depth_prok2_nrps$richness_new ~ gcf_rarefied_depth_prok2_nrps$depth_layer)

gcf_rarefied_depth_prok2_other = gcf_rarefied_depth_prok2 %>% filter(class == "Other")
wilcox.test(gcf_rarefied_depth_prok2_other$richness_new ~ gcf_rarefied_depth_prok2_other$depth_layer)

gcf_rarefied_depth_prok2_ripps = gcf_rarefied_depth_prok2 %>% filter(class == "RiPPs (Ribosomal Natural Products)")
wilcox.test(gcf_rarefied_depth_prok2_ripps$richness_new ~ gcf_rarefied_depth_prok2_ripps$depth_layer)

gcf_rarefied_depth_prok2_terpenes = gcf_rarefied_depth_prok2 %>% filter(class == "Terpenes")
wilcox.test(gcf_rarefied_depth_prok2_terpenes$richness_new ~ gcf_rarefied_depth_prok2_terpenes$depth_layer)

gcf_rarefied_depth_prok2_t1pks = gcf_rarefied_depth_prok2 %>% filter(class == "Type I Polyketide Synthases")
wilcox.test(gcf_rarefied_depth_prok2_t1pks$richness_new ~ gcf_rarefied_depth_prok2_t1pks$depth_layer)

gcf_rarefied_depth_prok2_t23pks = gcf_rarefied_depth_prok2 %>% filter(class == "Type II/III Polyketide Synthases")
wilcox.test(gcf_rarefied_depth_prok2_t23pks$richness_new ~ gcf_rarefied_depth_prok2_t23pks$depth_layer)

# NRPS/RiPPs enriched in P-A -------------

gcf_rarefied_pa_bat = gcf_rarefied_to_plot %>%
  filter(fraction_facet %in% c("0.2-0.8", "0.8-20")) %>%
  filter(depth_layer %in% c("BAT")) 
plot_test_pa = gcf_rarefied_pa_bat %>%
  #filter(class %in% c("Non-Ribosomal Peptide Synthetases", "RiPPs (Ribosomal Natural Products)")) %>%
  ggplot() +
  geom_boxplot(aes(x = fraction_facet, y = richness_new, fill = class, color = class), alpha = .6) +
  facet_grid(. ~ depth_layer*class, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = bgc_colors_serina) +
  scale_fill_manual(values = bgc_colors_serina) +
  #scale_y_log10() +
  theme_bw() +
  theme(axis.title.x = element_blank())
plot_test_pa

gcf_rarefied_pa_bat_nrps = gcf_rarefied_pa_bat %>% filter(class == "Non-Ribosomal Peptide Synthetases")
wilcox.test(gcf_rarefied_pa_bat_nrps$richness_new ~ gcf_rarefied_pa_bat_nrps$fraction_facet)

gcf_rarefied_pa_bat_other = gcf_rarefied_pa_bat %>% filter(class == "Other")
wilcox.test(gcf_rarefied_pa_bat_other$richness_new ~ gcf_rarefied_pa_bat_other$fraction_facet)

gcf_rarefied_pa_bat_ripps = gcf_rarefied_pa_bat %>% filter(class == "RiPPs (Ribosomal Natural Products)")
wilcox.test(gcf_rarefied_pa_bat_ripps$richness_new ~ gcf_rarefied_pa_bat_ripps$fraction_facet)

gcf_rarefied_pa_bat_terpenes = gcf_rarefied_pa_bat %>% filter(class == "Terpenes")
wilcox.test(gcf_rarefied_pa_bat_terpenes$richness_new ~ gcf_rarefied_pa_bat_terpenes$fraction_facet)

gcf_rarefied_pa_bat_t1pks = gcf_rarefied_pa_bat %>% filter(class == "Type I Polyketide Synthases")
wilcox.test(gcf_rarefied_pa_bat_t1pks$richness_new ~ gcf_rarefied_pa_bat_t1pks$fraction_facet)

gcf_rarefied_pa_bat_t23pks = gcf_rarefied_pa_bat %>% filter(class == "Type II/III Polyketide Synthases")
wilcox.test(gcf_rarefied_pa_bat_t23pks$richness_new ~ gcf_rarefied_pa_bat_t23pks$fraction_facet)

# Maybe by product class is more relevant? -------------

polar_colors = c("non-polar" = "#D28B38", "polar" = "#4E68B0")

# Combined plot
gcf_rarefied_to_plot_combined = gcf_rarefied_to_plot %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2")) %>%
  mutate(class = gsub("Type.*Polyketide", "Polyketide", class))
gcf_rarefied_to_plot_combined %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6, position = position_dodge(preserve = "single"), size = 0.3, outlier.size = .5) +
  facet_grid(class ~ fraction_facet, scales = "free", space = "free_x") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'pt')), 
        text = element_text(size = 6),
        #rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = .5),
        axis.ticks.length = unit(.5, 'mm'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(1, 1, 1, 1, "mm")),
        strip.background = element_rect(color = NA, fill = "black"),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0))
ggsave("data/processed/figures/Figure-2/Eucl/Figure-2CX-discovery_NP.pdf", width = 113, height = 113, unit = 'mm')


# NRPS
# ---
gcf_rarefied_to_plot_nrps = gcf_rarefied_to_plot %>%
  filter(class == "Non-Ribosomal Peptide Synthetases") %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2"))
gcf_rarefied_to_plot_nrps %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6) +
  facet_grid(. ~ class*fraction_facet, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())
  
paste0(gcf_rarefied_to_plot_nrps$fraction_facet, ";", gcf_rarefied_to_plot_nrps$depth_layer, ";", gcf_rarefied_to_plot_nrps$polar) %>% table()
pairwise.wilcox.test(gcf_rarefied_to_plot_nrps$richness_new,
                     paste0(gcf_rarefied_to_plot_nrps$fraction_facet, ";", gcf_rarefied_to_plot_nrps$depth_layer, ";", gcf_rarefied_to_plot_nrps$polar),
                     p.adjust.method = "fdr")

# RiPPS
# ---
gcf_rarefied_to_plot_ripps = gcf_rarefied_to_plot %>%
  filter(class == "RiPPs (Ribosomal Natural Products)") %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2"))
gcf_rarefied_to_plot_ripps %>% ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6) +
  facet_grid(. ~ class*fraction_facet, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())

pairwise.wilcox.test(gcf_rarefied_to_plot_ripps$richness_new,
                     paste0(gcf_rarefied_to_plot_ripps$fraction_facet, ";", gcf_rarefied_to_plot_ripps$depth_layer, ";", gcf_rarefied_to_plot_ripps$polar),
                     p.adjust.method = "fdr")

# T1PKS
# ---
gcf_rarefied_to_plot_t1pks = gcf_rarefied_to_plot %>%
  filter(class == "Type I Polyketide Synthases") %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2"))
gcf_rarefied_to_plot_t1pks %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6) +
  facet_grid(. ~ class*fraction_facet, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())

pairwise.wilcox.test(gcf_rarefied_to_plot_t1pks$richness_new,
                     paste0(gcf_rarefied_to_plot_t1pks$fraction_facet, ";", gcf_rarefied_to_plot_t1pks$depth_layer, ";", gcf_rarefied_to_plot_t1pks$polar),
                     p.adjust.method = "fdr")

# T2/3PKS
# ---
gcf_rarefied_to_plot_t23pks = gcf_rarefied_to_plot %>%
  filter(class == "Type II/III Polyketide Synthases") %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2"))
gcf_rarefied_to_plot_t23pks %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6) +
  facet_grid(. ~ class*fraction_facet, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())

pairwise.wilcox.test(gcf_rarefied_to_plot_t23pks$richness_new,
                     paste0(gcf_rarefied_to_plot_t23pks$fraction_facet, ";", gcf_rarefied_to_plot_t23pks$depth_layer, ";", gcf_rarefied_to_plot_t23pks$polar),
                     p.adjust.method = "fdr")

# PKS
# ---
gcf_rarefied_to_plot_pks = gcf_rarefied_to_plot %>%
  filter(class %in% c("Type II/III Polyketide Synthases", "Type I Polyketide Synthases")) %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2"))
gcf_rarefied_to_plot_pks %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6) +
  facet_grid(. ~ "PKS"*fraction_facet, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())

pairwise.wilcox.test(gcf_rarefied_to_plot_pks$richness_new,
                     paste0(gcf_rarefied_to_plot_pks$fraction_facet, ";", gcf_rarefied_to_plot_pks$depth_layer, ";", gcf_rarefied_to_plot_pks$polar),
                     p.adjust.method = "fdr")

# Terpenes
# ---
gcf_rarefied_to_plot_terpenes = gcf_rarefied_to_plot %>%
  filter(class == "Terpenes") %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2"))

gcf_rarefied_to_plot_terpenes %>% ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6) +
  facet_grid(. ~ class*fraction_facet, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())

pairwise.wilcox.test(gcf_rarefied_to_plot_terpenes$richness_new,
                     paste0(gcf_rarefied_to_plot_terpenes$fraction_facet, ";", gcf_rarefied_to_plot_terpenes$depth_layer, ";", gcf_rarefied_to_plot_terpenes$polar),
                     p.adjust.method = "fdr")
# Other
# ---
gcf_rarefied_to_plot_other = gcf_rarefied_to_plot %>%
  filter(class == "Other") %>%
  filter(!(depth_layer %in% c("BAT", "ABY") & fraction_facet == ">0.2"))
gcf_rarefied_to_plot_other %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = richness_new, fill = polar, color = polar), alpha = .6) +
  facet_grid(. ~ class*fraction_facet, scales = "free", space = "free") +
  ylab("Number of new GCFs for 2,000 cells") +
  scale_color_manual(values = polar_colors) +
  scale_fill_manual(values = polar_colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())

pairwise.wilcox.test(gcf_rarefied_to_plot_other$richness_new,
                     paste0(gcf_rarefied_to_plot_other$fraction_facet, ";", gcf_rarefied_to_plot_other$depth_layer, ";", gcf_rarefied_to_plot_other$polar),
                     p.adjust.method = "fdr")
