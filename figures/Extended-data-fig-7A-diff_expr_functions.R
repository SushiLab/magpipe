# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 7A - diff expr functions ============================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(patchwork)
library(viridis)
library(treeio) # YuLab-SMU/treeio
library(ggtree) # YuLab-SMU/treeio
library(tidytree)
#source('code/figures/figure4D-dataprep.R') # expects it to have been run

source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Data preparation =======================================================================

# functional annotations & metadata ------------------------------------------------------

formatted_metadata = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")
MGs = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-fetchMGs")
TXSSCAN = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-txsscan")
prokka = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-prokka")
kegg = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-kegg")
eggnog = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-eggnog")
BGCs = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-antismash")
pred_index = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-predatory")
metat_featurecounts_summary_processed = read_tsv("data/raw/marine_eremios/marine_eremios-featurecounts-summary_PIAMPJPB_metat-processed.tsv.gz", guess_max = 200000)
samples_clustering = read_tsv("data/raw/marine_eremios/PIAMPJPB-transcriptome_clusters-sample_order_membership.tsv")
cluster_tree = read.newick("data/raw/marine_eremios/PIAMPJPB-transcriptome_clusters.newick")

# Metat abundances featurecounts ---------------------------------------------------------

metat_featurecounts_summary_processed %>%
  filter(Gene %in% (MGs %>% filter(`mOTUs-10`) %>% pull(bestMGs))) %>%
  ggplot() +
  geom_density(aes(x = value))

metat_featurecounts_summary_spread = metat_featurecounts_summary_processed %>%
  select(Gene, material, value) %>% 
  filter(!is.na(material)) %>%
  spread(Gene, value)

# Functional analysis ====================================================================

metat_featurecounts_summary_clustered = metat_featurecounts_summary_processed %>%
  left_join(samples_clustering) %>%
  mutate(material = factor(material, levels = unique(samples_clustering$material))) %>%
  arrange(metat_featurecounts_summary_clustered)

statistics_table = NULL

# BGCs - Biosynthetic --------------------------------------------------------------------

antismash_biosynthetic = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (BGCs %>% filter(biosynthetic_gene | biosynthetic_add_gene) %>% pull(gene))) %>%
  left_join(BGCs, by = c("Gene" = "gene")) %>%
  group_by(biosynthetic_region, material) %>%
  summarize(n_biosynthetic_genes = n(),
            value = median(value),
            cluster = unique(cluster)) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster),
            n = n())

antismash_biosynthetic %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

antismash_biosynthetic_kruskal = kruskal.test(antismash_biosynthetic %>% pull(value), antismash_biosynthetic %>% pull(cluster))
antismash_biosynthetic_kruskal
pairwise.wilcox.test(antismash_biosynthetic %>% pull(value), antismash_biosynthetic %>% pull(cluster))

statistics_table = rbind(statistics_table,
                         tibble(category = "BGCs",
                                n = unique(antismash_biosynthetic$n),
                                statistic = antismash_biosynthetic_kruskal$statistic,
                                pvalue = antismash_biosynthetic_kruskal$p.value))

# BGCs - Transport -----------------------------------------------------------------------

antismash_transport = metat_featurecounts_summary_clustered %>%
  left_join(BGCs, by = c("Gene" = "gene")) %>%
  filter(transport_gene) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster))

antismash_transport %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

kruskal.test(bgc_t %>% pull(value), bgc_t %>% pull(cluster))
pairwise.wilcox.test(bgc_t %>% pull(value), bgc_t %>% pull(cluster))

# TXSSCAN - Secretion --------------------------------------------------------------------

txss_secretion = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (TXSSCAN %>% filter(Profile_coverage > .7 & as.numeric(`i-evalue`) < 10**-5) %>% filter(grepl("T.SS", Predicted_system) & grepl("T.SS", Gene)) %>% pull(Hit_Id))) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster),
            n = n())

txss_secretion %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

txss_secretion_kruskal = kruskal.test(txss_secretion %>% pull(value), txss_secretion %>% pull(cluster))
txss_secretion_kruskal
pairwise.wilcox.test(txss_secretion %>% pull(value), txss_secretion %>% pull(cluster))

statistics_table = rbind(statistics_table,
                         tibble(category = "Secretion",
                                n = unique(txss_secretion$n),
                                statistic = txss_secretion_kruskal$statistic,
                                pvalue = txss_secretion_kruskal$p.value))

# TXSSCAN - Flagellum --------------------------------------------------------------------

txss_flagellum = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (TXSSCAN %>% filter(Profile_coverage > .7 & as.numeric(`i-evalue`) < 10**-5) %>% filter(Predicted_system == "Flagellum" & grepl("Flg_", Gene)) %>% pull(Hit_Id))) %>%
  filter(grepl("scaffold_18|scaffold_78", Gene)) %>% # Where you have co-localization
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster),
            n = n())

txss_flagellum %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

txss_flagellum_kruskal = kruskal.test(txss_flagellum %>% pull(value), txss_flagellum %>% pull(cluster))
txss_flagellum_kruskal
pairwise.wilcox.test(txss_flagellum %>% pull(value), txss_flagellum %>% pull(cluster))

statistics_table = rbind(statistics_table,
                         tibble(category = "Flagellum",
                                n = unique(txss_flagellum$n),
                                statistic = txss_flagellum_kruskal$statistic,
                                pvalue = txss_flagellum_kruskal$p.value))

# Prokka - Motility ----------------------------------------------------------------------

prok_motility = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (prokka %>% filter(grepl('[M|m]otil|Flagel', product)) %>% pull(locus_tag))) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster))

prok_motility %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

kruskal.test(prok_motility %>% pull(value), prok_motility %>% pull(cluster))
pairwise.wilcox.test(prok_motility %>% pull(value), prok_motility %>% pull(cluster))

# Integrated - Motility ----------------------------------------------------------------------

prok_motility_genes = prokka %>% filter(grepl('[M|m]otil|Flagel', product)) %>% pull(locus_tag)
txsscan_flagellum_genes = TXSSCAN %>% filter(Profile_coverage > .7 & as.numeric(`i-evalue`) < 10**-5) %>% filter(Predicted_system == "Flagellum" & grepl("Flg_", Gene)) %>% pull(Hit_Id)
prok_motility_genes[!prok_motility_genes %in% txsscan_flagellum_genes]

integrated_motility = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% c(prok_motility_genes, txsscan_flagellum_genes)) %>%
  filter(grepl("scaffold_18|scaffold_78", Gene)) %>% # Where you have co-localization
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster))

integrated_motility %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

kruskal.test(integrated_motility %>% pull(value), integrated_motility %>% pull(cluster))
pairwise.wilcox.test(integrated_motility %>% pull(value), integrated_motility %>% pull(cluster))

# Prokka - Degradation -------------------------------------------------------------------

prokka_degrad = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (prokka %>% filter(grepl('[P|p]rotease|[P|p]eptidase', product)) %>% pull(locus_tag))) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster),
            n = n())

prokka_degrad %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

prokka_degrad_kruskal = kruskal.test(prokka_degrad %>% pull(value), prokka_degrad %>% pull(cluster))
prokka_degrad_kruskal
pairwise.wilcox.test(prokka_degrad %>% pull(value), prokka_degrad %>% pull(cluster))

statistics_table = rbind(statistics_table,
                         tibble(category = "Degradation",
                                n = unique(prokka_degrad$n),
                                statistic = prokka_degrad_kruskal$statistic,
                                pvalue = prokka_degrad_kruskal$p.value))

# KEGG - Degradation ---------------------------------------------------------------------

# From KO to Pathways
ko_to_pathway = read_tsv("data/processed/marine_eremios/transcriptomes/ko_pathway.list", col_names = c("KO", "Pathway"))
pathway_description = read_tsv("data/processed/marine_eremios/transcriptomes/pathway.list", col_names = c("Pathway", "Pathway description"), comment = "#")
ko_to_pathway_description = ko_to_pathway %>%
  filter(grepl('map', Pathway)) %>%
  mutate(KO = gsub('ko:', '', KO),
         Pathway = gsub('path:map', '', Pathway)) %>%
  left_join(pathway_description)

# Pathways to BRITE hierarchy
brite <- fread("grep -e '^A' -e '^B ' -e '^C ' data/processed/marine_eremios/transcriptomes/ko00001.keg", sep=" ", data.table=F, fill=T) %>%
  as_tibble() %>%
  mutate(V3 = ifelse(grepl("^A", V1), V2, V3),
         V2 = ifelse(grepl("^A", V1), gsub("^A", "", V1), V2),
         V1 = ifelse(grepl("^A", V1), gsub("[0-9]+$", "", V1), V1)) %>%
  mutate(V2 = str_pad(V2,5,pad="0")) %>%
  unite("V3", V3:V11,sep=" ") %>%
  mutate(V3 = gsub(" +$", "", V3)) %>%
  dplyr::rename(Level = V1, Pathway = V2, Description = V3)

brite_hierarchy <- NULL
for (i in 1:nrow(brite)){
  if (brite$Level[i] == "A") A_level <- brite[i,2:3] %>% dplyr::rename(Pathway_A = Pathway, Description_A = Description)
  if (brite$Level[i] == "B") B_level <- brite[i,2:3] %>% dplyr::rename(Pathway_B = Pathway, Description_B = Description)
  if (brite$Level[i] == "C") {
    C_level = brite[i,2:3] %>% dplyr::rename(Pathway_C = Pathway, Description_C = Description)
    brite_hierarchy <- rbind(brite_hierarchy, cbind(A_level, cbind(B_level, C_level))) %>% as_tibble()
  }
}

kegg_kos_levC <- kegg_kos %>%
  left_join(brite_hierarchy, by = c("Pathway" = "Pathway_C")) %>%
  select(gene, Description_C) %>%
  unique() %>%
  group_by(gene) %>%
  mutate(n_pathway = n()) %>%
  filter(!(n_pathway > 1 & is.na(Description_C)))

kegg_kos_levC %>% filter(grepl('[D|d]egradation', Description_C)) %>% pull(Description_C) %>% table() %>% sort()

kegg_degrad = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (kegg_kos_levC %>% filter(grepl('[D|d]egradation', Description_C)) %>% pull(gene))) %>% 
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster))

kegg_degrad %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

kruskal.test(kegg_degrad %>% pull(value), kegg_degrad %>% pull(cluster))
pairwise.wilcox.test(kegg_degrad %>% pull(value), kegg_degrad %>% pull(cluster))

# Pred. Index - Predatory ----------------------------------------------------------------

predatory = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (pred_index %>% filter(Type == "Predatory") %>% pull(Gene))) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster),
            n = n())

predatory %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

predatory_kruskal = kruskal.test(predatory %>% pull(value), predatory %>% pull(cluster))
predatory_kruskal
pairwise.wilcox.test(predatory %>% pull(value), predatory %>% pull(cluster))


statistics_table = rbind(statistics_table,
                         tibble(category = "Predatory markers",
                                n = unique(predatory$n),
                                statistic = predatory_kruskal$statistic,
                                pvalue = predatory_kruskal$p.value))

# Pred. Index - Non-predatory ------------------------------------------------------------

nonpredatory = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (pred_index %>% filter(Type == "Non-predatory") %>% pull(Gene))) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster),
            n = n())

nonpredatory %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

nonpredatory_kruskal = kruskal.test(nonpredatory %>% pull(value), nonpredatory %>% pull(cluster))
nonpredatory_kruskal
pairwise.wilcox.test(nonpredatory %>% pull(value), nonpredatory %>% pull(cluster))

statistics_table = rbind(statistics_table,
                         tibble(category = "Non-predatory markers",
                                n = unique(nonpredatory$n),
                                statistic = nonpredatory_kruskal$statistic,
                                pvalue = nonpredatory_kruskal$p.value))

# Prokka - Ribosomal ---------------------------------------------------------------------

prokka_ribosomal = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (prokka %>% filter(grepl('ribosomal protein', product)) %>% pull(locus_tag))) %>%
  group_by(material) %>%
  summarize(value = median(value),
            cluster = unique(cluster))

prokka_ribosomal %>%
  ggplot() +
  geom_boxplot(aes(x = cluster, y = value)) +
  geom_point(aes(x = cluster, y = value))

kruskal.test(prokka_ribosomal %>% pull(value), prokka_ribosomal %>% pull(cluster))
pairwise.wilcox.test(prokka_ribosomal %>% pull(value), prokka_ribosomal %>% pull(cluster))

# Summary figure =========================================================================

summary_table = antismash_biosynthetic %>% mutate(category = "BGCs") %>%
  rbind(txss_secretion %>% mutate(category = "Secretion")) %>%
  rbind(txss_flagellum %>% mutate(category = "Flagellum")) %>%
  rbind(prokka_degrad %>% mutate(category = "Degradation")) %>%
  rbind(predatory %>% mutate(category = "Pred.")) %>%
  rbind(nonpredatory %>% mutate(category = "Non-pred.")) %>%
  mutate(material = factor(material, levels = unique(material))) %>%
  arrange(material) %>%
  group_by(category) %>%
  mutate(scaled_value = scales::rescale(value, to = c(0,1))) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = rev(c("BGCs", "Secretion", "Degradation", "Pred.", "Non-pred.", "Flagellum"))))

tree = ggtree(cluster_tree) +
  #geom_tiplab(angle = 90) +
  coord_flip() +
  scale_x_reverse() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0,0,0,0),
        line = element_line(size = unit(0.1, 'mm')))

metadata_heatmap = summary_table %>%
  left_join(formatted_metadata) %>%
  select(material, `Ocean province` = ocean_province, `Size fraction` = size_fraction, cluster) %>%
  unique() %>%
  mutate(`Size fraction` = c("0.22-1.6" = "prokaryote-enriched",
                           "0.22-3" = "prokaryote-enriched",
                           "0.8->" = "particle-enriched",
                           "0.8-5" = "particle-enriched",
                           "5-20" = "particle-enriched")[`Size fraction`]) %>%
  gather(key = category, value = value, -material, -cluster) %>%
  mutate(material = factor(material, levels = unique(material))) %>%
  arrange(material) %>%
  ggplot() +
  geom_tile(aes(x = material, y = category, fill = value)) +
  facet_grid(.~cluster, scales = "free_x", space = "free") +
  scale_fill_manual(values = c(ocean_material2_colors, "prokaryote-enriched" = "#DD4968", "particle-enriched" = "#721F81FF")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size = 6),
        axis.text.y = element_blank(),
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, .5, .5, .5, "mm")),
        strip.background.x = element_rect(color = NA, fill = "black"),
        strip.background.y = element_rect(color = NA, fill = "black"),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.height = unit(.5, 'mm'))
metadata_heatmap

fun_heatmap = summary_table %>%
  ggplot() +
  geom_tile(aes(x = material, y = category, fill = scaled_value)) +
  ylab('Functional pathways') +
  xlab('Sample') +
  facet_grid(.~cluster, scales = "free_x", space = "free") +
  scale_fill_viridis(option = "D", begin = .05, end = .95,
                     guide = guide_colorbar(title = "Scaled median expression",
                                            title.position = "top",
                                            title.vjust = .9,
                                            barwidth = 4, barheight = 0.5)) +
  theme_minimal() +
  theme(#axis.text.x = element_text(angle = 90, hjust = 0, size = 6),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6),
        strip.text = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(.5, 'mm'),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6))
fun_heatmap

tree + metadata_heatmap + fun_heatmap + plot_layout(ncol = 1, heights = c(.4, .15, .45)) + plot_annotation(theme = theme(plot.margin = margin(0,0,0,0)))

ggsave("data/processed/figures/Figure-4/Figure-4D.raw.pdf", width = 70, height = 50, unit = "mm")

googlesheets4::write_sheet(summary_table, "1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "scaled-expr-values-supervised-groups")

# Statistics =============================================================================

vegan::adonis(metat_featurecounts_summary_spread %>% select(-material) ~ 
                as.character(metat_featurecounts_summary_hdbscan$cluster),
              permutations = 10000, by = 'margin',method = "euclidean")

vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                antismash_biosynthetic %>% arrange(material) %>% pull(value),
              permutations = 10000, by = 'margin',method = "euclidean")

vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                txss_flagellum %>% arrange(material) %>% pull(value),
              permutations = 10000, by = 'margin',method = "euclidean")
vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                prokka_degrad %>% arrange(material) %>% pull(value),
              permutations = 10000, by = 'margin',method = "euclidean")
vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                predatory %>% arrange(material) %>% pull(value),
              permutations = 10000, by = 'margin',method = "euclidean")
vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                txss_secretion %>% arrange(material) %>% pull(value),
              permutations = 10000, by = 'margin',method = "euclidean")

vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                antismash_biosynthetic %>% arrange(material) %>% pull(value) + 
                txss_flagellum %>% arrange(material) %>% pull(value),
              permutations = 10000, by = 'margin',method = "euclidean")


# some tests
antismash_bgcs = metat_featurecounts_summary_clustered %>%
  filter(Gene %in% (BGCs %>% filter(biosynthetic_gene | biosynthetic_add_gene) %>% pull(gene))) %>%
  left_join(BGCs, by = c("Gene" = "gene")) %>%
  group_by(biosynthetic_region, material) %>%
  summarize(n_bgc_b_genes = n(),
            value = median(value),
            cluster = unique(cluster))

antismash_bgcs_indiv = antismash_bgcs %>% select(-n_bgc_b_genes, -cluster) %>% spread(biosynthetic_region, value) %>% arrange(material) %>% select(-material)

vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                ,
              permutations = 10000, by = 'margin', method = "euclidean")

for (i in unique(antismash_bgcs$biosynthetic_region)){
  print(i)
  print(vegan::adonis(metat_featurecounts_summary_spread %>% arrange(material) %>% select(-material) ~ 
                        antismash_bgcs %>% filter(biosynthetic_region == i) %>% arrange(material) %>% pull(value),
                      permutations = 10000, by = 'margin', method = "euclidean"))
  print("")
}
