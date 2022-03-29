# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 2B - GCF Novelty ================================================================

# Libraries ==============================================================================

library(tidyverse)
library(patchwork)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# load data ==============================================================================

blacklisted_genomes = read_lines("data/raw/go_microbiomics/summaries/go_microbiomics-blacklisted_genomes.list")

scaffolds_length = read_tsv("https://microbiomics.io/ocean/db/1.0/antismash-output/genomes_bgcs-scaffolds_length.tsv")

genomes_summary = load_prettified_summary()

bgc_clustering = read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_metadata.tsv") %>%
  filter(dataset == "go_micro_genomes") %>% #pull(contig_edge) %>% table
  filter(!genome %in% blacklisted_genomes) %>%
  left_join(scaffolds_length %>% dplyr::rename(scaffold_length = length)) %>%
  filter(scaffold_length >= 5000) %>% #pull(contig_edge) %>% table
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_to_gcf_gcc.tsv")) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_classes.tsv")) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/gcc_abundance_prevalence.tsv")) %>%
  mutate(gcc = paste0("gcc_", gcc), gcf = paste0("gcf_", gcf)) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_dists_to_refseq_mibig.tsv")) %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, species = `dRep Dereplication Cluster`, phylum = `GTDB Taxonomy`) %>% mutate(phylum = gsub(".*;p__|;c__.*", "", phylum))) %>%
  group_by(gcf, species) %>% 
  filter(bgc_id == bgc_id[length == max(length)][1]) %>%
  #filter(bgc_id == bgc_id[contig_edge == min(contig_edge)][1]) %>%
  ungroup() #%>% pull(contig_edge) %>% table

# Plots ==================================================================================

gcf_dists = bgc_clustering %>%
  group_by(gcf) %>%
  summarize(n_bgcs = n(),
            RefSeq = mean(refseq),
            MIBiG = mean(mibig),
            `Non-Ribosomal Peptide Synthetases` = sum(`Non-Ribosomal Peptide Synthetases`),
            `Type I Polyketide Synthases` = sum(`Type I Polyketide Synthases`),
            `Type II/III Polyketide Synthases` = sum(`Type II/III Polyketide Synthases`),
            `RiPPs (Ribosomal Natural Products)` = sum(`RiPPs (Ribosomal Natural Products)`),
            Terpenes = sum(Terpenes),
            Other = sum(Other),
            phylum = paste0(phylum, collapse = ";"))

nrow(gcf_dists)
gcf_dists %>% filter(RefSeq > 0.2) %>% nrow()
gcf_dists %>% filter(RefSeq > 0.2) %>% nrow() / nrow(gcf_dists)
gcf_dists %>% filter(MIBiG > 0.2) %>% nrow()
gcf_dists %>% filter(MIBiG > 0.2) %>% nrow() / nrow(gcf_dists)


# BGC Class ------------------------------------------------------------------------------

gcf_class = gcf_dists %>%
  gather(key = class, value = distrib, -gcf, -RefSeq, -MIBiG, -n_bgcs, -phylum) %>%
  group_by(gcf, RefSeq, MIBiG, class) %>%
  summarize(n = sum(distrib)) %>%
  group_by(gcf) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  filter(r != 0) %>%
  gather(key = db, value = distance, -gcf, -class, -n, -r) %>%
  mutate(db = factor(db, levels = c("RefSeq", "MIBiG")))

gcf_class %>% filter(db == "RefSeq" & distance > 0.2) %>% pull(gcf) %>% unique %>% length()
gcf_class %>% filter(db == "RefSeq" & distance > 0.2) %>% filter(class %in% c("RiPPs (Ribosomal Natural Products)", "Terpenes", "Other")) %>% pull(gcf) %>% unique %>% length()
3000/3861

# GCF Phyla distrib ----------------------------------------------------------------------

phyla_bgcs = c("Actinobacteriota", "Proteobacteria", "Firmicutes", "Cyanobacteria")
phyla_repr = c("Bacteroidota", "Thermoplasmatota", "Marinisomatota", "Chloroflexota", "Verrucomicrobiota", "Planctomycetota")

phyla_colors = c(
  "Actinobacteriota" = "#f6bd82",
  "Proteobacteria" = "#b3c7e5",
  "Firmicutes" = "#c2b2d3",
  "Cyanobacteria" = "#a8dc93",
  "Bacteroidota" = "#f29d99",
  "Thermoplasmatota" = "#be9e96",
  "Marinisomatota" = "#eeb9d1",
  "Chloroflexota" = "#dcda96",
  "Verrucomicrobiota" = "#aad9e3",
  "Planctomycetota" = "#f9e8b9",
  "Other phyla" = "#7f7f7f"#"#c7c7c7"
)

gcf_phyla = gcf_dists %>%
  separate_rows(phylum, sep = ";") %>%
  mutate(phylum = ifelse(phylum %in% c(phyla_bgcs, phyla_repr), phylum, "Other phyla")) %>%
  group_by(gcf, RefSeq, MIBiG, phylum) %>%
  summarize(n = n()) %>%
  group_by(gcf) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  filter(r != 0) %>%
  gather(key = db, value = distance, -gcf, -phylum, -n, -r) %>%
  mutate(db = factor(db, levels = c("RefSeq", "MIBiG"))) %>%
  mutate(phylum = factor(phylum, levels = names(phyla_colors)))

gcf_phyla %>% filter(db == "RefSeq" & distance > 0.2) %>% pull(gcf) %>% unique %>% length()
gcf_phyla %>% filter(db == "RefSeq" & distance > 0.2) %>% filter(!phylum %in% phyla_bgcs) %>% pull(gcf) %>% unique %>% length()
1816/3861


# combined for plot ----------------------------------------------------------------------

assertthat::assert_that(all(names(bgc_colors_serina) %in% unique(gcf_class$class)))

gcf_phyla %>% rename(type = phylum) %>% mutate(facet = "Phyla") %>%
  rbind(gcf_class %>% rename(type = class) %>% mutate(facet = "BGC class")) %>%
  mutate(type = factor(type, levels = c(names(phyla_colors),
                                        c("Non-Ribosomal Peptide Synthetases",
                                          "Type I Polyketide Synthases", 
                                          "Type II/III Polyketide Synthases",
                                          "RiPPs (Ribosomal Natural Products)",
                                          "Terpenes",
                                          "Other")))) %>%
  ggplot() +
  #geom_density(aes(x = distance, y = ..scaled.., fill = type, color = type), color = "white", position = "stack", alpha = .8, bw = 0.02) +
  geom_histogram(aes(x = distance, fill = type, color = type), color = "white", bins = 30, size =.1) +
  #geom_histogram(aes(x = distance, fill = type, color = type), bins = 40, aplha = .8) +
  geom_vline(xintercept = 0.2, size = 0.3) +
  xlab("Avg. best cosine distance (d)")+
  ylab("GCF count") +
  scale_fill_manual(values = c(bgc_colors_serina, phyla_colors)) +
  scale_color_manual(values = c(bgc_colors_serina, phyla_colors)) +
  facet_grid(facet~db) +
  theme_bw() +
  theme(line = element_line(size = .25), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = .25),
        axis.ticks.length = unit(.5, 'mm'),
        legend.position = 'none',
        strip.placement = "outside",
        strip.text = element_text(color = "white", face = "bold", size = 6, margin = margin(.5, .5, .5, .5, "mm")),
        strip.background = element_rect(color = NA, fill = "black"),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(1, 'mm'),
        panel.grid.major = element_line(size = unit(0.25, 'pt')))
ggsave("data/processed/Figures/Figure-2/Figure-2B.pdf", width = 89, height = 52, unit = "mm")

# Number of New GCF in euclidean clustering ----------------------------------------------

gcf_eucl = read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/bgc_to_gcf_euclidean.tsv") %>% rename(gcf_eucl = gcf)
gcf_eucl %>% pull(gcf_eucl) %>% unique() %>% length()

bgc_clustering_w_eucl = bgc_clustering %>% left_join(gcf_eucl)
bgc_clustering_w_eucl %>% pull(gcf_eucl) %>% unique %>% length

bgc_clustering_w_eucl %>%
  filter(gcf %in% (gcf_dists %>% filter(RefSeq > 0.2) %>% pull(gcf))) %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length()

bigfam_gcf_eucl = read_tsv("data/raw/go_microbiomics/bgcs/bigfam/all_bigfam_bgc_gcf_membership.tsv") %>% rename(gcf_eucl = gcf)
bigfam_gcf_eucl %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length()

2877/29955

bgc_clustering_w_eucl %>%
  filter(`RiPPs (Ribosomal Natural Products)`) %>%
  filter(gcf %in% (gcf_dists %>% filter(RefSeq > 0.2) %>% pull(gcf))) %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length() /
read_tsv("data/raw/go_microbiomics/bgcs/bigfam/all_bigfam_bgc_classes.tsv") %>%
  left_join(bigfam_gcf_eucl) %>%
  filter(`RiPPs (Ribosomal Natural Products)`) %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length()

bgc_clustering_w_eucl %>%
  filter(`Terpenes`) %>%
  filter(gcf %in% (gcf_dists %>% filter(RefSeq > 0.2) %>% pull(gcf))) %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length() /
  read_tsv("data/raw/go_microbiomics/bgcs/bigfam/all_bigfam_bgc_classes.tsv") %>%
  left_join(bigfam_gcf_eucl) %>%
  filter(`Terpenes`) %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length()

bgc_clustering_w_eucl %>%
  filter(`Other`) %>%
  filter(gcf %in% (gcf_dists %>% filter(RefSeq > 0.2) %>% pull(gcf))) %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length() /
  read_tsv("data/raw/go_microbiomics/bgcs/bigfam/all_bigfam_bgc_classes.tsv") %>%
  left_join(bigfam_gcf_eucl) %>%
  filter(`Other`) %>%
  pull(gcf_eucl) %>%
  unique() %>%
  length()

# Class enrichment compared to bigFAM ?? ----------

bigfam_gcf_eucl_class = read_tsv("data/raw/go_microbiomics/bgcs/bigfam/all_bigfam_bgc_classes.tsv") %>%
  left_join(bigfam_gcf_eucl) %>%
  group_by(gcf_eucl) %>%
  summarize(n_bgcs = n(),
            `Non-Ribosomal Peptide Synthetases` = sum(`Non-Ribosomal Peptide Synthetases`),
            `Type I Polyketide Synthases` = sum(`Type I Polyketide Synthases`),
            `Type II/III Polyketide Synthases` = sum(`Type II/III Polyketide Synthases`),
            `RiPPs (Ribosomal Natural Products)` = sum(`RiPPs (Ribosomal Natural Products)`),
            Terpenes = sum(Terpenes),
            Other = sum(Other)) %>%
  gather(key = class, value = distrib, -gcf_eucl, -n_bgcs) %>%
  group_by(gcf_eucl, class) %>%
  summarize(n = sum(distrib)) %>%
  group_by(gcf_eucl) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  group_by(class) %>%
  summarize(proportion = sum(r)) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  arrange(desc(proportion))

gom_eucl_class = bgc_clustering_w_eucl %>%
  group_by(gcf_eucl) %>%
  summarize(n_bgcs = n(),
            `Non-Ribosomal Peptide Synthetases` = sum(`Non-Ribosomal Peptide Synthetases`),
            `Type I Polyketide Synthases` = sum(`Type I Polyketide Synthases`),
            `Type II/III Polyketide Synthases` = sum(`Type II/III Polyketide Synthases`),
            `RiPPs (Ribosomal Natural Products)` = sum(`RiPPs (Ribosomal Natural Products)`),
            Terpenes = sum(Terpenes),
            Other = sum(Other)) %>%
  gather(key = class, value = distrib, -gcf_eucl, -n_bgcs) %>%
  group_by(gcf_eucl, class) %>%
  summarize(n = sum(distrib)) %>%
  group_by(gcf_eucl) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  group_by(class) %>%
  summarize(proportion = sum(r)) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  arrange(desc(proportion))

gom_class = bgc_clustering_w_eucl %>%
  group_by(gcf) %>%
  summarize(n_bgcs = n(),
            `Non-Ribosomal Peptide Synthetases` = sum(`Non-Ribosomal Peptide Synthetases`),
            `Type I Polyketide Synthases` = sum(`Type I Polyketide Synthases`),
            `Type II/III Polyketide Synthases` = sum(`Type II/III Polyketide Synthases`),
            `RiPPs (Ribosomal Natural Products)` = sum(`RiPPs (Ribosomal Natural Products)`),
            Terpenes = sum(Terpenes),
            Other = sum(Other)) %>%
  gather(key = class, value = distrib, -gcf, -n_bgcs) %>%
  group_by(gcf, class) %>%
  summarize(n = sum(distrib)) %>%
  group_by(gcf) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  group_by(class) %>%
  summarize(proportion = sum(r)) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  arrange(desc(proportion))

bigfam_gcf_eucl_class
gom_eucl_class
gom_class

diagnostic_plots = bigfam_gcf_eucl_class %>% rename(prop_bigfam = proportion) %>%
  left_join(gom_eucl_class %>% rename(prop_eucl = proportion)) %>%
  left_join(gom_class %>% rename(prop_cosine = proportion)) %>%
  mutate(class = factor(class, levels = c(c("Non-Ribosomal Peptide Synthetases",
                                            "Type I Polyketide Synthases", 
                                            "Type II/III Polyketide Synthases",
                                            "RiPPs (Ribosomal Natural Products)",
                                            "Terpenes",
                                            "Other"))))

p1 = diagnostic_plots %>%
  mutate(y = prop_eucl - prop_bigfam) %>%
  ggplot() +
  geom_col(aes(x = class, y = 100*y, fill = class)) +
  scale_fill_manual(values = bgc_colors_serina) +
  ggtitle("A - OMD vs BiGFAM BGC classes") +
  ylab("Difference (%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

p1_2 = diagnostic_plots %>%
  mutate(y = prop_eucl / prop_bigfam) %>%
  ggplot() +
  geom_col(aes(x = class, y = y, fill = class)) +
  scale_fill_manual(values = bgc_colors_serina) +
  ggtitle("A - OMD vs BiGFAM BGC classes") +
  ylab("Fold increase") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
p1_2

p2 = diagnostic_plots %>%
  mutate(y = prop_cosine - prop_eucl) %>%
  ggplot() +
  geom_col(aes(x = class, y = 100*y, fill = class)) +
  scale_fill_manual(values = bgc_colors_serina) +
  ggtitle("B - Impact of eucl. vs cosine clustering BGC classes") +
  ylab("Difference (%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

p1 | p2

# gom classes gcc level

gom_class_gcc = bgc_clustering_w_eucl %>%
  group_by(gcc) %>%
  summarize(n_bgcs = n(),
            `Non-Ribosomal Peptide Synthetases` = sum(`Non-Ribosomal Peptide Synthetases`),
            `Type I Polyketide Synthases` = sum(`Type I Polyketide Synthases`),
            `Type II/III Polyketide Synthases` = sum(`Type II/III Polyketide Synthases`),
            `RiPPs (Ribosomal Natural Products)` = sum(`RiPPs (Ribosomal Natural Products)`),
            Terpenes = sum(Terpenes),
            Other = sum(Other)) %>%
  gather(key = class, value = distrib, -gcc, -n_bgcs) %>%
  group_by(gcc, class) %>%
  summarize(n = sum(distrib)) %>%
  group_by(gcc) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  group_by(class) %>%
  summarize(proportion = sum(r)) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  arrange(desc(proportion))

gom_class_gcc
