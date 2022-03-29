# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 2A - GCC overview ===============================================================

# Libraries ==============================================================================

library(tidyverse)
library(patchwork)
library(treeio) # YuLab-SMU/treeio
library(ggtree) # YuLab-SMU/ggtree
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
  filter(dataset == "go_micro_genomes") %>%
  filter(!genome %in% blacklisted_genomes) %>%
  left_join(scaffolds_length %>% dplyr::rename(scaffold_length = length)) %>%
  filter(scaffold_length >= 5000) %>% 
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_to_gcf_gcc.tsv")) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_classes.tsv")) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/gcc_abundance_prevalence.tsv")) %>%
  mutate(gcc = paste0("gcc_", gcc)) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_dists_to_refseq_mibig.tsv")) %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, species = `dRep Dereplication Cluster`, phylum = `GTDB Taxonomy`) %>% mutate(phylum = gsub(".*;p__|;c__.*", "", phylum))) %>%
  group_by(gcf, species) %>% 
  filter(bgc_id == bgc_id[length == max(length)][1]) %>%
  ungroup()

tree = "data/processed/Figures/Figure-2/2A/gcc_tree_rerooted.newick"
#tree = "data/processed/Figures/Figure-2/2A/for_itol/gcc_tree.newick"

bgc_clustering %>% select(-dataset, -folder) %>% mutate(gcf = paste0("gcf_", gcf)) %>% googlesheets4::write_sheet(., ss = "1VIETOjN8M9bjUMftr_2vJ7SnOLPe1yMEqI093FLNAW4", "BGCs clustering")

# prepare plot  ==========================================================================

gcc_tree <- drop.tip(read.newick(tree), "origin")

layer_1_6 = bgc_clustering %>%
  group_by(gcc) %>%
  summarize(n = n(),
            mean_d = mean(refseq),
            median_d = median(refseq),
            mibig_mean = mean(mibig),
            mibig_min = min(mibig),
            no_mag = all(!grepl("_METAG_", genome)),
            any_mag = any(grepl("_METAG_", genome)),
            any_new_mag = any(grepl("_METAG_[A-Z]{8}$", genome)),
            only_mag = all(grepl("_METAG_", genome)),
            only_new_mag = all(grepl("_METAG_[A-Z]{8}$", genome)))

gcc_metadata = 
  as_tibble(gcc_tree) %>%
  mutate(gcc = gsub("GCC0+", "gcc_", label)) %>%
  mutate(gcc = ifelse(gcc == "gcc_", "gcc_0", gcc)) %>%
  left_join(layer_1_6)

assertthat::are_equal(sum(grepl('gcc', gcc_metadata$gcc)), nrow(layer_1_6))

# Tree
gcc_plot = ggtree(gcc_tree, layout="fan", open.angle = 5, size = 0.15, alpha = 1) %<+%
  gcc_metadata + # add metadata to the tree data
  geom_tiplab(align=TRUE, linetype='dotted', linesize=.15, mapping = aes(label = NA)) +
  scale_color_viridis_c() +
  coord_polar(theta = 'y', start = 0, direction = -1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0,0,0,0))
gcc_plot
ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-tree.pdf", width = 60, height = 60, unit = "mm")

backbone = gcc_plot$data %>% filter(isTip)
circle_ratio = 355/360
total_elt = nrow(backbone) * (1/circle_ratio)

# Outer layers  ==========================================================================

# Layer 1: MAGs representatives ----------------------------------------------------------

novelty_colors_bluegreen = c(
  "no_mag" = "#DFF3DB",
  "any_mag"   = "#A8DCB5",
  "only_mag"     = "#43A2CA"
)

gcc_plot$data %>%
  filter(isTip) %>%
  select(y, no_mag, any_mag, only_mag) %>%
  mutate(any_mag = ifelse(only_mag, FALSE, any_mag)) %>%
  gather(key = type, value = filt, -y) %>%
  filter(filt) %>%
  ggplot() +
  geom_tile(aes(x = y, y = 2, fill = type), color = "white") +
  scale_fill_manual(values = novelty_colors_bluegreen) +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-20,20))
ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-genomes_type.pdf", width = 113, height = 113, unit = "mm")



# Layer 2: Taxonomy ----------------------------------------------------------------------

phyla_bgcs = c("Actinobacteriota", "Proteobacteria", "Firmicutes", "Cyanobacteria")
phyla_repr = c("Bacteroidota", "Thermoplasmatota", "Marinisomatota", "Chloroflexota", "Verrucomicrobiota", "Planctomycetota")

layer_taxo = bgc_clustering %>%
  mutate(phylum = ifelse(phylum %in% c(phyla_bgcs, phyla_repr), phylum, "Other")) %>%
  group_by(gcc, phylum) %>%
  #summarize(n = n()) %>%
  summarize(n = log(n()+1)) %>%
  group_by(gcc) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  left_join(gcc_plot$data %>% select(gcc, y))

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
  "Other" = "#7f7f7f"#"#c7c7c7"
)

layer_taxo %>%
  mutate(phylum = factor(phylum, levels = names(phyla_colors))) %>%
  ggplot() +
  geom_col(aes(x = y, y = 10*r, fill = phylum),color = "white", width = .75, size = .3) +
  scale_fill_manual(values = phyla_colors) +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position = "none",
        legend.position = "bottom",
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(5, "mm"),
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-20,15))
ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-taxo.pdf", width = 113, height = 113, unit = "mm")

# Layer 3: Product Classes ---------------------------------------------------------------

names(bgc_colors_serina) = c("Type I Polyketide Synthases", 
                             "RiPPs (Ribosomal Natural Products)",
                             "Terpenes",
                             "Type II/III Polyketide Synthases",
                             "Non-Ribosomal Peptide Synthetases",
                             "Other")

layer_classes = bgc_clustering %>%
  select(gcc, `Non-Ribosomal Peptide Synthetases`, `Type I Polyketide Synthases`, `Type II/III Polyketide Synthases`, `RiPPs (Ribosomal Natural Products)`, Terpenes, Other) %>%
  gather(key = class, value = distrib, -gcc) %>%
  group_by(gcc, class) %>%
  summarize(n = sum(distrib)) %>%
  group_by(gcc) %>%
  mutate(r = n/sum(n)) %>%
  ungroup() %>%
  left_join(gcc_plot$data %>% select(gcc, y)) %>%
  mutate(r = ifelse(r == 0, NA, r))

layer_classes %>%
  mutate(class = factor(class, levels =c("Non-Ribosomal Peptide Synthetases",
                                         "Type I Polyketide Synthases", 
                                         "Type II/III Polyketide Synthases",
                                         "RiPPs (Ribosomal Natural Products)",
                                         "Terpenes",
                                         "Other"))) %>%
  ggplot() +
  geom_point(aes(x = y, y = as.numeric(class), size = 10*r, fill = class, color = class), shape = 21) +
  scale_fill_manual(values = bgc_colors_serina) +
  scale_color_manual(values = bgc_colors_serina) +
  scale_size_continuous(range = c(.2, 1)) +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = .15),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,0,0,0),
        line = element_line(size = .15)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5, 6), limits = c(-39,6))
ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-classes.pdf", width = 113, height = 113, unit = "mm")


layer_classes %>%
  group_by(class) %>%
  summarize(s = sum(r, na.rm = T),
            n = sum(r > 0, na.rm = T)) %>%
  arrange(desc(n))

# Layer 4: Number of BGCs  ---------------------------------------------------------------

layer_bgcs = bgc_clustering %>%
  group_by(gcc, prevalence) %>%
  summarize(n = log(n() + 1),
            number_raw = n()) %>%
  left_join(gcc_plot$data %>% select(gcc, y))

layer_bgcs %>%
  ggplot() +
  geom_tile(aes(x = y, y = 15, fill = n), color = "white") +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position = "none",
        legend.position = "bottom",
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(5, "mm"),
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-20,16))
ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-n_bgcs.pdf", width = 123, height = 123, unit = "mm")

# Layer 5: GCC prevalence  ---------------------------------------------------------------

layer_prevalence = bgc_clustering %>%
  select(gcc, prevalence) %>%
  unique() %>%
  left_join(gcc_plot$data %>% select(gcc, y))

layer_prevalence %>%
  ggplot() +
  geom_tile(aes(x = y, y = 10, fill = prevalence), color = "white") +
  scale_fill_gradient(low = "#ece7f2", high = "#0570b0") +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position = "none",
        legend.position = "bottom",
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(5, "mm"),
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-13.3,10.6))

ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-prevalence.pdf", width = 125, height = 125, unit = "mm")


# Layer 6: Distance to RefSeq ------------------------------------------------------------

gcc_plot$data %>%
  filter(isTip) %>%
  ggplot() +
  geom_tile(aes(x = y, y = 7, fill = mean_d), color = "white") +
  scale_fill_gradient(low = "#FEE5D9", high = "#CB181D") +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position = "none",
        legend.position = "bottom",
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(5, "mm"),
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-10,8))

ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-distance.pdf", width = 130, height = 130, unit = "mm")

layer_1_6 %>% filter(mean_d > 0.4) %>% nrow()

# Layer 7: MIBIG -------------------------------------------------------------------------

layer_1_6 %>% pull(mibig_mean) %>% sort
layer_1_6 %>% pull(mibig_min) %>% sort
layer_1_6 %>% pull(mibig_min) %>% sort < 0.01
layer_1_6 %>% pull(mibig_min) %>% sort < 0.1
layer_1_6 %>% arrange(mibig_min) %>% View()

# MIBiG matches: 
# gcc_22: BGC0000939 (aerobactin), BGC0001499 (aerobactin), BGC0001870 (putrebactin, avaroferrin), BGC0001408 (xanthoferrin) --> Siderophores
# gcc_123: BGC0000858 (Ectoine)
# gcc_48: BGC0000577 (sphingonodin II), BGC0000573 (caulonodin I/II)
# gcc_135: BGC0001601 (gassericin-S)

go_micro_antismash = read_tsv(go_micro_antismash_file)
bgc_clustering_w_products = bgc_clustering %>%
  mutate(region = paste0(scaffold, "-biosynth_", gsub(".*region0+|.gbk", "", file))) %>%
  left_join(go_micro_antismash %>% select(region, products))

# siderophores
bgc_clustering_w_products %>% filter(grepl("siderophore", products)) %>% pull(gcc) %>% table %>% sort
layer_bgcs %>% filter(gcc == "gcc_22")

# ectoine
bgc_clustering_w_products %>% filter(grepl("ectoine", products)) %>% pull(gcc) %>% table %>% sort
layer_bgcs %>% filter(gcc == "gcc_123")

# hglE
bgc_clustering_w_products %>% filter(grepl("hglE", products)) %>% pull(gcc) %>% table %>% sort
layer_bgcs %>% filter(gcc == "gcc_92")

bgc_clustering_w_products %>% filter(gcc == "gcc_92") %>% pull(products) %>% table %>% sort

# arylpolyene
bgc_clustering_w_products %>% filter(grepl("arylpolyene", products)) %>% pull(gcc) %>% table %>% sort
layer_bgcs %>% filter(gcc == "gcc_41")

# PUFA
bgc_clustering_w_products %>% filter(grepl("PUFA", products)) %>% pull(gcc) %>% table %>% sort

# Prevalent ones
layer_prevalence %>% arrange(desc(prevalence))
prev_gcc = layer_prevalence %>% filter(prevalence > 0.9) %>% pull(gcc)
layer_taxo %>% group_by(gcc) %>% summarize(n = sum(n > 0)) %>% arrange(desc(n))
wide_gcc = layer_taxo %>% group_by(gcc) %>% summarize(n = sum(n > 0)) %>% filter(n >= 6) %>% pull(gcc)
both_gcc = prev_gcc[prev_gcc %in% wide_gcc]
layer_1_6 %>% arrange(mibig_min) %>% filter(gcc %in% both_gcc) %>% View()

bgc_clustering_w_products %>% filter(gcc == "gcc_14") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_92") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_13") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_20") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_77") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_55") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_24") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_84") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_21") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_9") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_17") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_8") %>% pull(products) %>% table %>% sort
bgc_clustering_w_products %>% filter(gcc == "gcc_29") %>% pull(products) %>% table %>% sort

gcc_plot$data %>%
  filter(isTip) %>%
  mutate(plot = ifelse(mibig_min < 0.01, 15, NA)) %>%
  ggplot() +
  geom_tile(aes(x = y, y = 14, fill = mean_d), color = "white") +
  geom_point(aes(x = y, y = plot), shape = 21, fill = "black", size = 1) +
  geom_text(aes(x = y, y = plot, label = gcc), size = 3) +
  scale_fill_gradient(low = "#FEE5D9", high = "#CB181D") +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-40,16))

gcc_plot$data %>%
  filter(isTip) %>%
  mutate(plot = ifelse(mibig_min < 0.1 & mibig_min > 0.01, 15, NA)) %>%
  ggplot() +
  geom_tile(aes(x = y, y = 14, fill = mean_d), color = "white") +
  geom_point(aes(x = y, y = plot), shape = 21, fill = "black", size = 1) +
  geom_text(aes(x = y, y = plot, label = gcc), size = 3) +
  scale_fill_gradient(low = "#FEE5D9", high = "#CB181D") +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-40,16))

gcc_plot$data %>%
  filter(isTip) %>%
  mutate(plot = ifelse(gcc %in% c("gcc_92", "gcc_41"), 15, NA)) %>%
  ggplot() +
  geom_tile(aes(x = y, y = 14, fill = mean_d), color = "white") +
  geom_point(aes(x = y, y = plot), shape = 21, fill = "black", size = 1) +
  geom_text(aes(x = y, y = plot, label = gcc), size = 3) +
  scale_fill_gradient(low = "#FEE5D9", high = "#CB181D") +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.5*size_converter),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, total_elt) +
  scale_y_continuous(breaks = NULL, limits = c(-40,16))

ggsave("data/processed/Figures/Figure-2/2A-full_circle/Figure-2A-mibig.pdf", width = 130, height = 130, unit = "mm")

# Check a bit more

# Previous ones
GCC00022 = c("BGC0000939", "BGC0000947", "BGC0001499", "BGC0000942", "BGC0001870", "BGC0001572", "BGC0000295", "BGC0001408", "BGC0000946", "BGC0001478") # Siderophores
GCC00048 = c("BGC0001786", "BGC0000577", "BGC0000573", "BGC0000570", "BGC0000580", "BGC0000574") # Ripp
GCC00135 = c("BGC0001601") # RiPP
GCC00123 = c("BGC0000855", "BGC0000857", "BGC0000853", "BGC0000858", "BGC0000856", "BGC0000859", "BGC0000854", "BGC0000860") # Ectoine

# New ones, interesting
GCC00014 = c("BGC0000635", "BGC0000640", "BGC0000656", "BGC0000647", "BGC0000650") # Carotenoid
GCC00092 = c("BGC0000837", "BGC0000284", "BGC0001760", "BGC0000067", "BGC0000041", "BGC0001273", "BGC0000056", "BGC0001047", "BGC0001275", "BGC0000838", "BGC0001631", "BGC0000467", "BGC0001050", "BGC0001124", "BGC0001163", "BGC0001276", "BGC0001089", "BGC0000231", "BGC0002008", "BGC0001560", "BGC0000299")
# GCC00092 --> APE Vf (Arylpolyene), phenolic lipids (T3/hglE-KS), flexirubin (arylpolyene/resorcinol), aryl polyenes
# GCC00092 --> oronofacic acid (T1KS), asperlactone (T1PKS) + 5
# GCC00092 --> 10*NRPS or NRPS hybrids

# New ones, Nothing too interesting, they are because of the REFs
GCC00017 = c("BGC0001844", "BGC0000351", "BGC0000400", "BGC0001479", "BGC0000352", "BGC0001135", "BGC0001890") # paenibacterin, amphi-enterobactin, ... (NRPS)
GCC00170 = c("BGC0000528", "BGC0000554", "BGC0001579", "BGC0001660") # Lanthipeptide
GCC00052 = c("BGC0001181", "BGC0001642") # geosmin, koraiol (Terpenes)
GCC00104 = c("BGC0000616") # amylocyclicin (RiPP)
GCC00088 = c("BGC0000538") # nisin Z (lantibiotic)
GCC00057 = c("BGC0000602") # subtilosin A RiPP

gcc_to_look_at = GCC00123
for (i in gcc_to_look_at){
  browseURL(paste0("https://mibig.secondarymetabolites.org/repository/", i, "/index.html#r1c1"))
  Sys.sleep(1)
}

