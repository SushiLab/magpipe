# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 2D - Biosynthetic div in BGC-rich lineages ======================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Load Data ------------------------------------------------------------------------------

# Assumes prep data has been run --> should be good now

# Prepare Data ---------------------------------------------------------------------------

go_micro_antismash_file = "data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-antismash_summary.tsv.gz"
eremio_antismash_file = "data/raw/marine_eremios/eremiobacterota_superproducer-integrated-cpl50_ctn10-antismash.filtered.tsv"

antismash_table = read_tsv(go_micro_antismash_file) %>%
  left_join(scaffolds_length %>% dplyr::rename(scaffold_length = length)) %>%
  filter(scaffold_length >= 5000) %>%
  filter(!genome %in% blacklisted_genomes) %>%
  select(-scaffold_length) %>%
  select(-contains("CDS"))

# Then we want to replace Eremio with the representative values:
eremio_antismash = read_tsv(eremio_antismash_file) %>%
  filter(genome == "MALA_SAMN05422137_METAG_HLLJDLBE") %>% # select representative
  mutate(genome = "MALA_SAMN05422189_METAG_HFLHJDGN") %>% # mock name
  select(-contains("CDS"))

antismash_table_fixed = antismash_table %>%
  filter(!genome %in% eremio_antismash$genome) %>%
  rbind(eremio_antismash) %>%
  separate_rows(`protoclusters products`, sep = ';') %>%
  group_by(genome, region) %>%
  summarize(products = paste0(sort(`protoclusters products`), collapse = ";"))

antismash_summary = antismash_table_fixed %>%
  group_by(genome) %>%
  summarize(`# Biosynthetic Regions` = length(unique(region)),
            `# Biosynthetic Products` = n())

genomes_of_interest = antismash_summary %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, `dRep Dereplication Cluster`, `Representative Genome`, `GTDB Taxonomy`)) %>%
  filter(`# Biosynthetic Regions` > 15) %>%
  group_by(`dRep Dereplication Cluster`) %>%
  filter(`# Biosynthetic Regions` == max(`# Biosynthetic Regions`)) %>%
  pull(genome)

data_to_plot = antismash_table_fixed  %>%
  group_by(genome, products) %>%
  summarize(n = n()) %>%
  filter(genome %in% genomes_of_interest) %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, `dRep Dereplication Cluster`, `Representative Genome`, `GTDB Taxonomy`)) %>%
  left_join(antismash_summary %>% select(genome, `# Biosynthetic Regions`, `# Biosynthetic Products`)) %>%
  group_by(genome, `dRep Dereplication Cluster`, `GTDB Taxonomy`, `# Biosynthetic Regions`, `# Biosynthetic Products`) %>%
  summarize(diversity = vegan::diversity(n, index = 'shannon'),
            richness = sum(n > 0)) %>%
  group_by(`dRep Dereplication Cluster`) %>%
  filter(diversity == max(diversity)) %>%
  ungroup() %>%
  arrange(desc(diversity)) %>%
  filter(!duplicated(`dRep Dereplication Cluster`))

View(data_to_plot)

googlesheets4::write_sheet(data_to_plot, "1AuxYbYrFQpVsxnvUPTs_M-Q9LSwZH32LhbpaQalP_wE", sheet = "BGC-rich lineages")

data_to_plot_filt = data_to_plot %>%
  filter(diversity >= 1.5)

data_to_plot_filt %>% pull(`GTDB Taxonomy`)

species_name_dict = c(
  "d__Bacteria;p__Eremiobacterota;c__UBP9;o__UBA4705;f__;g__;s__" = "Unknown Eremiobacterota",
  "d__Bacteria;p__Myxococcota;c__Polyangia;o__Polyangiales;f__Sandaracinaceae;g__;s__" = "Unknown Sandaracinaceae",
  "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus velezensis" = "Bacillus velezensis",
  "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus;s__Rhodococcus sp000813105" = "Rhodococcus sp000813105",
  "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus;s__Rhodococcus sp001942265" = "Rhodococcus sp001942265",
  "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus;s__Rhodococcus sp000333955" = "Rhodococcus sp000333955",
  "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Rhodococcus;s__Rhodococcus qingshengii" = "Rhodococcus qingshengii",
  "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis" = "Bacillus subtilis",
  "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Tistrellales;f__Tistrellaceae;g__Tistrella;s__Tistrella mobilis" = "Tistrella mobilis",
  "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacteroides;s__Mycobacteroides chelonae" = "Mycobacteroides chelonae",
  "d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Gordonia;s__Gordonia sp002009645" = "Gordonia sp002009645",
  "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Microcystaceae;g__Crocosphaera;s__Crocosphaera watsonii" = "Crocosphaera watsonii",
  "d__Bacteria;p__Planctomycetota;c__UBA8742;o__;f__;g__;s__" = "Unknown Planctomycetota",
  "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Alteromonadaceae;g__Pseudoalteromonas;s__Pseudoalteromonas elyakovii|d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Alteromonadaceae;g__Pseudoalteromonas;s__Pseudoalteromonas sp002850255" = "Pseudoalteromonas elyakovii",
  "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Ketobacteraceae;g__Ketobacter;s__Ketobacter sp002471665" = "Ketobacter sp002471665",
  "d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Roseovarius;s__Roseovarius indicus" = "Roseovarius indicus"
)

data_to_plot_filt %>%
  mutate(species = species_name_dict[`GTDB Taxonomy`]) %>%
  mutate(species = factor(species, levels = species)) %>%
  ggplot() +
  geom_col(aes(x = species, y = diversity), fill = "#737373") +
  ylab("Biosynthetic diversity") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(line = element_line(size = unit(.2, 'mm')), 
        text = element_text(size = 6),
        rect = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(5, 0, 0, 0)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = unit(.2, 'mm')), 
        axis.ticks.length = unit(.3, 'mm'),
        panel.spacing = unit(0, 'mm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(size = unit(.2, 'mm')),
        plot.margin = margin(5,1,1,1, 'mm'),
        plot.background = element_rect(color = "black", size = unit(0.2, 'mm'), fill = "white"))

ggsave(paste0(figures_path_proj, "Figure-3/Figure-3.inset.pdf"), width = 45, height = 55, units = 'mm')
