# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 7B - Eudoremicrobium gene content ===================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(DESeq2)
#source('code/figures/figure4D-dataprep.R')

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

# Featurecounts-based analysis ===========================================================

# MetaG detections featurecounts ---------------------------------------------------------

metag_featurecounts_summary_raw = read_tsv("data/raw/marine_eremios/marine_eremios-featurecounts-summary_PIAMPJPB.tsv", skip = 1) %>%
  select(!contains("ETHSEQ")) %>% 
  select(!contains("METAT")) %>%
  select(!contains("_T")) %>%
  gather(key = Sample, value = inserts, -c(Geneid, Chr, Start, End, Strand, Length)) %>%
  mutate(Feature = Geneid,
         Gene = gsub("^[^-]*-", "", Geneid),
         Genome = gsub("-.*", "", Geneid),
         Sample = gsub(".filter.*|.*/", "", Sample)) %>%
  select(Sample, Genome, Chr, Feature, Gene, Start, End, Strand, Length, inserts)

all_metag_samples = googlesheets4::read_sheet("1JDgQoi8bY5N8L_HokD2R2zWay3ACtbwj4YMbyQzdyMw") %>%
  filter(grepl("METAG$", representative_barcode))

pangaea_stream = tara_metadata_from_pangea(type = "Env_meso") %>%
  select(barcode = `Sample ID (TARA_barcode#, registered at ...)`,
         ocean_province = `OS region ([abbreviation] full name (MRG...)`,
         station = `Station (TARA_station#, registered at ...)`,
         material = `Sample material (TARA_station#_environmental-f...)`,
         depth = `Depth, nominal (from which this sample was co...)`) %>%
  mutate(size_fraction = gsub(".*_", "", material))

metag_metadata = all_metag_samples %>%
  select(representative_barcode) %>%
  mutate(representative_barcode = gsub("_METAG", "", representative_barcode)) %>%
  left_join(pangaea_stream, by = c("representative_barcode" = "barcode")) %>% 
  filter(!duplicated(representative_barcode)) %>%
  mutate(Sample = paste0(representative_barcode, "_METAG"))

metag_featurecounts_summary = metag_featurecounts_summary_raw %>%
  mutate(MGs = Gene %in% (MGs %>% filter(`mOTUs-10`) %>% pull(bestMGs))) %>%
  left_join(metag_metadata) %>%
  group_by(Sample) %>%
  summarize(n_genes = n(),
            n_genes_detected = sum(inserts > 0),
            n_mgs_detected = sum(inserts[MGs] > 0),
            mgs_median = median(inserts[MGs]),
            material = unique(material),
            ocean_province = unique(ocean_province),
            size_fraction = unique(size_fraction),
            depth = unique(depth))

metag_featurecounts_summary %>%
  filter(n_mgs_detected >= 10 & mgs_median >= 10) %>% #View()
  mutate(ocean_province = gsub(".*\\] | \\(.*", "", ocean_province),
         size_fraction = c("0.22-0.45" = "prokaryote-enriched",
                           "0.45-0.8" = "prokaryote-enriched",
                           "0.22-3" = "prokaryote-enriched",
                           "0.8->" = "particle-enriched",
                           "0.8-3" = "particle-enriched",
                           "0.8-5" = "particle-enriched",
                           "5-20" = "particle-enriched")[size_fraction]) %>%
  ggplot() +
  geom_hline(yintercept = unique(metag_featurecounts_summary$n_genes)) +
  geom_point(aes(x = mgs_median, y = n_genes_detected, color = ocean_province, fill = ocean_province, shape = size_fraction), size = 3, alpha = .8) +
  scale_color_manual(values = ocean_material2_colors, name = "Ocean Provinces") +
  scale_fill_manual(values = ocean_material2_colors, name = "Ocean Provinces") +
  scale_shape_manual(values = c("prokaryote-enriched" = 21, "particle-enriched" = 25), name = "Size Fraction (um)") +
  scale_x_log10() +
  xlab("Median insert count across marker genes (log)") +
  ylab("Number of genes detected") +
  theme_bw()

ggsave("data/processed/figures/Figure-SX-Eudore_gene_content.pdf", width = 183, height = 120, unit = 'mm')  
