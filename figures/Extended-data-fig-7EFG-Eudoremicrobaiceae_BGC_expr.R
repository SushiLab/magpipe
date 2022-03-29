# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 7EFG - BGC expression Ca. Eudoremicrobiaceae ========================

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(viridis)
#source('code/figures/figure4D-dataprep.R') # expects it to have been run

source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Data preparation =======================================================================

get_bgc_expression <- function(
  genome_id = "HLLJDLBE",
  raw_prefix = "data/raw/marine_eremios/Eudoremicrobiaceae-",
  raw_suffix = "-profiling-raw.tsv", # make sure files are not compressed
  MGs_genes = "data/processed/marine_eremios/distributions/fetchmgs_10_genes.list",
  formatted_metadata = metat_metadata,
  antismash_summary = "data/raw/marine_eremios/eremiobacterota_superproducer-integrated-cpl50_ctn10-antismash.filtered.tsv"
) {
  
  # Load convenient gene lists
  MGs = read_lines(MGs_genes) %>% gsub("^[^-]*-", "", .)
  
  # Nicely format antismash table
  antismash_table = read_tsv(antismash_summary) %>% filter(grepl(genome_id, genome))
  
  cds_biosynthetic = antismash_table %>% select(`CDS biosynthetic`) %>% separate_rows(`CDS biosynthetic`, sep = ";")
  cds_biosynthetic_add = antismash_table %>% select(`CDS biosynthetic additional`) %>% separate_rows(`CDS biosynthetic additional`, sep = ";")
  cds_transport = antismash_table %>% select(`CDS transport`) %>% separate_rows(`CDS transport`, sep = ";")
  cds_regulatory = antismash_table %>% select(`CDS regulatory`) %>% separate_rows(`CDS regulatory`, sep = ";")
  cds_resistance = antismash_table %>% select(`CDS resistance`) %>% separate_rows(`CDS resistance`, sep = ";")
  
  annotations_antismash = antismash_table %>%
    select(scaffold, biosynthetic_region = region, biosynthetic_products = products, gene = `CDS list`) %>%
    separate_rows(gene, sep = ";") %>%
    select(gene, scaffold, biosynthetic_region, biosynthetic_products) %>%
    mutate(biosynthetic_gene = gene %in% cds_biosynthetic$`CDS biosynthetic`,
           biosynthetic_add_gene = gene %in% cds_biosynthetic_add$`CDS biosynthetic additional`,
           transport_gene = gene %in% cds_transport$`CDS transport`,
           regulatory_gene = gene %in% cds_regulatory$`CDS regulatory`,
           resistance_gene = gene %in% cds_resistance$`CDS resistance`)
  
  # Load and prepare raw profile
  metat_featurecounts_summary_raw = read_tsv(paste0(raw_prefix, genome_id, raw_suffix), skip = 1) %>%
    select(!contains("ETHSEQ")) %>% 
    select(!contains("METAG")) %>%
    gather(key = Sample, value = inserts, -c(Geneid, Chr, Start, End, Strand, Length)) %>%
    mutate(Feature = Geneid,
           Gene = gsub("^[^-]*-", "", Geneid),
           Genome = gsub("-.*", "", Geneid),
           Sample = gsub(".filter.*|.*/", "", Sample)) %>%
    select(Sample, Genome, Chr, Feature, Gene, Start, End, Strand, Length, inserts)
  
  # Length normalisation
  metat_featurecounts_summary = metat_featurecounts_summary_raw %>%
    left_join(formatted_metadata, by = c("Sample" = "internal_sample_name")) %>%
    mutate(norm_cov = inserts / (Length/1000) / (as.numeric(n_inserts)/10**6)) # Get RPKM values
  
  # Marker gene normalise
  featurecounts_pseudo_count = 10**(floor(log10(metat_featurecounts_summary %>% filter(norm_cov > 0) %>% pull(norm_cov) %>% min())))/10
  metat_featurecounts_mgs_median = metat_featurecounts_summary %>%
    filter(Gene %in% MGs) %>%
    group_by(Sample) %>%
    summarize(mgs_median = median(norm_cov[inserts > 0] + featurecounts_pseudo_count),
              n_mgs_detected = sum(inserts > 0))
  metat_featurecounts_summary_processed = metat_featurecounts_summary %>%
    left_join(metat_featurecounts_mgs_median) %>%
    filter(n_mgs_detected >= 6) %>%
    mutate(value = log2((norm_cov + featurecounts_pseudo_count) / (mgs_median))) %>% # Use an additive log ratio to the marker genes / Log 2 is more standard
    group_by(Sample) %>%
    mutate(n_genes_detected = sum(inserts > 0)) %>%
    ungroup()
  
  # Summarize for BGCs:
  metat_featurecounts_summary_processed_biosynthetic = metat_featurecounts_summary_processed %>%
    filter(Gene %in% (annotations_antismash %>% filter(biosynthetic_gene | biosynthetic_add_gene) %>% pull(gene))) %>%
    left_join(annotations_antismash, by = c("Gene" = "gene")) %>%
    group_by(biosynthetic_region, material) %>%
    summarize(n_biosynthetic_genes = n(),
              n_inserts = sum(inserts),
              median_inserts = median(inserts),
              value = median(value),
              products = unique(biosynthetic_products))
  
  return(metat_featurecounts_summary_processed_biosynthetic)
}


# BGC Analysis ===========================================================================

metat_metadata = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")

genome_id = "PIAMPJPB"
genome_id = "LGBFILLL"
genome_id = "OCMKBGHM"
genome_id = "OLPPLKCL"
#genome_id = "HLLJDLBE" # Nothing, as expected because bathy species

metat_feature_counts_summary_processed_bgcs = get_bgc_expression(genome_id = genome_id)
length(unique(metat_feature_counts_summary_processed_bgcs$material))
summary(metat_feature_counts_summary_processed_bgcs$value)

# Plot:
metat_feature_counts_summary_processed_bgcs %>%
  mutate(shape = ifelse(median_inserts > 0, "Evidence for expression", "No evidence")) %>%# View()
  mutate(color = ifelse(abs(value) >= 2, "#f28e2c", "#4e79a7")) %>%
  mutate(x = paste0(products, " (", biosynthetic_region, ")")) %>%
  ggplot() +
  geom_boxplot(aes(x = x, y = value), outlier.shape = NA, size = 0.2) +
  geom_jitter(aes(x = x, y = value, shape = shape, color = color), width = 0.2, height = 0, alpha = .7, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(-14, 8) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_identity() +
  ylab("Median expression rate per sample (log2)\n(0 Indicates values similar to housekeeping genes)") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    axis.ticks = element_line(size = unit(.3, 'mm')),
    legend.title = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    plot.margin = margin(0,0,0,0),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.title.y = element_blank(),
    line = element_line(size = unit(0.3, 'mm')),
    #aspect.ratio = 1/2,
  ) +
  coord_flip()


ggsave(filename = paste0("data/submission/editorial-nature/Figures/Material/", genome_id, "-expr.raw.pdf"),
       width = 110, height = 3*length(unique(metat_feature_counts_summary_processed_bgcs$biosynthetic_region)), units = "mm")

ggsave(filename = paste0("data/processed/figures/Figure-S5-eremios_bgcs_expr/Figure-S5.", genome_id, ".raw.pdf"),
       width = 170, height = 4*length(unique(metat_feature_counts_summary_processed_bgcs$biosynthetic_region)), units = "mm")
       