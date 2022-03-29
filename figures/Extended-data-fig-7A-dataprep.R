# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 7A - Data prep ======================================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(stringi)
library(googlesheets4)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))
}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('/Users/paolil/polybox/PhD/Exploratorium/sushilab-colors/palettes-paoli.R')

# Data preparation =======================================================================

# Metadata -------------------------------------------------------------------------------

all_samples = googlesheets4::read_sheet("1JDgQoi8bY5N8L_HokD2R2zWay3ACtbwj4YMbyQzdyMw")
all_metat_samples = all_samples %>% filter(grepl("METAT$", dataset))

seqdepth = read_tsv("data/processed/marine_eremios/distributions/TOPC-METAG.stats", col_names = F) %>%
  rbind(read_tsv("data/processed/marine_eremios/distributions/TOPC-METAT.stats", col_names = F)) %>%
  rbind(read_tsv("data/processed/marine_eremios/distributions/TOPC-lsf.stats", col_names = F)) %>%
  rbind(read_tsv("data/processed/marine_eremios/distributions/GEOTRACES.stats", col_names = F)) %>%
  rbind(read_tsv("data/processed/marine_eremios/distributions/ALOHA.stats", col_names = F)) %>%
  rbind(read_tsv("data/processed/marine_eremios/distributions/MALASPINA.stats", col_names = F)) %>%
  mutate(sample = gsub("\\.[m|pe|s].*", "", X1),
         type = gsub(".readstats.*", "", X1)) %>%
  mutate(type = gsub(".*\\.", "", type)) %>%
  mutate(n_reads = as.numeric(gsub(".*:", "", X1))) %>%
  select(-X1) %>%
  spread(key = type, value = n_reads) %>%
  mutate(n_inserts = m + pe/2 + s) %>%
  select(sample, n_inserts)

pangaea_stream = tara_metadata_from_pangea(type = "Env_meso") %>%
  select(barcode = `Sample ID (TARA_barcode#, registered at ...)`,
         ocean_province = `OS region ([abbreviation] full name (MRG...)`,
         station = `Station (TARA_station#, registered at ...)`,
         material = `Sample material (TARA_station#_environmental-f...)`,
         depth = `Depth, nominal (from which this sample was co...)`) %>%
  mutate(size_fraction = gsub(".*_", "", material))

metat_dict = all_metat_samples %>%
  filter(grepl("_T$", representative_barcode)) %>%
  select(representative_barcode) %>%
  mutate(material = gsub("_T", "", representative_barcode)) %>%
  mutate(material = gsub("TARA_180_ZZZ_0.22-3", "TARA_180_DCM_0.22-3", material)) %>% 
  mutate(material = gsub("TARA_168_IZZ_0.22-3", "TARA_168_DCM_0.22-3", material)) %>% 
  mutate(material = gsub("TARA_189_ZZZ_0.22-3", "TARA_189_DCM_0.22-3", material)) %>%
  mutate(material = gsub("MXL", "MIX", material)) %>%
  mutate(material = gsub("IZZ", "ZZZ", material)) %>%
  mutate(material = gsub("TARA_100_MES", "TARA_100_DCM", material)) %>% # Deep DCM sample
  left_join(select(pangaea_stream, material, barcode)) %>% 
  select(representative_barcode, barcode) %>%
  filter(!duplicated(representative_barcode))

formatted_metadata = all_metat_samples %>% 
  left_join(seqdepth, by = c("representative_barcode" = "sample")) %>%
  select(representative_barcode, dataset, n_inserts) %>%
  left_join(metat_dict) %>%
  mutate(barcode = ifelse(grepl("TARA_.*_META.$", representative_barcode), gsub("_META.$", "", representative_barcode), barcode)) %>%
  mutate(barcode = gsub("SUB", "", barcode)) %>%
  left_join(pangaea_stream) %>%
  mutate(ocean_province = gsub(".*\\] | \\(.*", "", ocean_province),
         station = ifelse(grepl("TARA|Cruise", station), station, paste0("MALA_", station))) %>%
  rename(sample = representative_barcode)

write_tsv(formatted_metadata, "data/raw/marine_eremios/marine_eremios-metat_analysis-metadata.tsv")
googlesheets4::write_sheet(formatted_metadata, "1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")

# Fix sample naming:
tmp_metat_metadata = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")
tmp_metat_metadata = tmp_metat_metadata %>% rename(internal_sample_name = sample)
metat_dict = rbind(
  googlesheets4::read_sheet("1UVr124uMLvBY2IlYevgJXDzC2n8nDYh_xklgDRtETA4", sheet = "Tara prokaryote enriched metatranscriptomes") %>% select(sample_name, internal_sample_name),
  googlesheets4::read_sheet("1UVr124uMLvBY2IlYevgJXDzC2n8nDYh_xklgDRtETA4", sheet = "Tara eukarote enriched metatranscriptomes") %>% select(sample_name, internal_sample_name)) %>%
  select(sample = sample_name, internal_sample_name)
tmp_metat_metadata = metat_dict %>%
  left_join(tmp_metat_metadata)
googlesheets4::write_sheet(tmp_metat_metadata, "1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")

# Simplify featurecounts data ------------------------------------------------------------

metat_featurecounts_summary_raw = read_tsv("data/raw/marine_eremios/marine_eremios-featurecounts-summary_PIAMPJPB.tsv", skip = 1) %>%
  select(!contains("ETHSEQ")) %>% 
  select(!contains("METAG")) %>%
  gather(key = Sample, value = inserts, -c(Geneid, Chr, Start, End, Strand, Length)) %>%
  mutate(Feature = Geneid,
         Gene = gsub("^[^-]*-", "", Geneid),
         Genome = gsub("-.*", "", Geneid),
         Sample = gsub(".filter.*|.*/", "", Sample)) %>%
  select(Sample, Genome, Chr, Feature, Gene, Start, End, Strand, Length, inserts)

# Save raw data
write_tsv(metat_featurecounts_summary_raw, "data/raw/marine_eremios/marine_eremios-featurecounts-summary_PIAMPJPB_metat.tsv.gz")

# Process data
MGs = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-fetchMGs")
formatted_metadata = googlesheets4::read_sheet("1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "metat-metadata")

metat_featurecounts_summary = metat_featurecounts_summary_raw %>%
  left_join(formatted_metadata, by = c("Sample" = "sample")) %>%
  mutate(norm_cov = inserts / (Length/1000) / (as.numeric(n_inserts)/10**6)) # Get RPKM values

featurecounts_pseudo_count = 10**(floor(log10(metat_featurecounts_summary %>% filter(norm_cov > 0) %>% pull(norm_cov) %>% min())))/10

metat_featurecounts_mgs_median = metat_featurecounts_summary %>%
  filter(Gene %in% (MGs %>% filter(`mOTUs-10`) %>% pull(bestMGs))) %>%
  group_by(Sample) %>%
  summarize(mgs_median = median(norm_cov[inserts > 0] + featurecounts_pseudo_count),
            n_mgs_detected = sum(inserts > 0))

metat_featurecounts_mgs_median %>% filter(n_mgs_detected >= 6) %>% nrow()

metat_featurecounts_summary_processed = metat_featurecounts_summary %>%
  left_join(metat_featurecounts_mgs_median) %>%
  filter(n_mgs_detected >= 6) %>%
  mutate(value = log2((norm_cov + featurecounts_pseudo_count) / (mgs_median))) %>% # Use an additive log ratio to the marker genes / Log 2 is more standard
  group_by(Sample) %>%
  mutate(n_genes_detected = sum(inserts > 0)) %>%
  ungroup()

# We want it to be well centered in 0
metat_featurecounts_summary_processed %>%
  ggplot() +
  geom_density(aes(x = value))

metat_featurecounts_summary_processed %>%
  filter(Gene %in% (MGs %>% filter(`mOTUs-10`) %>% pull(bestMGs))) %>%
  ggplot() +
  geom_density(aes(x = value))

metat_featurecounts_summary_processed %>%
  filter(Gene %in% (MGs %>% filter(`mOTUs-10`) %>% pull(bestMGs))) %>%
  filter(inserts > 0) %>%
  ggplot() +
  geom_density(aes(x = value))

write_tsv(metat_featurecounts_summary_processed, "data/raw/marine_eremios/marine_eremios-featurecounts-summary_PIAMPJPB_metat-processed.tsv.gz")

# Get antismash annotations in nice format -----------------------------------------------

annotations_antismash_raw = read_tsv("data/raw/marine_eremios/eremiobacterota_superproducer-integrated-cpl50_ctn10-antismash.filtered.tsv") %>%
  filter(grepl("PIAM", genome)) %>%
  select(-genome)

cds_biosynthetic = annotations_antismash_raw %>% select(`CDS biosynthetic`) %>% separate_rows(`CDS biosynthetic`, sep = ";")
cds_biosynthetic_add = annotations_antismash_raw %>% select(`CDS biosynthetic additional`) %>% separate_rows(`CDS biosynthetic additional`, sep = ";")
cds_transport = annotations_antismash_raw %>% select(`CDS transport`) %>% separate_rows(`CDS transport`, sep = ";")
cds_regulatory = annotations_antismash_raw %>% select(`CDS regulatory`) %>% separate_rows(`CDS regulatory`, sep = ";")
cds_resistance = annotations_antismash_raw %>% select(`CDS resistance`) %>% separate_rows(`CDS resistance`, sep = ";")

annotations_antismash = annotations_antismash_raw %>%
  select(scaffold, biosynthetic_region = region, biosynthetic_products = products, gene = `CDS list`) %>%
  separate_rows(gene, sep = ";") %>%
  select(gene, scaffold, biosynthetic_region, biosynthetic_products) %>%
  mutate(biosynthetic_gene = gene %in% cds_biosynthetic$`CDS biosynthetic`,
         biosynthetic_add_gene = gene %in% cds_biosynthetic_add$`CDS biosynthetic additional`,
         transport_gene = gene %in% cds_transport$`CDS transport`,
         regulatory_gene = gene %in% cds_regulatory$`CDS regulatory`,
         resistance_gene = gene %in% cds_resistance$`CDS resistance`)

write_tsv(annotations_antismash, "data/raw/marine_eremios/marine_eremios-PIAMPJPB-antismash-processed.tsv")
googlesheets4::write_sheet(annotations_antismash, "1aDBI-gAB1k1ob7uSDXnOvg4mxHXgRK0NdFASD18m6_0", sheet = "annotations-antismash")
