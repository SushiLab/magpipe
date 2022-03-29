# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 2C - Data prep for trees ========================================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Prepare data ===========================================================================

genomes_summary = load_prettified_summary() %>%
  mutate(`Representative Genome` = as.logical(`Representative Genome`)) %>%
  left_join(read_tsv("data/raw/go_microbiomics/motus/v1.0/gom_mOTUsv2.mag-memberships.tsv", col_names = c("Genome Id", "motus_cluster"))) %>%
  mutate(motus_cluster = gsub("ext_mOTU_v26_", "gom_", motus_cluster))

samples_metadata = load_general_metadata() %>%
  mutate(`Internal Sample Name` = ifelse(grepl("TARA_", `Internal Sample Name`), paste0(`Internal Sample Name`, "_METAG"), `Internal Sample Name`))

motus_abundances = read_tsv("data/raw/go_microbiomics/motus/v1.0/gom.motus", skip = 2) %>%
  select(all_of(c("#consensus_taxonomy", samples_metadata$`Internal Sample Name`))) %>%
  gather(key = sample, value = motus_count, -`#consensus_taxonomy`) %>%
  mutate(taxonomy = gsub(" \\[.*", "", `#consensus_taxonomy`),
         motus_cluster = gsub(".*\\[|\\]", "", `#consensus_taxonomy`))

mardb_metadata = read_tsv("../MARINOMICS_EAN/scratch/processed/metadata/mardb/mardb_current.tsv.gz", guess_max = 10000)

go_micro_antismash_file = "data/raw/go_microbiomics/summaries/go_microbiomics-integrated-cpl50_ctn10-antismash_summary.tsv.gz"
eremio_antismash_file = "data/raw/marine_eremios/eremiobacterota_superproducer-integrated-cpl50_ctn10-antismash.filtered.tsv"

scaffolds_length = read_tsv("data/raw/go_microbiomics/summaries/genomes_bgcs-scaffolds_length.tsv")

# Investigate genomes not detected in the dataset ========================================

# Among the genomes with a motus cluster, which ones are not detected in the dataset?

genomes_not_detected = motus_abundances %>%
  filter(motus_cluster %in% unique(genomes_summary$motus_cluster)) %>%
  group_by(motus_cluster) %>%
  summarize(motu_cum_abd = sum(motus_count)) %>%
  filter(motu_cum_abd == 0) %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, motus_cluster))
assertthat::are_equal(nrow(genomes_not_detected), nrow(genomes_summary %>% filter(motus_cluster %in% genomes_not_detected$motus_cluster)))

gsub("[^_]*_[^_]*_|_[^_]*$", "", genomes_not_detected$genome) %>% table
gsub("_.*", "", genomes_not_detected$genome) %>% table
980/(980+9)

genomes_not_detected %>%
  group_by(motus_cluster) %>%
  summarize(n=n(),
            n_ref = sum(grepl("_REFG_", genome)),
            n_sags = sum(grepl("_SAGS_", genome))) %>%
  mutate(r_ref = n_ref/n,
         r_sags = n_sags/n) %>%
  mutate(type = ifelse(r_ref == 1, "REFs", "SAGs")) %>% #View() # Everything is homogeneous
  pull(type) %>% table()
637/(637+7)

# What about the genomes that don't make it to motus?

genomes_summary %>%
  filter(is.na(motus_cluster)) %>%
  pull(`Genome Id`) %>%
  gsub("_[^_]*$", "", .) %>%
  gsub(".*_", "", .) %>%
  table()

genomes_summary %>%
  filter(is.na(motus_cluster)) %>%
  pull(`Genome Id`) %>%
  gsub("_.*", "", .) %>%
  table()

# Inspect isolation source ===============================================================

source_not_detected = genomes_not_detected %>%
  mutate(mmp_ID = gsub(".*REFG_", "", genome)) %>%
  left_join(mardb_metadata %>% select(mmp_ID, isolation_source))

source_not_detected %>%
  pull(isolation_source) %>%
  table() %>%
  sort() %>% 
  View()

# Unknown sources ------------------------------------------------------------------------

source_not_detected_missing = source_not_detected %>%
  filter(grepl("[Mm]issing|[Nn]ot [Kk]nown", isolation_source) | is.na(isolation_source))
source_not_detected_remaining = source_not_detected %>% filter(! genome %in% source_not_detected_missing$genome)

# Sediment-like sources ------------------------------------------------------------------

source_not_detected_sediment = source_not_detected %>%
  filter(grepl("[Ss]ediment|[Mm]ud|[Ss]and|[Ss]oil", isolation_source))
source_not_detected_remaining = source_not_detected_remaining %>% filter(! genome %in% source_not_detected_sediment$genome)

# Hydrothermal region sources ------------------------------------------------------------

source_not_detected_smokers = source_not_detected %>%
  filter(grepl("[Hh]ydrothermal|[Ss]moker|[Bb]asaltic", isolation_source))
source_not_detected_remaining = source_not_detected_remaining %>% filter(! genome %in% source_not_detected_smokers$genome)

# Saltern/brine sources ------------------------------------------------------------------

source_not_detected_salterns = source_not_detected %>%
  filter(grepl("[Ss]alt|[Bb]rine", isolation_source))
source_not_detected_remaining = source_not_detected_remaining %>% filter(! genome %in% source_not_detected_salterns$genome)

# Coastal sources ------------------------------------------------------------------------

source_not_detected_coastal = source_not_detected %>%
  filter(grepl("[Cc]oast|[Tt]idal|[Pp]ool|[Ss]hore|[Ee]stuary|[Bb]rackish|[Mm]angrove", isolation_source))
source_not_detected_remaining = source_not_detected_remaining %>% filter(! genome %in% source_not_detected_coastal$genome)

# Host-associated sources ----------------------------------------------------------------

source_not_detected_host = source_not_detected %>%
  filter(grepl(paste(c("[Ff]ish,[Ss]quirt,[Cc]oral,[Ss]ponge,[Aa]lga,[Ll]obster,[Hh]omarus,[Oo]yster,[Cc]lam,[Ss]eagull,[Ss]almon,
                       [Bb]ass,[Tt]rout,[Ss]eaweed,[Ss]eagrass,[Gg]ill,[Gg]ut,[Ii]ntestine,[Aa]nemone,[Ff]ucus,[Hh]epatopancreas",
                       "[Aa]quaculture", "[Aa]quarium", "[Pp]enguin", "[Mm]at", "[Ee]lephant", "[Cc]ulture", "[Ff]eces"),
                     collapse = "|"), isolation_source))
source_not_detected_remaining = source_not_detected_remaining %>% filter(! genome %in% source_not_detected_host$genome)

# Iterate here until only ocean left -----------------------------------------------------

source_not_detected_remaining %>%
  pull(isolation_source) %>%
  table() %>%
  sort() %>% 
  View()

# Actual open ocean ----------------------------------------------------------------------

source_not_detected_water = source_not_detected_host = source_not_detected %>%
  filter(grepl(paste(c("[Ss]ea.*[Ww]ater", "[Dd]eep", "[Ww]ater"),
                     collapse = "|"), isolation_source))
source_not_detected_water %>%
  pull(isolation_source) %>%
  table() %>%
  sort() %>% 
  View()

# Conclusion =============================================================================

# On that basis we want to:
# Remove MARD genomes that did not get a motus cluster
# And remove MARD genomes not detected in motus
# Because isolation source may be a bit messy. (See bottom for actual evaluation of that)

genomes_summary %>%
  filter(grepl("MARD_", `Genome Id`)) %>%
  group_by(`dRep Dereplication Cluster`) %>%
  mutate(missing = all(is.na(motus_cluster))) %>%
  filter(missing) %>% View()

blacklisted_genomes = c(
  genomes_not_detected %>% filter(grepl("^MARD_", genome)) %>% pull(genome),
  genomes_summary %>% filter(is.na(motus_cluster) & grepl("^MARD_", `Genome Id`)) %>% pull(`Genome Id`)
)
assertthat::are_equal(length(blacklisted_genomes), 988 + 180)

write_lines(blacklisted_genomes, "data/raw/go_microbiomics/summaries/go_microbiomics-blacklisted_genomes.list")

genomes_summary %>%
  filter(`Genome Id` %in% blacklisted_genomes) %>%
  pull(`dRep Dereplication Cluster`) %>%
  table() %>% 
  sort()

# Setup summary for Figure 3 =============================================================

max_bgcs <- function(x, y){
  if (!all(is.na(y))){
    m = max(y, na.rm = T)
    first_max = which(y == m)[1]
    res = x[first_max]
  } else {
    res = NA
  }
  return(res)
}

antismash_table = read_tsv(go_micro_antismash_file) %>%
  separate_rows(`protoclusters products`, sep = ';') %>%
  mutate(`Antismash Product` = get_bcg_category(`protoclusters products`)) %>%
  left_join(scaffolds_length %>% dplyr::rename(scaffold_length = length)) %>%
  filter(scaffold_length >= 5000) %>%
  filter(!genome %in% blacklisted_genomes) %>%
  select(-scaffold_length) %>%
  select(-contains("CDS"))

# Then we want to replace Eremio with the representative values:
eremio_antismash = read_tsv(eremio_antismash_file) %>%
  separate_rows(`protoclusters products`, sep = ';') %>%
  mutate(`Antismash Product` = get_bcg_category(`protoclusters products`)) %>%
  filter(genome == "MALA_SAMN05422137_METAG_HLLJDLBE") %>% # select representative
  mutate(genome = "MALA_SAMN05422189_METAG_HFLHJDGN") %>% # mock name
  select(-contains("CDS"))

antismash_table_fixed = antismash_table %>%
  filter(!genome %in% eremio_antismash$genome) %>%
  rbind(eremio_antismash)

# create summary
antismash_summary = antismash_table_fixed %>%
  group_by(genome) %>%
  summarize(`# Biosynthetic Regions` = length(unique(region)),
            `# Biosynthetic Products` = n())

# Add numbers per category
for (product in unique(antismash_table_fixed$`Antismash Product`)){
  antismash_summary[, product] =
    antismash_table_fixed %>%
    mutate(category_bool = `Antismash Product` == product) %>%
    group_by(genome) %>%
    summarise(n = sum(category_bool)) %>%
    pull(n)
}
# Add all the information you may want
for (product in unique(antismash_table_fixed$`protoclusters products`)){
  antismash_summary[, paste0('product:', product)] =
    antismash_table_fixed %>%
    mutate(category_bool = `protoclusters products` == product) %>%
    group_by(genome) %>%
    summarise(n = sum(category_bool)) %>%
    pull(n)
}

Species_biosynth = genomes_summary %>%
  filter(!`Genome Id` %in% blacklisted_genomes) %>%
  select(`Genome Id`, `dRep Dereplication Cluster`, `Representative Genome`, `GTDB Taxonomy`) %>%
  left_join(antismash_summary %>% dplyr::rename(`Genome Id` = genome)) %>%
  group_by(`dRep Dereplication Cluster`) %>%
  summarize(`Representative Genome` = `Genome Id`[`Representative Genome`],
            `RiPPs (Ribosomal Natural Products)` = max_bgcs(RiPPs, `# Biosynthetic Regions`),
            `Non-Ribosomal Peptide Synthases` = max_bgcs(NRPS, `# Biosynthetic Regions`),
            `Type I Polyketide Synthases` = max_bgcs(PKSI, `# Biosynthetic Regions`),
            `Type II & III Polyketide Synthases` = max_bgcs(`Other PKS`, `# Biosynthetic Regions`),
            Terpenes = max_bgcs(Terpenes, `# Biosynthetic Regions`),
            Other = max_bgcs(Other, `# Biosynthetic Regions`) + max_bgcs(Saccharides, `# Biosynthetic Regions`),
            `# Biosynthetic Products` = max_bgcs(`# Biosynthetic Products`, `# Biosynthetic Regions`),
            `# Biosynthetic Regions` = max_bgcs(`# Biosynthetic Regions`, `# Biosynthetic Regions`), # MAKE SURE PRODUCTS ARE CHANGED LAST
            Taxonomy = paste(`GTDB Taxonomy`, collapse = "|"),
            Genomes = paste(`Genome Id`, collapse = "|"),
            `# Genomes` = n(),
            `# Archaea` = sum(grepl("d__Archaea", `GTDB Taxonomy`)),
            `# Bacteria` = sum(grepl("d__Bacteria", `GTDB Taxonomy`)),
            `# MAGS` = sum(grepl("_METAG_", `Genome Id`)),
            `# SAGS` = sum(grepl("_SAGS_", `Genome Id`)),
            `# REFG` = sum(grepl("_REFG_", `Genome Id`))) %>%
  arrange(desc(`# Biosynthetic Products`))

Species_biosynth %>% filter(`# Biosynthetic Regions` > 15) %>% View()


genomes_summary %>%
  #filter(grepl("_METAG_", `Genome Id`)) %>%
  #filter(grepl("_SAGS_", `Genome Id`)) %>%
  filter(grepl("_REFG_", `Genome Id`)) %>%
  filter(`GTDB Taxonomy` != "N/A") %>%
  pull(`GTDB Taxonomy`) %>%
  gsub(";c__.*", "", .) %>%
  table() %>%
  length()

genomes_summary %>%
  filter(! `Genome Id` %in% blacklisted_genomes) %>%
  #filter(grepl("_METAG_", `Genome Id`)) %>%
  #filter(grepl("_SAGS_", `Genome Id`)) %>%
  filter(grepl("_REFG_", `Genome Id`)) %>%
  filter(`GTDB Taxonomy` != "N/A") %>%
  pull(`GTDB Taxonomy`) %>%
  gsub(";c__.*", "", .) %>%
  table() %>%
  length()

mags_phyla = genomes_summary %>% filter(! `Genome Id` %in% blacklisted_genomes) %>% filter(`GTDB Taxonomy` != "N/A") %>% filter(grepl("_METAG_", `Genome Id`)) %>%   pull(`GTDB Taxonomy`) %>% gsub(";c__.*", "", .) %>% unique
sags_phyla = genomes_summary %>% filter(! `Genome Id` %in% blacklisted_genomes) %>% filter(`GTDB Taxonomy` != "N/A") %>% filter(grepl("_SAGS_", `Genome Id`)) %>%   pull(`GTDB Taxonomy`) %>% gsub(";c__.*", "", .) %>% unique
refs_phyla = genomes_summary %>% filter(! `Genome Id` %in% blacklisted_genomes) %>% filter(`GTDB Taxonomy` != "N/A") %>% filter(grepl("_REFG_", `Genome Id`)) %>%   pull(`GTDB Taxonomy`) %>% gsub(";c__.*", "", .) %>% unique
sum(!sags_phyla %in% mags_phyla)
sum(!refs_phyla %in% mags_phyla)

read_tsv(go_micro_antismash_file) %>%
  left_join(scaffolds_length %>% dplyr::rename(scaffold_length = length)) %>%
  filter(scaffold_length >= 5000) %>%
  filter(!genome %in% blacklisted_genomes) %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, `dRep Dereplication Cluster`, `Representative Genome`, `GTDB Taxonomy`)) %>%
  filter(grepl(";s__$", `GTDB Taxonomy`)) %>%
  nrow()
15659/39055

# Additional BGC info -------

bgc_clustering = read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_metadata.tsv") %>%
  filter(dataset == "go_micro_genomes") %>%
  filter(!genome %in% blacklisted_genomes) %>%
  left_join(scaffolds_length %>% dplyr::rename(scaffold_length = length)) %>%
  filter(scaffold_length >= 5000) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_to_gcf_gcc.tsv")) %>%
  left_join(read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_dists_to_refseq_mibig.tsv")) %>%
  left_join(genomes_summary %>% select(genome = `Genome Id`, species = `dRep Dereplication Cluster`))

# produce new selector
bgc_selector = read_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_metadata.tsv") %>%
  select(bgc_id) %>%
  mutate(selected = bgc_id %in% bgc_clustering$bgc_id)
nrow(bgc_selector)
sum(bgc_selector$selected)
bgc_selector %>% write_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/all_bgc_selector_lucas.tsv")
  
# Numbers of GCFs, GCCs
bgc_clustering %>% pull(gcf) %>% unique %>% length
bgc_clustering %>% pull(gcc) %>% unique %>% length

# New GCFs
gcf_dist = bgc_clustering %>%
  # group_by(gcf, species) %>% # Uncomment these 3 lines for dereplication
  # filter(bgc_id == bgc_id[length == max(length)][1]) %>%
  # ungroup() %>%
  group_by(gcf) %>%
  summarize(n = n(),
            mean_d = mean(refseq),
            median_d = median(refseq),
            mean_d_mibig = mean(mibig),
            median_d_mibig = median(mibig))
gcf_dist %>% filter(mean_d >= 0.2) %>% nrow
gcf_dist %>% filter(median_d >= 0.2) %>% nrow
gcf_dist %>% filter(mean_d_mibig >= 0.2) %>% nrow / 6907
gcf_dist %>% filter(median_d_mibig >= 0.2) %>% nrow / 6907
# size distrib
gcf_dist %>% pull(n) %>% table %>% sort 
gcf_dist %>% filter(median_d >= 0.2) %>% pull(n) %>% table %>% sort 

gcf_dist %>% write_tsv("data/raw/go_microbiomics/bgcs/tables-v0.1/curated_gcfs_dist_nofilter.tsv")

# New GCCs
gcc_dist = bgc_clustering %>%
  # group_by(gcf, species) %>% # Uncomment these 3 lines for dereplication
  # filter(bgc_id == bgc_id[length == max(length)][1]) %>%
  # ungroup() %>%
  group_by(gcc) %>%
  summarize(n = n(),
            mean_d = mean(refseq),
            median_d = median(refseq),
            no_mag = all(!grepl("_METAG_", genome)),
            any_mag = any(grepl("_METAG_", genome)),
            only_mag = all(grepl("_METAG_", genome)),
            only_new_mag = all(grepl("_METAG_[A-Z]{8}$", genome)))
gcc_dist %>% filter(mean_d >= 0.4) %>% nrow
gcc_dist %>% filter(median_d >= 0.4) %>% nrow

# GCC genome composition
nrow(gcc_dist)
sum(gcc_dist$no_mag)
sum(gcc_dist$any_mag)
sum(gcc_dist$any_new_mag)
sum(gcc_dist$only_mag)
sum(gcc_dist$only_new_mag)
gcc_dist %>% filter(median_d >= 0.4) %>% nrow
gcc_dist %>% filter(median_d >= 0.4) %>% pull(any_mag) %>% sum
gcc_dist %>% filter(median_d >= 0.4) %>% pull(any_new_mag) %>% sum
gcc_dist %>% filter(median_d >= 0.4) %>% pull(only_mag) %>% sum
gcc_dist %>% filter(median_d >= 0.4) %>% pull(only_new_mag) %>% sum

