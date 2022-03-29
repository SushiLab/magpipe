# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 2A - GCC overview prep ==========================================================

library(tidyverse)

# Based on BGCs

gtdb = read_tsv("data/raw/gtdb-metadata/ar122_metadata_r89.tsv") %>%
  rbind(read_tsv("data/raw/gtdb-metadata/bac120_metadata_r89.tsv"))

mibig = read_tsv("data/processed/go_microbiomics/mibig/mibig_genera.tsv", col_names = c("genus", "n"))

gtdb_tax = gtdb %>%
  select(gtdb_taxonomy, ncbi_taxonomy) %>%
  unique()

gtdb_genus = gtdb_tax %>%
  select(gtdb_taxonomy) %>%
  unique() %>%
  mutate(genus = gsub(".*;g__|;s__.*", "", gtdb_taxonomy),
         phylum = gsub(".*;p__|;c__.*", "", gtdb_taxonomy)) %>%
  select(-gtdb_taxonomy)


mibig %>%
  filter(!genus %in% gtdb_genus$genus) %>%
  View()

l = NULL
for (g in mibig %>% filter(!genus %in% gtdb_genus$genus) %>% pull(genus)){
  if (any(grepl(g, c(gtdb_tax$ncbi_taxonomy, gtdb_tax$gtdb_taxonomy)))){
    message(g)
    tmp = gtdb_tax %>%
      filter(grepl(g, ncbi_taxonomy))
    print(tmp)
    l = c(l, g)
  }
}

l

ncbi_genus = gtdb %>%
  select(ncbi_taxonomy) %>%
  unique() %>%
  mutate(genus = gsub(".*;g__|;s__.*", "", ncbi_taxonomy),
         phylum = gsub(".*;p__|;c__.*", "", ncbi_taxonomy)) %>%
  select(-ncbi_taxonomy) %>%
  rbind(tibble(genus = "Polyangium", phylum = "Myxococcota"))

expected_phyla = mibig %>%
  left_join(rbind(gtdb_genus, ncbi_genus) %>% unique) %>%
  filter(!is.na(phylum)) %>%
  filter(!duplicated(genus)) %>%
  group_by(phylum) %>%
  summarize(n = sum(n)) %>%
  mutate(r = round(n / sum(n)*100)) %>%
  arrange(desc(r))

bgc_based = expected_phyla %>%
  filter(r > 5)
bgc_based

# Based on abundance in the db:

highly_repr_phyla = genomes_summary %>%
  filter(!`Genome Id` %in% blacklisted_genomes) %>%
  filter(as.logical(`Representative Genome`)) %>%
  select(taxo = `GTDB Taxonomy`) %>%
  mutate(phylum = gsub(".*;p__|;c__.*", "", taxo)) %>%
  group_by(phylum) %>%
  summarize(n = n()) %>%
  mutate(r = round(n / sum(n)*100)) %>%
  arrange(desc(r))

repr_based = highly_repr_phyla %>%
  filter(r > 1) %>%
  filter(!phylum %in% bgc_based$phylum)
repr_based
