# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 2CD - Compare MAG quality ===========================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')

# Functions ------------------------------------------------------------------------------

compute_drep_score <- function(table, use_gtdb = F){
  # From https://drep.readthedocs.io/en/latest/module_descriptions.html#choose
  # 𝐴(𝑐𝑜𝑚𝑝𝑙𝑒𝑡𝑒𝑛𝑒𝑠𝑠)–𝐵(𝑐𝑜𝑛𝑡𝑎𝑚𝑖𝑛𝑎𝑡𝑖𝑜𝑛)+𝐶(𝐶𝑜𝑛𝑡𝑎𝑚𝑖𝑛𝑎𝑡𝑖𝑜𝑛∗(𝑠𝑡𝑟𝑎𝑖𝑛ℎ𝑒𝑡𝑒𝑟𝑜𝑔𝑒𝑛𝑒𝑖𝑡𝑦/100))+𝐷(𝑙𝑜𝑔(𝑁50))+𝐸(𝑙𝑜𝑔(𝑠𝑖𝑧𝑒))
  # defaults: A=1, B=5, C=1, D=0.5, E=0
  table = table %>%
    mutate(
      drep_checking_score = `Mean Completeness` -
        5*`Mean Contamination` +
        (`Mean Contamination`*(`CheckM Strain heterogeneity`/100)) +
        0.5*(log10(`CheckM N50 (scaffolds)`)),
      drep_checkm_score = `CheckM Completeness` -
        5*`CheckM Contamination` +
        (`CheckM Contamination`*(`CheckM Strain heterogeneity`/100)) +
        0.5*(log10(`CheckM N50 (scaffolds)`))
    )
  # some small difference...
  table %>%
    filter(round(table$drep_score) != round(table$drep_checking_score)) %>%
    select(drep_score, drep_checking_score) %>%
    View("dRep score discrepancies")
  if (use_gtdb == T) {
    table = table %>%
      mutate(
        gtdb_drep_score = checkm_completeness -
          5*checkm_contamination +
          (checkm_contamination*(checkm_strain_heterogeneity/100)) +
          0.5*(log10(n50_scaffolds)),
        gtdb_quality = checkm_completeness - 5*checkm_contamination
      )
  }
  return(table)
}

compare_scores <- function(table, pattern1 = "TARA_.*METAG_[A-Z]{8}$", pattern2 = "TARA_.*_[A-Z]{3}[0-9]"){
  clusters_info = table %>%
    group_by(drep_cluster) %>%
    summarize(cat1 = any(grepl(pattern1, `Bin Id`)),
              cat2 = any(grepl(pattern2, `Bin Id`)))
  print(paste("Categories share", clusters_info %>% filter(cat1 & cat2) %>% nrow(), " clusters."))
  print(paste("Category 1 has", clusters_info %>% filter(cat1 & !cat2) %>% nrow(), "unique clusters."))
  print(paste("Category 1 has", clusters_info %>% filter(!cat1 & cat2) %>% nrow(), "unique clusters."))
  select_clusters = clusters_info %>%
    filter(cat1 & cat2) %>%
    pull(drep_cluster)
  filtered = table %>%
    filter(drep_cluster %in% select_clusters) %>%
    filter(grepl(paste(pattern1, pattern2, sep = "|"), `Bin Id`))
  processed = filtered %>%
    group_by(drep_cluster) %>%
    summarize(
      cat1_rep = grepl(pattern1, `Bin Id`)[`drep_rep_genome`],
      cat2_rep = grepl(pattern2, `Bin Id`)[`drep_rep_genome`],
      cat1_best_score = max(`drep_score`[grepl(pattern1,  `Bin Id`)]),
      cat2_best_score = max(`drep_score`[grepl(pattern2,  `Bin Id`)])    
    ) %>%
    mutate(cat1_is_better = cat1_best_score > cat2_best_score,
           cat1_cat2_diff = cat1_best_score - cat2_best_score,
           cat1_cat2_perc = cat1_cat2_diff/abs(cat2_best_score)*100)
}

# Load data ------------------------------------------------------------------------------

summary_table = load_raw_summary() %>%
  left_join(load_general_metadata()) %>%
  mutate(drep_rep_genome = as.logical(drep_rep_genome))

gtdb_metadata = rbind(
  read_tsv("data/raw/gtdb-metadata/bac120_metadata_r89.tsv"),
  read_tsv("data/raw/gtdb-metadata/ar122_metadata_r89.tsv")
)

# Format data ----------------------------------------------------------------------------

mags_with_species_annot = summary_table %>%
  filter(grepl("_METAG_", `Bin Id`)) %>%
  filter(!grepl(";s__$|N/A", classification))

gtdb_relevant_metadata = gtdb_metadata %>%
  filter(gtdb_taxonomy %in% (mags_with_species_annot %>% pull(classification) %>% unique()))

summary_table_mags_w_gtdb = mags_with_species_annot %>%
  full_join(gtdb_relevant_metadata, by = c("classification" = "gtdb_taxonomy"))

# Comparison with Delmont ----------------------------------------------------------------

Figure_2B_1_table = summary_table %>%
  # We want to compare with the TARA samples, prok, EPI, non-arctic
  filter(dataset %in% c("Tara Oceans", NA)) %>% 
  filter(depth_layer == "EPI" | is.na(depth_layer)) %>%
  filter(size_fraction %in% c("0.22-1.6", "0.22-3") | is.na(size_fraction)) %>%
  filter(as.numeric(gsub("TARA_", "", station)) <= 152 | is.na(station)) %>%
  compare_scores(.) %>% 
  mutate(comparison = "Rec. Genomes vs Delmont") 

Figure_2B_1_table %>% summary()
Figure_2B_1_table %>% pull(cat1_cat2_perc) %>% t.test()
Figure_2B_1_table %>% pull(cat1_cat2_perc) %>% wilcox.test()
sum(Figure_2B_1_table$cat1_cat2_perc>0)/nrow(Figure_2B_1_table)

Figure_2B_1 = Figure_2B_1_table %>%
  mutate(cat1_cat2_perc = ifelse(cat1_cat2_perc > 100, 100, cat1_cat2_perc)) %>%
  mutate(facet = "Vs Delmont") %>%
  ggplot() +
  geom_violin(aes(x = "dummy", y = cat1_cat2_perc)) +
  geom_boxplot(aes(x = "dummy", y = cat1_cat2_perc), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(position = "right", limits = c(-100, 100)) +
  coord_flip() +
  ylab("Quality score differences (%)") +
  facet_grid(rows = vars(facet), switch = "y") +
  theme_bw() +
  theme_cell +
  theme(rect = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "mm"),
        axis.ticks.length = unit(0.5, "mm"),
        axis.title.x =  element_text(size = 6, colour = "black"),
        axis.text.x =  element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background.y = element_rect(colour = NA, fill = "black"),
        strip.text.y = element_text(size = 6, colour = "white", face = "bold"))

Figure_2B_1

ggsave("data/processed/figures/Figure-S8-benchmark/Figure-S8-1.pdf", Figure_2B_1, width = 113, height = 40, units = "mm")


# Comparison with Tully ------------------------------------------------------------------

# Selecting the Tully genomes in GTDB
gtdb_metadata %>%
  filter(ncbi_submitter == "Tara Oceans Consortium") %>%
  pull(ncbi_bioproject) %>%
  table()

# These are TMED and the second Tully paper

tully_table = gtdb_metadata %>% filter(ncbi_submitter == "Tara Oceans Consortium")
nrow(tully_table)
tully_species = gtdb_metadata %>% filter(ncbi_submitter == "Tara Oceans Consortium") %>% 
  pull(gtdb_taxonomy) %>% unique()
length(tully_species)
tully_table %>% pull(gtdb_representative) %>% sum()

tully_comparison_table = summary_table %>%
  # We want to compare with the TARA samples, only non-arctic
  filter(grepl("_METAG_", `Bin Id`)) %>%
  filter(dataset == "Tara Oceans") %>% 
  filter(as.numeric(gsub("TARA_", "", station)) <= 152 | is.na(size_fraction)) %>%
  filter(!grepl(";s__$|N/A", classification)) %>%
  inner_join(tully_table, by = c("classification" = "gtdb_taxonomy")) %>%
  compute_drep_score(., use_gtdb = T)

tully_comparison_table %>% pull(accession) %>% unique() %>% length()
tully_comparison_table %>% pull(classification) %>% unique() %>% length()

Figure_2B_2_table = tully_comparison_table %>%
  group_by(classification) %>%
  summarize(mag_score = max(drep_checkm_score),
            mag_quality = max(Quality),
            gtdb_score = max(gtdb_drep_score),
            gtdb_quality = max(gtdb_quality)) %>%
  mutate(score_diff = mag_score - gtdb_score,
         score_perc = score_diff/abs(gtdb_score)*100,
         qual_diff = mag_quality - gtdb_quality,
         mag_score_is_better = mag_score > gtdb_score) %>%
  select(score_perc, qual_diff, mag_score_is_better)

Figure_2B_2_table %>% summary()
sum(Figure_2B_2_table$mag_score_is_better) / nrow(Figure_2B_2_table)
t.test(Figure_2B_2_table$score_perc)
wilcox.test(Figure_2B_2_table$score_perc)

Figure_2B_2_table %>%
  mutate(score_perc = ifelse(score_perc > 100, 100, score_perc)) %>%
  ggplot() +
  geom_violin(aes(x = "dummy", y = score_perc)) +
  geom_boxplot(aes(x = "dummy", y = score_perc), width = 0.5) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Quality score differences") +
  ylim(-100, 100) +
  theme_bw() +
  theme_cell +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

sum(!tully_species %in% (tully_comparison_table %>% pull(classification) %>% unique()))

gtdb_metadata %>%
  filter(ncbi_submitter == "Tara Oceans Consortium") %>%
  mutate(matched = ifelse(gtdb_taxonomy %in% tully_species[!tully_species %in% (tully_comparison_table %>% pull(classification) %>% unique())], "missing", "matched")) %>%
  ggplot() +
  geom_density(aes(x = checkm_completeness - 5*checkm_contamination + (checkm_contamination*(checkm_strain_heterogeneity/100)) + 0.5*(log10(n50_scaffolds)), color = matched))

# Parks comparison =======================================================================

# Selecting the Parks genomes in GTDB
gtdb_metadata %>%
  filter(ncbi_submitter == "University of Queensland") %>%
  pull(ncbi_bioproject) %>%
  table()

parks_table = gtdb_metadata %>% filter(ncbi_bioproject == "PRJNA348753")
nrow(parks_table)
parks_species = gtdb_metadata %>% filter(ncbi_bioproject == "PRJNA348753") %>% pull(gtdb_taxonomy) %>% unique()
length(parks_species)
parks_table %>% pull(gtdb_representative) %>% sum()

# estimate number of genomes from tara, download suppl from the paper
tara_projects_all = c("ERP001736", "ERP001737", "ERP003628", "ERP003708", "ERP006155", "ERP006156", "ERP010880", "ERP010881", "ERP010882")
tara_projects_ofinterest = c("ERP001736", "ERP001737", "ERP003708")
ts2 = readxl::read_xls("~/Downloads/41564_2017_12_MOESM2_ESM.xls")
ts3 = readxl::read_xls("~/Downloads/41564_2017_12_MOESM3_ESM.xls")
ts2 %>% filter(`Study Accession` %in% tara_projects_ofinterest) %>% pull(`Study Title`) %>% table
ts2_srx = ts2 %>% filter(`Study Accession` %in% tara_projects_all) %>% pull(`SRA Experiment Accession`)
ts3 %>% filter(`SRA Bin ID` %in% ts2_srx)

parks_comparison_table = summary_table %>%
  # We want to compare with the TARA samples, only non-arctic
  filter(grepl("_METAG_", `Bin Id`)) %>%
  filter(dataset == "Tara Oceans") %>% 
  filter(as.numeric(gsub("TARA_", "", station)) <= 152 | is.na(size_fraction)) %>%
  filter(!grepl(";s__$|N/A", classification)) %>%
  inner_join(parks_table, by = c("classification" = "gtdb_taxonomy")) %>%
  compute_drep_score(., use_gtdb = T)

parks_comparison_table %>% pull(accession) %>% unique() %>% length()
parks_comparison_table %>% pull(classification) %>% unique() %>% length()

Figure_2B_3_table = parks_comparison_table %>%
  group_by(classification) %>%
  summarize(mag_score = max(drep_checkm_score),
            mag_quality = max(Quality),
            gtdb_score = max(gtdb_drep_score),
            gtdb_quality = max(gtdb_quality)) %>%
  mutate(score_diff = mag_score - gtdb_score,
         score_perc = score_diff/abs(gtdb_score)*100,
         qual_diff = mag_quality - gtdb_quality,
         mag_score_is_better = mag_score > gtdb_score) %>%
  select(score_perc, qual_diff, mag_score_is_better)

Figure_2B_3_table %>% summary()
sum(Figure_2B_3_table$mag_score_is_better) / nrow(Figure_2B_3_table)
t.test(Figure_2B_3_table$score_perc)
wilcox.test(Figure_2B_3_table$score_perc)

Figure_2B_3_table %>%
  mutate(score_perc = ifelse(score_perc > 100, 100, score_perc)) %>%
  ggplot() +
  geom_violin(aes(x = "dummy", y = score_perc)) +
  geom_boxplot(aes(x = "dummy", y = score_perc), width = 0.5) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Quality score differences") +
  ylim(-80, 80) +
  theme_bw() +
  theme_cell +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Combined Parks and Tully ===============================================================




Figure_2B_2 = rbind(
  Figure_2B_2_table %>% mutate(facet = "Vs Tully"),
  Figure_2B_3_table %>% mutate(facet = "Vs Parks")
) %>%
  mutate(score_perc = ifelse(score_perc > 100, 100, score_perc)) %>%
  ggplot() +
  geom_violin(aes(x = "dummy", y = score_perc)) +
  geom_boxplot(aes(x = "dummy", y = score_perc), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(position = "right", limits = c(-100, 100)) +
  coord_flip() +
  facet_grid(rows = vars(facet), switch = "y") +
  theme_bw() +
  theme_cell +
  theme(rect = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "mm"),
        axis.ticks.length = unit(0.5, "mm"),
        axis.title.x =  element_text(size = 6, colour = "black"),
        axis.text.x =  element_text(size = 6, colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background.y = element_rect(colour = NA, fill = "black"),
        strip.text.y = element_text(size = 6, colour = "white", face = "bold"))

Figure_2B_2

ggsave("data/processed/figures/Figure-S8-benchmark/Figure-S8-2.pdf", Figure_2B_2, width = 113, height = 70, units = "mm")

