# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 2A - Abundance correlation binning ==================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('/Users/paolil/polybox/PhD/Exploratorium/sushilab-colors/palettes-paoli.R')

# Load data ------------------------------------------------------------------------------

general_metadata = load_general_metadata()

samples_list = rbind(read_tsv("data/raw/go_microbiomics/samples/GEOTRACES.rand10-5-5.samples", col_names = "Internal Sample Name") %>% mutate(backmapping = "n=610", description = "Biogeotraces-HOTS-BATS"),
                     read_tsv("data/raw/go_microbiomics/samples/TARA_OCEANS_prok.rand10.samples", col_names = "Internal Sample Name") %>% mutate(backmapping = "n=180", description = "Tara Oceans\n(Prok)"),
                     read_tsv("data/raw/go_microbiomics/samples/TARA_OCEANS_vir.rand10.samples", col_names = "Internal Sample Name") %>% mutate(backmapping = "n=190", description = "Tara Oceans\n(Viral)"),
                     read_tsv("data/raw/go_microbiomics/samples/MALASPINA.rand5.samples", col_names = "Internal Sample Name") %>% mutate(backmapping = "n=58", description = "Malaspina")) %>%
  mutate(`Internal Sample Name` = gsub("_METAG", "", `Internal Sample Name`)) %>%
  left_join(general_metadata %>% select(`Internal Sample Name`, size_fraction))

summary_table = load_raw_summary()

without_diffcov = read_tsv("data/raw/go_microbiomics/summaries/diffcov_go_micro_metabat2.c2k.e500_cpl50_ctn10-evaluate_summary.tsv") %>%
  mutate(Sample = gsub("_[^_]*$", "", `Bin Id`),
         `Internal Sample Name` = gsub('_METAG', '', basename(dirname(`Folder Location`))),
         Method = gsub('_a[0-9]+\\.', '', basename(`Folder Location`))) %>%
  add_genomes_quality(weight = 5)

# Prepare data ------------------------------------------------------------------------------

sum(samples_list$`Internal Sample Name` %in% summary_table$`Internal Sample Name`)
samples_list$`Internal Sample Name`[!(samples_list$`Internal Sample Name` %in% summary_table$`Internal Sample Name`)]

col_subset = c("Internal Sample Name", "Folder Location", "Bin Id", "Mean Completeness", "Mean Contamination", "CheckM Completeness", "CheckM Contamination", "CheckM Strain heterogeneity", "CheckM N50 (scaffolds)", "Anvio Domain", "Anvio Domain Confidence", "Anvio Completion", "Anvio Redundancy", "Anvio # Scaffolds", "Anvio Length", "Folder Location", "Method", "Quality", "HighQ", "GoodQ", "MediumQ", "LowQ")

figure_2B_table = 
  rbind(
    left_join(samples_list, without_diffcov %>% select(all_of(col_subset))) %>% mutate(`Differential Coverage` = FALSE),
    left_join(samples_list, summary_table %>% select(all_of(col_subset))) %>% mutate(`Differential Coverage` = TRUE)
  ) %>%
  group_by(`Internal Sample Name`, `Differential Coverage`) %>%
  summarize(backmapping = unique(backmapping),
            description = unique(description),
            `# MAGs` = n(),
            `# HighQ` = sum(HighQ),
            `# GoodQ` = sum(GoodQ),
            `# MediumQ` = sum(MediumQ),
            `# LowQ` = sum(LowQ),
            `Cumulative Q-score` = sum(Quality),
            `Cumulative Q'-score` = sum(`Mean Completeness` - 5*`Mean Contamination` + (`Mean Contamination`*(`CheckM Strain heterogeneity`/100)) + 0.5*(log10(`CheckM N50 (scaffolds)`))),
            `Cumulative Q'-checkm-score` = sum(`CheckM Completeness` - 5*`CheckM Contamination` + (`CheckM Contamination`*(`CheckM Strain heterogeneity`/100)) + 0.5*(log10(`CheckM N50 (scaffolds)`))),
            `Average Q-score` = mean(Quality),
            `Average Q'-score` = mean(`Mean Completeness` - 5*`Mean Contamination` + (`Mean Contamination`*(`CheckM Strain heterogeneity`/100)) + 0.5*(log10(`CheckM N50 (scaffolds)`))),
            `Average Q'-checkm-score` = mean(`CheckM Completeness` - 5*`CheckM Contamination` + (`CheckM Contamination`*(`CheckM Strain heterogeneity`/100)) + 0.5*(log10(`CheckM N50 (scaffolds)`)))
  ) %>%
  replace(is.na(.), 0) %>%
  left_join(load_general_metadata()) %>%
  group_by(Sample) %>%
  summarize(dataset = unique(dataset),
            backmapping = unique(backmapping),
            description = unique(description),
            size_fraction = unique(size_fraction),
            `# MAGs ratio` = `# MAGs`[`Differential Coverage`]/`# MAGs`[!`Differential Coverage`],
            `Cumulative Q-score ratio` = `Cumulative Q-score`[`Differential Coverage`]/`Cumulative Q-score`[!`Differential Coverage`],
            `Cumulative Q'-score ratio` = `Cumulative Q'-score`[`Differential Coverage`]/`Cumulative Q'-score`[!`Differential Coverage`],
            `Cumulative Q'-checkm-score ratio` = `Cumulative Q'-checkm-score`[`Differential Coverage`]/`Cumulative Q'-checkm-score`[!`Differential Coverage`],
            `Average Q-score ratio` = `Average Q-score`[`Differential Coverage`]/`Average Q-score`[!`Differential Coverage`],
            `Average Q'-score ratio` = `Average Q'-score`[`Differential Coverage`]/`Average Q'-score`[!`Differential Coverage`],
            `Average Q'-checkm-score ratio` = `Average Q'-checkm-score`[`Differential Coverage`]/`Average Q'-checkm-score`[!`Differential Coverage`]) %>%
  replace(is.na(.), 1)

# Plot figure ----------------------------------------------------------------------------

fraction_dict = c("<-0.22" = " (Viral)",
                  "0.1-0.22" = " (Viral)",
                  "0.22-0.45" = " (F.L.1)",
                  "0.45-0.8" = " (F.L.1)",
                  "0.22-1.6" = " (F.L.2)",
                  "0.22-3" = " (F.L.2)",
                  "0.2-0.8" = " (F.L.1)",
                  "0.8-20" = " (P.A.)")

plot_order = c("n=190, Tara Oceans (Viral)", "n=190, Tara Oceans (F.L.1)", "n=180, Tara Oceans (F.L.2)", "n=58, Malaspina (F.L.1)", "n=58, Malaspina (P.A.)", "n=610, Biogeotraces", "n=610, Hawaiian Ocean Time-Series", "n=610, Bermuda-Atlantic Time-Series")

figure_2B_table %>%
  #filter(`Cumulative Q'-score ratio` != Inf) %>%
  summary()
figure_2B_table %>%
  filter(`Cumulative Q'-score ratio` != Inf) %>%
  summary()

figure_2B = figure_2B_table %>%
  mutate(fraction = ifelse(grepl("Tara|Mala", dataset), fraction_dict[size_fraction], "")) %>%
  mutate(facet = factor(paste0(backmapping, ", ", dataset, fraction), levels = plot_order)) %>%
  mutate(dot_type = ifelse(`Cumulative Q'-score ratio` %in% c(Inf, 1), "Special Case", "Normal")) %>%
  ggplot() +
  geom_boxplot(aes(x = facet, y = `Cumulative Q'-score ratio`, color = dataset), size = 0.3, outlier.shape = NA, fill = NA, position = position_dodge(width = 1)) +
  geom_point(aes(x = facet, y = `Cumulative Q'-score ratio`, color = dataset, fill = dataset), size = .5, shape = 21, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0)) +
  #geom_point(aes(x = facet, y = mean(`# MAGs ratio`)), shape = 3, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.2) +
  scale_color_manual(values = dataset_colors[unique(figure_2B_table$dataset)]) +
  scale_fill_manual(values = alpha(dataset_colors[unique(figure_2B_table$dataset)], alpha = 0.6)) +
  facet_grid(.~facet, scales = "free_x", space = "free_x") +
  coord_cartesian(clip = "off") +
  scale_y_log10() +
  theme_bw() +
  theme_cell +
  theme(plot.margin = margin(0,0,0,0, 'mm'),
        legend.position = "none",
        rect = element_rect(size = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        panel.spacing.x = unit(0, 'mm'),
        strip.background = element_rect(color = "lightgrey", fill = "lightgrey"),
        strip.text = element_text(size = 6, colour = "white", face = "bold", margin = margin(0.5, 1, 0.5, 1, "mm")))

figure_2B

ggsave(paste0(figures_path_proj, "Figure-2/Figure-2B.pdf"), figure_2B, width = 48, height = 25, units = "mm")




