# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 1C - Mapping rates ==============================================================

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('/Users/paolil/polybox/PhD/Exploratorium/sushilab-colors/palettes-paoli.R')

# Load data ------------------------------------------------------------------------------

mapping_go = read_tsv("data/raw/go_microbiomics/summaries/go_microbiomics-mapping_rates.tsv", col_names = c("dataset", "readfile", "total_reads", "mapped_reads", "perc"))
mapping_gorg = read_tsv("data/raw/go_microbiomics/summaries/gorg_sags-mapping_rates.tsv", col_names = c("dataset", "readfile", "total_reads", "mapped_reads", "perc"))
mapping_gem = read_tsv("data/raw/go_microbiomics/summaries/gem_mags-mapping_rates.tsv", col_names = c("dataset", "readfile", "total_reads", "mapped_reads", "perc"))
samples_metadata = load_general_metadata()

fraction_categories = c("<-0.22"    = "<0.2",
                        "0.1-0.22"  = "<0.2",
                        "0.22-0.45" = "0.2-0.8",
                        "0.45-0.8"  = "0.2-0.8",
                        "0.2-0.8"   = "0.2-0.8",
                        "0.22-1.6"  = "0.2-3",
                        "0.22-3"    = "0.2-3",
                        "0.8-20"    = "0.8-20",
                        "0.2<-"     = ">0.2"
)

formatted_table = rbind(
  mapping_go %>%
    filter(!grepl(".s$", readfile)) %>% # remove single reads
    filter(!grepl("_T_", readfile)) %>% # remove metat
    mutate(`Internal Sample Name` = gsub("_.RR.*", "", readfile)) %>%
    mutate(`Internal Sample Name` = gsub("_METAG.*", "", `Internal Sample Name`)) %>%
    group_by(`Internal Sample Name`) %>%
    summarize(mapping_rates = mean(perc)) %>%
    left_join(samples_metadata) %>%
    mutate(Database = "GO Micro"),
  mapping_gorg %>%
    filter(dataset != "TARA_LSF") %>% # remove LSF samples
    mutate(`Internal Sample Name` = gsub("_.RR.*", "", readfile)) %>%
    mutate(`Internal Sample Name` = gsub("_METAG.*", "", `Internal Sample Name`)) %>%
    group_by(`Internal Sample Name`) %>%
    summarize(mapping_rates = mean(perc)) %>%
    left_join(samples_metadata) %>%
    mutate(Database = "GORG (filtered)"),
  mapping_gem %>%
    filter(dataset != "TARA_LSF") %>% # remove LSF samples
    mutate(`Internal Sample Name` = gsub("_.RR.*", "", readfile)) %>%
    mutate(`Internal Sample Name` = gsub("_METAG.*", "", `Internal Sample Name`)) %>%
    group_by(`Internal Sample Name`) %>%
    summarize(mapping_rates = mean(perc)) %>%
    left_join(samples_metadata) %>%
    mutate(Database = "GEM (MAGs)")
) %>%
  mutate(depth_layer = factor(depth_layer, levels = c("EPI", "MES", "BAT", "ABY"))) %>%
  mutate(size_fraction_binned = factor(fraction_categories[size_fraction], levels = unique(fraction_categories))) %>%
  mutate(latitude_binned = ifelse(abs(latitude) <= 30, "<30º", "30-60º")) %>%
  mutate(latitude_binned = ifelse(abs(latitude) > 60, ">60º", latitude_binned)) %>%
  mutate(latitude_binned = factor(latitude_binned, levels = c("<30º", "30-60º", ">60º"))) %>%
  mutate(Database = factor(Database, levels = c("GO Micro", "GEM (MAGs)", "GORG (filtered)")))

  
p1 = formatted_table %>%
  ggplot() +
  geom_boxplot(aes(x = size_fraction_binned, y = mapping_rates, fill = Database, color = Database), alpha = .6, outlier.colour = NA, size = .3) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0,20,40,60,80,100), minor_breaks = c(10,30,50,70,90), expand = c(.01,.01)) +
  facet_grid(.~size_fraction_binned, scales = "free_x") +
  scale_color_manual(values = database_colors_bluegreen) +
  scale_fill_manual(values = database_colors_bluegreen) +
  ylab("% Reads mapped") +
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'mm')), 
        text = element_text(size = 6),
        axis.title.y = element_text(size = 6, margin = margin(0,1,0,0, 'mm')),
        axis.text.y = element_text(size = 6, hjust = 1, margin = margin(0,.5,0,0, 'mm')),
        legend.position = 'none',
        rect = element_blank(),
        strip.background = element_rect(color = NA, fill = "black"),
        strip.text = element_text(color = "white", face = "bold", size = 6, angle = 90, margin = margin(.5, 1, .5, 1, "mm")),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(0, 'mm'))
p1

p2 = formatted_table %>%
  filter(size_fraction_binned %in% c("0.2-0.8", "0.2-3")) %>%
  ggplot() +
  geom_boxplot(aes(x = latitude_binned, y = mapping_rates, fill = Database, color = Database), alpha = .6, outlier.colour = NA, size = .3) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0,20,40,60,80,100), minor_breaks = c(10,30,50,70,90), expand = c(.01,.01)) +
  facet_grid(.~latitude_binned, scales = "free_x") +
  scale_color_manual(values = database_colors_bluegreen) +
  scale_fill_manual(values = database_colors_bluegreen) + 
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'mm')), 
        text = element_text(size = 6),
        legend.position = 'none',
        rect = element_blank(),
        strip.background = element_rect(color = NA, fill = "black"),
        strip.text = element_text(color = "white", face = "bold", size = 6, angle = 90, margin = margin(.5, 1, .5, 1, "mm")),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(0, 'mm'))
p2

p3 = formatted_table %>%
  filter(size_fraction_binned %in% c("0.2-0.8", "0.2-3")) %>%
  ggplot() +
  geom_boxplot(aes(x = depth_layer, y = mapping_rates, fill = Database, color = Database), alpha = .6, outlier.colour = NA, size = .3) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0,20,40,60,80,100), minor_breaks = c(10,30,50,70,90), expand = c(.01,.01)) +
  facet_grid(.~depth_layer, scales = "free_x") +
  scale_color_manual(values = database_colors_bluegreen) +
  scale_fill_manual(values = database_colors_bluegreen) + 
  theme_bw() +
  theme(line = element_line(size = unit(.3, 'mm')), 
        text = element_text(size = 6),
        rect = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(color = NA, fill = "black"),
        strip.text = element_text(color = "white", face = "bold", size = 6, angle = 90, margin = margin(.5, 1, .5, 1, "mm")),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(0, 'mm'))
p3

p = p1 + p2 + p3 + plot_layout(nrow = 1, widths = c(5/11, 3/11, 3/11)) + plot_annotation(theme = theme(plot.margin = margin()))
p

ggsave(paste0(figures_path_proj, "Figure-1/Figure-1D.test.pdf"), p, width = 69, height = 30, units = 'mm')


formatted_table %>%
  filter(Database == "GO Micro") %>%
  filter(size_fraction %in% c("0.8-20", "0.2<-")) %>%
  pull(mapping_rates) %>%
  median()
formatted_table %>%
  filter(Database == "GORG (filtered)") %>%
  pull(mapping_rates) %>%
  median()
formatted_table %>%
  filter(Database == "GEM (MAGs)") %>%
  pull(mapping_rates) %>%
  median()
40.8/19
40.8/12.6

formatted_table %>%
  filter(Database == "GO Micro") %>%
  filter(size_fraction %in% c("0.22-0.45", "0.45-0.8", "0.2-0.8", "0.22-1.6", "0.22-3")) %>%
  pull(mapping_rates) %>%
  median()
formatted_table %>%
  filter(Database == "GORG (filtered)") %>%
  filter(size_fraction %in% c("0.22-0.45", "0.45-0.8", "0.2-0.8", "0.22-1.6", "0.22-3")) %>%
  pull(mapping_rates) %>%
  median()
formatted_table %>%
  filter(Database == "GEM (MAGs)") %>%
  filter(size_fraction %in% c("0.22-0.45", "0.45-0.8", "0.2-0.8", "0.22-1.6", "0.22-3")) %>%
  pull(mapping_rates) %>%
  median()
60/22
60/21

formatted_table %>%
  filter(Database == "GO Micro") %>%
  filter(size_fraction %in% c("0.8-20", "0.2<-")) %>%
  pull(mapping_rates) %>%
  median()
formatted_table %>%
  filter(Database == "GORG (filtered)") %>%
  filter(size_fraction %in% c("0.8-20", "0.2<-")) %>%
  pull(mapping_rates) %>%
  median()
formatted_table %>%
  filter(Database == "GEM (MAGs)") %>%
  filter(size_fraction %in% c("0.8-20", "0.2<-")) %>%
  pull(mapping_rates) %>%
  median()
41/26
41/12
