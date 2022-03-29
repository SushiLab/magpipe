# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 3D - Ca. Eudoremicrobium distribution ===========================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Variables ------------------------------------------------------------------------------

ordered_fractions = c("<-0.22", "0.1-0.22", "0.22-0.45", "0.2-0.8", "0.45-0.8", "0.22-1.6", "0.22-3", "0.2<-", "0.8-3", "0.8-5", "0.8-20", "0.8->", "3->", "5-20", "20-180", "180-2000")

# Load abundance data --------------------------------------------------------------------

eudore_clusters = read_tsv("data/raw/go_microbiomics/motus/v1.0/gom_mOTUsv2.mag-memberships.tsv", col_names = c("Genome Id", "motu_cluster")) %>%
  mutate(motu_cluster = gsub("ext_mOTU_v26_", "gom_", motu_cluster)) %>%
  filter(grepl("PIAMPJPB|HLLJDLBE|LGBFILLL|OLPPLKCL|OCMKBGHM", `Genome Id`))
eudore_clusters_dict = gsub(".*_", "", eudore_clusters$`Genome Id`)
names(eudore_clusters_dict) = eudore_clusters$motu_cluster

motu_profile = read_tsv("data/raw/go_microbiomics/motus/v1.0/gom.motus", skip = 2) %>%
  gather(key = sample, value = motu_abd, -`#consensus_taxonomy`) %>%
  mutate(taxonomy = gsub(" \\[.*", "", `#consensus_taxonomy`),
         motu_cluster = gsub(".*\\[|\\]", "", `#consensus_taxonomy`))

motu_profile_rel = motu_profile %>%
  group_by(sample) %>%
  mutate(rel_abd = motu_abd/sum(motu_abd)) %>%
  ungroup()

motu_profile_rel_eudore = motu_profile_rel %>%
  filter(motu_cluster %in% eudore_clusters$motu_cluster)

# Load metadata --------------------------------------------------------------------------

pangaea_stream = tara_metadata_from_pangea(type = "Env_meso") %>%
  select(barcode = `Sample ID (TARA_barcode#, registered at ...)`,
         ocean_province = `OS region ([abbreviation] full name (MRG...)`,
         station = `Station (TARA_station#, registered at ...)`,
         material = `Sample material (TARA_station#_environmental-f...)`,
         depth = `Depth, nominal (from which this sample was co...)`) %>%
  left_join(tara_metadata_from_pangea(type = "Env_sensor") %>% # add this for E taraoceanii specific analysis
               select(barcode = `Sample ID (TARA_barcode#, registered at ...)`,
                      temperature = `Temp [°C] (median value (50th percentile...)`,
                      salinity = `Sal (median value (50th percentile...)`)) %>%
  mutate(size_fraction = gsub(".*_", "", material))


metadata_ethseq = load_general_metadata() %>% 
  select(barcode = `Internal Sample Name`, station, ocean_province, depth, depth_layer, size_fraction) %>%
  mutate(material = paste(station, depth_layer, size_fraction, sep = "_")) %>%
  select(-depth_layer) %>%
  filter(grepl("ETHSEQ", barcode))

metat_dict = motu_profile_rel_eudore %>%
  filter(grepl("_T$", sample)) %>%
  select(sample) %>%
  mutate(material = gsub("_T", "", sample)) %>%
  left_join(select(pangaea_stream, material, barcode)) %>% 
  select(sample, barcode) %>%
  filter(!duplicated(sample))

motu_profile_rel_eudore_w_metadata = motu_profile_rel_eudore %>% 
  mutate(Genome = eudore_clusters_dict[motu_cluster]) %>%
  left_join(metat_dict) %>%
  mutate(barcode = ifelse(grepl("TARA_.*_META.$", sample), gsub("_META.$", "", sample), barcode)) %>%
  mutate(barcode = ifelse(is.na(barcode), sample, barcode)) %>%
  mutate(barcode = gsub("SUB", "", barcode)) %>%
  left_join(rbind(pangaea_stream, metadata_ethseq)) %>%
  mutate(depth_layer = add_depth_layers(depth),
         ocean_province = gsub(".*\\] | \\(.*", "", ocean_province),
         station = ifelse(grepl("TARA|Cruise", station), station, paste0("MALA_", station))) %>%
  mutate(rel_abd = ifelse(rel_abd == 0, NA, rel_abd))

# quick test:
motu_profile_rel_eudore %>%
  filter(motu_cluster == "gom_001516") %>%
  filter(grepl("TARA", sample)) %>%
  left_join(metat_dict) %>%
  mutate(barcode = ifelse(grepl("TARA_.*_META.$", sample), gsub("_META.$|SUB", "", sample), barcode)) %>%
  left_join(pangaea_stream) %>% 
  filter(rel_abd > 0) %>% googlesheets4::write_sheet(ss = "1TGnPmlikvzhNGBa4fFU0ftcFlMQCGgFEcf3BvM2lXKA", sheet = "E. taraoceanii distrib w. metadata")
  #View()

# Marker genes Distribution Size fractions -----------------------------------------------

# Check where it's not
motu_profile_rel_eudore_w_metadata %>%
  filter(!grepl("T$", sample)) %>%
  mutate(size_fraction = factor(size_fraction, levels = ordered_fractions)) %>%
  filter(ocean_province == "Arctic Ocean" | grepl("Cruise", station)) %>%
  ggplot() +
  geom_tile(aes(x=station,y=size_fraction,fill=rel_abd)) +
  facet_grid(factor(gsub(".*_", "", Genome), levels=c("HLLJDLBE", "PIAMPJPB", "LGBFILLL", "OCMKBGHM", "OLPPLKCL"))*depth_layer~ocean_province, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(direction = -1, end = .9, na.value = "lightgrey", limits = c(0,0.08), name = "Relative Abundance") +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

motu_profile_rel_eudore_w_metadata %>%
  filter(!grepl("T$", sample)) %>%
  filter(ocean_province != "Arctic Ocean" & !grepl("Cruise", station)) %>%
  mutate(size_fraction = factor(size_fraction, levels = ordered_fractions)) %>%
  filter(depth_layer %in% c("EPI")) %>%
  ggplot() +
  geom_tile(aes(x=station,y=size_fraction,fill=rel_abd)) +
  facet_grid(factor(gsub(".*_", "", Genome), levels=c("HLLJDLBE", "PIAMPJPB", "LGBFILLL", "OCMKBGHM", "OLPPLKCL"))*depth_layer~ocean_province, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(direction = -1, end = .9, na.value = "lightgrey", limits = c(0,0.08), name = "Relative Abundance") +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# ggsave(plot_distrib,
#        filename = paste0(figures_path_proj, "Figure-5/Marine_eremio_distrib_EPI.pdf"),
#        width = 184, height = 60, units = col_unit)

motu_profile_rel_eudore_w_metadata %>%
  filter(!grepl("T$", sample)) %>%
  filter(ocean_province != "Arctic Ocean" & !grepl("Cruise", station)) %>%
  mutate(size_fraction = factor(size_fraction, levels = ordered_fractions)) %>%
  filter(depth_layer %in% c("MES")) %>%
  ggplot() +
  geom_tile(aes(x=station,y=size_fraction,fill=rel_abd)) +
  facet_grid(factor(gsub(".*_", "", Genome), levels=c("HLLJDLBE", "PIAMPJPB", "LGBFILLL", "OCMKBGHM", "OLPPLKCL"))*depth_layer~ocean_province, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(direction = -1, end = .9, na.value = "lightgrey", limits = c(0,0.87), name = "Relative Abundance") +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

motu_profile_rel_eudore_w_metadata %>%
  filter(!grepl("T$", sample)) %>%
  filter(ocean_province != "Arctic Ocean" & !grepl("Cruise", station)) %>%
  mutate(size_fraction = factor(size_fraction, levels = ordered_fractions)) %>%
  filter(depth_layer %in% c("BAT")) %>%
  ggplot() +
  geom_tile(aes(x=station,y=size_fraction,fill=rel_abd)) +
  facet_grid(factor(gsub(".*_", "", Genome), levels=c("HLLJDLBE", "PIAMPJPB", "LGBFILLL", "OCMKBGHM", "OLPPLKCL"))*depth_layer~ocean_province, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(direction = -1, end = .9, na.value = "lightgrey", limits = c(0,0.08), name = "Relative Abundance") +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Marker genes Distribution Overview -----------------------------------------------------

motu_profile_rel_eudore_w_metadata_overview = motu_profile_rel_eudore_w_metadata %>%
  filter(!grepl("T$", sample)) %>%
  group_by(Genome, ocean_province, station, depth_layer) %>%
  summarize(max_rel_abd = max(rel_abd, na.rm = T)) %>%
  ungroup() %>%
  mutate(max_rel_abd = ifelse(max_rel_abd == -Inf, NA, max_rel_abd))
         
ocean_province_dict = 
  c("North Pacific Ocean" = "NPO",
    "South Pacific Ocean" = "SPO",
    "North Atlantic Ocean" = "NAO",
    "South Atlantic Ocean" = "SAO",
    "Mediterranean Sea" = "MS",
    "Red Sea" = "RS",
    "Indian Ocean" = "IO",
    "Southern Ocean" = "SO")

genome_dict = c(
  "HLLJDLBE" = "C1",
  "PIAMPJPB" = "C2",
  "LGBFILLL" = "B1",
  "OLPPLKCL" = "A1",
  "OCMKBGHM" = "A2"
)

plot_distrib = motu_profile_rel_eudore_w_metadata_overview %>%
  filter(ocean_province != "Arctic Ocean" & !grepl("Cruise", station)) %>%
  filter(!grepl('Mediterranean', ocean_province)) %>%
  mutate(depth_layer = factor(depth_layer, levels = rev(c("EPI", "MES", "BAT")))) %>%
  mutate(ocean_province = factor(ocean_province_dict[ocean_province],
                                 levels = ocean_province_dict[c("North Pacific Ocean", "South Pacific Ocean", "North Atlantic Ocean", "South Atlantic Ocean", "Mediterranean Sea", "Red Sea", "Indian Ocean", "Southern Ocean")])) %>%
  mutate(Genome = factor(genome_dict[gsub(".*_", "", Genome)],
                         levels=genome_dict[c("OLPPLKCL", "OCMKBGHM", "LGBFILLL", "HLLJDLBE", "PIAMPJPB")])) %>% #View()
  ggplot() +
  geom_tile(aes(x=station,y=depth_layer,fill=max_rel_abd)) +
  facet_grid(Genome ~ ocean_province,
             scales = "free_x", space = "free_x", switch = "y") +
  scale_fill_viridis_c(direction = -1, end = .9, limits = c(-0.0001, 0.08), breaks = c(0, 0.02, 0.04, 0.06, 0.08), 
                      label = c(0,2,4,6,8), na.value = "lightgrey",
                       guide = guide_colorbar(title = "Relative Abundance (%)", #title.position = "top",
                                              title.vjust = .9,
                                              barwidth = 10, barheight = 0.5)) +
  scale_y_discrete(expand=c(0,0), position = 'right') +
  scale_x_discrete(expand=c(0,0)) +
  xlab("Sampling Stations") +
  theme_bw() +
  theme_cell +
  theme(text = element_text(size = 6),
        panel.border = element_rect(size = unit(.2, 'mm')),
        plot.margin = margin(0,0,0,0, 'mm'),
        panel.spacing = unit(0.3, "mm"),
        strip.background = element_rect(color = NA, fill = "black"),
        strip.text = element_text(size = 6, colour = "white", face = "bold", margin = margin(.5, 0, .5, 0, "mm")),
        #strip.text = element_text(size = 6, colour = "black", margin = margin(.5, 0, .5, 0, "mm")),
        strip.placement = "outside",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.title.y =  element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0, "mm"))
plot_distrib

ggsave(plot_distrib,
        filename = "data/processed/figures/Figure-S6-eremios_distributions/Figure-S6-eremios_distributions.raw.pdf",
        width = 183, height = 40, units = 'mm')


