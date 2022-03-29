# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 1A - Depth distrib ==================================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(ggimage)
library(facetscales)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Prepare data ---------------------------------------------------------------------------

# Samples data
figure_1B_table =
  load_general_metadata() %>%
  select(-`Internal Sample Name`) %>%
  # Remove GORG for now
  # rbind(gorg_metadata_from_local() %>%
  #         mutate(dataset = "GORG SAGs",
  #                Sample = paste("GORG", Sample),
  #                station = paste(station, gsub("HOT[0-9]*|BATS[0-9]*", "", cruise)),
  #                depth_num = as.numeric(depth),
  #                depth_layer = add_depth_layers(as.character(depth)),
  #                `size_fraction` = "N/A",
  #                `temperature [°C]` = "N/A",
  #                `oxygen [µmol/kg]` = "N/A") %>%
  #        select(-cruise, -plate, -links)) %>%
mutate(
  size = ifelse(grepl("Time-Series", dataset), "Time-Series", "Survey"),
  image = c(
    "Biogeotraces" = "code/analysis/lib/icon_circle_biogeotraces.png",
    "Malaspina" = "code/analysis/lib/icon_circle_malaspina.png",
    "Tara Oceans" = "code/analysis/lib/icon_circle_tara.png",
    "Hawaiian Ocean Time-Series" = "code/analysis/lib/icon_clock_hots.png",
    "Bermuda-Atlantic Time-Series" = "code/analysis/lib/icon_clock_bats.png",
    "GORG SAGs" = "code/analysis/lib/icon_cross_gorg.png"
  )[dataset]) %>%
  mutate(depth_layer = factor(depth_layer, levels=c("EPI", "MES", "BAT", "ABY")))

# Plot samples through depths using transformed scale ------------------------------------

trans <- function(x) -sqrt(x)
inv <- function(x) x**2

figure_1B = ggplot(figure_1B_table) +
  geom_image(aes(x = longitude, y = depth_num, image = image, size = size), asp = 1.85) +
  scale_size_manual(values = c(0.01, 0.0166)) +
  scale_y_continuous(trans = scales::trans_new("revsqrt_trans", trans, inv, minor_breaks = scales::regular_minor_breaks(reverse = TRUE)),
                     limits = c(6000, 0),
                     breaks = c(0, 200, 1000, 2000, 4000, 6000),
                     expand=c(0,0)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
  ylab('Depth (m)') +
  xlab('Longitude') +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme_cell +
  theme(rect = element_blank(),
        text = element_text(size = unit(6, "pt")),
        plot.margin = margin(1, 1, 1, 1, "mm"),
        legend.position = "none",
        axis.text.y = element_text(size = unit(6, "pt")),
        axis.text.x = element_text(size = unit(6, "pt"), hjust = 1),
        axis.title.y = element_text(size = unit(6, "pt")),
        axis.title.x = element_text(size = unit(6, "pt")),
        axis.ticks.length = unit(0.5, "mm"))

figure_1B

# Plot samples through depths using facets ------------------------------------

scales_y <- list(
  `EPI` = scale_y_continuous(trans = "reverse", limits = c(200, 0), breaks = c(0, 100, 200), expand = c(0,0)),
  `MES` = scale_y_continuous(trans = "reverse", limits = c(1000, 200), breaks = c(500, 1000), minor_breaks = c(750), expand = c(0,0)),
  `BAT` = scale_y_continuous(trans = "reverse", limits = c(4000, 1000), breaks = c(2000, 4000), expand = c(0,0)),
  `ABY` = scale_y_continuous(trans = "reverse", limits = c(6000, 4000), breaks = c(5000, 6000), expand = c(0,0))
)

figure_1B = ggplot(figure_1B_table) +
  geom_image(aes(x = longitude, y = depth_num, image = image, size = size), asp = 12.3) +
  facet_grid_sc(depth_layer~., scales = list(y = scales_y)) +
  scale_size_manual(values = c(0.01, 0.0166)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
  ylab('Depth (m)') +
  xlab('Longitude (°)') +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme_cell +
  theme(rect = element_blank(),
        text = element_text(size = unit(6, "pt")),
        plot.margin = margin(1, 0, 1, 0, "mm"),
        legend.position = "none",
        axis.text.y = element_text(size = unit(6, "pt"), hjust = 1),
        axis.text.x = element_text(size = unit(6, "pt")),
        axis.title.y = element_text(size = unit(6, "pt")),
        axis.title.x = element_text(size = unit(6, "pt")),
        axis.ticks.length = unit(0.5, "mm"),
        strip.background = element_rect(colour = NA, fill = "lightgrey"),
        panel.spacing.y = unit(0, 'mm'),
        strip.text = element_text(size = 6, colour = "white", face = "bold", margin = margin(0.5, 0.5, 0.5, 0.5, "mm")))

figure_1B

# Save figure -------

ggsave(paste0(figures_path_proj, "Figure-1/Figure-1B.pdf"), figure_1B, width = 113, height = 41.5, units = col_unit)

