# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 1A - Map ========================================================================

# Libraries ------------------------------------------------------------------------------

rm(list = ls())
library(leaflet)
library(tidyverse)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R') # dataset colors

# Define some helpful functions ----------------------------------------------------------

format_links <-function(Sample, depth_layer, size_fraction, station, links = NULL){
  # Helpful debug statement:
  print(paste("processing a group of", length(Sample), "samples"))
  print(tibble(station, Sample, depth_layer, size_fraction))
  if (length(unique(station)) > 1){message(paste("Careful here, you could be pulling together things that are from different stations.",
                                                 paste(Sample, collapse = ", ")))}
  # Define helper functions
  row_tag <- function(row){return(paste0('<tr>', paste(row, collapse=''), '</tr>'))}
  header_tag <- function(header){return(paste0('<th>', header, '</th>'))}
  cell_tag <- function(cell){return(paste0('<td>', cell, '</td>'))}
  # Specify the order of the fraction possibilities and init the table
  fraction_levels = c("<-0.22", "0.1-0.22", "0.22-0.45", "0.2-0.8", "0.45-0.8", "0.22-1.6", "0.22-3", "0.8-20","0.2<-", "0.22<-")
  assertthat::assert_that(all(size_fraction %in% fraction_levels))
  ordered_size_fractions = sort(factor(unique(size_fraction), levels = fraction_levels, labels = fraction_levels))
  content = row_tag(header_tag(c('Fraction (&mu;m)', as.character(ordered_size_fractions))))
  # Build up the table row by row
  # Specify the order and loop through the depth possibilities
  depth_levels = c("EPI", "MES", "BAT", "ABY")
  assertthat::assert_that(all(depth_layer %in% depth_levels))
  for (i in sort(factor(unique(depth_layer), levels = depth_levels))){
    row = cell_tag(paste0("<b>", i, "</b>"))
    # Loop through the ordered size_fractions
    for (j in ordered_size_fractions){
      flag = sum(depth_layer == i & size_fraction == j)
      if (flag == 0) {
        row = paste0(row, cell_tag('None'))
      } else {
        target = Sample[which(depth_layer == i & size_fraction == j)]
        assertthat::are_equal(flag, length(target))
        if (length(target) == 1){
          link_tag = ifelse(is.null(links), paste0("<a href='../table/?filter=", target, "' target='_blank'>", gsub("^[^_]*_|_[^_]*$", "", target),"</a>"), links)
          row = paste0(row, cell_tag(link_tag))
        } else {
          link_tag = ifelse(is.null(links),
                            paste0("<a href='../table/?filter=", target, "' target='_blank'>", gsub("^[^_]*_|_[^_]*$", "", target), "</a>", collapse = ", "),
                            paste0(links, collapse = ", "))
          row = paste0(row, cell_tag(link_tag))
        }
      }
    }
    # debug
    print(row)
    # update the content
    content = paste0(content, row_tag(row))
  }
  # add some basic styling
  style = 
    "
table {
  display: block;
  overflow: auto;
}
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
}
tr:first-child td, tr:first-child th {
  border-top: none;
}
tr:last-child td {
  border-bottom: none;
}
tr td:first-child, tr th:first-child {
  border-left: none;
}
tr td:last-child, tr th:last-child {
  border-right: none;
}
th, td {
  text-align: left;
}
"
styled_table = paste0("<table>", content, "</table>", "<br><style>", style, "</style>")
return(styled_table)
}

# Prepare data ---------------------------------------------------------------------------

# Load metagenomes metadata
general_metadata = load_general_metadata() %>%
  select(-`Internal Sample Name`) %>%
  mutate(`temperature [°C]` = ifelse(is.na(`temperature [°C]`), "N/A", `temperature [°C]`),
         `oxygen [µmol/kg]` = ifelse(is.na(`oxygen [µmol/kg]`), "N/A", `oxygen [µmol/kg]`))

# Prepare spatial surveys
# SANITY CHECK: Are all stations  closely located ?
general_metadata %>%
  filter(!dataset %in% c('Hawaiian Ocean Time-Series', 'Bermuda-Atlantic Time-Series')) %>%
  mutate(position = paste(DescTools::RoundTo(latitude, 0.1), DescTools::RoundTo(longitude, 0.1), sep = ", "),
         rounded = paste(DescTools::RoundTo(latitude, 1), DescTools::RoundTo(longitude, 1), sep = ", ")) %>%
  group_by(station) %>%
  summarize(positions = paste(unique(position), collapse = " | "),
            rounded = paste(unique(rounded), collapse = " | "),
            var_lat = var(latitude),
            var_lng = var(longitude),
            n_positions = length(unique(position))) %>%
  arrange(desc(n_positions)) #%>% View()
# Except a few cases where the station are a bit spread out, this seems reasonable
spatial_surveys = general_metadata %>%
  filter(!dataset %in% c('Hawaiian Ocean Time-Series', 'Bermuda-Atlantic Time-Series')) %>%
  group_by(station) %>%
  # Prepare text popup
  mutate(text = paste0("<b>", Sample, "</b><br>",
                       "Station: ", station, "<br>", 
                       "Dataset: ", dataset, "<br><br>", 
                       "<u>Links to genomes from the station:</u>",
                       format_links(Sample, depth_layer, size_fraction, station))) %>%
  mutate(text = paste0(text,
                       "<u>Sample metadata:</u><br>",
                       "Latitude: ", latitude, "<br>",
                       "Longitude: ", longitude, "<br>",
                       "Size Fraction: ", size_fraction, " [µm]<br>",
                       "Depth Layer: ", depth_layer, "<br>",
                       "Depth: ", depth, " [m]<br>",
                       "Temperature: ", `temperature [°C]`, " [°C]<br>",
                       "Oxygen: ", `oxygen [µmol/kg]`, " [µmol/kg]"))  %>%
  ungroup()

# Prepare time series
# Some safeety checks
assertthat::assert_that(all(
  length(general_metadata %>% filter(dataset == 'Hawaiian Ocean Time-Series') %>% pull(latitude) %>% unique()) == 1,
  length(general_metadata %>% filter(dataset == 'Hawaiian Ocean Time-Series') %>% pull(longitude) %>% unique()) == 1,
  length(general_metadata %>% filter(dataset == 'Hawaiian Ocean Time-Series') %>% pull(depth_layer) %>% unique()) == 1
),
msg = "You can't really collapse all the data can you?")
assertthat::assert_that(all(
  length(general_metadata %>% filter(dataset == 'Bermuda-Atlantic Time-Series') %>% pull(latitude) %>% unique()) == 1,
  length(general_metadata %>% filter(dataset == 'Bermuda-Atlantic Time-Series') %>% pull(longitude) %>% unique()) == 1,
  length(general_metadata %>% filter(dataset == 'Bermuda-Atlantic Time-Series') %>% pull(depth_layer) %>% unique()) == 1
),
msg = "You can't really collapse all the data can you?")
time_series = general_metadata %>%
  filter(dataset %in% c('Hawaiian Ocean Time-Series', 'Bermuda-Atlantic Time-Series')) %>%
  # Remove duplicates
  filter(!duplicated(dataset)) %>%
  # Prepare text popup
  mutate(Sample = dataset,
         text = paste0("<b>", dataset, "</b><br>",
                       "<a href='../table/?filter=", gsub("_.*", "_", Sample), "' target='_blank'>Link to genomes from the time-series</a>", "<br><br>",
                       "Latitude: ", latitude, "<br>",
                       "Longitude: ", longitude, "<br>",
                       "Size Fraction: ", size_fraction, " [µm]<br>",
                       "Depth Layer: ", depth_layer))

# Finalise the data
sample_markers = spatial_surveys %>%
  rbind(time_series) %>%
  # Add groups for layered display
  mutate(groups = paste0(dataset, ' - ', depth_layer))

# Leaflet --------------------------------------------------------------------------------

# Define custom icons
custom_icons <- iconList(
  "Biogeotraces" = makeIcon("/Users/paolil/polybox/PhD/Projects/GENOME_RESOLVED_OCEAN_EAN/code/analysis/lib/icon_circle_biogeotraces.png", iconWidth = 12, iconHeight = 12),
  "Malaspina" = makeIcon("/Users/paolil/polybox/PhD/Projects/GENOME_RESOLVED_OCEAN_EAN/code/analysis/lib/icon_circle_malaspina.png", iconWidth = 12, iconHeight = 12),
  "Tara Oceans" = makeIcon("/Users/paolil/polybox/PhD/Projects/GENOME_RESOLVED_OCEAN_EAN/code/analysis/lib/icon_circle_tara.png", iconWidth = 12, iconHeight = 12),
  "Hawaiian Ocean Time-Series" = makeIcon("/Users/paolil/polybox/PhD/Projects/GENOME_RESOLVED_OCEAN_EAN/code/analysis/lib/icon_clock_hots.png", iconWidth = 20, iconHeight = 20),
  "Bermuda-Atlantic Time-Series" = makeIcon("/Users/paolil/polybox/PhD/Projects/GENOME_RESOLVED_OCEAN_EAN/code/analysis/lib/icon_clock_bats.png", iconWidth = 20, iconHeight = 20)
)

# Generate legend -> Need to save the file and open html to see the legends
html_legend <- paste(
  "<img src='", custom_icons$`Tara Oceans`$iconUrl, "'style='vertical-align:middle;width:20px;height:20px;'>", "<div style='display:inline-block;vertical-align:middle'><span style='font-style:italic'>Tara</span> Oceans</div>", "<br/>",
  "<img src='", custom_icons$Malaspina$iconUrl, "' style='vertical-align:middle;width:20px;height:20px;'>", "<div style='display:inline-block;vertical-align:middle'>Malaspina</div>", "<br/>",
  "<img src='", custom_icons$Biogeotraces$iconUrl, "' style='vertical-align:middle;width:20px;height:20px;'>", "<div style='display:inline-block;vertical-align:middle'>BioGEOTRACES</div>", "<br/>",
  "<img src='", custom_icons$`Hawaiian Ocean Time-Series`$iconUrl, "' style='vertical-align:middle;width:20px;height:20px;'>", "<div style='display:inline-block;vertical-align:middle'>HOT</div>", "<br/>",
  "<img src='", custom_icons$`Bermuda-Atlantic Time-Series`$iconUrl, "' style='vertical-align:middle;width:20px;height:20px;'>", "<div style='display:inline-block;vertical-align:middle'>BATS</div>", "<br/>"
)

# On the real dataset
leaflet_map = sample_markers %>%
  leaflet(data = ., options = leafletOptions(worldCopyJump = T, zoomControl = F, dragging = T)) %>% # , 
  addProviderTiles(providers$Esri.OceanBasemap) %>% # add `, options = providerTileOptions(noWrap = TRUE)` to remove duplicate worlds. Pb, it shows a ugly data not available. Fix here, but couldn't implement so far. https://stackoverflow.com/questions/47477956/nowrap-option-on-tilelayer-is-only-partially-working
  #addGraticule(style = list(color = "lightgrey", weight = 0.5)) %>%
  addMapPane("bgeo", zIndex = 605) %>% 
  addMapPane("mala", zIndex = 610) %>% 
  addMapPane("tara", zIndex = 615) %>% 
  addMapPane("hots", zIndex = 625) %>%
  addMapPane("bats", zIndex = 630) %>%
  addMarkers(data = sample_markers %>% filter(dataset == "Biogeotraces"), ~longitude, ~latitude, popup = ~text, label = ~Sample,
             icon = ~custom_icons[dataset], group = ~groups, options = pathOptions(pane = "bgeo")) %>%
  addMarkers(data = sample_markers %>% filter(dataset == "Malaspina"), ~longitude, ~latitude, popup = ~text, label = ~Sample,
             icon = ~custom_icons[dataset], group = ~groups, options = pathOptions(pane = "mala")) %>%
  addMarkers(data = sample_markers %>% filter(dataset == "Tara Oceans"), ~longitude, ~latitude, popup = ~text, label = ~Sample,
             icon = ~custom_icons[dataset], group = ~groups, options = pathOptions(pane = "tara")) %>% # , popupOptions = popupOptions(keepInView = TRUE) tested as an option to keep maxBounds, but both options are incompatible
  addMarkers(data = sample_markers %>% filter(dataset == "Hawaiian Ocean Time-Series"), ~longitude, ~latitude, popup = ~text, label = ~Sample,
             icon = ~custom_icons[dataset], group = ~groups, options = pathOptions(pane = "hots")) %>%
  addMarkers(data = sample_markers %>% filter(dataset == "Bermuda-Atlantic Time-Series"), ~longitude, ~latitude, popup = ~text, label = ~Sample,
             icon = ~custom_icons[dataset], group = ~groups, options = pathOptions(pane = "bats")) %>%
  #setMaxBounds(lng1 = -180, lat1 = -90, lng2 = 180, lat2 = 90) %>% # Nicer without larger bounds for user experience
  setView(lng = 0, lat = 20, zoom = 2) %>%
  #addLayersControl(overlayGroups = ~unique(groups)) %>%
  addControl(html = html_legend, position = "bottomleft")

leaflet_map

wrapper_save_html(leaflet_map, "data/processed/figures/Figure-1/material/Map-interactive.legend.html")
