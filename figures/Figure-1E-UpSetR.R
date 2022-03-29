# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Figure 1E - UpSetR =====================================================================

# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(UpSetR)
library(VennDiagram)
for (src_file in list.files('code/analysis/R')){
  print(paste('Sourcing', src_file, '...'))
  source(paste('code/analysis/R', src_file, sep = '/'))}
source('code/analysis/lib/Cell_Press_guidelines.R')
source('../../Exploratorium/sushilab-colors/palettes-paoli.R')

# Load data ------------------------------------------------------------------------------

summary_table = load_raw_summary()

# UpSetR version -------------------------------------------------------------------------

p = upset(
  fromList(
    list(
      "Internal_MAGs" = summary_table %>% filter(grepl("METAG_[A-Z]{8}$", `Bin Id`)) %>% pull(drep_cluster) %>% unique(),
      "Delmont_MAGs" = summary_table %>% filter(grepl("TARA_.*_METAG_[A-Z]{3}[0-9]", `Bin Id`)) %>% pull(drep_cluster) %>% unique(),
      "SAGs" = summary_table %>% filter(grepl("_SAGS_", `Bin Id`)) %>% pull(drep_cluster) %>% unique(),
      "REFs" = summary_table %>% filter(grepl("_REFG_", `Bin Id`)) %>% pull(drep_cluster) %>% unique()
    )),
  intersections = list(
    list("REFs"),
    list("SAGs"),
    list("Delmont_MAGs"),
    list("Internal_MAGs"),
    list("Internal_MAGs", "Delmont_MAGs"),
    list("Internal_MAGs", "REFs"),
    list("Internal_MAGs", "SAGs"),
    list("Delmont_MAGs", "REFs"),
    list("Delmont_MAGs", "SAGs"),
    list("REFs", "SAGs"),
    list("Internal_MAGs", "Delmont_MAGs", "REFs"),
    list("Internal_MAGs", "Delmont_MAGs", "SAGs"),
    list("Internal_MAGs", "REFs", "SAGs"),
    list("Delmont_MAGs", "REFs", "SAGs"),
    list("Internal_MAGs", "Delmont_MAGs", "REFs", "SAGs")
  ),
  #order.by = "degree",
  keep.order = T,
  empty.intersections = "on",
  text.scale = 0.8,
  line.size = NA,
  #sets.bar.color = c("red", "blue", "green", "yellow"),
  mb.ratio = c(.6, .4),
  sets.x.label = "# Species",
  mainbar.y.label = "# Species", 
)

pdf(file=paste0(figures_path_proj, "Figure-1/Figure-1F.pdf"), width = 6.4, height = 2)
p
dev.off()

# Venn Diagram version -------------------------------------------------------------------

diagram_cols = c(viridis::magma(begin = .2, end = .8, n=4)[1],
                 dataset_colors["MarDB"],
                 dataset_colors["Delmont 2018 MAGs"],
                 dataset_colors["GORG SAGs"])

#Make the plot
venn.diagram(
  x = list(
    summary_table %>% filter(grepl("METAG_[A-Z]{8}$", `Bin Id`)) %>% pull(drep_cluster) %>% unique(),
    summary_table %>% filter(grepl("_REFG_", `Bin Id`)) %>% pull(drep_cluster) %>% unique(),
    summary_table %>% filter(grepl("TARA_.*_METAG_[A-Z]{3}[0-9]", `Bin Id`)) %>% pull(drep_cluster) %>% unique(),
    summary_table %>% filter(grepl("_SAGS_", `Bin Id`)) %>% pull(drep_cluster) %>% unique()
  ),
  category.names = c(
    paste(c("Reconstr. MAGs",
            paste("# Genomes:", summary_table %>% filter(grepl("METAG_[A-Z]{8}$", `Bin Id`)) %>% nrow() %>% prettyNum(big.mark = ",")),
            paste("# Clusters", summary_table %>% filter(grepl("METAG_[A-Z]{8}$", `Bin Id`)) %>% pull(drep_cluster) %>% unique() %>% length()) %>% prettyNum(big.mark = ",")),
          collapse = "\n"),
    paste(c("Ref. Genomes",
            paste("# Genomes:", summary_table %>% filter(grepl("_REFG_", `Bin Id`)) %>% nrow() %>% prettyNum(big.mark = ",")),
            paste("# Clusters", summary_table %>% filter(grepl("_REFG_", `Bin Id`)) %>% pull(drep_cluster) %>% unique() %>% length()) %>% prettyNum(big.mark = ",")),
          collapse = "\n"),
    paste(c("Delmont MAGs",
            paste("# Genomes:", summary_table %>% filter(grepl("TARA_.*_METAG_[A-Z]{3}[0-9]", `Bin Id`)) %>% nrow() %>% prettyNum(big.mark = ",")),
            paste("# Clusters", summary_table %>% filter(grepl("TARA_.*_METAG_[A-Z]{3}[0-9]", `Bin Id`)) %>% pull(drep_cluster) %>% unique() %>% length()) %>% prettyNum(big.mark = ",")),
          collapse = "\n"),
    paste(c("SAGs",
            paste("# Genomes:", summary_table %>% filter(grepl("_SAGS_", `Bin Id`)) %>% nrow() %>% prettyNum(big.mark = ",")),
            paste("# Clusters", summary_table %>% filter(grepl("_SAGS_", `Bin Id`)) %>% pull(drep_cluster) %>% unique() %>% length()) %>% prettyNum(big.mark = ",")),
          collapse = "\n")
  ),
  filename = 'data/processed/figures/Figure-2/Figure-2D.png',
  output = TRUE,
  imagetype="png",
  height = 1560, 
  width = 1920, 
  resolution = 500,
  lwd = 1,
  col = diagram_cols,
  fill = c(alpha(diagram_cols[1],0.3), alpha(diagram_cols[2],0.3),alpha(diagram_cols[3],0.3),  alpha(diagram_cols[4],0.3)),
  cex = 1,
  fontfamily = "sans",
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = diagram_cols,
  cat.just = list(c(.4,.3), c(.6,.3), c(.5,.3), c(.5,.3)))
