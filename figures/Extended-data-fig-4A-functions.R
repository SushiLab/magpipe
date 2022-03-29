# Biosynthetic potential of the global ocean microbiome ==================================
# Code used to produce the figures and analysis of the paper
# Extended Data Fig. 4A - Helper functions ===============================================

# Libraries ==============================================================================

library(tidyverse)

# Functions ==============================================================================

plot_umap_hdbscan_helper <- function(umap_layout, hdbscan, spread_table, metadata){
  res = NULL
  message(paste0("Adding background info for ", length(unique(hdbscan$cluster)), " clusters."))
  for (i in unique(hdbscan$cluster)){
    res = rbind(res,
                as_tibble(umap_layout) %>% 
                  mutate(sample = spread_table$sample) %>%
                  left_join(metadata %>% rename(sample = `Internal Sample Name`)) %>%
                  mutate(ocean_province = "Background") %>%
                  mutate(clusters = paste0("Cluster ", i)))
  }
  res = rbind(res,
              as_tibble(umap_layout) %>% 
                mutate(clusters = paste0("Cluster ", hdbscan$cluster)) %>%
                mutate(sample = spread_table$sample) %>%
                left_join(metadata %>% rename(sample = `Internal Sample Name`)))
  return(res)
}

hdbscan_auto_optimum <- function(umap_layout, hard_min = NULL, hard_max = 50){
  scores = tibble()
  if (is.null(hard_min)){hard_min = 2}
  for (k in hard_min:min(round(nrow(umap_layout)/2), hard_max)){
    tmp_hdbscan = dbscan::hdbscan(umap_layout, minPts = k)
    res = tibble(k = k,
                 score = sum(tmp_hdbscan$membership_prob),
                 n_unassign = sum(tmp_hdbscan$cluster == 0),
                 n_clusters = length(unique(tmp_hdbscan$cluster)))
    scores = rbind(scores, res)
    message(paste('For minPts =', res$k, '|', res$n_unassign, 'unassigned,',  'score of', res$score, 'and', res$n_clusters, 'clusters.'))
  }
  
  scores = scores %>% 
    mutate(adj_score = score - n_unassign) %>%
    arrange(desc(adj_score))
  
  message(paste('\nBest score of', scores$adj_score[1], 'with minPts =', scores$k[1], 'leads to', scores$n_clusters[1], 'clusters.'))
  
  hdbscan_res = dbscan::hdbscan(umap_layout, minPts = scores$k[1])
  return(hdbscan_res)
}

balanced_cluster_subselection <- function(table, size){
  res = NULL
  for (i in (table %>% filter(clusters != "Cluster 0") %>% pull(clusters) %>% unique())){
    res = c(res,
            table %>% filter(clusters == i) %>% pull(sample) %>% sample(., size = size))
  }
  return(res)
}

balanced_anova_R2_estimate <- function(cluster_tbl, dist_tbl, subsample_size, rep = 30){
  res_tbl = NULL
  for (i in 1:rep){
    message(paste0("Running iteration number ", i, "..."))
    subselection = balanced_cluster_subselection(cluster_tbl, subsample_size)
    subselected_distance = as.dist(dist_tbl[which(cluster_tbl$sample %in% subselection), which(cluster_tbl$sample %in% subselection)])
    subselected_clusters = cluster_tbl %>% filter(sample %in% subselection) %>% pull(clusters)
    permanova = vegan::adonis2(subselected_distance ~ as.character(subselected_clusters), by = 'margin')
    res_tbl = rbind(res_tbl, tibble(iteration = i, R2 = permanova$R2[1]))
  }
  return(res_tbl)
}