
## ---- Load libraries ---- ##
library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(igraph)
library(ggraph)

## ---- Load results ---- ##
# Update path if needed
results_file <- file.path(Sys.getenv("HOME"), "r_tmp", "pairwise_results.csv")

pairwise_long <- read_csv(results_file, show_col_types = FALSE)

## ---- Helper: convert to matrix ---- ##
to_matrix <- function(df, value_col) {
  
  mat <- df %>%
    select(raster_i, raster_j, !!sym(value_col)) %>%
    pivot_wider(names_from = raster_j, values_from = !!sym(value_col))
  
  mat_df <- as.data.frame(mat)
  rownames(mat_df) <- mat_df$raster_i
  mat_df$raster_i <- NULL
  
  mat_mat <- as.matrix(mat_df)
  
  # Make symmetric
  mat_mat[lower.tri(mat_mat)] <- t(mat_mat)[lower.tri(mat_mat)]
  diag(mat_mat) <- 0
  
  return(mat_mat)
}

## ---- Build matrices ---- ##
mat_jaccard  <- to_matrix(pairwise_long, "d_jaccard")
mat_sorensen <- to_matrix(pairwise_long, "d_sorensen")
mat_kappa    <- to_matrix(pairwise_long, "kappa")

## ---- Heatmaps ---- ##
pheatmap(mat_jaccard,
         main = "Jaccard Dissimilarity",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

pheatmap(mat_sorensen,
         main = "Sorensen Dissimilarity")

pheatmap(mat_kappa,
         main = "Cohen's Kappa",
         color = colorRampPalette(c("blue","white","red"))(50))

## ---- Long format for ggplot ---- ##
pairwise_plot <- pairwise_long %>%
  pivot_longer(cols = c(d_jaccard, d_sorensen, kappa),
               names_to = "metric",
               values_to = "value")

## ---- Distribution plots ---- ##
ggplot(pairwise_plot, aes(x = value, fill = metric)) +
  geom_histogram(bins = 40, alpha = 0.6) +
  facet_wrap(~metric, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Pairwise Agreement Metrics")

## ---- Scatter relationships ---- ##
ggplot(pairwise_long, aes(d_jaccard, kappa)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Jaccard vs Kappa")

ggplot(pairwise_long, aes(d_sorensen, kappa)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Sorensen vs Kappa")

## ---- Clustering dendrogram ---- ##
dist_mat <- as.dist(mat_jaccard)
hc <- hclust(dist_mat)

plot(hc, main = "Raster Clustering (Jaccard)")

## ---- Network visualization (threshold) ---- ##
# Define similarity threshold (adjust!)
threshold <- 0.7

edges <- pairwise_long %>%
  filter(kappa > threshold)

g <- graph_from_data_frame(edges, directed = FALSE)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle(paste("Raster Similarity Network (Kappa >", threshold, ")"))