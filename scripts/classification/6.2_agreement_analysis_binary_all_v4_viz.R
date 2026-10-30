
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (BINARY GRASS VS NON-GRASS)              ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# step 6: calculate pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent.

#### VISUALIZE AGREEMENT RESULTS FROM V4 SCRIPT ####

## ---- Load libraries ---- ##
library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(igraph)
library(ggraph)
library(ape) # for PCoA
library(ggrepel)

## ---- Load results ---- ##

# Set directories
results_dir <- file.path("R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis/outputs_binary_pairwise")

# Results file
results_file <- file.path(results_dir, "results_1Jun2026/pairwise_binary_spatial_agreement_01Jun2026.csv")

# Load results file
pairwise_long <- read_csv(results_file, show_col_types = FALSE)

# Range of values for MCC
range(pairwise_long$mcc, na.rm = TRUE)

# Inspect pairs with no overlapping valid data (provincial datasets)
pairwise_long %>% filter(!is.finite(d_jaccard) | !is.finite(kappa))

# Filter out pairs with no overlapping valid data
pairwise_clean <- pairwise_long %>%
  filter(N > 0)

## ---- Helper: convert to matrix ---- ##
#to_matrix <- function(df, value_col) {
  
#  mat <- df %>%
#    select(raster_i, raster_j, !!sym(value_col)) %>%
#    pivot_wider(names_from = raster_j, values_from = !!sym(value_col))
  
#  mat_df <- as.data.frame(mat)
#  rownames(mat_df) <- mat_df$raster_i
#  mat_df$raster_i <- NULL
  
#  mat_mat <- as.matrix(mat_df)
  
  # Make symmetric
#  mat_mat[lower.tri(mat_mat)] <- t(mat_mat)[lower.tri(mat_mat)]
#  diag(mat_mat) <- 0
  
#  return(mat_mat)
#}

to_matrix <- function(df, value_col) {
  
  # Get full set of raster names
  rasters <- sort(unique(c(df$raster_i, df$raster_j)))
  
  # Create full grid of all combinations
  full_grid <- expand.grid(
    raster_i = rasters,
    raster_j = rasters,
    stringsAsFactors = FALSE
  )
  
  # Join values onto full grid
  df_full <- full_grid %>%
    left_join(df, by = c("raster_i", "raster_j"))
  
  # Pivot to wide format with enforced column order
  mat <- df_full %>%
    select(raster_i, raster_j, !!sym(value_col)) %>%
    pivot_wider(
      names_from = raster_j,
      values_from = !!sym(value_col)
    )
  
  # Convert to matrix
  mat_df <- as.data.frame(mat)
  rownames(mat_df) <- mat_df$raster_i
  mat_df$raster_i <- NULL
  
  mat_mat <- as.matrix(mat_df)
  
  # Enforce same ordering on rows & columns
  mat_mat <- mat_mat[rasters, rasters]
  
  # Symmetrize safely
  for (i in seq_len(nrow(mat_mat))) {
    for (j in seq_len(ncol(mat_mat))) {
      if (is.na(mat_mat[i, j]) && !is.na(mat_mat[j, i])) {
        mat_mat[i, j] <- mat_mat[j, i]
      } else if (!is.na(mat_mat[i, j]) && is.na(mat_mat[j, i])) {
        mat_mat[j, i] <- mat_mat[i, j]
      }
    }
  }
  
  diag(mat_mat) <- 0
  
  return(mat_mat)
}

## ---- Build matrices ---- ##
mat_jaccard  <- to_matrix(pairwise_long, "d_jaccard")
mat_sorensen <- to_matrix(pairwise_long, "d_sorensen")
mat_kappa    <- to_matrix(pairwise_long, "kappa")
mat_f1       <- to_matrix(pairwise_long, "f1")
mat_mcc      <- to_matrix(pairwise_long, "mcc")

# Make cells with value = 1 NA (no valid overlap)
mat_jaccard[mat_jaccard == 1]   <- NA
mat_sorensen[mat_sorensen == 1] <- NA
mat_kappa[mat_kappa == 1]       <- NA
mat_f1[mat_f1 == 1]             <- NA
mat_mcc[mat_mcc == 1]           <- NA

## ---- PCoA ---- ##
# Remove datasets with incomplete comparisons
keep <- complete.cases(mat_jaccard)
mat_jaccard_clean <- mat_jaccard[keep, keep]

pcoa_res <- pcoa(as.dist(mat_jaccard_clean))

# Extract coordinates
scores <- as.data.frame(pcoa_res$vectors[,1:2])
scores$dataset <- rownames(scores)

# Plot
ggplot(scores,
       aes(Axis.1, Axis.2, label = dataset)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5) +
  theme_minimal() +
  labs(
    x = paste0(
      "PCoA1 (",
      round(pcoa_res$values$Relative_eig[1] * 100, 1),
      "%)"),
    y = paste0(
      "PCoA2 (",
      round(pcoa_res$values$Relative_eig[2] * 100, 1),
      "%)"))

# Or replace missing comparisons with the maximum observed distance
mat_jaccard[is.na(mat_jaccard)] <- max(mat_jaccard, na.rm = TRUE)

pcoa_res2 <- pcoa(as.dist(mat_jaccard))

# Extract coordinates
scores2 <- as.data.frame(pcoa_res2$vectors[,1:3])
scores2$dataset <- rownames(scores2)

# Plot
ggplot(scores2,
       aes(Axis.1, Axis.2, label = dataset)) +
  geom_point(size = 3) +
  #geom_text(vjust = -0.5) +
  geom_text_repel(
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0,
    segment.alpha = 0.5) +
  theme_minimal() +
  labs(
    x = paste0(
      "PCoA1 (",
      round(pcoa_res2$values$Relative_eig[1] * 100, 1),
      "%)"),
    y = paste0(
      "PCoA2 (",
      round(pcoa_res2$values$Relative_eig[2] * 100, 1),
      "%)"))

# amount of variance explained
pcoa_res2$values$Relative_eig[1:10]

# look at 3rd axis
ggplot(scores2,
       aes(Axis.1, Axis.3, label = dataset)) +
  geom_point() +
  geom_text_repel()

ggplot(scores2,
       aes(Axis.2, Axis.3, label = dataset)) +
  geom_point() +
  geom_text_repel()

# evaluate ordination quality
# If the relationship is reasonably strong, the 2D plot is still useful
ord_dist <- dist(scores2[,1:2])

plot(
  as.vector(as.dist(mat_jaccard)),
  as.vector(ord_dist)
)
abline(0,1,col="red")

# scree plot
eig <- 100 * pcoa_res2$values$Relative_eig

qplot(
  seq_along(eig),
  eig,
  geom = c("line","point")
) +
  labs(
    x = "PCoA Axis",
    y = "% Variation Explained"
  )


## ---- Heatmaps ---- ##
pheatmap(mat_jaccard,
         main = "Jaccard Dissimilarity",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         na_col = "black",
         filename = file.path(results_dir, "heatmap_jaccard.png"))

pheatmap(mat_sorensen,
         main = "Sorensen Dissimilarity",
         na_col = "black",
         filename = file.path(results_dir, "heatmap_sorensen.png"))

pheatmap(mat_kappa,
         main = "Cohen's Kappa",
         color = colorRampPalette(c("blue","white","red"))(50),
         na_col = "black",
         filename = file.path(results_dir, "heatmap_kappa.png"))

pheatmap(mat_f1,
         main = "F1 Score (Macro)",
         color = colorRampPalette(c("white", "orange", "red"))(50),
         na_col = "black",
         filename = file.path(results_dir, "heatmap_f1.png"))

#pheatmap(mat_mcc,
#         main = "Matthews Correlation Coefficient (MCC)",
#         color = colorRampPalette(c("blue","white","red"))(50),
#         breaks = seq(-1, 1, length.out = 51),
#         na_col = "black",
#         filename = file.path(results_dir, "heatmap_mcc.png"))

pheatmap(mat_mcc,
         main = "Matthews Correlation Coefficient (MCC)",
         color = colorRampPalette(c("blue","white","red"))(50),
         breaks = seq(-1, 1, length.out = 51),
         na_col = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         filename = file.path(results_dir, "heatmap_mcc.png"))


## ---- Long format for ggplot ---- ##
pairwise_plot <- pairwise_clean %>%
  pivot_longer(cols = c(d_jaccard, d_sorensen, kappa, f1, mcc),
               names_to = "metric",
               values_to = "value")

## ---- Distribution plots ---- ##
dist <- ggplot(pairwise_plot, aes(x = value, fill = metric)) +
  geom_histogram(bins = 40, alpha = 0.6) +
  facet_wrap(~metric, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Pairwise Agreement Metrics")

# Save as PNG
ggsave(filename = file.path(results_dir, "distribution_metrics.png"),
       plot = dist,
       width = 8,
       height = 6,
       dpi = 300)

## ---- Scatter relationships ---- ##
# Does overlap translate into real agreement?
scatt.fm <- ggplot(pairwise_clean, aes(f1, mcc, color = N)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "F1 vs MCC")

ggsave(file.path(results_dir, "scatter_f1_vs_mcc_N.png"),
       scatt.fm, width = 8, height = 6, dpi = 300)

# Direct comparison of overlap (Jaccard) and balanced agreement (MCC)
scatt.jm <- ggplot(pairwise_clean, aes(1 - d_jaccard, mcc, color = N)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Jaccard Similarity vs MCC")

ggsave(file.path(results_dir, "scatter_jaccard_vs_mcc_N.png"),
       scatt.jm, width = 8, height = 6, dpi = 300)

# How distorted is Kappa relative to a better metric?
scatt.km <- ggplot(pairwise_clean, aes(kappa, mcc, color = N)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Kappa vs MCC")

ggsave(file.path(results_dir, "scatter_kappa_vs_mcc_N.png"),
       scatt.km, width = 8, height = 6, dpi = 300)

# Consistency check - should be almost perfectly correlated
scatt.fj <- ggplot(pairwise_clean, aes(1 - d_jaccard, f1)) +
  geom_point() +
  theme_minimal() +
  labs(title = "F1 vs Jaccard Similarity")

ggsave(file.path(results_dir, "scatter_f1_vs_jaccard.png"),
       scatt.fj, width = 8, height = 6, dpi = 300)

# Does agreement depend on how much area overlaps?
scatt.mn <- ggplot(pairwise_clean, aes(N, mcc)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "MCC vs Overlap (N)")

ggsave(file.path(results_dir, "scatter_mcc_vs_N.png"),
       scatt.mn, width = 8, height = 6, dpi = 300)

# As Jaccard dissimilarity increases, Kappa should decrease 
# Strong negative relationship means metrics are consistent, both capture the same structure of agreement
# Kappa is frequently near zero even when Jaccard varies widely - some raster pairs have little better-than-random agreement,
# even when their overlap (Jaccard) differs. Kappa is very sensitive to marginal totals while Jaccard is not -
# Kappa penalizes imbalance, chance agreement, absence (zeroes)
# Kappa can also drop due to spatial mismatch structure - some overlap but disagreement in WHERE values occur
# Kappa may not be a reliable summary metric for land cover classification
scatt.jk <- ggplot(pairwise_clean, aes(d_jaccard, kappa, color = N)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Jaccard vs Kappa")

# Save as PNG
ggsave(filename = file.path(results_dir, "scatter_jaccard_vs_kappa_N.png"),
       plot = scatt.jk,
       width = 8,
       height = 6,
       dpi = 300)

scatt.sk <- ggplot(pairwise_clean, aes(d_sorensen, kappa, color = N)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Sorensen vs Kappa")

# Save as PNG
ggsave(filename = file.path(results_dir, "scatter_sorensen_vs_kappa_N.png"),
       plot = scatt.sk,
       width = 8,
       height = 6,
       dpi = 300)


scatt.js <- ggplot(pairwise_clean, aes(d_jaccard, d_sorensen)) +
  geom_point() +
  theme_minimal()

# Save as PNG
ggsave(filename = file.path(results_dir, "scatter_sorensen_vs_jaccard.png"),
       plot = scatt.js,
       width = 8,
       height = 6,
       dpi = 300)

# Correlation between jaccard and kappa
cor(pairwise_clean$d_jaccard, pairwise_clean$kappa)

# Calculate proportion of presence in each raster to see if this explains pattern

## ---- Clustering dendrogram ---- ##
# Rebuild matrices but don't replace 1 with NA
mat_jaccard  <- to_matrix(pairwise_long, "d_jaccard")
dist_mat <- as.dist(mat_jaccard)
hc <- hclust(dist_mat)

png(file.path(results_dir, "dendrogram_jaccard.png"),
    width = 1200, height = 1000, res = 150)

plot(hc, main = "Raster Clustering (Jaccard)")

dev.off()

# Exclude pairs with no overlap
# Keep only rows/cols with no NA
mat_jaccard[mat_jaccard == 1]   <- NA
keep <- complete.cases(mat_jaccard)
mat_jaccard_clean <- mat_jaccard[keep, keep]

dist_mat <- as.dist(mat_jaccard_clean)
hc <- hclust(dist_mat)

png(file.path(results_dir, "dendrogram_jaccard_validpairsonly.png"),
    width = 1200, height = 1000, res = 150)

plot(hc, main = "Raster Clustering (Jaccard)")

dev.off()

## ---- MCC clustering (valid pairs only) ---- ##
# Raster clustering should be based on balanced, global similarity measure
# MCC is best:
#     Uses full confusion matrix
#     Robust to class imbalance
#     Not overly influenced by dominant classes
#     Handles multi-class structure correctly
#     Produces true correlation-like structure

# Remove HFI
mat_mcc_nohfi <- mat_mcc[rownames(mat_mcc) != "hfi",
                         colnames(mat_mcc) != "hfi"]
# Remove invalid MCC values
mat_mcc_nohfi[!is.finite(mat_mcc_nohfi)] <- NA

# Keep only rasters with full pairwise comparisons
keep <- complete.cases(mat_mcc_nohfi)
mat_mcc_clean <- mat_mcc_nohfi[keep, keep]

# Convert to distance
dist_mat <- as.dist(1 - mat_mcc_clean)

hc <- hclust(dist_mat)

png(file.path(results_dir, "dendrogram_mcc_validpairsonly.png"),
    width = 1200, height = 1000, res = 150)

plot(hc, main = "Raster Clustering (MCC, valid pairs only)")

dev.off()

# Keep all datasets - NAs become maximum observed distance
keep <- rowSums(is.finite(mat_mcc)) > 1
mat_mcc_clean <- mat_mcc[keep, keep]

# Convert similarity → distance
dist_mat <- as.matrix(1 - mat_mcc_clean)

# Replace NA distances with max distance
dist_mat[!is.finite(dist_mat)] <- max(dist_mat, na.rm = TRUE)

# Convert back to dist
dist_mat <- as.dist(dist_mat)
dist_mat

# Cluster
hc <- hclust(dist_mat)

png(file.path(results_dir, "dendrogram_mcc_allpairs.png"),
    width = 1200, height = 1000, res = 150)

plot(hc, main = "Raster Clustering (MCC, NAs as max distance)")

dev.off()

## ---- Network visualization (threshold) ---- ##
# Define similarity threshold (adjust!)
threshold <- 0.6

# Edges (aka lines) will only be retained where Kappa is above threshold
# Lines in plot mean high agreement beyond chance
edges <- pairwise_clean %>%
  filter(kappa > threshold)

g <- graph_from_data_frame(edges, directed = FALSE)

# Highly similar rasters cluster close together
# Nodes with many edges agree strongly with many others and may represent consensus product
net.plot <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle(paste("Raster Similarity Network (Kappa >", threshold, ")"))

# Save as PNG
ggsave(filename = file.path(results_dir, "network_kappa_threshold_060.png"),
       plot = net.plot,
       width = 8,
       height = 6,
       dpi = 300)

#---- Same but for MCC ----#
# Define similarity threshold (adjust!)
threshold <- 0.7

# Edges (aka lines) will only be retained where Kappa is above threshold
# Lines in plot mean high agreement beyond chance
edges <- pairwise_clean %>%
  filter(mcc > threshold)

g <- graph_from_data_frame(edges, directed = FALSE)

# Highly similar rasters cluster close together
# Nodes with many edges agree strongly with many others and may represent consensus product
net.plot <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle(paste("Raster Similarity Network (MCC >", threshold, ")"))

# Save as PNG
ggsave(filename = file.path(results_dir, "network_mcc_threshold_070.png"),
       plot = net.plot,
       width = 8,
       height = 6,
       dpi = 300)
