
###################################################################################
##          ANALYSIS: SPATIAL AGREEMENT (NATIVE VS TAME VS NON-GRASS)            ##
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

## ---- Load results ---- ##

# Set directories
results_dir <- file.path("R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis/outputs_natvtam_pairwise")

# Results file
results_file <- file.path(results_dir, "pairwise_nattam_spatial_agreement_01Jun2026.csv")

# Load results file
pairwise_long <- read_csv(results_file, show_col_types = FALSE)

# Inspect pairs with no overlapping valid data (provincial datasets)
#pairwise_long %>% filter(!is.finite(d_jaccard) | !is.finite(kappa))

# Filter out pairs with no overlapping valid data
pairwise_clean <- pairwise_long %>%
  filter(N > 0)

## ---- Helper: convert to matrix ---- ##
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
mat_kappa      <- to_matrix(pairwise_clean, "kappa")
mat_jacc_macro <- to_matrix(pairwise_clean, "jaccard_macro")
mat_f1_macro   <- to_matrix(pairwise_clean, "f1_macro")
mat_mcc        <- to_matrix(pairwise_clean, "mcc_multi")
mat_acc        <- to_matrix(pairwise_clean, "accuracy")

# Make cells with value = 1 NA (no valid overlap)
#mat_jaccard[mat_jaccard == 1]   <- NA
#mat_sorensen[mat_sorensen == 1] <- NA
#mat_kappa[mat_kappa == 1]       <- NA

## ---- Heatmaps ---- ##
pheatmap(mat_kappa,
         main = "Cohen's Kappa",
         color = colorRampPalette(c("blue","white","red"))(50),
         na_col = "black",
         filename = file.path(results_dir, "heatmap_kappa.png"))

pheatmap(mat_mcc,
         main = "Multi-class MCC",
         color = colorRampPalette(c("blue","white","red"))(50),
         na_col = "black",
         filename = file.path(results_dir, "heatmap_mcc.png"))

pheatmap(mat_jacc_macro,
         main = "Macro Jaccard",
         na_col = "black",
         filename = file.path(results_dir, "heatmap_jaccard_macro.png"))

pheatmap(mat_f1_macro,
         main = "Macro F1 Score",
         na_col = "black",
         filename = file.path(results_dir, "heatmap_f1_macro.png"))

## ---- Class-specific heatmaps ---- ##
for (cls in c("0","1","2")) {
  
  mat <- to_matrix(pairwise_clean, paste0("jaccard_class", cls))
  
  pheatmap(mat,
           main = paste("Jaccard - Class", cls),
           na_col = "black",
           filename = file.path(results_dir,
                                paste0("heatmap_jaccard_class_", cls, ".png")))
}

## ---- Distribution plots ---- ##
pairwise_plot <- pairwise_clean %>%
  select(kappa, jaccard_macro, f1_macro, mcc_multi, accuracy) %>%
  pivot_longer(everything(),
               names_to = "metric",
               values_to = "value")

ggplot(pairwise_plot, aes(x = value, fill = metric)) +
  geom_histogram(bins = 40, alpha = 0.6) +
  facet_wrap(~metric, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Agreement Metrics")

ggsave(file.path(results_dir, "distribution_metrics.png"),
       width = 8, height = 6, dpi = 300)

## ---- Metric relationships ---- ##
ggplot(pairwise_clean, aes(jaccard_macro, kappa, color = N)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Jaccard (macro) vs Kappa")

ggsave(file.path(results_dir, "scatter_jaccard_vs_kappa.png"),
       width = 8, height = 6, dpi = 300)


ggplot(pairwise_clean, aes(mcc_multi, kappa, color = N)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "MCC vs Kappa")

ggsave(file.path(results_dir, "scatter_mcc_vs_kappa.png"),
       width = 8, height = 6, dpi = 300)


ggplot(pairwise_clean, aes(jaccard_macro, mcc_multi)) +
  geom_point() +
  theme_minimal()

ggsave(file.path(results_dir, "scatter_jaccard_vs_mcc.png"),
       width = 8, height = 6, dpi = 300)

## ---- Clustering ---- ##
dist_mat <- as.dist(1 - mat_mcc)
hc <- hclust(dist_mat)

png(file.path(results_dir, "dendrogram_mcc.png"),
    width = 1200, height = 1000, res = 150)

plot(hc, main = "Raster Clustering (MCC)")
dev.off()

## ---- Network visualization ---- ##
threshold <- 0.6

edges <- pairwise_clean %>%
  filter(mcc_multi > threshold)

g <- graph_from_data_frame(edges, directed = FALSE)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle(paste("Network (MCC >", threshold, ")"))

ggsave(file.path(results_dir, "network_mcc.png"),
       width = 8, height = 6, dpi = 300)

# For Kappa
threshold <- 0.6

edges <- pairwise_clean %>%
  filter(kappa > threshold)

g <- graph_from_data_frame(edges, directed = FALSE)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle(paste("Network (Kappa >", threshold, ")"))

ggsave(file.path(results_dir, "network_kappa.png"),
       width = 8, height = 6, dpi = 300)

## ---- Class driving disagreement ---- ##
pairwise_clean %>%
  mutate(
    worst_class = case_when(
      jaccard_class0 <= jaccard_class1 & jaccard_class0 <= jaccard_class2 ~ "class0",
      jaccard_class1 <= jaccard_class2 ~ "class1",
      TRUE ~ "class2"
    )
  ) %>%
  ggplot(aes(worst_class)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Class Driving Disagreement")
