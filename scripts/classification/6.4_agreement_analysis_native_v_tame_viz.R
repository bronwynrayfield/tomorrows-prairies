
###################################################################################
##      VISUALIZATION: SPATIAL AGREEMENT (NATIVE VS TAME VS NON-GRASS)           ##
###################################################################################

##----------------------- Alison Long - May 2026 -------------------------

# STEP 6: Visualize pairwise spatial agreement metrics to assess similarities and 
#         differences among existing grassland datasets. This version handles all 
#         grassland data products regardless of extent and calculates metrics at
#         the ecodistrict and CCS level as well as globally.

# original file name: agreement_analysis_native_v_tame_v4_viz

#### VISUALIZE AGREEMENT RESULTS FROM V4 SCRIPT ####

## ---- Load libraries ---- ##
library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(igraph)
library(ggraph)
library(ape)
library(vegan)
library(ggrepel)

## ---- Setup ---- ##

# Set directories
results_dir <- file.path("R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis/outputs_nattam_pairwise")

# Results file
#results_file <- file.path(results_dir, "pairwise_nattam_spatial_agreement_01Jun2026.csv")
results_file <- file.path(results_dir, "pairwise_nattam_spatial_agreement_02Jun2026_plicombined.csv")

# Metadata file
meta_file <- 
  
  ## ---- Load results ---- ##

# Load results file
pairwise_long <- read_csv(results_file, show_col_types = FALSE)

# Inspect pairs with no overlapping valid data (provincial datasets)
#pairwise_long %>% filter(!is.finite(d_jaccard) | !is.finite(kappa))

# Filter out pairs with no overlapping valid data
pairwise_clean <- pairwise_long %>%
  filter(N > 0)

# Filter out pairs with very little overlapping data (PLI subsections)
pairwise_clean <- pairwise_clean %>%
  filter(N > 3400)

# Filter out HFI
pairwise_clean <- pairwise_clean %>%
  filter(raster_i != "hfi",
         raster_j != "hfi")

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
  
  diag(mat_mat) <- 1 # 1 for similarity, 0 for dissimilarity
  
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

# Macro Jaccard = "How similar are the products across all classes, giving each class equal importance?"
# Macro F1 = "Across all classes, how balanced is agreement once precision and recall are both considered?"
# Multiclass MCC = "What is the overall correlation between the two classification maps?"
# Kappa = "How much more agreement exists than would be expected given the class frequencies of the two datasets?"

## ---- PCoA with explanatory variables ---- ##
# dist_mat <- as.matrix(1 - mat_mcc)
# dist_mat[!is.finite(dist_mat)] <-max(dist_mat, na.rm = TRUE)
# pcoa_res <- pcoa(as.dist(dist_mat))

D <- 1 - mat_jacc_macro

diag(D) <- 0

D_dist <- as.dist(D)

# run pcoa
pcoa_res <- pcoa(D_dist)

# extract first two axes
scores <- as.data.frame(pcoa_res$vectors[, 1:2])

scores$raster <- rownames(scores)

# percent variance explained
var_exp <- round(100 * pcoa_res$values$Relative_eig[1:2], 1)
var_exp

# join metadata
dataset_meta <- dataset_meta %>%
  filter(raster %in% rownames(mat_jacc_macro))

scores <- left_join(scores, dataset_meta, by = "raster")

# reorder metadata to match distance matrix
meta_ord <- dataset_meta %>%
  slice(match(
    rownames(mat_jacc_macro),
    raster
  ))

identical(
  meta_ord$raster,
  rownames(mat_jacc_macro)
)

# basic pcoa plot
ggplot(
  scores,
  aes(Axis.1, Axis.2)
) +
  geom_point(
    aes(color = source_type),
    size = 4
  ) +
  geom_text(
    aes(label = raster),
    nudge_y = 0.01
  ) +
  theme_bw() +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)")
  )

# continuous variables
fit_cont <- envfit(
  pcoa_res$vectors[, 1:2],
  meta_ord %>%
    select(resolution, year),
  permutations = 999
)

fit_cont

# categorical variables
fit_cat <- envfit(
  pcoa_res$vectors[, 1:2],
  meta_ord %>%
    select(source_type, extent, nattam),
  permutations = 999
)

fit_cat

# plot envfit vectors
plot(
  pcoa_res$vectors[,1:2],
  pch = 19,
  xlab = paste0("PCoA1 (", var_exp[1], "%)"),
  ylab = paste0("PCoA2 (", var_exp[2], "%)")
)

text(
  pcoa_res$vectors[,1:2],
  labels = rownames(pcoa_res$vectors),
  pos = 3
)

plot(fit_cont, add = TRUE, col = "red")

plot(fit_cat, add = TRUE, col = "blue")

# PERMANOVA using adonis2

# categorical variables
adonis2(
  D_dist ~ source_type + extent + nattam,
  data = meta_ord,
  permutations = 999
)

# mixed variables
adonis2(
  D_dist ~ source_type +
    extent +
    nattam +
    resolution +
    year,
  data = meta_ord,
  permutations = 999
)

# sequential tests by term
adonis2(
  D_dist ~ source_type +
    extent +
    nattam +
    resolution +
    year,
  data = meta_ord,
  permutations = 999,
  by = "margin" # since attributes are likely correlated
)

# before interpreting significant factor effects,
# run for any categorical variables that come out
# as sig in PERMANOVA to determine whether sig 
# reflects centroid separation vs differences in 
# dispersion among groups:
betadisper(D_dist, meta_ord$source_type) |> anova()


## ---- PCoA ---- ##
dist_mat <- as.matrix(1 - mat_mcc)

dist_mat[!is.finite(dist_mat)] <-
  max(dist_mat, na.rm = TRUE)

pcoa_res <- pcoa(as.dist(dist_mat))

# Keep only datasets with complete pairwise comparisons
# keep <- complete.cases(mat_mcc_pcoa)
# mat_mcc_pcoa <- mat_mcc_pcoa[keep, keep]

# Extract scores
scores <- as.data.frame(pcoa_res$vectors[, 1:2])
scores$dataset <- rownames(scores)

# Plot
p <- ggplot(
  scores,
  aes(Axis.1, Axis.2, label = dataset)
) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0,
    segment.alpha = 0.5,
    seed = 123
  ) +
  theme_bw() +
  labs(
    title = "PCoA of Triclass Dataset Agreement (MCC)",
    x = paste0(
      "PCoA 1 (",
      round(100 * pcoa_res$values$Relative_eig[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA 2 (",
      round(100 * pcoa_res$values$Relative_eig[2], 1),
      "%)"
    )
  )

p


## ---- Heatmaps ---- ##
pheatmap(mat_kappa,
         main = "Cohen's Kappa",
         color = colorRampPalette(c("blue","white","red"))(50),
         na_col = "black",
         filename = file.path(results_dir, "heatmap_kappa.png"))

# Remove HFI
mat_mcc_nohfi <- mat_mcc[rownames(mat_mcc) != "hfi",
                         colnames(mat_mcc) != "hfi"]

pheatmap(mat_mcc_nohfi,
         main = "Multi-class MCC",
         color = colorRampPalette(c("red","white","blue"))(50),
         na_col = "black",
         filename = file.path(results_dir, "heatmap_mcc.png"))

pheatmap(mat_jacc_macro,
         main = "Macro Jaccard Similarity",
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
                                paste0("heatmap_jaccard_similarity_class_", cls, ".png")))
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
  labs(title = "Jaccard similarity (macro) vs Kappa")

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
  theme_minimal() +
  labs(title = "Jaccard similarity (macro) vs MCC")

ggsave(file.path(results_dir, "scatter_jaccard_vs_mcc.png"),
       width = 8, height = 6, dpi = 300)

## ---- Clustering ---- ##

# Keep only rasters with full pairwise comparisons
keep <- complete.cases(mat_mcc_nohfi)
mat_mcc_clean <- mat_mcc_nohfi[keep, keep]

dist_mat <- as.dist(1 - mat_mcc_clean)
hc <- hclust(dist_mat)

png(file.path(results_dir, "dendrogram_mcc.png"),
    width = 1200, height = 1000, res = 150)

plot(hc, main = "Raster Clustering (MCC)")
dev.off()

## ---- Network visualization ---- ##
threshold <- 0.5

edges <- pairwise_clean %>%
  filter(mcc_multi > threshold)

g <- graph_from_data_frame(edges, directed = FALSE)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle(paste("Network (MCC >", threshold, ")"))

ggsave(file.path(results_dir, "network_mcc_050.png"),
       width = 8, height = 6, dpi = 300)

# For macro Jaccard
threshold <- 0.5

edges <- pairwise_clean %>%
  filter(jaccard_macro > threshold)

g <- graph_from_data_frame(edges, directed = FALSE)

ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle(paste("Network (Macro Jaccard >", threshold, ")"))

ggsave(file.path(results_dir, "network_macrojaccard_050.png"),
       width = 8, height = 6, dpi = 300)

## ---- Class driving disagreement ---- ##
# Across all raster pairs, which class most often has the poorest spatial agreement?
# Assign each raster pair to the class with the lowest Jaccard similarity
# Low Jaccard = low overlap = high disagreement
# The “worst class” = the class where the two rasters agree the least
# Limitations: ignores effect size, sensitive to ties/small differences, doesn't account for class importance
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

ggsave(file.path(results_dir, "class_driving_disagreement.png"),
       width = 8, height = 6, dpi = 300)

# Add magnitude (how much worse the worst class is)
pairwise_clean %>%
  mutate(
    min_jacc = pmin(jaccard_class0, jaccard_class1, jaccard_class2, na.rm = TRUE)
  )

# Evaluate how strongly one class drives disagreement
# Compare max and min Jaccard
pairwise_clean %>%
  mutate(
    max_jacc = pmax(jaccard_class0, jaccard_class1, jaccard_class2, na.rm = TRUE),
    min_jacc = pmin(jaccard_class0, jaccard_class1, jaccard_class2, na.rm = TRUE),
    disagreement_gap = max_jacc - min_jacc
  )

# Focus on only strong disagreement


## ---- Compute disagreement diagnostics ---- ##

pairwise_diag <- pairwise_clean %>%
  rowwise() %>%
  mutate(
    jacc_vec = list(c(jaccard_class0, jaccard_class1, jaccard_class2)),
    class_names = list(c("class0", "class1", "class2")),
    
    # Handle NA safely
    valid_idx = list(which(!is.na(jacc_vec))),
    
    worst_class = if (length(valid_idx) == 0) {
      NA_character_
    } else {
      class_names[valid_idx[which.min(jacc_vec[valid_idx])]]
    },
    
    best_class = if (length(valid_idx) == 0) {
      NA_character_
    } else {
      class_names[valid_idx[which.max(jacc_vec[valid_idx])]]
    },
    
    min_jacc = if (length(valid_idx) == 0) NA else min(jacc_vec, na.rm = TRUE),
    max_jacc = if (length(valid_idx) == 0) NA else max(jacc_vec, na.rm = TRUE),
    
    disagreement_gap = max_jacc - min_jacc
  ) %>%
  ungroup()

# Plot 1: Which class most often drives disagreement
# Shows which class is weakest most often
plot_worst_class <- ggplot(pairwise_diag, aes(worst_class)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Class Most Frequently Driving Disagreement",
       x = "Class",
       y = "Number of Raster Pairs")

ggsave(file.path(results_dir, "worst_class_counts.png"),
       plot_worst_class, width = 7, height = 5, dpi = 300)

# Plot 2: Strength of disagreement (distribution)
# Near 0 = classes behave similarly
# Large values = one class is much worse
plot_gap <- ggplot(pairwise_diag, aes(disagreement_gap)) +
  geom_histogram(bins = 30, fill = "tomato", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Distribution of Disagreement Strength",
       x = "Jaccard Difference (best - worst)",
       y = "Count")

ggsave(file.path(results_dir, "disagreement_gap_hist.png"),
       plot_gap, width = 7, height = 5, dpi = 300)

# Plot 3: Worst class colored by magnitude
# Shows not just frequency, but severity
# Identifies which class causes strongest disagreement
plot_worst_vs_gap <- ggplot(pairwise_diag,
                            aes(worst_class, disagreement_gap, fill = worst_class)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Strength of Disagreement by Class",
       x = "Worst Class",
       y = "Disagreement Gap")

ggsave(file.path(results_dir, "worst_class_vs_gap.png"),
       plot_worst_vs_gap, width = 8, height = 6, dpi = 300)

# Plot 4: Focus on strong disagreements only
plot_strong <- pairwise_diag %>%
  filter(min_jacc < 0.3) %>%
  ggplot(aes(worst_class)) +
  geom_bar(fill = "darkred") +
  theme_minimal() +
  labs(title = "Class Driving Strong Disagreement (Jaccard < 0.3)")

ggsave(file.path(results_dir, "worst_class_strong_only.png"),
       plot_strong, width = 7, height = 5, dpi = 300)

# Plot 5: Which class is most variable overall?
# Shows which class is most consistent (tight spread) vs most variable (wide spread)
pairwise_long_format <- pairwise_clean %>%
  pivot_longer(
    cols = starts_with("jaccard_class"),
    names_to = "class",
    values_to = "jaccard"
  )

plot_variability <- ggplot(pairwise_long_format, aes(class, jaccard, fill = class)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Jaccard Similarity by Class",
       x = "Class",
       y = "Jaccard Similarity")

ggsave(file.path(results_dir, "jaccard_by_class_boxplot.png"),
       plot_variability, width = 8, height = 6, dpi = 300)

# Plot 6: Link disagreement to MCC
# Shows whether class-specific disagreement explains overall disagreement
plot_gap_mcc <- ggplot(pairwise_diag,
                       aes(disagreement_gap, mcc_multi, color = worst_class)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Disagreement Strength vs MCC")

ggsave(file.path(results_dir, "gap_vs_mcc.png"),
       plot_gap_mcc, width = 8, height = 6, dpi = 300)

#---------- Disagreement drivers by class -----------#
# Which raster datasets consistently produce low agreement for each class?

# Reshape to long format
pairwise_long_class <- pairwise_clean %>%
  pivot_longer(
    cols = starts_with("jaccard_class"),
    names_to = "class",
    values_to = "jaccard"
  )

# Identify low-agreement pairs
threshold <- 0.3

low_agreement <- pairwise_long_class %>%
  filter(jaccard < threshold)

# Count datasets contributing to disagreement
dataset_disagreement <- low_agreement %>%
  select(raster_i, raster_j, class) %>%
  pivot_longer(cols = c(raster_i, raster_j),
               values_to = "dataset") %>%
  group_by(dataset, class) %>%
  summarise(n = n(), .groups = "drop")

# Plot
plot_dataset_disagreement <- ggplot(dataset_disagreement,
                                    aes(x = dataset, y = n, fill = class)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(title = "Datasets Driving Disagreement by Class",
       x = "Dataset",
       y = "Number of Low-Agreement Pairs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "datasets_driving_disagreement.png"),
       plot_dataset_disagreement,
       width = 10, height = 6, dpi = 300)

# Normalize by number of comparisons

# Total number of comparisons per dataset
total_pairs <- pairwise_long_class %>%
  pivot_longer(cols = c(raster_i, raster_j),
               values_to = "dataset") %>%
  count(dataset, name = "total_pairs")

# Number of low-agreement cases per dataset
dataset_disagreement <- dataset_disagreement %>%
  rename(low_pairs = n) %>%
  left_join(total_pairs, by = "dataset") %>%
  mutate(prop = low_pairs / total_pairs)

# Plot
plot_dataset_disagreement <- ggplot(dataset_disagreement,
                                    aes(x = dataset, y = prop, fill = class)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(title = "Datasets Driving Disagreement by Class (Normalized)",
       x = "Dataset",
       y = "Proportion of Low-Agreement Pairs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "datasets_driving_disagreement_norm.png"),
       plot_dataset_disagreement,
       width = 10, height = 6, dpi = 300)

#------------- Map spatial disagreement hotspots -----------#
# Stack rasters
r_stack <- rast(rfiles_local)

# Compute disagreement per pixel

# Class-specific disagreement

# Variance based disagreement
# Variance = disagreement intensity

class_disagreement <- function(target_class) {
  
  binary <- r_stack == target_class
  
  app(binary, fun = function(x) {
    if (all(is.na(x))) return(NA)
    
    # variance of binary vector
    v <- var(x, na.rm = TRUE)
    
    return(v/0.25) # normalize to 0-1 where 0=full agreement, 1=max disagreement
  })
}

# Easier to interpret: mean-based disagreement
class_disagreement <- function(target_class) {
  
  binary <- r_stack == target_class
  
  app(binary, fun = function(x) {
    if (all(is.na(x))) return(NA)
    
    p <- mean(x, na.rm = TRUE)
    
    # Disagreement highest at p = 0.5
    return(1 - abs(p - 0.5) * 2)
  })
}

# Compute maps
class0_map <- class_disagreement(0)
class1_map <- class_disagreement(1)
class2_map <- class_disagreement(2)

# Save rasters 
writeRaster(class0_map,
            file.path(local_dir, "disagreement_class0.tif"),
            overwrite = TRUE)

writeRaster(class1_map,
            file.path(local_dir, "disagreement_class1.tif"),
            overwrite = TRUE)

writeRaster(class2_map,
            file.path(local_dir, "disagreement_class2.tif"),
            overwrite = TRUE)

# Quick plot
plot(class1_map, main = "Spatial Disagreement (Class 1)")


