
###################################################################################
##                   ANALYSIS: SPATIAL AGREEMENT SUMMARIES                       ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets. This version includes only the grassland
#        data products that cover the full Prairies Ecozone extent.

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/class_inputs"
proc_grass_data_dir <- "./processeddatasets/grass"
out_grass_dir <- "./outputs/grass"

#-------- PIXELWISE AGREEMENT ------------

# List your raster files (modify pattern/path as needed)
files <- list.files("path/to/rasters", pattern = "\\.tif$", full.names = TRUE)

# Load as a SpatRaster stack
rstack <- rast(files)

# Count how many rasters have *non-NA* data for each pixel
count_nonNA <- app(rstack, fun = function(x) sum(!is.na(x)))

# Save result 
#writeRaster(count_nonNA, "pixel_data_count.tif", overwrite = TRUE)

#-------- AREA PER ECODISTRICT -----------------

# 1) Raster with data/NA
r <- rast("path/to/your_raster.tif")

# 2) Ecodistrict polygons
ed <- vect("path/to/ecodistricts.shp")

# Make sure CRS matches (reproject polygons to raster CRS if needed)
if (!same.crs(r, ed)) {
  ed <- project(ed, crs(r))
}

# (Optional but recommended) Create a stable ID field if none exists
if (is.null(ed$ecod_id)) {
  ed$ecod_id <- seq_len(nrow(ed))
}

# Create a raster of cell areas in hectares
area_ha <- cellSize(r, unit = "ha")   # SpatRaster of per-cell area in hectares

# Create a 1/0 mask of "has data"
has_data <- ifel(!is.na(r), 1, 0)

# Multiply mask by area to get hectares-with-data per cell
ha_with_data <- has_data * area_ha

# Zonal sum by ecodistrict; exact=TRUE weights partial cell overlaps by coverage
zsum <- extract(ha_with_data, ed, fun = sum, na.rm = TRUE, exact = TRUE)

# Join results back to polygons or keep as a table
out <- merge(as.data.frame(ed)[, c("ecod_id")], zsum, by.x = "ecod_id", by.y = "ID", all.x = TRUE)

# Clean up column names
names(out)[names(out) == names(out)[ncol(out)]] <- "hectares_with_data"

# Save as CSV (optional)
write.csv(out, "ecodistrict_hectares_with_data.csv", row.names = FALSE)

# If you also want the polygon layer with the attribute joined:
ed_out <- ed
ed_out$hectares_with_data <- out$hectares_with_data[match(ed_out$ecod_id, out$ecod_id)]
writeVector(ed_out, "ecodistricts_with_hectares.gpkg", filetype = "GPKG", overwrite = TRUE)

############################

#!/usr/bin/env Rscript

# ============================================================
# Compare multiple aligned binary (0/1) grassland rasters:
# 1) Pixelwise similarity (Jaccard/Sørensen) + PCoA + clustering
# 2) Optional agreement (Cohen's kappa)
# 3) Consensus + uncertainty maps
# 4) Landscape configuration metrics + PCA
# 5) PERMANOVA to test model/source/time effects
#
# Assumptions:
# - All rasters are perfectly aligned (same grid/extent/resolution)
# - Values are 0/1 with NA for NoData
#
# Datasets to include (cover full Prairies Ecozone extent):
# ESA, WRI, GLC, CGR, CEC, PUDL, Plowprint, PUDL, Comer, LUTS, ACI, CNVC, LCC, LCR, ECCC, NGI
#
# Outputs written to OUT_DIR
# ============================================================

# -----------------------
# User settings (EDIT ME)
# -----------------------
RASTER_DIR   <- "path/to/rasters"       # folder with rasters
RASTER_GLOB  <- "\\.tif$"               # pattern to match files
META_CSV     <- "path/to/metadata.csv"  # metadata file (see format below)
OUT_DIR      <- "outputs_compare_grass"
SAMPLE_N     <- 200000                  # set NULL to use all pixels (may be large). If rasters are large, keep between 50k-500k. For stable Jaccard patterns across 20ish rasters, 200k is often plenty.
SEED         <- 42

# Choose a pair of rasters for a detailed FP/FN/TP difference map (optional)
# Set to NULL to skip. Use raster names (without extension) or full filename match.
PAIR_MAP_A   <- NULL  # e.g., "raster_01"
PAIR_MAP_B   <- NULL  # e.g., "raster_02"

# Landscape metrics to compute (landscapemetrics naming)
LSM_WHAT <- c(
  "lsm_c_ca",    # total class area
  "lsm_c_pland", # % of landscape (grass proportion)
  "lsm_c_np",    # number of patches
  "lsm_c_lpi",   # largest patch index
  "lsm_c_ed",    # edge density
  "lsm_c_split", # splitting index
  "lsm_c_ai",    # aggregation index
  "lsm_c_cohesion" # cohesion (connectivity proxy)
)

# -----------------------
# Packages
# -----------------------
req_pkgs <- c(
  "terra", "vegan", "ggplot2", "pheatmap", "cluster",
  "landscapemetrics", "irr", "dplyr", "tidyr", "readr"
)

to_install <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

library(terra)
library(vegan)
library(ggplot2)
library(pheatmap)
library(cluster)
library(landscapemetrics)
library(irr)
library(dplyr)
library(tidyr)
library(readr)

set.seed(SEED)

# -----------------------
# Create output folder
# -----------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------
# Metadata format
# -----------------------
# metadata.csv must contain at least:
#   raster   : filename (with extension) or basename used to match
#   timestamp: parseable date or datetime string (e.g., "2020-07-15")
#   model    : model id/name
#   source   : data source id/name
#
# Example:
# raster,timestamp,model,source
# map1.tif,2019-06-01,RF,Sentinel2
# map2.tif,2020-06-01,UNet,NAIP

meta <- read_csv(META_CSV, show_col_types = FALSE) %>%
  mutate(
    timestamp = as.POSIXct(timestamp, tz = "UTC", tryFormats = c(
      "%Y-%m-%d %H:%M:%S", "%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"
    )),
    model  = as.factor(model),
    source = as.factor(source)
  )

# -----------------------
# Load rasters
# -----------------------
files <- list.files(RASTER_DIR, pattern = RASTER_GLOB, full.names = TRUE)
if (length(files) == 0) stop("No rasters found in RASTER_DIR with RASTER_GLOB pattern.")

# Match metadata to files
# Try matching by exact filename first; if not, match by basename.
files_df <- tibble(
  file = files,
  raster = basename(files),
  raster_nox = tools::file_path_sans_ext(basename(files))
)

meta2 <- meta %>%
  left_join(files_df, by = "raster") %>%
  mutate(file = ifelse(is.na(file),
                       files_df$file[match(raster, files_df$raster_nox)],
                       file)) %>%
  mutate(file = as.character(file))

if (any(is.na(meta2$file))) {
  missing <- meta2$raster[is.na(meta2$file)]
  stop("Could not match these metadata rasters to files: ", paste(missing, collapse = ", "))
}

# Keep in metadata order
files_use <- meta2$file

rasters <- lapply(files_use, terra::rast)

# -----------------------
# Alignment checks (basic)
# -----------------------
ref <- rasters[[1]]
ok <- sapply(rasters[-1], function(r) {
  all(terra::ext(r) == terra::ext(ref)) &&
    all(terra::res(r) == terra::res(ref)) &&
    terra::crs(r) == terra::crs(ref) &&
    terra::nrow(r) == terra::nrow(ref) &&
    terra::ncol(r) == terra::ncol(ref)
})
if (!all(ok)) stop("Not all rasters appear perfectly aligned (extent/res/crs/dim mismatch).")

# -----------------------
# Ensure 0/1 values (best-effort)
# -----------------------
# Convert to integer and clamp to {0,1} for safety
rasters <- lapply(rasters, function(r) {
  r <- round(r)
  r[r < 0] <- 0
  r[r > 1] <- 1
  r
})

names(rasters) <- meta2$raster

# ============================================================
# PART 1: Pixelwise similarity / distance matrices
# ============================================================

# Create a common NA mask across all rasters
# Only pixels that are non-NA for ALL rasters are used in comparisons.
stack_all <- rast(rasters)
common_ok <- app(stack_all, fun = function(...) {
  v <- c(...)
  as.integer(all(!is.na(v)))
})

# Sample pixels or take all
if (!is.null(SAMPLE_N)) {
  message("Sampling pixels for distance computations: ", SAMPLE_N)
  samp <- spatSample(stack_all, size = SAMPLE_N, method = "random",
                     na.rm = TRUE, as.df = TRUE, values = TRUE)
  # Remove rows with any NA (should already be removed, but keep safe)
  samp <- samp[complete.cases(samp), , drop = FALSE]
  X <- as.matrix(samp)
} else {
  message("Using all pixels (may be memory heavy).")
  vals <- values(stack_all)
  vals <- vals[complete.cases(vals), , drop = FALSE]
  X <- as.matrix(vals)
}

# X is n_pixels x n_rasters; we want objects=rasters as rows for vegdist
X_t <- t(X)

# Jaccard distance (presence/absence-focused)
dist_jac <- vegdist(X_t, method = "jaccard", binary = TRUE)

# Sørensen/Dice distance (implemented via Bray-Curtis on binary)
dist_sor <- vegdist(X_t, method = "bray", binary = TRUE)

# Save distance matrices
saveRDS(dist_jac, file.path(OUT_DIR, "dist_jaccard.rds"))
saveRDS(dist_sor, file.path(OUT_DIR, "dist_sorensen_dice.rds"))

# Heatmap of Jaccard distances
jac_mat <- as.matrix(dist_jac)
pheatmap(jac_mat,
         main = "Pairwise Jaccard Dissimilarity (rasters)",
         filename = file.path(OUT_DIR, "heatmap_jaccard.png"),
         width = 9, height = 8)

# Hierarchical clustering dendrogram
hc <- hclust(dist_jac, method = "average")
png(file.path(OUT_DIR, "dendrogram_jaccard.png"), width = 1200, height = 800)
plot(hc, main = "Hierarchical clustering (average) on Jaccard distances",
     xlab = "", sub = "")
dev.off()

# ============================================================
# PART 2: PCoA on distance matrix
# ============================================================
# cmdscale = classical MDS / PCoA
pcoa <- cmdscale(dist_jac, eig = TRUE, k = 4)
scores <- as.data.frame(pcoa$points)
colnames(scores) <- paste0("PCoA", 1:ncol(scores))

scores <- bind_cols(meta2, scores)

# Plot PCoA (color by source, shape by model)
p <- ggplot(scores, aes(PCoA1, PCoA2, color = source, shape = model)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 13) +
  labs(title = "PCoA of rasters using Jaccard dissimilarity",
       subtitle = "Color: data source | Shape: model",
       x = "PCoA1", y = "PCoA2")

ggsave(file.path(OUT_DIR, "pcoa_jaccard.png"), p, width = 9, height = 6, dpi = 200)

# Save PCoA scores
write_csv(scores, file.path(OUT_DIR, "pcoa_scores.csv"))

# ============================================================
# PART 3: Optional Cohen's kappa matrix (agreement-style)
# ============================================================

# Kappa can be slow on large pixel counts; use sample X (already possibly sampled).
# Compute pairwise kappa on sampled pixels.
nR <- nrow(X_t)
kappa_mat <- matrix(NA_real_, nR, nR, dimnames = list(meta2$raster, meta2$raster))

message("Computing Cohen's kappa matrix (pairwise)...")
for (i in seq_len(nR)) {
  for (j in seq_len(nR)) {
    # irr::kappa2 expects two columns with categorical ratings
    # If a raster has only one class in the sample, kappa2 may error; handle gracefully.
    v1 <- X_t[i, ]
    v2 <- X_t[j, ]
    out <- tryCatch({
      irr::kappa2(cbind(v1, v2))$value
    }, error = function(e) NA_real_)
    kappa_mat[i, j] <- out
  }
}

saveRDS(kappa_mat, file.path(OUT_DIR, "kappa_matrix.rds"))
write_csv(as.data.frame(kappa_mat), file.path(OUT_DIR, "kappa_matrix.csv"))

# Heatmap of kappa (agreement)
pheatmap(kappa_mat,
         main = "Pairwise Cohen's Kappa (agreement)",
         filename = file.path(OUT_DIR, "heatmap_kappa.png"),
         width = 9, height = 8)

# ============================================================
# PART 4: Consensus and uncertainty maps
# ============================================================

message("Building consensus and uncertainty maps...")
# consensus count: number of rasters classifying grass at each pixel
consensus_count <- app(stack_all, sum, na.rm = FALSE)

# proportion of rasters classifying grass
consensus_prop <- consensus_count / nlyr(stack_all)

# uncertainty via Bernoulli entropy: -p log p - (1-p) log(1-p)
# handle p=0 or p=1 safely
entropy_fun <- function(p) {
  eps <- 1e-12
  p <- pmin(pmax(p, eps), 1 - eps)
  -(p * log(p) + (1 - p) * log(1 - p))
}
uncertainty <- app(consensus_prop, entropy_fun)

writeRaster(consensus_count, file.path(OUT_DIR, "consensus_count.tif"), overwrite = TRUE)
writeRaster(consensus_prop,  file.path(OUT_DIR, "consensus_proportion.tif"), overwrite = TRUE)
writeRaster(uncertainty,     file.path(OUT_DIR, "uncertainty_entropy.tif"), overwrite = TRUE)

# Quick PNG maps
png(file.path(OUT_DIR, "map_consensus_proportion.png"), width = 1200, height = 900)
plot(consensus_prop, main = "Consensus proportion (fraction of rasters = grass)")
dev.off()

png(file.path(OUT_DIR, "map_uncertainty_entropy.png"), width = 1200, height = 900)
plot(uncertainty, main = "Uncertainty (entropy; high near p=0.5)")
dev.off()

# ============================================================
# PART 5: Pairwise difference map for a selected pair (optional)
# ============================================================

if (!is.null(PAIR_MAP_A) && !is.null(PAIR_MAP_B)) {
  # find layers by name match (either exact raster filename or no-extension)
  idxA <- which(meta2$raster == PAIR_MAP_A | tools::file_path_sans_ext(meta2$raster) == PAIR_MAP_A)
  idxB <- which(meta2$raster == PAIR_MAP_B | tools::file_path_sans_ext(meta2$raster) == PAIR_MAP_B)
  
  if (length(idxA) == 0 || length(idxB) == 0) {
    warning("PAIR_MAP_A/B did not match metadata rasters. Skipping pairwise map.")
  } else {
    rA <- rasters[[idxA[1]]]
    rB <- rasters[[idxB[1]]]
    
    # Encode pixelwise outcomes:
    # 0 = TN (0/0), 1 = FP (1/0), 2 = FN (0/1), 3 = TP (1/1)
    diff_map <- app(c(rA, rB), fun = function(a, b) {
      ifelse(is.na(a) | is.na(b), NA,
             ifelse(a == 0 & b == 0, 0,
                    ifelse(a == 1 & b == 0, 1,
                           ifelse(a == 0 & b == 1, 2, 3))))
    })
    
    out_name <- paste0("pairwise_diff_", PAIR_MAP_A, "_vs_", PAIR_MAP_B, ".tif")
    writeRaster(diff_map, file.path(OUT_DIR, out_name), overwrite = TRUE)
    
    png(file.path(OUT_DIR, paste0("pairwise_diff_", PAIR_MAP_A, "_vs_", PAIR_MAP_B, ".png")),
        width = 1200, height = 900)
    plot(diff_map, main = paste("Pairwise difference:", PAIR_MAP_A, "vs", PAIR_MAP_B))
    legend("bottomleft",
           legend = c("TN (0/0)", "FP (1/0)", "FN (0/1)", "TP (1/1)"),
           fill = terrain.colors(4), bty = "n")
    dev.off()
  }
}

# ============================================================
# PART 6: Landscape configuration metrics + PCA
# ============================================================

message("Computing landscape metrics per raster...")
# landscapemetrics expects categorical values; our rasters are 0/1 already.
# We'll compute CLASS-level metrics (class = 1), and use values.
# Some metrics need patches; ensure 8-neighbor connectivity if desired via directions=8.

metric_list <- lapply(seq_along(rasters), function(i) {
  r <- rasters[[i]]
  
  # landscapemetrics can compute a set of metrics; if some fail (rare class),
  # we keep NA.
  m <- tryCatch({
    calculate_lsm(r, what = LSM_WHAT, directions = 8)
  }, error = function(e) NULL)
  
  if (is.null(m)) {
    return(tibble(raster = meta2$raster[i], metric = LSM_WHAT, value = NA_real_))
  }
  
  # Keep class=1 entries when applicable; if metric is landscape-level, class may be NA
  m <- m %>%
    mutate(raster = meta2$raster[i]) %>%
    # Prefer class==1 where class exists
    group_by(metric) %>%
    arrange(is.na(class), desc(class == 1)) %>%
    slice(1) %>%
    ungroup() %>%
    select(raster, metric, value)
  
  m
})

metric_long <- bind_rows(metric_list)

metric_wide <- metric_long %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(match(raster, meta2$raster))

write_csv(metric_wide, file.path(OUT_DIR, "landscape_metrics.csv"))

# PCA on metrics (standardize)
metric_mat <- metric_wide %>%
  select(-raster) %>%
  as.data.frame()

# Handle all-NA or near-constant columns
keep_cols <- sapply(metric_mat, function(v) {
  ok <- sum(!is.na(v)) >= 3
  var_ok <- suppressWarnings(sd(v, na.rm = TRUE)) > 0
  ok && var_ok
})
metric_mat2 <- metric_mat[, keep_cols, drop = FALSE]

# Impute remaining NAs with column median (simple and robust)
for (j in seq_len(ncol(metric_mat2))) {
  v <- metric_mat2[[j]]
  if (anyNA(v)) metric_mat2[[j]][is.na(v)] <- median(v, na.rm = TRUE)
}

pca <- prcomp(metric_mat2, scale. = TRUE)

# Save PCA scores
pca_scores <- as.data.frame(pca$x)
pca_scores$raster <- metric_wide$raster
pca_scores <- left_join(meta2, pca_scores, by = "raster")
write_csv(pca_scores, file.path(OUT_DIR, "pca_landscape_metric_scores.csv"))

# PCA plot
p2 <- ggplot(pca_scores, aes(PC1, PC2, color = source, shape = model)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 13) +
  labs(title = "PCA of landscape configuration metrics",
       subtitle = "Color: data source | Shape: model",
       x = "PC1", y = "PC2")

ggsave(file.path(OUT_DIR, "pca_landscape_metrics.png"), p2, width = 9, height = 6, dpi = 200)

# ============================================================
# PART 7: PERMANOVA (adonis2) to test effects of model/source/time
# ============================================================

message("Running PERMANOVA on Jaccard distances...")

# PERMANOVA expects a data frame with matching rownames/order
meta_perm <- meta2 %>%
  mutate(
    time_num = as.numeric(timestamp) # numeric time for model term
  )

# Ensure distance order matches meta order
# vegdist order uses row order of X_t (rasters), which is meta2 order by construction.
perm <- adonis2(dist_jac ~ model + source + time_num, data = meta_perm)

# Save PERMANOVA result
sink(file.path(OUT_DIR, "permanova_jaccard.txt"))
print(perm)
sink()

message("Done. Outputs written to: ", normalizePath(OUT_DIR))

