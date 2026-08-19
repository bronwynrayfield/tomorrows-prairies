# =============================================================================
# 5a-visualize-spatial-predictors.R
# =============================================================================
#
# Purpose:
#   Produces maps of key spatial predictor variables for presentation
#   and publication. Reads processed predictor rasters and applies
#   consistent styling across all maps.
#
# Inputs:
#   - Processed predictor rasters from 5-spatial-predictors.R
#   - Study area boundary (sa)
#   - Ecoregion boundaries (optional overlay)
#   - Water binary layer (for masking lakes)
#
# Outputs:
#   - ./output/presentation/fig_pred_crop_suitability.png
#   - ./output/presentation/fig_pred_[others].png
#
# Notes:
#   - All maps use sa_r extent (Prairie ecozone only)
#   - Lakes masked using AAFC water class (code 31)
#   - Colour schemes chosen for accessibility and visual clarity
#   - Maps designed for 16:9 presentation format
#
# Dependencies:
#   - 0-setup.R
#   - 5-spatial-predictors.R outputs
#
# =============================================================================

library(ggplot2)
library(tidyterra)
library(scales)


source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "functions-general.R"))

# ---- Resolution setting ----
# Change to 30 for final run
target_res <- 300
agg_fact   <- target_res / 30

# ---- Directories ----
predictors_filepath <- "predictors"
indir <- file.path("./processeddatasets", predictors_filepath, paste0(target_res, "m"))

# ---- Load files ----
sa   <- st_read(file.path(sa_filepath, sa_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

# ---- Resample to target resolution (skip if target_res == 30) ----
if (agg_fact > 1) {
  sa_r <- aggregate(sa_r, fact = agg_fact, fun = "modal")
}  
  
# Water mask — for masking lakes
water2010_bin <- rast("./processeddatasets/predictors/300m/water2010_prop_3200m_300m.tif")

# Crop suitability
crop_suit <- rast("./processeddatasets/predictors/300m/gp_tillage_300m_filled.tif")

# Proportion annual cropland
prop_cropland <- rast("./processeddatasets/predictors/300m/crop2015_prop_3200m_300m.tif")

# Distance to roads
dist_to_road <- rast("./processeddatasets/predictors/30m/dist_to_road_30m.tif")


# ---- Prepare crop suitability ----

# Align & crop to study area
water_sa <- resample(water2010_bin, sa_r, method = "max")

# Crop and mask to study area, then mask out water
crop_suit_sa <- resample(crop_suit, sa_r, method = "bilinear")
crop_suit_sa <- mask(crop_suit_sa, sa_r)
crop_suit_sa <- ifel(water_sa == 1, NA, crop_suit_sa)

prop_cropland_sa <- resample(prop_cropland, sa_r, method = "bilinear")
prop_cropland_sa <- mask(prop_cropland_sa, sa_r)
prop_cropland_sa <- ifel(water_sa == 1, NA, prop_cropland_sa)

dist_to_road_sa <- resample(dist_to_road, sa_r, method = "bilinear")
dist_to_road_sa <- mask(dist_to_road_sa, sa_r)
dist_to_road_sa <- ifel(water_sa == 1, NA, dist_to_road_sa)


# ---- Plot ----

# Define major Prairie cities
cities <- data.frame(
  city = c("Calgary", "Edmonton", "Saskatoon", "Regina", "Winnipeg", 
           "Lethbridge", "Brandon"),
  lon  = c(-114.066, -113.323, -106.670, -104.618, -97.138, 
           -112.843, -99.952),
  lat  = c(51.045, 53.631, 52.132, 50.445, 49.899, 
           49.694, 49.840)
)

# Convert to sf in same CRS as your raster
cities_sf <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs(sa_r))


# Cropland suitability
# Convert to data frame for ggplot
crop_suit_df <- as.data.frame(crop_suit_sa, xy = TRUE) |>
  filter(!is.na(conv_prob_GP))  # adjust column name if needed

# Plot
map_crop_suit <- ggplot() +
  geom_raster(data = crop_suit_df, aes(x = x, y = y, fill = conv_prob_GP)) +
  geom_sf(data = sa, fill = NA, colour = "grey40", linewidth = 0.4) +
  geom_sf(data = cities_sf, colour = "white", size = 2, shape = 21,
          fill = "white") +
  geom_sf_label(data = cities_sf, aes(label = city),
                size = 2.5, fill = "white", alpha = 0.7,
                label.size = 0.1, nudge_y = 30000) +
  scale_fill_viridis_c(name = "Cropland\nsuitability",
                       option = "magma", direction = 1) +
  theme_void() +
  theme(legend.position = "right",
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

ggsave("./output/presentation/fig_pred_crop_suitability.png",
       plot = map_crop_suit, width = 10, height = 7, dpi = 300, bg = "white")


# Proportion cropland
# Convert to data frame
prop_cropland_df <- as.data.frame(prop_cropland_sa, xy = TRUE) |>
  filter(!is.na(focal_mean))

# Map
map_prop_cropland <- ggplot() +
  geom_raster(data = prop_cropland_df, aes(x = x, y = y, fill = focal_mean)) +
  geom_sf(data = sa, fill = NA, colour = "grey40", linewidth = 0.4) +
  geom_sf(data = cities_sf, colour = "white", size = 2, shape = 21,
          fill = "white") +
  geom_sf_label(data = cities_sf, aes(label = city),
                size = 2.5, fill = "white", alpha = 0.7,
                label.size = 0.1, nudge_y = 30000) +
  scale_fill_viridis_c(name = "Proportion\nannual cropland",
                       option = "magma", direction = 1) +
  theme_void() +
  theme(legend.position = "right",
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

ggsave("./output/presentation/fig_pred_prop_cropland.png",
       plot = map_prop_cropland, width = 10, height = 7, dpi = 300, bg = "white")


# Distance to road
# Convert to data frame
dist_to_road_df <- as.data.frame(dist_to_road_sa, xy = TRUE) |>
  filter(!is.na(layer))

# Map
map_dist_road <- ggplot() +
  geom_raster(data = dist_to_road_df, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = sa, fill = NA, colour = "grey40", linewidth = 0.4) +
  geom_sf(data = cities_sf, colour = "white", size = 2, shape = 21,
          fill = "white") +
  geom_sf_label(data = cities_sf, aes(label = city),
                size = 2.5, fill = "white", alpha = 0.7,
                label.size = 0.1, nudge_y = 30000) +
  scale_fill_viridis_c(
    name = "Distance\nto road",
    option = "magma",
    direction = -1,
    trans  = "sqrt",
    limits = c(0, quantile(dist_to_road_df$focal_mean, 0.95, na.rm = TRUE)),
    oob    = squish,                          # tail values → max colour, not NA
    breaks = c(0, 500, 1500, 10000),      # human-readable breaks in real units
    labels = scales::comma
  ) +
  theme_void() +
  theme(legend.position    = "right",
        plot.background    = element_rect(fill = "white", colour = NA),
        panel.background   = element_rect(fill = "white", colour = NA))

ggsave("./output/presentation/fig_pred_dist_to_road.png",
       plot = map_dist_road, width = 10, height = 7, dpi = 300, bg = "white")
