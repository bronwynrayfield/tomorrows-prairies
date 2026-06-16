# =============================================================================
# 6-response.R
# =============================================================================
#
# Purpose:
#   Derives the binary response variable for the Random Forest conversion
#   risk model: grassland pixels in 2010 that converted to annual cropland
#   by 2020. Also summarises conversion by ecoregion and produces a map.
#
# Inputs:
#   - AAFC LUTS 2010 and 2020 (cropped to modelling extent)
#   - Modelling extent raster (me_r) : defines final mask
#   - Study area shapefile (sa) : for map intersection
#   - Ecozones and ecoregions (cropped to modelling extent)
#
# Outputs:
#   - response_30m.tif : binary raster (1 = converted, 0 = not converted)
#   - response_{target_res}m.tif : binary raster (1 = converted, 0 = not converted)
#   - grassland_conversion_to_cropland_map.png : map of conversion events by ecoregion
#
# Notes:
#   - Response defined as: native grassland in 2010 AND annual cropland in 2020
#   - Native grassland: codes 63-68
#   - Annual cropland (2020 only): codes 52, 56
#   - Total conversion events at 30m: ~200,946 pixels (~18,085 ha)
#   - target_res must match setting in 5-spatial-predictors.R
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 3-crop-data.R outputs (AAFC LUTS 2010, 2020)
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "functions-general.R"))
source(file.path("scripts/landcover/functions", "reclass-aafc-luts.R"))

# ---- Resolution setting ----
# Must match setting in 5-spatial-predictors.R
target_res <- 300
agg_fact   <- target_res / 30

# ---- Load files ----
me_r <- rast(file.path(sa_filepath, me_raster_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
sa   <- st_read(file.path(sa_filepath, sa_filename))

lulc2010 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent"))
lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))

ecozones   <- st_read(tagged_path("./processeddatasets", ecozones_filepath,   ecozones_filename,   "modelling-extent"))
ecoregions <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "modelling-extent"))

# ---- Mask to Prairie ecozone study area ----
lulc2010 <- mask(crop(lulc2010, sa_r), sa_r)
lulc2020 <- mask(crop(lulc2020, sa_r), sa_r)

# ---- Response variable ----
was_native <- reclass_aafc(lulc2010, scheme = "binary", binary_class = "native", check = TRUE)
is_crop    <- reclass_aafc(lulc2020, scheme = "binary", binary_class = "crop",   check = TRUE)

response <- ifel(was_native == 1 & is_crop == 1, 1, 0)

# QA/QC
freq(response)
#   layer value     count
# 1     1     0 518363681
# 2     1     1    200946

# ---- Summarise conversion by ecoregion ----
# Filter to Prairie ecozone only (ECOZONE_ID = 10)
ecoregions_prairie <- ecoregions |>
  filter(ECOZONE_ID == 10)

ecoregions_r <- rasterize(vect(ecoregions_prairie), sa_r, field = "ECOREGION1")

conversion_by_ecoregion <- zonal(response, ecoregions_r, fun = "sum", na.rm = TRUE)
names(conversion_by_ecoregion) <- c("ecoregion", "n_conversions")
conversion_by_ecoregion <- conversion_by_ecoregion[order(conversion_by_ecoregion$n_conversions), ]

# Summary
summary(conversion_by_ecoregion$n_conversions)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 107    1927    9468   28705   39413  108679 
print(conversion_by_ecoregion)
# ecoregion n_conversions
# 7 Southwest Manitoba Uplands           107
# 4        Lake Manitoba Plain          1423
# 2             Cypress Upland          2431
# 3           Fescue Grassland          9468
# 6      Moist Mixed Grassland         35095
# 5            Mixed Grassland         43730
# 1             Aspen Parkland        108679

# How many ecoregions have fewer than 100 conversion events?
sum(conversion_by_ecoregion$n_conversions < 100, na.rm = TRUE)
# 0
# How many have zero?
sum(conversion_by_ecoregion$n_conversions == 0, na.rm = TRUE)
# 0

# ---- Map ----
# Aggregate to 3x target resolution for plotting efficiency
response_plot <- aggregate(response, fact = 20, fun = "max")
response_df   <- as.data.frame(response_plot, xy = TRUE) |>
  filter(LU2010 == 1)

ecoregions_sf <- st_as_sf(ecoregions) |> st_intersection(sa)
ecozone_sf    <- st_as_sf(ecozones)   |> st_intersection(sa)

ecoregions_label <- ecoregions_sf |>
  left_join(conversion_by_ecoregion, by = c("ECOREGION1" = "ecoregion")) |>
  filter(n_conversions >= 1000)

map <- ggplot() +
  geom_tile(data = response_df, aes(x = x, y = y),
            fill = "darkred", alpha = 0.8) +
  geom_sf(data = ecoregions_sf, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_sf(data = ecozone_sf,    fill = NA, colour = "black",  linewidth = 0.7) +
  geom_sf_label(data = ecoregions_label, aes(label = ECOREGION1),
                size = 2.5, fill = "white", alpha = 0.7, label.size = 0.1) +
  theme_void()

ggsave(file.path("./output/landcover/conversion-response", "grassland_conversion_to_cropland_map.png"),
       plot = map, width = 10, height = 7, dpi = 300, bg = "white")


# ---- Resample to target resolution and save ----
# Always save 30m response
write_rast_tagged(response, "./processeddatasets", "response", "response.tif",
                  "30m", datatype = "INT1U")

# Aggregate and save at target resolution if different from 30m
if (agg_fact > 1) {
  sa_r_target     <- aggregate(sa_r, fact = agg_fact, fun = "modal")
  response_target <- resample(response, sa_r_target, method = "max")
  write_rast_tagged(response_target, "./processeddatasets", "response", "response.tif",
                    paste0(target_res, "m"), datatype = "INT1U")
}
