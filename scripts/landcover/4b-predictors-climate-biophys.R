# =============================================================================
# 5-spatial-predictors.R
# =============================================================================
#
# Purpose:
#   Derives spatial predictor variables at a target resolution for use in
#   the Random Forest grassland conversion risk model. Includes landscape context
#   variables (moving window proportions) and cropland suitability imputation.
#   Change target_res to 30 for final model run.
#
# Inputs:
#   - Modelling extent raster template (me_r)
#   - AAFC LUTS binary rasters (2010, 2015) from 4-prepare-spatial-inputs.R
#   - Protected areas binary rasters from 4-prepare-spatial-inputs.R
#   - Olimb & Robinson 2019 cropland suitability
#   - Ecodistricts (for cropland suitability imputation)
#
# Outputs:
#   - Landscape context proportion rasters at target resolution
#     saved to processeddatasets/predictors/{target_res}m/
#   - Cropland suitability (gap-filled) at target resolution
#
# Notes:
#   - Window sizes defined in win_sizes_m
#   - target_res = 300 for prototyping, change to 30 for final run
#   - At 300m, AAFC binary layers become proportions after resampling;
#     results are for prototyping only and not directly comparable to
#     30m final run
#   - pa_hard threshold at >0.5 : cell treated as protected if majority
#     of underlying 30m pixels are protected
#   - Cropland suitability NA imputation uses ecodistrict means
#   - crop class not present in 2010 AAFC LUTS : crop2010 will be
#     all zeros, included for completeness
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 4-prepare-spatial-inputs.R outputs
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "reclass-aafc-luts.R"))

# ---- Resolution setting ----
# Change to 30 for final run
target_res <- 300
agg_fact   <- target_res / 30

# ---- Window sizes ----
# Comment out sizes to exclude from this run
win_sizes_m <- c(
  3200#,
  # 6000,
  # 10000
)

# ---- Output directory ----
predictors_filepath <- "predictors"
outdir <- file.path("./processeddatasets", predictors_filepath, paste0(target_res, "m"))

# ---- Load files ----

# Modelling extent
me_r <- rast(file.path(sa_filepath, me_raster_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

# Ecodistricts (for cropland suitability imputation)
ecodistricts <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath,
                                    ecodistricts_filename, "modelling-extent"))

# Binary AAFC filenames = cleaner than using the merged filename
aafc2015_binary_filename <- "AAFC2015.tif"
aafc2010_binary_filename <- "AAFC2010.tif"

# Binary filepath
aafc_binary_filepath <- file.path(aafc_filepath, "binary")

# AAFC LUTS binary rasters — 2015
native2015_bin     <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "native"))
tame2015_bin       <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "tame"))
crop2015_bin       <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "crop"))
settlement2015_bin <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "settlement"))
forest2015_bin     <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "forest"))
water2015_bin      <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2015_binary_filename, "water"))

# AAFC LUTS binary rasters — 2010
native2010_bin     <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "native"))
tame2010_bin       <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "tame"))
settlement2010_bin <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "settlement"))
forest2010_bin     <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "forest"))
water2010_bin      <- rast(tagged_path("./processeddatasets", aafc_binary_filepath, aafc2010_binary_filename, "water"))

# Protected areas binary rasters
pa_hard_bin        <- rast(tagged_path("./processeddatasets", pa_filepath, pa_filename, "hard", "binary"))
pa_provpasture_bin <- rast(tagged_path("./processeddatasets", pa_filepath, pa_filename, "provpasture", "binary"))

# Cropland suitability
crop_suit <- rast(tagged_path("./processeddatasets", crop_suit_2019_filepath,
                              crop_suit_2019_filename, "modelling-extent"))

# Roads
roads <- rast(tagged_path("./processeddatasets", road_filepath,
                              road_filename, "modelling-extent"))
roads <- rast("./processeddatasets/landcover/Roads - Poley et al 2021/roads_modelling-extent.tif")
roads <- resample(roads, me_r, method = "bilinear")

# ---- Resample to target resolution (skip if target_res == 30) ----
if (agg_fact > 1) {
  sa_r <- aggregate(sa_r, fact = agg_fact, fun = "modal")
  me_r <- aggregate(me_r, fact = agg_fact, fun = "modal")
  
  # AAFC binary layers — average only, no thresholding
  # Note: at 300m values are proportions (0-1), not binary
  # Results at 300m are for prototyping only and not directly
  # comparable to 30m final run
  native2015_bin     <- resample(native2015_bin,     me_r, method = "average")
  tame2015_bin       <- resample(tame2015_bin,       me_r, method = "average")
  crop2015_bin       <- resample(crop2015_bin,       me_r, method = "average")
  settlement2015_bin <- resample(settlement2015_bin, me_r, method = "average")
  forest2015_bin     <- resample(forest2015_bin,     me_r, method = "average")
  water2015_bin      <- resample(water2015_bin,      me_r, method = "average")
  native2010_bin     <- resample(native2010_bin,     me_r, method = "average")
  tame2010_bin       <- resample(tame2010_bin,       me_r, method = "average")
  settlement2010_bin <- resample(settlement2010_bin, me_r, method = "average")
  forest2010_bin     <- resample(forest2010_bin,     me_r, method = "average")
  water2010_bin      <- resample(water2010_bin,      me_r, method = "average")
  crop_suit          <- resample(crop_suit,          me_r, method = "bilinear")
  
  # PA rasters — average then threshold at >0.5
  # Cell is treated as protected if majority of underlying 30m pixels are protected
  pa_hard_bin        <- resample(pa_hard_bin,        me_r, method = "average")
  pa_hard_bin        <- ifel(pa_hard_bin > 0.5, 1, 0)
  pa_provpasture_bin <- resample(pa_provpasture_bin, me_r, method = "average")
  pa_provpasture_bin <- ifel(pa_provpasture_bin > 0.5, 1, 0)
}

# ---- Apply pa_hard mask to native and tame grass ----
# Required step to map unprotected grasslands
native2015_unprot <- ifel(pa_hard_bin == 1, NA, native2015_bin)
tame2015_unprot   <- ifel(pa_hard_bin == 1, NA, tame2015_bin)
native2010_unprot <- ifel(pa_hard_bin == 1, NA, native2010_bin)

# ---- Cropland suitability imputation ----
# Fill NAs in crop_suit with ecodistrict means
ecodistrict_zones <- rasterize(ecodistricts, me_r, field = "ECODISTRIC")
eco_means         <- zonal(crop_suit, ecodistrict_zones, fun = "mean", na.rm = TRUE)
names(eco_means)  <- c("ECODISTRIC", "mean_suit")
ecodistrict_means <- classify(ecodistrict_zones, eco_means)
crop_suit_filled  <- ifel(is.na(crop_suit) & !is.na(me_r), ecodistrict_means, crop_suit)

# QA/QC
freq(is.na(crop_suit_filled))

# Save
write_rast_tagged(crop_suit_filled, outdir, crop_suit_2019_filename, paste0(target_res, "m"), "filled")

# Distance to roads at 30m
# ---- Distance to cropland edge ----
roads_for_dist <- ifel(roads == 1, 1, NA)
dist_to_road  <- distance(roads_for_dist)
dist_to_road  <- mask(dist_to_road, me_r)

# Save
writeRaster(dist_to_road, "./processeddatasets/predictors/30m/dist_to_road_30m.tif")

# ---- Moving window focal statistics ----
# Define layers and names for focal loop
focal_layers <- list(
  native2015_unprot  = native2015_unprot,
  tame2015_unprot    = tame2015_unprot,
  crop2015           = crop2015_bin,
  settlement2015     = settlement2015_bin,
  forest2015         = forest2015_bin,
  water2015          = water2015_bin,
  native2010_unprot  = native2010_unprot,
  tame2010           = tame2010_bin,
  settlement2010     = settlement2010_bin,
  forest2010         = forest2010_bin,
  water2010          = water2010_bin,
  provpasture        = pa_provpasture_bin
)

# Loop over window sizes and layers
terraOptions(progress = 1, memfrac = 0.8)

for (win_m in win_sizes_m) {
  message("Processing window size: ", win_m, "m")
  win_cells <- round(win_m / res(me_r)[1])
  w         <- matrix(1, nrow = win_cells, ncol = win_cells)
  
  for (layer_name in names(focal_layers)) {
    message("  Layer: ", layer_name)
    out_filename <- paste0(layer_name, "_prop_", win_m, "m_", target_res, "m.tif")
    out_path     <- file.path(outdir, out_filename)
    
    focal(focal_layers[[layer_name]],
          w         = w,
          fun       = "mean",
          na.policy = "omit",
          na.rm     = TRUE,
          filename  = out_path,
          overwrite = TRUE)
  }
}

# ---- Distance to cropland edge ----
crop2015_for_dist <- ifel(crop2015_bin == 1, 1, NA)
dist_to_crop2015  <- distance(crop2015_for_dist)
dist_to_crop2015  <- mask(dist_to_crop2015, me_r)

# Always save 30m version
writeRaster(dist_to_crop2015,
            file.path("./processeddatasets", predictors_filepath, "30m", "dist_to_crop2015_30m.tif"),
            overwrite = TRUE,
            wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES")))

# Resample and save at target resolution if different from 30m
if (agg_fact > 1) {
  dist_to_crop2015_target <- resample(dist_to_crop2015, me_r, method = "bilinear")
  writeRaster(dist_to_crop2015_target,
              file.path(outdir, paste0("dist_to_crop2015_", target_res, "m.tif")),
              overwrite = TRUE,
              wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES")))
} 


# ---- Resample predictors to study area grid and save ----
outdir_sa <- file.path("./processeddatasets", predictors_filepath,
                       paste0(target_res, "m"), "sa")

if (agg_fact > 1) {
  sa_r_target <- aggregate(sa_r, fact = agg_fact, fun = "modal")
} else {
  sa_r_target <- sa_r
}

pred_files <- list.files(outdir, pattern = "\\.tif$", full.names = TRUE)
for (f in pred_files) {
  r_sa <- resample(rast(f), sa_r_target, method = "average")
  r_sa <- mask(r_sa, sa_r_target)
  writeRaster(r_sa, file.path(outdir_sa, basename(f)),
              overwrite = TRUE,
              wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES")))
}

# ---- QA/QC visualization ----
# ---- Helper function for QA/QC visualization ----
plot_predictors <- function(dir, title = "") {
  files <- list.files(dir, pattern = "\\.tif$", full.names = TRUE)
  stack <- rast(files)
  names(stack) <- tools::file_path_sans_ext(basename(files))
  message(title, ": ", nlyr(stack), " layers")
  global(stack, fun = c("min", "max", "mean"), na.rm = TRUE) |> print()
  plot(stack,
       col  = hcl.colors(100, "viridis"),
       axes = FALSE,
       mar  = c(0.5, 0.5, 1.5, 2),
       main = names(stack))
  invisible(stack)
}

pred_stack_me <- plot_predictors(outdir,    title = "Modelling extent")
pred_stack_sa <- plot_predictors(outdir_sa, title = "Study area")


















#--- Setup ----
source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "0-functions-general.R"))

# Resolution setting
target_res <- 30
agg_fact   <- target_res / 30

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Land cover (AAFC LUTS)
lulc2010 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent"))
lulc2015 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent"))
lulc2020 <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))

# Protected areas
pa_hard <- rast(tagged_path("./processeddatasets", pa_filepath, pa_filename, "modelling-extent", "hard", "binary"))
pa_provpasture <- rast(tagged_path("./processeddatasets", pa_filepath, pa_filename, "modelling-extent", "provpasture", "binary"))

# Olimb & Robinson 2019 cropland suitability
crop_suit <- rast(tagged_path("./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "modelling-extent"))

# ---- Resample to target resolution (skip if target_res == 30) ----
if (agg_fact > 1) {
  me_r      <- aggregate(me_r, fact = agg_fact, fun = "modal")
  lulc2010  <- resample(lulc2010, me_r, method = "near")
  lulc2015  <- resample(lulc2015, me_r, method = "near")
  lulc2020  <- resample(lulc2020, me_r, method = "near")
  crop_suit <- resample(crop_suit, me_r, method = "bilinear")
  # For binary pa rasters, average then threshold; mask if >50% protected
  pa_hard   <- resample(pa_hard, me_r, method = "average")
  pa_hard   <- ifel(pa_hard > 0.5, 1, 0)
  pa_provpasture <- resample(pa_provpasture, me_r, method = "average")
  pa_provpasture <- ifel(pa_provpasture > 0.5, 1, 0)
}



