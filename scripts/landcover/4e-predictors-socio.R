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

source(file.path("scripts/landcover", "0a-setup.R"))

# ---- Load files ----
# Modelling extent
me_r <- rast(file.path(sa_filepath, me_raster_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))

# Ecodistricts
ecodistricts <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath,
                                    ecodistricts_filename, "modelling-extent"))
# Canfax results
canfax2009 <- st_read(tagged_path("./processeddatasets", canfax_filepath, 
                                  canfax2009_filename, "modelling-extent"))
canfax2024 <- st_read(tagged_path("./processeddatasets", canfax_filepath, 
                                  canfax2024_filename, "modelling-extent"))

# Process elevators
elevator <- st_read(tagged_path("./processeddatasets", elevator_filepath, 
                              elevator_filename, "modelling-extent"))
# Roads
roads <- rast(tagged_path("./processeddatasets", road_filepath,
                              paste0(road_filename, ".tif"), "modelling-extent"))
# Rails
rails <- rast(tagged_path("./processeddatasets", rail_filepath,
                          paste0(rail_filename, ".tif"), "modelling-extent"))

# Population centres
popctr2011 <- st_read(tagged_path("./processeddatasets", popctr_filepath, 
                              paste0(popctr2011_filename, ".shp"), "modelling-extent"))
popctr2021 <- st_read(tagged_path("./processeddatasets", popctr_filepath, 
                                  paste0(popctr2021_filename, ".shp"), "modelling-extent"))

# ---- Rasterize Canfax results ----
# 2005 - 2009
canfax2009_v <- vect(canfax2009)
canfax2009_r <- rasterize(canfax2009_v, me_r, field = "nt_b_25")
write_rast_tagged(canfax2009_r, "./processeddatasets", pred_dir, "canfax_nt_b_25_2005-2009.tif", paste0(res_m, "m"))
# 2020 - 2024
canfax2024_v <- vect(canfax2024)
canfax2024_r <- rasterize(canfax2024_v, me_r, field = "nt_b_25")
write_rast_tagged(canfax2024_r, "./processeddatasets", pred_dir, "canfax_nt_b_25_2020-2024.tif", paste0(res_m, "m"))


# ---- Distance to ----
# Process elevators
elev_d <- distance(me_r, vect(elevator))   # metres, one value per cell
elev_d <- mask(elev_d, me_r)
plot(elev_d / 1000, main = "Distance to process elevator (km)")
points(vect(elevator), pch = 21, bg = "white", col = "black", lwd = 0.8,
       cex = 0.4 + 2 * sqrt(elevator$CAP_T / max(elevator$CAP_T, na.rm = TRUE)))

# Roads
roads_for_dist <- ifel(roads == 1, 1, NA)
dist_to_road  <- distance(roads_for_dist)
dist_to_road  <- mask(dist_to_road, me_r)
plot(dist_to_road / 1000, main = "Distance to roads (km)")

# Rails
rails_for_dist <- ifel(rails == 1, 1, NA)
dist_to_rail  <- distance(rails_for_dist)
dist_to_rail  <- mask(dist_to_rail, me_r)
plot(dist_to_rail / 1000, main = "Distance to railways (km)")

# ---- Population centre distance predictors -------------------------------
# POPCTR (>=1,000 people, >=400/km2), 2011: 175 polygons in extent.
# PCTYPE: 1 core (28), 2 fringe (30), 4 outside any CMA/CA (112), 6 secondary core (5).
# PCCLASS by population: 2 small (1,000-29,999), 3 medium (30,000-99,999), 4 large (100,000+).

# Subset population centres that serve as the core or secondary core of a 
# census metropolitan area (CMA) or census agglomeration (CA).
# A CMA needs a total population of at least 100,000 with 50,000 or more in the core
# A CA needs a core of at least 10,000.
pc_cma_ca_cores2011 <- popctr2011[popctr2011$PCTYPE %in% c(1, 6), ]
pc_cma_ca_cores2021 <- popctr2021[popctr2021$PCTYPE %in% c(1, 6), ]

# Distance to nearest population centre and distance to nearest regional centre
# 2011
pc_d_all_2011 <- distance(me_r, vect(popctr2011), rasterize = TRUE, unit = "m")
pc_d_all_2011 <- mask(pc_d_all_2011, me_r)
plot(pc_d_all_2011/1000, main = "Distance to nearest population centre (km)")
plot(st_geometry(popctr2011), border = "white", lwd = 0.4, add = TRUE)

pc_d_core_2011 <- distance(me_r, vect(pc_cma_ca_cores2011), rasterize = TRUE, unit = "m")
pc_d_core_2011 <- mask(pc_d_core_2011, me_r)
plot(pc_d_core_2011/1000, main = "Distance to nearest regional centre (km)")
plot(st_geometry(pc_cma_ca_cores2011), border = "white", lwd = 0.4, add = TRUE)

# 2021
pc_d_all_2021 <- distance(me_r, vect(popctr2021), rasterize = TRUE, unit = "m")
pc_d_all_2021 <- mask(pc_d_all_2021, me_r)
plot(pc_d_all_2021/1000, main = "Distance to nearest population centre (km)")
plot(st_geometry(popctr2021), border = "white", lwd = 0.4, add = TRUE)

pc_d_core_2021 <- distance(me_r, vect(pc_cma_ca_cores2021), rasterize = TRUE, unit = "m")
pc_d_core_2021 <- mask(pc_d_core_2021, me_r)
plot(pc_d_core_2021/1000, main = "Distance to nearest regional centre (km)")
plot(st_geometry(pc_cma_ca_cores2021), border = "white", lwd = 0.4, add = TRUE)


# Save
write_rast_tagged(elev_d, "./processeddatasets", elevator_filepath, "dist_to_elevator.tif", "modelling-extent")
write_rast_tagged(dist_to_road,  "./processeddatasets", road_filepath, "dist_to_road.tif", "modelling-extent")
write_rast_tagged(dist_to_rail,  "./processeddatasets", rail_filepath, "dist_to_rail.tif", "modelling-extent")
write_rast_tagged(pc_d_all_2011,  "./processeddatasets", popctr_filepath, "dist_to_popctr2011.tif", "modelling-extent")
write_rast_tagged(pc_d_all_2021,  "./processeddatasets", popctr_filepath, "dist_to_popctr2021.tif", "modelling-extent")
write_rast_tagged(pc_d_core_2011,  "./processeddatasets", popctr_filepath, "dist_to_regctr2011.tif", "modelling-extent")
write_rast_tagged(pc_d_core_2021,  "./processeddatasets", popctr_filepath, "dist_to_regctr2021.tif", "modelling-extent")


# Read files back in and check for correlation and ranges
f <- file.path(pred_dir, paste0("dist_to_", c("popctr", "regctr"), 
                                rep(c(2011, 2021), each = 2), "_30m.tif"))
pc <- rast(f)
names(pc) <- c("popctr2011", "regctr2011", "popctr2021", "regctr2021")

# 100k regular sample - full vectors are 1.1e9 elements
s <- spatSample(pc, 1e5, method = "regular", na.rm = TRUE, as.df = TRUE)
round(cor(s, method = "spearman"), 3)
#            popctr2011 regctr2011 popctr2021 regctr2021
# popctr2011      1.000      0.537      0.944      0.536
# regctr2011      0.537      1.000      0.538      0.972
# popctr2021      0.944      0.538      1.000      0.541
# regctr2021      0.536      0.972      0.541      1.000

# range: 2021 values outside the 2011 range could be a problem for RF model
sapply(s, range)
# popctr2011 regctr2011 popctr2021 regctr2021
# [1,]        0.0        0.0        0.0        0.0
# [2,]   175968.3   290293.9   175215.1   291176.6







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



