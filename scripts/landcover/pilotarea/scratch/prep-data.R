##############################################
# Data pre-processing including cropping, imputing, derived calculations

#--- Setup ----
source(file.path("scripts", "0-setup.R"))

#--- Directories ---
# Pilot area directory
pilot_dir <- "./processeddatasets/pilot_area"

#--- Load files ---
pilot <- st_read(file.path(pilot_dir, "pilot_area_test_aligned_cnvc.shp"))
pilot_rast_prairies <- rast(file.path(pilot_dir, "pilot_area_test_aligned_cnvc.tif"))

# Olimb & Robinson 2019 cropland suitability
# or read in the already reproected raster (see below)
#crop_suit_reproj<- rast(file.path("./processeddatasets/agriculture", "OlimbRobinson_GrassToGrain_EcologicalIndicators2019_datashare", "gp_tillage_reprojected.tif"))
crop_suit_raw <- rast(crop_suit_2019_file)

# Census subdivisions
census_raw <- st_read(census_subdivision_file)

# Roads
roads_raw <-st_read(road_file, layer = "NRN_SK_14_0_ROADSEG")

# LULC
# or read in the already reprojected raster (see below)
#lulc_prairies_reproj <- rast("./processeddatasets/landcover/aafc/cropped/AAFC2020_40kmbuffer_noUS_reproj.tif")
lulc_prairies_raw <- rast("./processeddatasets/landcover/aafc/cropped/AAFC2020_40kmbuffer_noUS.tif")

# ECCC grass loss
grass_loss_class_raw <- rast(grass_loss_class_file)
grass_loss_year_raw <- rast(grass_loss_year_file)  
ipast_loss_class_raw <- rast(ipast_loss_class_file)
ipast_loss_year_raw <- rast(ipast_loss_year_file)
upast_loss_class_raw <- rast(upast_loss_class_file)
upast_loss_year_raw <- rast(upast_loss_year_file)

#--- Reproject to match pilot shapefile CRS ---
crop_suit_reproj <- project(crop_suit_raw, pilot_rast_prairies, method = "bilinear")
lulc_prairies_reproj <- project(lulc_prairies_raw, pilot_rast_prairies, method = "near") #takes a long time, so saved on disk
census_reproj <- st_transform(census_raw, st_crs(pilot))
roads_reproj <- st_transform(roads_raw, st_crs(pilot))
grass_loss_class_reproj <- project(grass_loss_class_raw, pilot_rast_prairies, method = "near")
grass_loss_year_reproj <- project(grass_loss_year_raw, pilot_rast_prairies, method = "near")
ipast_loss_class_reproj <- project(ipast_loss_class_raw, pilot_rast_prairies, method = "near")
ipast_loss_year_reproj <- project(ipast_loss_year_raw, pilot_rast_prairies, method = "near")
upast_loss_class_reproj <- project(upast_loss_class_raw, pilot_rast_prairies, method = "near")
upast_loss_year_reproj <- project(upast_loss_year_raw, pilot_rast_prairies, method = "near")

#--- Write to disk ---
writeRaster(crop_suit_reproj, file.path("./processeddatasets/agriculture", "OlimbRobinson_GrassToGrain_EcologicalIndicators2019_datashare", "gp_tillage_reprojected.tif"))
writeRaster(lulc_prairies_reproj, "./processeddatasets/landcover/aafc/cropped/AAFC2020_40kmbuffer_noUS_reproj.tif")
writeRaster(census_reproj)
writeRaster(roads_reproj)
writeRaster(grass_loss_class_reproj, "./processeddatasets/landcover/MostRecentChangeByClass/2022/loss/c05_grass/recentchange_2022_loss_c05_grass_b01_class_reproj.tif")
writeRaster(grass_loss_year_reproj, "./processeddatasets/landcover/MostRecentChangeByClass/2022/loss/c05_grass/recentchange_2022_loss_c05_grass_b01_year_reproj.tif")
writeRaster(ipast_loss_class_reproj, "./processeddatasets/landcover/MostRecentChangeByClass/2022/loss/c10_ipast/recentchange_2022_loss_c10_ipast_b01_class_reproj.tif")
writeRaster(ipast_loss_year_reproj, "./processeddatasets/landcover/MostRecentChangeByClass/2022/loss/c10_ipast/recentchange_2022_loss_c10_ipast_b01_year_reproj.tif")
writeRaster(upast_loss_class_reproj, "./processeddatasets/landcover/MostRecentChangeByClass/2022/loss/c04_upast/recentchange_2022_loss_c04_upast_b01_class_reproj.tif")
writeRaster(upast_loss_year_reproj, "./processeddatasets/landcover/MostRecentChangeByClass/2022/loss/c04_upast/recentchange_2022_loss_c04_upast_b01_year_reproj.tif")

#--- Crop and mask to pilot area or padded pilot area ----
# Pilot raster template
pilot_rast <- crop(pilot_rast_prairies, pilot)
pilot_rast <- mask(pilot_rast, pilot)

# Cropland suitability
crop_suit_pilot <- crop(crop_suit_reproj, pilot)
crop_suit_pilot <- mask(crop_suit_pilot, pilot)

# Census subdivisions
census <- st_intersection(census_reproj, pilot)
#plot(census["CSDNAME"])
plot(crop_suit_pilot, main = "Cropland suitability (2019)")
plot(st_geometry(census),
     add = TRUE, col = NA, border = "black", lwd = 0.6)

grass_loss_class <- crop(grass_loss_class_reproj, pilot)
grass_loss_class <- mask(grass_loss_class, pilot)

grass_loss_year <- crop(grass_loss_year_reproj, pilot)
grass_loss_year <- mask(grass_loss_year, pilot)

#--- Pad the pilot area for use with some spatial predictors ----
# Pad the pilot extent by 20 cells
pad <- res(pilot_rast) * 120
pilot_v <- vect(st_transform(pilot, crs(pilot)))
pilot_pad <- buffer(pilot_v, width = pad)
pilot_sf <- st_as_sf(pilot)
pilot_pad_sf <- st_as_sf(pilot_pad)
# Plot both padded and un-padded together
plot(st_geometry(pilot_pad_sf), col = NA, border = "red", lwd = 2)
plot(st_geometry(pilot_sf), col = NA, border = "blue", lwd = 2, add = TRUE)


# Roads
# Crop using padded pilot area to avoid edge effects
roads <- st_intersection(roads_reproj, pilot_pad_sf)
plot(roads["ROADCLASS"], main = "Roads by Class", key.pos = 4)

# LULC
# Crop using padded pilot area to avoid edge effects
lulc_pilot <- crop(lulc_prairies_reproj, pilot_pad_sf)
lulc_pilot <- mask(lulc_pilot, pilot_pad_sf)

rm(crop_suit_reproj, census_reproj, roads_reproj, lulc_prairies_reproj) #cleanup

#--- Crop and mask to focal census subdivision (csd)----
# Create mask with focal csd
focal_csd <- census %>%
  filter(CSDNAME %in% c("Val Marie", "Val Marie No. 17"))

csd_mask <- rasterize(vect(focal_csd), pilot_rast,
                      field = 1, background = NA, touches = TRUE)

# Cropland suitability
crop_suit <- mask(crop_suit_pilot, csd_mask)
crop_suit <- crop(crop_suit, csd_mask)
plot(crop_suit, main = "Cropland suitability (2019)")

# # LULC
# lulc <- mask(lulc_pilot, csd_mask)
# lulc <- crop(lulc, csd_mask)
# plot(lulc, main = "LULC (2020)")



#--- Impute NA cells ----
# Cropland suitability
# Mean of all non-NA cells
mean_val <- global(crop_suit, fun = "mean", na.rm = TRUE)[[1]]

# Fill only where (a) the raster is NA and (b) we're inside the CSD mask
inside <- !is.na(csd_mask)                 
crop_suit_filled <- ifel(is.na(crop_suit) & inside, mean_val, crop_suit)
plot(crop_suit_filled, main = "Cropland suitability Imputed (2019)")

#--- Calculate distance to roads ----
# buffer by ~ half a cell (in raster units)
half_cell <- max(res(pilot_rast)) * 0.51
roads_buf <- st_buffer(roads, dist = half_cell)

roads_mask <- rasterize(vect(roads_buf), pilot_rast,
                        field = 1, background = NA, touches = TRUE)

dist_to_road_pad <- distance(roads_mask)

dist_to_road <- crop(dist_to_road_pad, pilot)
dist_to_road <- mask(dist_to_road, pilot)

dist_to_road <- mask(dist_to_road_pad, csd_mask)
dist_to_road <- crop(dist_to_road, csd_mask)
plot(dist_to_road)

#--- Calculate proportion of cropland within moving window ----
# Binary cropland raster
cropland <- ifel((lulc_pilot == 51) | (lulc_pilot == 52), 1, 0)
plot(cropland)

# window size in metres
win_size_m <- 3600

# convert to number of cells
win_size_cells <- round(win_size_m / res(cropland)[1]) + 1

# create a square (3.6 km × 3.6 km) window
w <- matrix(1, nrow = win_size_cells, ncol = win_size_cells)

# focal mean gives proportion of 1s
cropland_prop_pad <- focal(cropland, w = w, fun = "mean", na.policy = "omit", na.rm = TRUE)

# visualize
plot(cropland_prop_pad, main = "Proportion of cropland within 3.6 km radius")


#Loop over to test sensitivity of window size
# your binary cropland raster (0/1)
r <- cropland  

# window sizes in metres
win_sizes <- c(800, 1600, 3200, 3600, 7200)

# prepare an empty list to store results
prop_maps <- list()

for (w_m in win_sizes) {
  cat("Running window:", w_m, "m\n")
  
  # convert to number of cells
  w_cells <- round(w_m / res(r)[1])
  
  # force odd size (needed by terra::focal)
  if (w_cells %% 2 == 0) w_cells <- w_cells + 1
  
  # build square window
  w <- matrix(1, nrow = w_cells, ncol = w_cells)
  
  # focal mean → proportion of 1s in window
  prop_maps[[as.character(w_m)]] <- focal(r, w = w, fun = "mean", na.policy = "omit", na.rm = TRUE)
}

# --- Compute pairwise Pearson correlations among window sizes ---
vals <- lapply(prop_maps, function(x) values(x, na.rm = TRUE))

common_na <- Reduce("|", lapply(vals, is.na))
vals_clean <- lapply(vals, function(v) v[!common_na])

corr_mat <- matrix(NA, nrow = length(win_sizes), ncol = length(win_sizes),
                   dimnames = list(win_sizes, win_sizes))

for (i in seq_along(win_sizes)) {
  for (j in seq_along(win_sizes)) {
    corr_mat[i,j] <- cor(vals_clean[[i]], vals_clean[[j]], method = "spearman")
  }
}

print(round(corr_mat, 3))

library(reshape2)
library(ggplot2)

corr_df <- melt(corr_mat, varnames = c("Window1", "Window2"), value.name = "Correlation")

mean_corr <- apply(corr_mat, 1, function(x) mean(x[-which(names(x) == names(x)[which.max(x)])]))
barplot(mean_corr,
        names.arg = win_sizes,
        col = "grey",
        main = "Average correlation with other window sizes",
        ylab = "Mean correlation",
        ylim = c(0, 1))



# Visualize with the moving window size shown
# find the center of the raster
r <- cropland_prop_pad    
e <- ext(r)
cx <- (e$xmin + e$xmax) / 2
cy <- (e$ymin + e$ymax) / 2

# half the side length (3.6 km window)
half_side <- 3600 / 2

# build a square extent and turn it into a polygon
rect_ext <- ext(cx - half_side, cx + half_side, cy - half_side, cy + half_side)
square_poly <- as.polygons(rect_ext, crs = crs(r))   # SpatVector polygon

# plot on your cropland (or cropland_prop)
plot(cropland, col = c("grey85","goldenrod"), legend = FALSE,
     main = "Cropland with 3.6 × 3.6 km Window")
lines(square_poly, col = "red", lwd = 2)


cropland_prop <- mask(cropland_prop_pad, pilot)
cropland_prop <- crop(cropland_prop, pilot)

cropland_prop <- mask(cropland_prop, csd_mask)
cropland_prop <- crop(cropland_prop, csd_mask)
plot(cropland_prop)

cropland_focal <- mask(cropland, pilot)
cropland_focal <- crop(cropland_focal, pilot)

cropland_focal <- mask(cropland_focal, csd_mask)
cropland_focal <- crop(cropland_focal, csd_mask)
plot(cropland_focal)

#--- Write rasters to disk ----
# Cropland suitability
writeRaster(crop_suit, file.path(pilot_dir, "cropland_suitability_2019_withNAs.tif"))
writeRaster(crop_suit_filled, file.path(pilot_dir, "cropland_suitability_2019_noNAs.tif"))

#Census subdivisions

#roads
writeRaster(dist_to_road, file.path(pilot_dir, ""))

#lulc
writeRaster(lulc_pilot, file.path(pilot_dir, "lulc_aafc_2020.tif"))

#protected areas
#binary crop x2
