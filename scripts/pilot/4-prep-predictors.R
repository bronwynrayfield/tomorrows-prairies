##############################################
# Prepare predictor rasters for cropland suitability modelling
# Fill in NA's, dist to roads, proportion of cropland in moving window

#--- Setup ----
source(file.path("./scripts", "0-setup.R"))
source(file.path("./scripts/pilot", "functions.R"))

# Load files 
pilot <- st_read(file.path(pilot_filepath, "GNP_pilot_area_crop_CSD.shp"))
pilot_r <- rast(file.path(pilot_filepath, "GNP_pilot_area_crop_CSD.tif"))
lulc_pilot_pad <- rast(file.path(pilot_filepath, "AAFC2020_40kmbuffer_noUS_pilot_pad.tif"))
gnp <- st_read(file.path(pilot_filepath, "gnp_park_pilot.shp"))
crop_suit <- rast(file.path(pilot_filepath, paste0(crop_suit_2019_filename, ".tif")))
roads <- st_read(file.path(pilot_filepath, "NRN_SK_14_0_GPKG_en_pilot_pad.gpkg"), layer = "NRN_SK_14_0_ROADSEG")

# For plotting later, get centre coordinates of gnp
gnp_cent <- sf::st_coordinates(sf::st_centroid(gnp))

#--- Cropland suitability - impute NA cells ----
# Cropland suitability
# Mean of all non-NA cells
mean_val <- global(crop_suit, fun = "mean", na.rm = TRUE)[[1]]

# Fill only where (a) the raster is NA and (b) we're inside the CSD mask
inside <- !is.na(pilot_r)                 
crop_suit_filled <- ifel(is.na(crop_suit) & inside, mean_val, crop_suit)
plot(crop_suit_filled)


# Save to disk 
png("./results/landcover/pilot/crop_suitability.png", width = 2000, height = 1600, res = 300)

cols <- hcl.colors(100, "Inferno")  # or "Viridis", "Plasma", "Heat2"

plot(crop_suit_filled,
     col = cols,
     axes = FALSE,
     legend = TRUE
)

plot(st_geometry(gnp), add = TRUE, border = "white", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex=0.5, col="white")

dev.off()


#--- Calculate distance to roads ----
# Buffer by ~ half a cell (in raster units) to avoid broken roads
# Not necessary for the pilot area but keeping code here in case it's useful in future
#half_cell <- max(res(pilot_r)) * 0.51
#roads <- st_buffer(roads, dist = half_cell)

roads_mask <- rasterize(vect(roads), pilot_r,
                        field = 1, background = NA, touches = TRUE)

dist_to_road_pad <- distance(roads_mask)

dist_to_road <- crop(dist_to_road_pad, pilot)
dist_to_road <- mask(dist_to_road, pilot)

# Save to disk 
png("./results/landcover/pilot/dist_to_road.png", width = 2000, height = 1600, res = 300)

cols <- rev(hcl.colors(100, "Inferno"))  # or "Viridis", "Plasma", "Heat2"

plot(dist_to_road,
     col = cols,
     axes = FALSE,
     legend = TRUE
)

plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex=0.5)

dev.off()


#--- Calculate landscape context within moving window ----
# window size in metres
win_size_m <- 3200

# convert to number of cells
win_size_cells <- round(win_size_m / res(lulc_pilot_pad)[1])

# create a square (3.6 km × 3.6 km) window
w <- matrix(1, nrow = win_size_cells, ncol = win_size_cells)

#--- Cropland ----
# Binary cropland raster
cropland_pad <- ifel((lulc_pilot_pad == 51) | (lulc_pilot_pad == 52), 1, 0)
plot(cropland_pad)

# focal mean gives proportion of 1s
cropland_prop_pad <- focal(cropland_pad, w = w, fun = "mean", na.policy = "omit", na.rm = TRUE)

# visualize
cropland_prop <- mask(cropland_prop_pad, pilot)
cropland_prop <- crop(cropland_prop, pilot)

# Save to disk 
png("./results/landcover/pilot/cropland_proportion.png", width = 2000, height = 1600, res = 300)

cols <- hcl.colors(100, "Inferno")  # or "Viridis", "Plasma", "Heat2"

plot(cropland_prop,
     col = cols,
     axes = FALSE,
     legend = TRUE
)

plot(st_geometry(gnp), add = TRUE, border = "white", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex=0.5, col="white")

dev.off()


# # Loop over to test sensitivity of window size
# # binary cropland raster (0/1)
# r <- cropland_pad  
# 
# # window sizes in metres
# win_sizes <- c(800, 1600, 3200, 3600, 7200)
# 
# # prepare an empty list to store results
# prop_maps <- list()
# 
# for (w_m in win_sizes) {
#   cat("Running window:", w_m, "m\n")
#   
#   # convert to number of cells
#   w_cells <- round(w_m / res(r)[1])
#   
#   # force odd size (needed by terra::focal)
#   if (w_cells %% 2 == 0) w_cells <- w_cells + 1
#   
#   # build square window
#   w <- matrix(1, nrow = w_cells, ncol = w_cells)
#   
#   # focal mean → proportion of 1s in window
#   prop_maps[[as.character(w_m)]] <- focal(r, w = w, fun = "mean", na.policy = "omit", na.rm = TRUE)
# }
# 
# # --- Compute pairwise Pearson correlations among window sizes ---
# vals <- lapply(prop_maps, function(x) values(x, na.rm = TRUE))
# 
# common_na <- Reduce("|", lapply(vals, is.na))
# vals_clean <- lapply(vals, function(v) v[!common_na])
# 
# corr_mat <- matrix(NA, nrow = length(win_sizes), ncol = length(win_sizes),
#                    dimnames = list(win_sizes, win_sizes))
# 
# for (i in seq_along(win_sizes)) {
#   for (j in seq_along(win_sizes)) {
#     corr_mat[i,j] <- cor(vals_clean[[i]], vals_clean[[j]], method = "spearman")
#   }
# }
# 
# print(round(corr_mat, 3))
# 
# 
# corr_df <- melt(corr_mat, varnames = c("Window1", "Window2"), value.name = "Correlation")
# 
# mean_corr <- apply(corr_mat, 1, function(x) mean(x[-which(names(x) == names(x)[which.max(x)])]))
# barplot(mean_corr,
#         names.arg = win_sizes,
#         col = "grey",
#         main = "Average correlation with other window sizes",
#         ylab = "Mean correlation",
#         ylim = c(0, 1))



# Visualize with the moving window size shown
# find the center of the raster
# Response variable for conversion risk (grassland -> cropland)
cropland <- mask(cropland_pad, pilot)
cropland <- crop(cropland, pilot)

e <- ext(cropland_prop_pad)
cx <- (e$xmin-8000 + e$xmax) / 2
cy <- (e$ymin-8000 + e$ymax) / 2

# half the side length (3.6 km window)
half_side <- 3200 / 2

# build a square extent and turn it into a polygon
rect_ext <- ext(cx - half_side, cx + half_side, cy - half_side, cy + half_side)
square_poly <- as.polygons(rect_ext, crs = crs(r))   # SpatVector polygon

# Save to disk 
png("./results/landcover/pilot/cropland_show_moving_window.png", width = 2000, height = 1600, res = 300)

# plot cropland
plot(cropland, col = c("grey85","#FEA500"),
     axes = FALSE, 
     legend = FALSE)
lines(square_poly, col = "red", lwd = 2)

plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex=0.5)

dev.off()

# Visualize without moving window shown 
png("./results/landcover/pilot/cropland.png", width = 2000, height = 1600, res = 300)

# plot cropland
plot(cropland, col = c("grey85","#FEA500"),
     axes = FALSE,
     legend = FALSE)
plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex=0.5)

dev.off()


#--- Settlement ----
# Binary developed raster
# Do not include roads
settle_pad <- ifel((lulc_pilot_pad == 21) | (lulc_pilot_pad == 22) | (lulc_pilot_pad == 24) | (lulc_pilot_pad == 28) | (lulc_pilot_pad == 29), 1, 0)
plot(settle_pad)

# focal mean gives proportion of 1s
settle_prop_pad <- focal(settle_pad, w = w, fun = "mean", na.policy = "omit", na.rm = TRUE)

# visualize
settle_prop <- mask(settle_prop_pad, pilot)
settle_prop <- crop(settle_prop, pilot)

# Save to disk 
png("./results/landcover/pilot/settle_proportion.png", width = 2000, height = 1600, res = 300)

cols <- hcl.colors(100, "Inferno")  # or "Viridis", "Plasma", "Heat2"

plot(settle_prop,
     col = cols,
     axes = FALSE,
     legend = TRUE
)

plot(st_geometry(gnp), add = TRUE, border = "white", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex=0.5, col="white")

dev.off()

# Response variable for conversion risk (grassland -> settlement)
settle <- mask(settle_pad, pilot)
settle <- crop(settle, pilot)

# Visualize
png("./results/landcover/pilot/settle.png", width = 2000, height = 1600, res = 300)

# plot settlement
plot(settle, col = c("grey85","#7B3294"),
     axes = FALSE,
     legend = FALSE)
plot(st_geometry(gnp), add = TRUE, border = "black", lwd = 0.8)
text(gnp_cent[1], gnp_cent[2], "Grasslands National Park", cex=0.5)

dev.off()


#--- Write rasters to disk ----
# Cropland suitability
writeRaster(crop_suit, file.path(pilot_filepath, "cropland_suitability_2019_withNAs.tif"))
writeRaster(crop_suit_filled, file.path(pilot_filepath, "cropland_suitability_2019_noNAs.tif"))
# Roads
writeRaster(dist_to_road, file.path(pilot_filepath, "dist_to_road.tif"))
# Proportion cropland
writeRaster(cropland_prop, file.path(pilot_filepath, "cropland_proportion_3-2.tif"))
# Proportion settlement
writeRaster(settle_prop, file.path(pilot_filepath, "settle_proportion_3-2.tif"))
# Cropland - binary
writeRaster(cropland, file.path(pilot_filepath, "cropland_binary.tif"))
# Settlement - binary
writeRaster(settle, file.path(pilot_filepath, "settle_binary.tif"))
