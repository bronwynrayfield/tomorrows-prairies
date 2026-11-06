##############################################
# Data pre-processing including cropping, imputing, derived calculations

#--- Setup ----
source(file.path("./scripts", "0-setup.R"))

#--- Load files ---
pilot <- st_read(file.path("./processeddatasets", pilot_filepath, pilot_s_filename))
pilot_r_prairies <- rast(file.path("./processeddatasets", pilot_filepath, pilot_r_filename))

# Olimb & Robinson 2019 cropland suitability
crop_suit_reproj<- rast(file.path("./processeddatasets", crop_suit_2019_filepath, paste0(crop_suit_2019_filename, ".tif")))

# Census subdivisions
census_reproj <- st_read(file.path("./processeddatasets", census_subdivision_filepath, census_subdivision_filename))

# Roads
roads_reproj <-st_read(file.path("./processeddatasets", road_filepath, road_filename), layer = "NRN_SK_14_0_ROADSEG")

# LULC
lulc_prairies_reproj <- rast(file.path("./processeddatasets", "landcover/aafc/reprojected", aafc2020_filename))

# ECCC loss
grass_loss_class_reproj <- rast(file.path("./processeddatasets", grass_loss_filepath, grass_loss_class_filename))
grass_loss_year_reproj <- rast(file.path("./processeddatasets", grass_loss_filepath, grass_loss_year_filename))  
ipast_loss_class_reproj <- rast(file.path("./processeddatasets", ipast_loss_filepath, ipast_loss_class_filename))
ipast_loss_year_reproj <- rast(file.path("./processeddatasets", ipast_loss_filepath, ipast_loss_year_filename))
upast_loss_class_reproj <- rast(file.path("./processeddatasets", upast_loss_filepath, upast_loss_class_filename))
upast_loss_year_reproj <- rast(file.path("./processeddatasets", upast_loss_filepath, upast_loss_year_filename))

# Protected areas
pa_reproj <- st_read(file.path("./processeddatasets", pa_filepath, paste0(pa_filename, ".gdb")), layer = pa_filename)

#--- Crop and mask to pilot area or padded pilot area ----
# Pilot raster template
pilot_r <- crop(pilot_r_prairies, pilot)
pilot_r <- mask(pilot_r, pilot)

# Census subdivisions
census <- st_intersection(census_reproj, pilot)
# #plot(census["CSDNAME"])

# Protected areas
# This didn't work so I had to run an intersection manually in QGIS
# pa <- st_intersection(pa_reproj, pilot)
# gnp <- pa[pa$NAME_E == "Grasslands National Park Of Canada", ]

#--- Pad the pilot area for use with some spatial predictors ----
# Pad the pilot extent by 20 cells
pad <- res(pilot_r) * 120
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

#--- Crop and mask to focal census subdivision (csd)----
# Create mask with focal csd
focal_csd <- census %>%
  filter(CSDNAME %in% c("Val Marie", "Val Marie No. 17"))

csd_mask <- rasterize(vect(focal_csd), pilot_r,
                      field = 1, background = NA, touches = TRUE)

# Create new pilot area boundary
# merge the two focal csd polygons & remove hole
pilot_new <- pilot %>%
  st_union() %>%        # dissolve into one multipart geometry
  st_make_valid()      # fix donut holes, ensure simple polygon

# Cropland suitability
crop_suit_pilot <- crop(crop_suit_reproj, pilot)
crop_suit_pilot <- mask(crop_suit_pilot, pilot)
# plot(crop_suit_pilot, main = "Cropland suitability (2019)")
# plot(st_geometry(census),
#      add = TRUE, col = NA, border = "black", lwd = 0.6)

# Cropland suitability
crop_suit <- mask(crop_suit_pilot, csd_mask)
crop_suit <- crop(crop_suit, csd_mask)
plot(crop_suit, main = "Cropland suitability (2019)")

# # LULC
# lulc <- mask(lulc_pilot, csd_mask)
# lulc <- crop(lulc, csd_mask)
# plot(lulc, main = "LULC (2020)")

# grass loss
grass_loss_class <- crop(grass_loss_class_reproj, csd_mask)
grass_loss_class <- mask(grass_loss_class, csd_mask)
grass_loss_year <- crop(grass_loss_year_reproj, csd_mask)
grass_loss_year <- mask(grass_loss_year, csd_mask)

# ipast loss
ipast_loss_class <- crop(ipast_loss_class_reproj, csd_mask)
ipast_loss_class <- mask(ipast_loss_class, csd_mask)
ipast_loss_year <- crop(ipast_loss_year_reproj, csd_mask)
ipast_loss_year <- mask(ipast_loss_year, csd_mask)

# upast loss
upast_loss_class <- crop(upast_loss_class_reproj, csd_mask)
upast_loss_class <- mask(upast_loss_class, csd_mask)
upast_loss_year <- crop(upast_loss_year_reproj, csd_mask)
upast_loss_year <- mask(upast_loss_year, csd_mask)


#--- Write to disk ----
writeRaster(pilot_r, file.path(pilot_filepath, "GNP_pilot_area_crop_CSD.tif"), overwrite = FALSE)
st_write(pilot_new, file.path(pilot_filepath, "GNP_pilot_area_crop_CSD.shp"), delete_dsn = FALSE)
st_write(focal_csd, file.path(pilot_filepath, "census_pilot.shp"), delete_dsn = FALSE)
writeRaster(csd_mask, file.path("./processeddatasets", pilot_filepath, "census_pilot.tif"), overwrite = FALSE)
writeRaster(crop_suit, file.path("./processeddatasets", pilot_filepath, paste0(crop_suit_2019_filename, ".tif")), overwrite = FALSE)
writeRaster(lulc_pilot, file.path("./processeddatasets", pilot_filepath, "AAFC2020_40kmbuffer_noUS_pilot_pad.tif"), overwrite = FALSE)
st_write(census, file.path("./processeddatasets", pilot_filepath, census_subdivision_filename), delete_dsn = FALSE)
st_write(roads, file.path("./processeddatasets", pilot_filepath, "NRN_SK_14_0_GPKG_en_pilot_pad.gpkg"), layer = "NRN_SK_14_0_ROADSEG", delete_dsn = FALSE)
writeRaster(grass_loss_class, file.path("./processeddatasets", pilot_filepath, grass_loss_class_filename), overwrite = FALSE)
writeRaster(grass_loss_year, file.path("./processeddatasets", pilot_filepath, grass_loss_year_filename), overwrite = FALSE)
writeRaster(ipast_loss_class, file.path("./processeddatasets", pilot_filepath, ipast_loss_class_filename), overwrite = FALSE)
writeRaster(ipast_loss_year, file.path("./processeddatasets", pilot_filepath, ipast_loss_year_filename), overwrite = FALSE)
writeRaster(upast_loss_class, file.path("./processeddatasets", pilot_filepath, upast_loss_class_filename), overwrite = FALSE)
writeRaster(upast_loss_year, file.path("./processeddatasets", pilot_filepath, upast_loss_year_filename), overwrite = FALSE)
#st_write(pa, file.path("./processeddatasets", pilot_filepath, "pa_pilot.shp"), delete_dsn = TRUE)
#st_write(gnp, file.path("./processeddatasets", pilot_filepath, "gnp_park_pilot.shp"))

