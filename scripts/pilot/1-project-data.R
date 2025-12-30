##############################################
# Data projection
# Executed separately because very slow
# Only run once

#--- Setup ---
source(file.path("scripts", "0-setup.R"))

#--- Directories ---
# Pilot area directory
pilot_dir <- "./processeddatasets/pilot_area"

#--- Load files ---
# Pilot area
pilot <- st_read(file.path(pilot_dir, "pilot_area_test_aligned_cnvc.shp"))
pilot_rast_prairies <- rast(file.path(pilot_dir, "pilot_area_test_aligned_cnvc.tif"))

# Olimb & Robinson 2019 cropland suitability
crop_suit_raw <- rast(file.path("./existingdatasets", crop_suit_2019_filepath, crop_suit_2019_filename))

# Census subdivisions
census_raw <- st_read(file.path("./existingdatasets", census_subdivision_filepath, census_subdivision_filename))

# Roads
#roads_raw <-st_read(road_file, layer = "NRN_SK_14_0_ROADSEG")
roads_raw <-st_read(file.path("./existingdatasets", road_filepath, road_filename), layer = "NRN_SK_14_0_ROADSEG")

# LULC
lulc_prairies_raw <- rast(file.path("./processeddatasets", aafc2020_filepath, aafc2020_filename))

# ECCC grass loss
grass_loss_class_raw <- rast(file.path("./existingdatasets", grass_loss_filepath, grass_loss_class_filename))
grass_loss_year_raw <- rast(file.path("./existingdatasets", grass_loss_filepath, grass_loss_year_filename))  
ipast_loss_class_raw <- rast(file.path("./existingdatasets", ipast_loss_filepath, ipast_loss_class_filename))
ipast_loss_year_raw <- rast(file.path("./existingdatasets", ipast_loss_filepath, ipast_loss_year_filename))
upast_loss_class_raw <- rast(file.path("./existingdatasets", upast_loss_filepath, upast_loss_class_filename))
upast_loss_year_raw <- rast(file.path("./existingdatasets", upast_loss_filepath, upast_loss_year_filename))

# Protected Areas
pa_raw <- st_read(file.path("./existingdatasets", pa_filepath, paste0(pa_filename, ".gdb")), layer = pa_filename)

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
pa_reproj <- st_transform(pa_raw, st_crs(pilot))

#--- Write to disk ---
writeRaster(crop_suit_reproj, file.path("./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename), overwrite = FALSE)
writeRaster(lulc_prairies_reproj, file.path("./processeddatasets", "landcover/aafc/reprojected", aafc2020_filename), overwrite = FALSE)
st_write(census_reproj, file.path("./processeddatasets", census_subdivision_filepath, census_subdivision_filename), delete_dsn = FALSE)
st_write(roads_reproj, file.path("./processeddatasets", road_filepath, road_filename), layer = "NRN_SK_14_0_ROADSEG", delete_dsn = FALSE)
writeRaster(grass_loss_class_reproj, file.path("./processeddatasets", grass_loss_filepath, grass_loss_class_filename), overwrite = FALSE)
writeRaster(grass_loss_year_reproj, file.path("./processeddatasets", grass_loss_filepath, grass_loss_year_filename), overwrite = FALSE)
writeRaster(ipast_loss_class_reproj, file.path("./processeddatasets", ipast_loss_filepath, ipast_loss_class_filename), overwrite = FALSE)
writeRaster(ipast_loss_year_reproj, file.path("./processeddatasets", ipast_loss_filepath, ipast_loss_year_filename), overwrite = FALSE)
writeRaster(upast_loss_class_reproj, file.path("./processeddatasets", upast_loss_filepath, upast_loss_class_filename), overwrite = FALSE)
writeRaster(upast_loss_year_reproj, file.path("./processeddatasets", upast_loss_filepath, upast_loss_year_filename), overwrite = FALSE)
# Note that this throws a write error: GDAL Error 6: Unsupported geometry type
#st_write(pa_reproj, file.path("./processeddatasets", pa_filepath, paste0(pa_filename, ".gdb")), layer = pa_filename, delete_dsn = FALSE)

         