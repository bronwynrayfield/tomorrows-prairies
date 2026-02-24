# Data projection
# Apply NA values and scaling factor if necessary
# Executed separately because very slow for rasters especially
# Only run once

# Setup
source(file.path("scripts", "0-setup.R"))
source(file.path("scripts", "0-functions-general.R"))

#--- Load files ----

# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))

# Olimb & Robinson 2019 cropland suitability
crop_suit_raw <- rast(file.path("./existingdatasets", capitalize_first(crop_suit_2019_filepath), crop_suit_2019_filename))
# Census subdivisions
census_raw <- st_read(file.path("./existingdatasets", capitalize_first(census_subdivision_filepath), census_subdivision_filename))
# National ecological framework
ecozones_raw <- st_read(file.path("./existingdatasets", capitalize_first(ecozones_filepath), ecozones_filename))
#Note that ecoregions and ecodistrics were preprocessed manually to deal with misalignment
ecoregions_raw <- st_read(file.path("./processeddatasets", ecoregions_filepath, ecoregions_filename))
ecodistricts_raw <- st_read(file.path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename))
# Roads
#roads_raw <-st_read(road_file, layer = "NRN_SK_14_0_ROADSEG")
roads_raw <-st_read(file.path("./existingdatasets", capitalize_first(road_filepath), road_filename), layer = "NRN_SK_14_0_ROADSEG")
# LULC
lulc2020_prairies_raw <- rast(file.path("./processeddatasets", capitalize_first(aafc_filepath), "merged", aafc2020_filename))
lulc2015_prairies_raw <- rast(file.path("./processeddatasets", capitalize_first(aafc_filepath), "merged", aafc2015_filename))
lulc2010_prairies_raw <- rast(file.path("./processeddatasets", capitalize_first(aafc_filepath), "merged", aafc2010_filename))
lulc2000_prairies_raw <- rast(file.path("./processeddatasets", capitalize_first(aafc_filepath), "merged", aafc2000_filename))
# ECCC most recent loss
grass_loss_class_raw <- rast(file.path("./existingdatasets", capitalize_first(grass_loss_filepath), grass_loss_class_filename))
grass_loss_year_raw <- rast(file.path("./existingdatasets", capitalize_first(grass_loss_filepath), grass_loss_year_filename))  
ipast_loss_class_raw <- rast(file.path("./existingdatasets", capitalize_first(ipast_loss_filepath), ipast_loss_class_filename))
ipast_loss_year_raw <- rast(file.path("./existingdatasets", capitalize_first(ipast_loss_filepath), ipast_loss_year_filename))
upast_loss_class_raw <- rast(file.path("./existingdatasets", capitalize_first(upast_loss_filepath), upast_loss_class_filename))
upast_loss_year_raw <- rast(file.path("./existingdatasets", capitalize_first(upast_loss_filepath), upast_loss_year_filename))
crop_loss_class_raw <- rast(file.path("./existingdatasets", capitalize_first(crop_loss_filepath), crop_loss_class_filename))
crop_loss_year_raw <- rast(file.path("./existingdatasets", capitalize_first(crop_loss_filepath), crop_loss_year_filename))  
# ECCC most recent gain
crop_gain_class_raw <- rast(file.path("./existingdatasets", capitalize_first(crop_gain_filepath), crop_gain_class_filename))
crop_gain_year_raw <- rast(file.path("./existingdatasets", capitalize_first(crop_gain_filepath), crop_gain_year_filename))  
builtup_gain_class_raw <- rast(file.path("./existingdatasets", capitalize_first(builtup_gain_filepath), builtup_gain_class_filename))
builtup_gain_year_raw <- rast(file.path("./existingdatasets", capitalize_first(builtup_gain_filepath), builtup_gain_year_filename))  
# Protected Areas
pa_raw <- st_read(file.path("./existingdatasets", capitalize_first(pa_filepath), paste0(pa_filename, ".gdb")), layer = pa_filename)
# Soil
oc_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), oc_filename))
bd_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), bd_filename))
clay_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), clay_filename))
sand_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), sand_filename))
silt_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), silt_filename))
ph_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), ph_filename))

#--- Unscale to real units and apply NAs ----
oc_raw[oc_raw == oc_nodata] <- NA
oc_real <- oc_raw / oc_scale

bd_raw[bd_raw == bd_nodata] <- NA
bd_real <- bd_raw / bd_scale

clay_raw[clay_raw == clay_nodata] <- NA
sand_raw[sand_raw == sand_nodata] <- NA
silt_raw[silt_raw == silt_nodata] <- NA

ph_raw[ph_raw == ph_nodata] <- NA
ph_real <- ph_raw / ph_scale

#--- Reproject to match study area shapefile CRS ----

crop_suit_reproj <- project(crop_suit_raw, sa, method = "bilinear")
lulc2020_prairies_reproj <- project(lulc2020_prairies_raw, sa, method = "near")
lulc2015_prairies_reproj <- project(lulc2015_prairies_raw, sa, method = "near")
lulc2010_prairies_reproj <- project(lulc2010_prairies_raw, sa, method = "near")
lulc2000_prairies_reproj <- project(lulc2000_prairies_raw, sa, method = "near")
census_reproj <- st_transform(census_raw, st_crs(sa))
ecozones_reproj <- st_transform(ecozones_raw, st_crs(sa))
ecoregions_reproj <- st_transform(ecoregions_raw, st_crs(sa))
ecodistricts_reproj <- st_transform(ecodistricts_raw, st_crs(sa))
grass_loss_class_reproj <- project(grass_loss_class_raw, sa, method = "near")
grass_loss_year_reproj <- project(grass_loss_year_raw, sa, method = "near")
ipast_loss_class_reproj <- project(ipast_loss_class_raw, sa, method = "near")
ipast_loss_year_reproj <- project(ipast_loss_year_raw, sa, method = "near")
upast_loss_class_reproj <- project(upast_loss_class_raw, sa, method = "near")
upast_loss_year_reproj <- project(upast_loss_year_raw, sa, method = "near")
crop_loss_class_reproj <- project(crop_loss_class_raw, sa, method = "near")
crop_loss_year_reproj <- project(crop_loss_year_raw, sa, method = "near")
crop_gain_class_reproj <- project(crop_gain_class_raw, sa, method = "near")
crop_gain_year_reproj <- project(crop_gain_year_raw, sa, method = "near")
builtup_gain_class_reproj <- project(builtup_gain_class_raw, sa, method = "near")
builtup_gain_year_reproj <- project(builtup_gain_year_raw, sa, method = "near")
oc_reproj <- project(oc_real, sa, method = "bilinear")
bd_reproj <- project(bd_real, sa, method = "bilinear")
clay_reproj <- project(clay_raw, sa, method = "bilinear")
sand_reproj <- project(sand_raw, sa, method = "bilinear")
silt_reproj <- project(silt_raw, sa, method = "bilinear")
ph_reproj <- project(ph_real, sa, method = "bilinear")

#roads_reproj <- st_transform(roads_raw, st_crs(sa))
#pa_reproj <- st_transform(pa_raw, st_crs(sa))

#--- Write to disk ----

write_rast_tagged(crop_suit_reproj, "./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "proj-naaea")
write_sf_tagged(census_reproj, "./processeddatasets", census_subdivision_filepath, census_subdivision_filename, "proj-naaea")
write_sf_tagged(ecozones_reproj, "./processeddatasets", ecozones_filepath, ecozones_filename, "proj-naaea")
write_sf_tagged(ecoregions_reproj, "./processeddatasets", ecoregions_filepath, ecoregions_filename, "proj-naaea")
write_sf_tagged(ecodistricts_reproj, "./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "proj-naaea")
write_rast_tagged(lulc2020_prairies_reproj, "./processeddatasets", aafc_filepath, "reprojected", aafc2020_filename, "proj-naaea")
write_rast_tagged(lulc2015_prairies_reproj, "./processeddatasets", aafc_filepath, "reprojected", aafc2015_filename, "proj-naaea")
write_rast_tagged(lulc2010_prairies_reproj, "./processeddatasets", aafc_filepath, "reprojected", aafc2010_filename, "proj-naaea")
write_rast_tagged(lulc2000_prairies_reproj, "./processeddatasets", aafc_filepath, "reprojected", aafc2000_filename, "proj-naaea")
write_rast_tagged(grass_loss_class_reproj, "./processeddatasets", grass_loss_filepath, grass_loss_class_filename, "proj-naaea")
write_rast_tagged(grass_loss_year_reproj, "./processeddatasets", grass_loss_filepath, grass_loss_year_filename, "proj-naaea")
write_rast_tagged(ipast_loss_class_reproj, "./processeddatasets", ipast_loss_filepath, ipast_loss_class_filename, "proj-naaea")
write_rast_tagged(ipast_loss_year_reproj, "./processeddatasets", ipast_loss_filepath, ipast_loss_year_filename, "proj-naaea")
write_rast_tagged(upast_loss_class_reproj, "./processeddatasets", upast_loss_filepath, upast_loss_class_filename, "proj-naaea")
write_rast_tagged(upast_loss_year_reproj, "./processeddatasets", upast_loss_filepath, upast_loss_year_filename, "proj-naaea")
write_rast_tagged(crop_loss_class_reproj, "./processeddatasets", crop_loss_filepath, crop_loss_class_filename, "proj-naaea")
write_rast_tagged(crop_loss_year_reproj, "./processeddatasets", crop_loss_filepath, crop_loss_year_filename, "proj-naaea")
write_rast_tagged(crop_gain_class_reproj, "./processeddatasets", crop_gain_filepath, crop_gain_class_filename, "proj-naaea")
write_rast_tagged(crop_gain_year_reproj, "./processeddatasets", crop_gain_filepath, crop_gain_year_filename, "proj-naaea")
write_rast_tagged(builtup_gain_class_reproj, "./processeddatasets", builtup_gain_filepath, builtup_gain_class_filename, "proj-naaea")
write_rast_tagged(builtup_gain_year_reproj, "./processeddatasets", builtup_gain_filepath, builtup_gain_year_filename, "proj-naaea")
write_rast_tagged(oc_reproj, "./processeddatasets", soil_filepath, oc_filename, "proj-naaea", datatype="FLT4S")
write_rast_tagged(bd_reproj, "./processeddatasets", soil_filepath, bd_filename, "proj-naaea", datatype="FLT4S")
write_rast_tagged(clay_reproj, "./processeddatasets", soil_filepath, clay_filename, "proj-naaea", datatype="FLT4S")
write_rast_tagged(sand_reproj, "./processeddatasets", soil_filepath, sand_filename, "proj-naaea", datatype="FLT4S")
write_rast_tagged(silt_reproj, "./processeddatasets", soil_filepath, silt_filename, "proj-naaea", datatype="FLT4S")
write_rast_tagged(ph_reproj, "./processeddatasets", soil_filepath, ph_filename, "proj-naaea", datatype="FLT4S")

#st_write(roads_reproj, file.path("./processeddatasets", road_filepath, road_filename), layer = "NRN_SK_14_0_ROADSEG", delete_dsn = FALSE)
# Note that this throws a write error: GDAL Error 6: Unsupported geometry type
#st_write(pa_reproj, file.path("./processeddatasets", pa_filepath, paste0(pa_filename, ".gdb")), layer = pa_filename, delete_dsn = FALSE)

