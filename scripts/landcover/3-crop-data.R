# =============================================================================
# 3-crop-data.R
# =============================================================================
#
# Purpose:
#   Crops and masks all datasets to the modelling extent.
#   Also derives topographic layers (slope, aspect) from the DEM.
#   All outputs are saved to processeddatasets/ with "modelling-extent" tag.
#
# Inputs:
#   - Modelling extent (shapefile and raster template)
#   - Pilot area shapefile
#   - All reprojected datasets from 2-project-data.R:
#       AAFC LUTS 2000, 2010, 2015, 2020
#       ECCC most recent change layers (grass, improved pasture, 
#         unimproved pasture, crop, built-up)
#       Olimb & Robinson 2019 cropland suitability
#       Census subdivisions
#       National ecological framework (ecozones, ecoregions, ecodistricts)
#       Protected areas (CPCAD 2024)
#       Soil (OC, BD, clay, sand, silt, pH)
#       DEM (NRCan MRDEM 30m)
#
# Outputs:
#   - All datasets cropped/masked to modelling extent
#   - All datasets cropped/masked to pilot area (AAFC LUTS only)
#   - Slope and aspect derived from DEM
#
# Notes:
#   - Roads processing incomplete — rasterization pending
#   - pa uses st_crop() to bounding box rather than st_intersection()
#     for speed; NA handling deferred to rasterization step
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 2-project-data.R outputs
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Cropland
crop_suit_reproj<- rast(tagged_path("./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "proj-naaea"))
plowprint_reproj<- rast(tagged_path("./processeddatasets", plowprint_filepath, plowprint_filename, "proj-naaea"))

# Census subdivisions
census_reproj <- st_read(tagged_path("./processeddatasets", census_subdivision_filepath, census_subdivision_filename, "proj-naaea"))

# National ecological framework
ecozones_reproj <- st_read(tagged_path("./processeddatasets", ecozones_filepath, ecozones_filename, "proj-naaea"))
ecoregions_reproj <- st_read(tagged_path("./processeddatasets", ecoregions_filepath, ecoregions_filename, "proj-naaea"))
ecodistricts_reproj <- st_read(tagged_path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "proj-naaea"))

# Roads
roads_reproj <-st_read(file.path("./processeddatasets", road_filepath, road_filename), layer = "NRN_SK_14_0_ROADSEG")

# LULC - AAFC LUTS
lulc2020_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2020_filename, "proj-naaea"))
lulc2015_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2015_filename, "proj-naaea"))
lulc2010_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2010_filename, "proj-naaea"))
lulc2000_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2000_filename, "proj-naaea"))

# ECCC most recent change
grass_loss_class_reproj <- rast(tagged_path("./processeddatasets", grass_loss_filepath, grass_loss_class_filename, "proj-naaea"))
grass_loss_year_reproj <- rast(tagged_path("./processeddatasets", grass_loss_filepath, grass_loss_year_filename, "proj-naaea"))  
ipast_loss_class_reproj <- rast(tagged_path("./processeddatasets", ipast_loss_filepath, ipast_loss_class_filename, "proj-naaea"))
ipast_loss_year_reproj <- rast(tagged_path("./processeddatasets", ipast_loss_filepath, ipast_loss_year_filename, "proj-naaea"))
upast_loss_class_reproj <- rast(tagged_path("./processeddatasets", upast_loss_filepath, upast_loss_class_filename, "proj-naaea"))
upast_loss_year_reproj <- rast(tagged_path("./processeddatasets", upast_loss_filepath, upast_loss_year_filename, "proj-naaea"))
crop_loss_class_reproj <- rast(tagged_path("./processeddatasets", crop_loss_filepath, crop_loss_class_filename, "proj-naaea"))
crop_loss_year_reproj <- rast(tagged_path("./processeddatasets", crop_loss_filepath, crop_loss_year_filename, "proj-naaea"))
crop_gain_class_reproj <- rast(tagged_path("./processeddatasets", crop_gain_filepath, crop_gain_class_filename, "proj-naaea"))
crop_gain_year_reproj <- rast(tagged_path("./processeddatasets", crop_gain_filepath, crop_gain_year_filename, "proj-naaea"))
builtup_gain_class_reproj <- rast(tagged_path("./processeddatasets", builtup_gain_filepath, builtup_gain_class_filename, "proj-naaea"))
builtup_gain_year_reproj <- rast(tagged_path("./processeddatasets", builtup_gain_filepath, builtup_gain_year_filename, "proj-naaea"))

# Protected Areas
pa_reproj <- st_read(
  tagged_path("./processeddatasets", pa_filepath, pa_filename, "proj-naaea", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(pa_filename), "_proj-naaea")
)
ncc_reproj <- st_read(
  tagged_path("./processeddatasets", ncc_easement_filepath, ncc_easement_filename,
              "proj-naaea", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(ncc_easement_filename), "_proj-naaea")
)

# Soil
oc_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, oc_filename, "proj-naaea"))
bd_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, bd_filename, "proj-naaea"))
clay_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, clay_filename, "proj-naaea"))
sand_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, sand_filename, "proj-naaea"))
silt_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, silt_filename, "proj-naaea"))
ph_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, ph_filename, "proj-naaea"))

# Soil Moisture
psssm_mean_reproj <- rast(tagged_path("./processeddatasets", psssm_filepath, psssm_mean_filename, "proj-naaea"))
psssm_sd_reproj <- rast(tagged_path("./processeddatasets", psssm_filepath, psssm_sd_filename, "proj-naaea"))
tsm_planting_reproj <- rast(tagged_path("./processeddatasets", tsm_filepath, tsm_planting_filename, "proj-naaea"))
tsm_growing_reproj <- rast(tagged_path("./processeddatasets", tsm_filepath, tsm_growing_filename, "proj-naaea"))

# Topography
dem_reproj <- rast(tagged_path("./processeddatasets", dem_filepath, dem_filename, "proj-naaea"))

# Climate
gddgrow_5_hist_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, gddgrow_5_hist_filename, "proj-naaea"))
gddgrow_5_fut_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, gddgrow_5_fut_filename, "proj-naaea"))
hddheat_18_hist_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, hddheat_18_hist_filename, "proj-naaea"))
hddheat_18_fut_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, hddheat_18_fut_filename, "proj-naaea"))
prtot_hist_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, prtot_hist_filename, "proj-naaea"))
prtot_fut_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, prtot_fut_filename, "proj-naaea"))
tn_mean_hist_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, tn_mean_hist_filename, "proj-naaea"))
tn_mean_fut_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, tn_mean_fut_filename, "proj-naaea"))
tx_mean_hist_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, tx_mean_hist_filename, "proj-naaea"))
tx_mean_fut_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, tx_mean_fut_filename, "proj-naaea"))
nr_cdd_hist_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, nr_cdd_hist_filename, "proj-naaea"))
nr_cdd_fut_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, nr_cdd_fut_filename, "proj-naaea"))

# Socio-economic
canfax_reproj <- st_read(tagged_path("./processeddatasets", canfax_filepath, canfax_filename, "proj-naaea"))

#--- Crop and mask to pilot area ----
# LULC - AAFC LUTS
lulc2020_pilot <- mask(crop(lulc2020_prairies_reproj, pilot), pilot)
lulc2015_pilot <- mask(crop(lulc2015_prairies_reproj, pilot), pilot)
crop_suit_pilot <- mask(crop(crop_suit_reproj, pilot), pilot)

#--- Crop and mask to modelling extent ----
# Census subdivisions
census <- st_intersection(census_reproj, me)
# plot(census["CSDNAME"])

# National ecological framework
ecozones <- st_intersection(ecozones_reproj, me)
ecoregions <- st_intersection(ecoregions_reproj, me)
ecodistricts <- st_intersection(ecodistricts_reproj, me)

# Protected Areas
# using st_crop() instead of st_intersection() for speed
# will rasterize later which will take care of NAs
pa <- st_crop(pa_reproj, st_bbox(me_r))
ncc <- st_crop(ncc_reproj, st_bbox(me_r))

# Roads
roads <- st_intersection(roads_reproj, me)
plot(roads["ROADCLASS"], main = "Roads by Class", key.pos = 4)
roads <- mask(crop(roads_r, me), me)

# Cropland
crop_suit <- mask(crop(crop_suit_reproj, me_r), me_r)
plot(crop_suit, main = "Cropland suitability (2019)")
plowprint <- mask(crop(plowprint_reproj, me_r), me_r)

# LULC - AAFC LUTS
lulc2020 <- mask(crop(lulc2020_prairies_reproj, me_r), me_r)
lulc2015 <- mask(crop(lulc2015_prairies_reproj, me_r), me_r)
lulc2010 <- mask(crop(lulc2010_prairies_reproj, me_r), me_r)
lulc2000 <- mask(crop(lulc2000_prairies_reproj, me_r), me_r)

# ECCC most recent change
grass_loss_class <- mask(crop(grass_loss_class_reproj, me_r), me_r)
grass_loss_year <- mask(crop(grass_loss_year_reproj, me_r), me_r)
ipast_loss_class <- mask(crop(ipast_loss_class_reproj, me_r), me_r)
ipast_loss_year <- mask(crop(ipast_loss_year_reproj, me_r), me_r)
upast_loss_class <- mask(crop(upast_loss_class_reproj, me_r), me_r)
upast_loss_year <- mask(crop(upast_loss_year_reproj, me_r), me_r)
crop_loss_class <- mask(crop(crop_loss_class_reproj, me_r), me_r)
crop_loss_year <- mask(crop(crop_loss_year_reproj, me_r), me_r)
crop_gain_class <- mask(crop(crop_gain_class_reproj, me_r), me_r)
crop_gain_year <- mask(crop(crop_gain_year_reproj, me_r), me_r)
builtup_gain_class <- mask(crop(builtup_gain_class_reproj, me_r), me_r)
builtup_gain_year <- mask(crop(builtup_gain_year_reproj, me_r), me_r)

# Soil
oc <- mask(crop(oc_reproj, me_r), me_r)
bd <- mask(crop(bd_reproj, me_r), me_r)
clay <- mask(crop(clay_reproj, me_r), me_r)
sand <- mask(crop(sand_reproj, me_r), me_r)
silt <- mask(crop(silt_reproj, me_r), me_r)
ph <- mask(crop(ph_reproj, me_r), me_r)

# Soil Moisture
psssm_mean <- mask(crop(psssm_mean_reproj, me_r), me_r)
psssm_sd <- mask(crop(psssm_sd_reproj, me_r), me_r)
tsm_planting <- mask(crop(tsm_planting_reproj, me_r), me_r)
tsm_growing <- mask(crop(tsm_growing_reproj, me_r), me_r)

# Topography
dem <- mask(crop(dem_reproj, me_r), me_r)

# Climate
gddgrow_5_hist <- mask(crop(gddgrow_5_hist_reproj, me_r), me_r)
gddgrow_5_fut <- mask(crop(gddgrow_5_fut_reproj, me_r), me_r)
hddheat_18_hist <- mask(crop(hddheat_18_hist_reproj, me_r), me_r)
hddheat_18_fut <- mask(crop(hddheat_18_fut_reproj, me_r), me_r)
prtot_hist <- mask(crop(prtot_hist_reproj, me_r), me_r)
prtot_fut <- mask(crop(prtot_fut_reproj, me_r), me_r)
tn_mean_hist <- mask(crop(tn_mean_hist_reproj, me_r), me_r)
tn_mean_fut <- mask(crop(tn_mean_fut_reproj, me_r), me_r)
tx_mean_hist <- mask(crop(tx_mean_hist_reproj, me_r), me_r)
tx_mean_fut <- mask(crop(tx_mean_fut_reproj, me_r), me_r)
nr_cdd_hist <- mask(crop(nr_cdd_hist_reproj, me_r), me_r)
nr_cdd_fut <- mask(crop(nr_cdd_fut_reproj, me_r), me_r)

# Socio-economic
canfax <- st_intersection(canfax_reproj, me)

#--- Write to disk ----
# Pilot area
write_rast_tagged(lulc2020_pilot, "./processeddatasets", "landcover/aafc/cropped", aafc2020_filename, "pilot")
write_rast_tagged(lulc2015_pilot, "./processeddatasets", "landcover/aafc/cropped", aafc2015_filename, "pilot")
write_rast_tagged(crop_suit_pilot, "./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "pilot")


# Modelling extent
write_sf_tagged(census, "./processeddatasets", census_subdivision_filepath, census_subdivision_filename, "modelling-extent")
write_sf_tagged(ecozones, "./processeddatasets", ecozones_filepath, ecozones_filename, "modelling-extent")
write_sf_tagged(ecoregions, "./processeddatasets", ecoregions_filepath, ecoregions_filename, "modelling-extent")
write_sf_tagged(ecodistricts, "./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "modelling-extent")
write_sf_tagged(pa, "./processeddatasets", pa_filepath, pa_filename, "modelling-extent", driver = "GPKG")
write_sf_tagged(ncc, "./processeddatasets", ncc_easement_filepath, ncc_easement_filename, "modelling-extent", driver = "GPKG")
#write_sf_tagged(roads, "./processeddatasets", road_filepath, road_filename, "modelling-extent", layer = "roads_merge")
#write_rast_tagged(roads, "./processeddatasets", road_filepath, "roads.tif", "modelling-extent")
write_rast_tagged(crop_suit, "./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "modelling-extent")
write_rast_tagged(plowprint, "./processeddatasets", plowprint_filepath, plowprint_filename, "modelling-extent")
write_rast_tagged(lulc2020, file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent")
write_rast_tagged(lulc2015, file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent")
write_rast_tagged(lulc2010, file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent")
write_rast_tagged(lulc2000, file.path("./processeddatasets", aafc_filepath), "cropped", aafc2000_filename, "modelling-extent")
write_rast_tagged(grass_loss_class, "./processeddatasets", grass_loss_filepath, grass_loss_class_filename, "modelling-extent")
write_rast_tagged(grass_loss_year,  "./processeddatasets", grass_loss_filepath, grass_loss_year_filename,  "modelling-extent")
write_rast_tagged(ipast_loss_class, "./processeddatasets", ipast_loss_filepath, ipast_loss_class_filename, "modelling-extent")
write_rast_tagged(ipast_loss_year,  "./processeddatasets", ipast_loss_filepath, ipast_loss_year_filename,  "modelling-extent")
write_rast_tagged(upast_loss_class, "./processeddatasets", upast_loss_filepath, upast_loss_class_filename, "modelling-extent")
write_rast_tagged(upast_loss_year,  "./processeddatasets", upast_loss_filepath, upast_loss_year_filename,  "modelling-extent")
write_rast_tagged(crop_loss_class, "./processeddatasets", crop_loss_filepath, crop_loss_class_filename, "modelling-extent")
write_rast_tagged(crop_loss_year,  "./processeddatasets", crop_loss_filepath, crop_loss_year_filename,  "modelling-extent")
write_rast_tagged(crop_gain_class, "./processeddatasets", crop_gain_filepath, crop_gain_class_filename, "modelling-extent")
write_rast_tagged(crop_gain_year,  "./processeddatasets", crop_gain_filepath, crop_gain_year_filename,  "modelling-extent")
write_rast_tagged(builtup_gain_class, "./processeddatasets", builtup_gain_filepath, builtup_gain_class_filename, "modelling-extent")
write_rast_tagged(builtup_gain_year,  "./processeddatasets", builtup_gain_filepath, builtup_gain_year_filename,  "modelling-extent")
write_rast_tagged(oc, "./processeddatasets", soil_filepath, oc_filename, "modelling-extent", datatype="FLT4S")
write_rast_tagged(bd, "./processeddatasets", soil_filepath, bd_filename, "modelling-extent", datatype="FLT4S")
write_rast_tagged(clay, "./processeddatasets", soil_filepath, clay_filename, "modelling-extent", datatype="FLT4S")
write_rast_tagged(sand, "./processeddatasets", soil_filepath, sand_filename, "modelling-extent", datatype="FLT4S")
write_rast_tagged(silt, "./processeddatasets", soil_filepath, silt_filename, "modelling-extent", datatype="FLT4S")
write_rast_tagged(ph, "./processeddatasets", soil_filepath, ph_filename, "modelling-extent", datatype="FLT4S")
write_rast_tagged(dem, "./processeddatasets", dem_filepath, dem_filename, "modelling-extent")
write_rast_tagged(psssm_mean, "./processeddatasets", psssm_filepath, psssm_mean_filename, "modelling-extent")
write_rast_tagged(psssm_sd, "./processeddatasets", psssm_filepath, psssm_sd_filename, "modelling-extent")
write_rast_tagged(tsm_planting, "./processeddatasets", tsm_filepath, tsm_planting_filename, "modelling-extent")
write_rast_tagged(tsm_growing, "./processeddatasets", tsm_filepath, tsm_growing_filename, "modelling-extent")
write_rast_tagged(gddgrow_5_hist, "./processeddatasets", climate_filepath, gddgrow_5_hist_filename, "modelling-extent")
write_rast_tagged(gddgrow_5_fut, "./processeddatasets", climate_filepath, gddgrow_5_fut_filename, "modelling-extent")
write_rast_tagged(hddheat_18_hist, "./processeddatasets", climate_filepath, hddheat_18_hist_filename, "modelling-extent")
write_rast_tagged(hddheat_18_fut, "./processeddatasets", climate_filepath, hddheat_18_fut_filename, "modelling-extent")
write_rast_tagged(prtot_hist, "./processeddatasets", climate_filepath, prtot_hist_filename, "modelling-extent")
write_rast_tagged(prtot_fut, "./processeddatasets", climate_filepath, prtot_fut_filename, "modelling-extent")
write_rast_tagged(tn_mean_hist, "./processeddatasets", climate_filepath, tn_mean_hist_filename, "modelling-extent")
write_rast_tagged(tn_mean_fut, "./processeddatasets", climate_filepath, tn_mean_fut_filename, "modelling-extent")
write_rast_tagged(tx_mean_hist, "./processeddatasets", climate_filepath, tx_mean_hist_filename, "modelling-extent")
write_rast_tagged(tx_mean_fut, "./processeddatasets", climate_filepath, tx_mean_fut_filename, "modelling-extent")
write_rast_tagged(nr_cdd_hist, "./processeddatasets", climate_filepath, nr_cdd_hist_filename, "modelling-extent")
write_rast_tagged(nr_cdd_fut, "./processeddatasets", climate_filepath, nr_cdd_fut_filename, "modelling-extent")
write_sf_tagged(canfax, "./processeddatasets", canfax_filepath, canfax_filename, "modelling-extent")


#st_write(roads, file.path("./processeddatasets", pilot_filepath, "NRN_SK_14_0_GPKG_en_pilot_pad.gpkg"), layer = "NRN_SK_14_0_ROADSEG", delete_dsn = FALSE)


dist_to_road_m <- distance(roads)


#--- Extra ----

# Pilot raster template
pilot_r <- crop(pilot_r_prairies, pilot)
pilot_r <- mask(pilot_r, pilot)


# Protected areas
# This didn't work so I had to run an intersection manually in QGIS
# pa <- st_intersection(pa_reproj, pilot)
# gnp <- pa[pa$NAME_E == "Grasslands National Park Of Canada", ]

#--- Pad the pilot area for use with some spatial predictors 
# Pad the pilot extent by 20 cells
pad <- res(pilot_r) * 120
pilot_v <- vect(st_transform(pilot, crs(pilot)))
pilot_pad <- buffer(pilot_v, width = pad)
pilot_sf <- st_as_sf(pilot)
pilot_pad_sf <- st_as_sf(pilot_pad)
# Plot both padded and un-padded together
plot(st_geometry(pilot_pad_sf), col = NA, border = "red", lwd = 2)
plot(st_geometry(pilot_sf), col = NA, border = "blue", lwd = 2, add = TRUE)
