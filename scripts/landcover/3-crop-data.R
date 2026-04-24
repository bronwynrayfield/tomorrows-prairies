# Data pre-processing -> cropping

# Setup
source(file.path("scripts/landcover", "0-setup.R"))
source(file.path("scripts/landcover/functions", "0-functions-general.R"))

#--- Load files ----
# Study area
pilot <- st_read(file.path("./processeddatasets", pilot_filepath, pilot_filename))
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Olimb & Robinson 2019 cropland suitability
crop_suit_reproj<- rast(tagged_path("./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "proj-naaea"))

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

# Soil
oc_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, oc_filename, "proj-naaea"))
bd_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, bd_filename, "proj-naaea"))
clay_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, clay_filename, "proj-naaea"))
sand_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, sand_filename, "proj-naaea"))
silt_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, silt_filename, "proj-naaea"))
ph_reproj <- rast(tagged_path("./processeddatasets", soil_filepath, ph_filename, "proj-naaea"))

# Topography
dem_reproj <- rast(tagged_path("./processeddatasets", dem_filepath, dem_filename, "proj-naaea"))

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

# Roads
roads <- st_intersection(roads_reproj, me)
plot(roads["ROADCLASS"], main = "Roads by Class", key.pos = 4)
roads <- mask(crop(roads_r, me), me)

# Cropland suitability
crop_suit <- mask(crop(crop_suit_reproj, me_r), me_r)
plot(crop_suit, main = "Cropland suitability (2019)")

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

# Topography
dem <- mask(crop(dem_reproj, me_r), me_r)
slope  <- terrain(dem, "slope",  unit="degrees")
aspect <- terrain(dem, "aspect", unit="degrees")
aspect[slope < 0.5] <- NA


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
#write_sf_tagged(roads, "./processeddatasets", road_filepath, road_filename, "modelling-extent", layer = "roads_merge")
#write_rast_tagged(roads, "./processeddatasets", road_filepath, "roads.tif", "modelling-extent")
write_rast_tagged(crop_suit, "./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "modelling-extent")
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
write_rast_tagged(slope, "./processeddatasets", dem_filepath, slope_filename, "modelling-extent")
write_rast_tagged(aspect, "./processeddatasets", dem_filepath, aspect_filename, "modelling-extent")


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
