# =============================================================================
# Purpose:
#   Crops and masks all datasets to the modelling extent.
#   All outputs are saved to processeddatasets/ with "modelling-extent" tag.
#   Also creates binary land cover rasters from AAFC LUTS
#
# Inputs:
#   - Modelling extent (shapefile and raster template)
#   - All reprojected datasets from 2-project-data.R
#
# Outputs:
#   - All datasets cropped/masked to modelling extent
#   - Binary land cover rasters derived from AAFC LUTS
#
# Notes:
#   - pa uses st_crop() to bounding box rather than st_intersection()
#     for speed; NA handling deferred to rasterization step
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

#--- Load files ----
# Study area
sa <- st_read(file.path(sa_filepath, sa_filename))
sa_r <- rast(file.path(sa_filepath, sa_raster_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Cropland
crop_suit_reproj<- rast(tagged_path("./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "proj-naaea"))
plowprint_reproj<- rast(tagged_path("./processeddatasets", plowprint_filepath, plowprint_filename, "proj-naaea"))

# Census subdivisions
census_reproj <- st_read(tagged_path("./processeddatasets", census_subdivision_filepath, census_subdivision_filename, "proj-naaea"))

# National ecological framework
#ecozones_reproj <- st_read(file.path("./processeddatasets", ecozones_filepath, ecozones_filename))
ecoregions_reproj <- st_read(file.path("./processeddatasets", ecoregions_filepath, ecoregions_filename))
ecodistricts_reproj <- st_read(file.path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename))

# Protected Areas
pa_reproj <- st_read(
  tagged_path("./processeddatasets", pa_filepath, pa_filename, "proj-naaea", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(pa_filename), "_proj-naaea")
)

sk_cce_reproj <- st_read(
  tagged_path("./processeddatasets", sk_cce_filepath, sk_cce_filename, "proj-naaea", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(sk_cce_filename), "_proj-naaea")
)

# Land tenure
ncc_reproj <- st_read(
  tagged_path("./processeddatasets", ncc_easement_filepath, ncc_easement_filename,
              "proj-naaea", driver = "GPKG"),
  layer = paste0(tools::file_path_sans_ext(ncc_easement_filename), "_proj-naaea")
)

# Roads
roads_reproj <-rast(tagged_path("./processeddatasets", road_filepath, paste0(road_filename, ".tif"), "proj-naaea"))

# LULC - AAFC LUTS
lulc2020_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2020_filename, "proj-naaea"))
lulc2015_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2015_filename, "proj-naaea"))
lulc2010_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2010_filename, "proj-naaea"))
lulc2000_prairies_reproj <- rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2000_filename, "proj-naaea"))

# Grassland
lulc2025 <- rast(file.path("./processeddatasets", grass_consensus_filepath, grass_consensus_filename))

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
# Climate (imputed fut, projection use only)
gddgrow_5_fut_imputed_reproj  <- rast(tagged_path("./processeddatasets", climate_filepath, gddgrow_5_fut_filename,  "imputed_proj-naaea"))
hddheat_18_fut_imputed_reproj <- rast(tagged_path("./processeddatasets", climate_filepath, hddheat_18_fut_filename, "imputed_proj-naaea"))
prtot_fut_imputed_reproj      <- rast(tagged_path("./processeddatasets", climate_filepath, prtot_fut_filename,      "imputed_proj-naaea"))
tn_mean_fut_imputed_reproj    <- rast(tagged_path("./processeddatasets", climate_filepath, tn_mean_fut_filename,    "imputed_proj-naaea"))
tx_mean_fut_imputed_reproj    <- rast(tagged_path("./processeddatasets", climate_filepath, tx_mean_fut_filename,    "imputed_proj-naaea"))
nr_cdd_fut_imputed_reproj     <- rast(tagged_path("./processeddatasets", climate_filepath, nr_cdd_fut_filename,     "imputed_proj-naaea"))

# Socio-economic
canfax2009_reproj <- st_read(tagged_path("./processeddatasets", canfax_filepath, canfax2009_filename, "proj-naaea"))
canfax2024_reproj <- st_read(tagged_path("./processeddatasets", canfax_filepath, canfax2024_filename, "proj-naaea"))
elevator_reproj <- st_read(tagged_path("./processeddatasets", elevator_filepath, elevator_filename, "proj-naaea"))
popctr2011_reproj <- st_read(tagged_path("./processeddatasets", popctr_filepath, paste0(popctr2011_filename, ".shp"), "proj-naaea"))
popctr2021_reproj <- st_read(tagged_path("./processeddatasets", popctr_filepath, paste0(popctr2021_filename, ".shp"), "proj-naaea"))
db2011_reproj <- st_read(tagged_path("./processeddatasets", popdens_filepath, paste0(db2011_filename, ".shp"), "proj-naaea"))
db2021_reproj <- st_read(tagged_path("./processeddatasets", popdens_filepath, paste0(db2021_filename, ".shp"), "proj-naaea"))

#--- Crop and mask to modelling extent ----
# Census subdivisions
census <- st_intersection(census_reproj, me)
# plot(census["CSDNAME"])

# National ecological framework
# ecozones <- st_intersection(ecozones_reproj, me)
ecoregions <- st_intersection(ecoregions_reproj, me)
ecoregions_r <- terra::rasterize(terra::vect(ecoregions), me_r, field = "ECOREGION_") # also create a raster version of ecoregions
ecodistricts <- st_intersection(ecodistricts_reproj, me)
ecodistricts_r <- terra::rasterize(terra::vect(ecodistricts), me_r, field = "ECODISTRIC") # also create a raster version of ecodistricts

# Protected Areas
# CPCAD
# using st_crop() and then st_intersection() for speed, will rasterize later 
# clean geometries
pa <- st_crop(pa_reproj, st_bbox(me_r))
pa <- st_intersection(pa, me)
pa <- st_collection_extract(pa, "POLYGON")
pa <- st_cast(pa, "MULTIPOLYGON")

# cropping to study area once here because it takes a while to run and the 
# file size on disk is small
pa_sa <- st_intersection(pa, sa)
pa_sa <- st_collection_extract(pa_sa, "POLYGON")
pa_sa <- st_cast(pa_sa, "MULTIPOLYGON")

# SK Crown Conservation Easements
sk_cce <- st_crop(sk_cce_reproj, st_bbox(me_r))
sk_cce <- st_intersection(sk_cce, me)
sk_cce <- st_cast(sk_cce, "MULTIPOLYGON")


# Land tenure
ncc <- st_crop(ncc_reproj, st_bbox(me_r))
ncc <- st_intersection(ncc, me)
ncc <- st_collection_extract(ncc, "POLYGON")
ncc <- st_cast(ncc, "MULTIPOLYGON")

# cropping to study area once here because it takes a while to run and the 
# file size on disk is small
ncc_sa <- st_intersection(ncc, sa)
ncc_sa <- st_collection_extract(ncc_sa, "POLYGON")
ncc_sa <- st_cast(ncc_sa, "MULTIPOLYGON")

# Roads
roads <- mask(crop(roads_reproj, me_r), me_r)

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
# Climate (imputed fut)
gddgrow_5_fut_imputed  <- mask(crop(gddgrow_5_fut_imputed_reproj,  me_r), me_r)
hddheat_18_fut_imputed <- mask(crop(hddheat_18_fut_imputed_reproj, me_r), me_r)
prtot_fut_imputed      <- mask(crop(prtot_fut_imputed_reproj,      me_r), me_r)
tn_mean_fut_imputed    <- mask(crop(tn_mean_fut_imputed_reproj,    me_r), me_r)
tx_mean_fut_imputed    <- mask(crop(tx_mean_fut_imputed_reproj,    me_r), me_r)
nr_cdd_fut_imputed     <- mask(crop(nr_cdd_fut_imputed_reproj,     me_r), me_r)

# check: imputed layers should have no NA inside the modelling extent (expect 0)
for (x in list(gddgrow_5_fut_imputed, hddheat_18_fut_imputed, prtot_fut_imputed,
               tn_mean_fut_imputed, tx_mean_fut_imputed, nr_cdd_fut_imputed))
  print(global(is.na(x) & !is.na(me_r), "sum", na.rm = TRUE))

# Socio-economic
canfax2009 <- st_intersection(canfax2009_reproj, me)
canfax2024 <- st_intersection(canfax2024_reproj, me)
elevator <- st_intersection(elevator_reproj, me)
popctr2011 <- st_intersection(popctr2011_reproj, me)
popctr2021 <- st_intersection(popctr2021_reproj, me)
# me's dense outline made [me, ] run for minutes, and rasterize/mask to me_r later drops anything outside the exact exten; 
# simplify at 1 km, then buffer by the same 1 km so the simplified shape still contains me (simplification moves the boundary at most dTolerance)
me_s <- st_buffer(st_simplify(me, dTolerance = 1000), 1000)
db2011 <- db2011_reproj[me_s, ]
db2021 <- db2021_reproj[me_s, ]


#--- Write to disk ----
write_sf_tagged(census, "./processeddatasets", census_subdivision_filepath, census_subdivision_filename, "modelling-extent")
write_sf_tagged(ecozones, "./processeddatasets", ecozones_filepath, ecozones_filename, "modelling-extent")
write_sf_tagged(ecoregions, "./processeddatasets", ecoregions_filepath, ecoregions_filename, "modelling-extent")
write_rast_tagged(ecoregions_r, "./processeddatasets", ecoregions_filepath, paste0(tools::file_path_sans_ext(ecoregions_filename), ".tif"), "modelling-extent", datatype = "INT2U")
write_sf_tagged(ecodistricts, "./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "modelling-extent")
write_rast_tagged(ecodistricts_r, "./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "modelling-extent", datatype = "INT2U")
write_sf_tagged(pa, "./processeddatasets", pa_filepath, pa_filename, "modelling-extent", driver = "GPKG")
write_sf_tagged(pa_sa, "./processeddatasets", pa_filepath, pa_filename, "study-area", driver = "GPKG")
write_sf_tagged(sk_cce, "./processeddatasets", sk_cce_filepath, sk_cce_filename, "modelling-extent", driver = "GPKG")
write_sf_tagged(ncc, "./processeddatasets", ncc_easement_filepath, ncc_easement_filename, "modelling-extent", driver = "GPKG")
write_sf_tagged(ncc_sa, "./processeddatasets", ncc_easement_filepath, ncc_easement_filename, "study-area", driver = "GPKG")
write_rast_tagged(roads, "./processeddatasets", road_filepath, paste0(road_filename, ".tif"), "modelling-extent")
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
write_rast_tagged(gddgrow_5_fut_imputed,  "./processeddatasets", climate_filepath, gddgrow_5_fut_filename,  "imputed_modelling-extent")
write_rast_tagged(hddheat_18_fut_imputed, "./processeddatasets", climate_filepath, hddheat_18_fut_filename, "imputed_modelling-extent")
write_rast_tagged(prtot_fut_imputed,      "./processeddatasets", climate_filepath, prtot_fut_filename,      "imputed_modelling-extent")
write_rast_tagged(tn_mean_fut_imputed,    "./processeddatasets", climate_filepath, tn_mean_fut_filename,    "imputed_modelling-extent")
write_rast_tagged(tx_mean_fut_imputed,    "./processeddatasets", climate_filepath, tx_mean_fut_filename,    "imputed_modelling-extent")
write_rast_tagged(nr_cdd_fut_imputed,     "./processeddatasets", climate_filepath, nr_cdd_fut_filename,     "imputed_modelling-extent")
write_sf_tagged(canfax2009, "./processeddatasets", canfax_filepath, canfax2009_filename, "modelling-extent")
write_sf_tagged(canfax2024, "./processeddatasets", canfax_filepath, canfax2024_filename, "modelling-extent")
write_sf_tagged(elevator, "./processeddatasets", elevator_filepath, elevator_filename, "modelling-extent")
write_sf_tagged(popctr2011, "./processeddatasets", popctr_filepath, paste0(popctr2011_filename, ".shp"), "modelling-extent")
write_sf_tagged(popctr2021, "./processeddatasets", popctr_filepath, paste0(popctr2021_filename, ".shp"), "modelling-extent")
write_sf_tagged(db2011, "./processeddatasets", popdens_filepath, paste0(db2011_filename, ".shp"), "modelling-extent")
write_sf_tagged(db2021, "./processeddatasets", popdens_filepath, paste0(db2021_filename, ".shp"), "modelling-extent")

#--- Create binary layers from AAFC LUTS rasters ----
lulc <- list(
  `2010` = rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "modelling-extent")),
  `2015` = rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2015_filename, "modelling-extent")),
  `2020` = rast(tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "modelling-extent"))
)

classes <- list(
  `2010` = c("native", "tame&crop", "settlement", "forest", "water"),
  `2015` = c("native", "tame", "crop", "settlement", "forest", "water"),
  `2020` = c("native", "tame", "crop", "settlement", "forest", "water")
)

# LULC binary layers
bin_layers <- list()

# Reclass to binary layer and save result to disk
for (yr in names(lulc)) {
  bin_layers[[yr]] <- list()
  for (cls in classes[[yr]]) {
    r <- reclass_aafc(lulc[[yr]], scheme = "binary", binary_class = cls, check = TRUE)
    writeRaster(
      r, file.path(out_dir, sprintf("AAFC%s_%s.tif", yr, cls)),
      overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
    )
    bin_layers[[yr]][[cls]] <- r
  }
}

#--- Crop select layers to study area ----
# Binary LUTS layers for 2010 and 2020
classes <- list(
  `2010` = c("native", "tame&crop", "settlement", "forest", "water"),
  `2020` = c("native", "tame", "crop", "settlement", "forest", "water")
)

for (yr in names(classes)) for (cl in classes[[yr]]) {
  f <- file.path("./processeddatasets", aafc_filepath, "binary", sprintf("AAFC%s_%s.tif", yr, cl))
  writeRaster(mask(crop(rast(f), sa_r), sa_r), sub("\\.tif$", "_studyarea.tif", f),
              overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))   # 0/1 fits in INT1U, NA flag 255
}

# LUTS
mask(crop(lulc2010_prairies_reproj, sa_r), sa_r, filename = tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2010_filename, "studyarea"))
mask(crop(lulc2020_prairies_reproj, sa_r), sa_r, filename = tagged_path(file.path("./processeddatasets", aafc_filepath), "cropped", aafc2020_filename, "studyarea"))

# Grassland
mask(crop(lulc2025, sa_r), sa_r, filename = tagged_path("./processeddatasets", grass_consensus_filepath, grass_consensus_filename, "studyarea"))

# Boundaries
ecoregions_sa <- ecoregions[ecoregions$ECOZONE_ID == 10,]   # Prairies Ecozone
write_sf_tagged(ecoregions_sa, "./processeddatasets", ecoregions_filepath, ecoregions_filename, "studyarea")
terra::rasterize(terra::vect(ecoregions_sa), me_r, field = "ECOREGION_", filename = tagged_path("./processeddatasets", ecoregions_filepath, paste0(tools::file_path_sans_ext(ecoregions_filename), ".tif"), "studyarea"), overwrite = TRUE)

ecodistricts_sa <- ecodistricts[ecodistricts$ECOZONE_ID == 10,]   # Prairies Ecozone
write_sf_tagged(ecodistricts_sa, "./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "studyarea")
terra::rasterize(terra::vect(ecodistricts_sa), me_r, field = "ECODISTRIC", filename = tagged_path("./processeddatasets", ecodistricts_filepath, paste0(tools::file_path_sans_ext(ecodistricts_filename), ".tif"), "studyarea"), overwrite = TRUE)



