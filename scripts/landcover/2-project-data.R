# =============================================================================
# 2-project-data.R
# =============================================================================
#
# Purpose:
#   Reprojects all raw datasets to the common CRS (NAE Albers Equal Area).
#   Also applies scaling factors and NoData values to soil rasters.
#   Slow to run — only needs to be re-run if raw inputs change.
#
# Inputs:
#   - Modelling extent raster template (me_r) defines target CRS and raster grid
#   - Raw datasets from existingdatasets/:
#       Olimb & Robinson 2019 cropland suitability
#       Census subdivisions
#       National ecological framework (ecozones, ecoregions, ecodistricts)
#       AAFC LUTS 2000, 2010, 2015, 2020 (merged)
#       ECCC most recent change layers (grass, improved pasture,
#         unimproved pasture, crop, built-up)
#       Protected areas (CPCAD 2024)
#       Soil (OC, BD, clay, sand, silt, pH)
#       DEM (NRCan MRDEM 30m)
#
# Outputs:
#   - All datasets reprojected to NAE Albers Equal Area
#   - Soil rasters unscaled to real units and NoData values applied
#   - Protected areas cast to MULTIPOLYGON and validated
#   - All outputs saved to processeddatasets/ with "proj-naaea" tag
#
# Notes:
#   - Raster reprojection uses method = "near" for categorical data (LULC,
#     change layers) and method = "bilinear" for continuous data (soil, DEM)
#   - Roads processing is incomplete and commented out pending resolution
#     of the merged road layer workflow
#   - Ecoregions and ecodistricts were preprocessed manually in QGIS to
#     resolve geometry misalignment issues before being read in here
#   - CPCAD geometry cleaning (st_cast, st_make_valid) is done here
#     before reprojection to avoid propagating invalid geometries
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
#   - 1-merge-tiled-data.R outputs (AAFC LUTS merged rasters)
# =============================================================================

source(file.path("scripts/landcover", "0a-setup.R"))

#--- Load files ----
# Study area & Modelling extent
sa <- st_read(file.path(sa_filepath, sa_filename))
me <- st_read(file.path(sa_filepath, me_filename))
me_r <- rast(file.path(sa_filepath, me_raster_filename))

# Cropland
crop_suit_raw <- rast(file.path("./existingdatasets", capitalize_first(crop_suit_2019_filepath), crop_suit_2019_filename))
plowprint_raw <- rast(file.path("./existingdatasets", capitalize_first(plowprint_filepath), plowprint_filename))

# Census subdivisions
census_raw <- st_read(file.path("./existingdatasets", capitalize_first(census_subdivision_filepath), census_subdivision_filename))

# National ecological framework
ecozones_raw <- st_read(file.path("./existingdatasets", capitalize_first(ecozones_filepath), ecozones_filename))
#Note that ecoregions and ecodistrics were preprocessed manually to deal with misalignment
ecoregions_raw <- st_read(file.path("./processeddatasets", ecoregions_filepath, ecoregions_filename))
ecodistricts_raw <- st_read(file.path("./processeddatasets", ecodistricts_filepath, ecodistricts_filename))

# Protected Areas
pa_raw <- st_read(file.path("./existingdatasets", capitalize_first(pa_filepath), paste0(pa_filename, ".gdb")), layer = pa_filename)

# Land tenure
ncc_gdb <- file.path("./existingdatasets", capitalize_first(ncc_easement_filepath), paste0(ncc_easement_filename, ".gdb"))
lyr_names <- st_layers(ncc_gdb)$name
ncc_list <- lapply(lyr_names, \(x) st_read(ncc_gdb, layer = x, quiet = TRUE))
# lapply(ncc_list, names)  # check names to make sure the layers all have the same schema (they do! but their crs are different)

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

# Soil
oc_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), oc_filename))
bd_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), bd_filename))
clay_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), clay_filename))
sand_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), sand_filename))
silt_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), silt_filename))
ph_raw <- rast(file.path("./existingdatasets", capitalize_first(soil_filepath), ph_filename))

# Soil Moisture
psssm_mean_raw <- rast(file.path("./existingdatasets", psssm_filepath, psssm_mean_filename))
psssm_sd_raw <- rast(file.path("./existingdatasets", psssm_filepath, psssm_sd_filename))
tsm_planting_raw <- rast(file.path("./existingdatasets", tsm_filepath, tsm_planting_filename))
tsm_growing_raw <- rast(file.path("./existingdatasets", tsm_filepath, tsm_growing_filename))

# Topography
dem_raw <- rast(file.path("./existingdatasets", capitalize_first(dem_filepath), dem_filename))

# Climate
gddgrow_5_hist_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), gddgrow_5_hist_filename))
gddgrow_5_fut_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), gddgrow_5_fut_filename))
hddheat_18_hist_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), hddheat_18_hist_filename))
hddheat_18_fut_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), hddheat_18_fut_filename))
prtot_hist_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), prtot_hist_filename))
prtot_fut_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), prtot_fut_filename))
tn_mean_hist_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), tn_mean_hist_filename))
tn_mean_fut_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), tn_mean_fut_filename))
tx_mean_hist_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), tx_mean_hist_filename))
tx_mean_fut_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), tx_mean_fut_filename))
nr_cdd_hist_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), nr_cdd_hist_filename))
nr_cdd_fut_raw <- rast(file.path("./existingdatasets", capitalize_first(climate_filepath), nr_cdd_fut_filename))

# Socio-economic
# Roads
# Merged road layer is huge so only read in the portion within the modelling extent bounding box
# read just the layer definition/CRS without loading everything
roads_crs <- st_crs(st_read(file.path("./existingdatasets", capitalize_first(road_filepath), paste0(road_filename, ".gdb")),
                            layer="roads_merge", quiet=TRUE, query="SELECT * FROM roads_merge LIMIT 1"))

bb_poly <- st_as_sfc(st_bbox(me))
bb_poly_roads <- st_transform(bb_poly, roads_crs)

roads_raw <- st_read(file.path("./existingdatasets", capitalize_first(road_filepath), paste0(road_filename, ".gdb")),
                             layer="roads_merge",
                             wkt_filter = st_as_text(bb_poly_roads),
                             quiet=TRUE)

# Canfax net returns
canfax_raw <- st_read(file.path("./existingdatasets", capitalize_first(canfax_filepath), canfax_filename))
elevator_raw <- st_read(file.path("./existingdatasets", capitalize_first(elevator_filepath), elevator_filename))
popctr2011_raw <- st_read(file.path("./existingdatasets", capitalize_first(popctr_filepath), paste0(popctr2011_filename, ".shp")))
popctr2021_raw <- st_read(file.path("./existingdatasets", capitalize_first(popctr_filepath), paste0(popctr2021_filename, ".shp")))

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

#--- clean PA data ----
# CPCAD
# Cast to multipolygon to eliminate exotic geometry types
pa_cast <- st_cast(pa_raw, "MULTIPOLYGON")

# Make valid
pa_valid <- st_make_valid(pa_cast)

# Verify
sum(!st_is_valid(pa_valid))

# Check geometry types are now all standard
table(st_geometry_type(pa_valid))

# Check how many are still invalid
sum(!st_is_valid(pa_valid))

# NCC easements
ncc_all <- Map(\(x, nm) x |> mutate(src_layer = nm) |> st_transform(crs(me_r)),
               ncc_list, lyr_names) |>
  bind_rows() |>
  st_cast("MULTIPOLYGON") |>
  st_make_valid()


#--- Reproject to match study area shapefile CRS ----
crop_suit_reproj <- project(crop_suit_raw, me_r, method = "bilinear")
plowprint_reproj <- project(plowprint_raw, me_r, method = "near")
lulc2020_prairies_reproj <- project(lulc2020_prairies_raw, me_r, method = "near")
lulc2015_prairies_reproj <- project(lulc2015_prairies_raw, me_r, method = "near")
lulc2010_prairies_reproj <- project(lulc2010_prairies_raw, me_r, method = "near")
lulc2000_prairies_reproj <- project(lulc2000_prairies_raw, me_r, method = "near")
census_reproj <- st_transform(census_raw, st_crs(sa))
ecozones_reproj <- st_transform(ecozones_raw, st_crs(sa))
ecoregions_reproj <- st_transform(ecoregions_raw, st_crs(sa))
ecodistricts_reproj <- st_transform(ecodistricts_raw, st_crs(sa))

grass_loss_class_reproj <- project(grass_loss_class_raw, me_r, method = "near")
grass_loss_year_reproj <- project(grass_loss_year_raw, me_r, method = "near")
ipast_loss_class_reproj <- project(ipast_loss_class_raw, me_r, method = "near")
ipast_loss_year_reproj <- project(ipast_loss_year_raw, me_r, method = "near")
upast_loss_class_reproj <- project(upast_loss_class_raw, me_r, method = "near")
upast_loss_year_reproj <- project(upast_loss_year_raw, me_r, method = "near")
crop_loss_class_reproj <- project(crop_loss_class_raw, me_r, method = "near")
crop_loss_year_reproj <- project(crop_loss_year_raw, me_r, method = "near")
crop_gain_class_reproj <- project(crop_gain_class_raw, me_r, method = "near")
crop_gain_year_reproj <- project(crop_gain_year_raw, me_r, method = "near")
builtup_gain_class_reproj <- project(builtup_gain_class_raw, me_r, method = "near")
builtup_gain_year_reproj <- project(builtup_gain_year_raw, me_r, method = "near")
oc_reproj <- project(oc_real, me_r, method = "bilinear")
bd_reproj <- project(bd_real, me_r, method = "bilinear")
clay_reproj <- project(clay_raw, me_r, method = "bilinear")
sand_reproj <- project(sand_raw, me_r, method = "bilinear")
silt_reproj <- project(silt_raw, me_r, method = "bilinear")
ph_reproj <- project(ph_real, me_r, method = "bilinear")
dem_reproj <- project(dem_raw, me_r, method = "bilinear")
psssm_mean_reproj <- project(psssm_mean_raw, me_r, method = "bilinear")
psssm_sd_reproj <- project(psssm_sd_raw, me_r, method = "bilinear")
tsm_planting_reproj <- project(tsm_planting_raw, me_r, method = "bilinear")
tsm_growing_reproj <- project(tsm_growing_raw, me_r, method = "bilinear")
gddgrow_5_hist_reproj <- project(gddgrow_5_hist_raw, me_r, method = "bilinear")
gddgrow_5_fut_reproj <- project(gddgrow_5_fut_raw, me_r, method = "bilinear")
hddheat_18_hist_reproj <- project(hddheat_18_hist_raw, me_r, method = "bilinear")
hddheat_18_fut_reproj <- project(hddheat_18_fut_raw, me_r, method = "bilinear")
prtot_hist_reproj <- project(prtot_hist_raw, me_r, method = "bilinear")
prtot_fut_reproj <- project(prtot_fut_raw, me_r, method = "bilinear")
tn_mean_hist_reproj <- project(tn_mean_hist_raw, me_r, method = "bilinear")
tn_mean_fut_reproj <- project(tn_mean_fut_raw, me_r, method = "bilinear")
tx_mean_hist_reproj <- project(tx_mean_hist_raw, me_r, method = "bilinear")
tx_mean_fut_reproj <- project(tx_mean_fut_raw, me_r, method = "bilinear")
nr_cdd_hist_reproj <- project(nr_cdd_hist_raw, me_r, method = "bilinear")
nr_cdd_fut_reproj <- project(nr_cdd_fut_raw, me_r, method = "bilinear")
pa_reproj <- st_transform(pa_valid, crs = crs(me_r))
canfax_reproj <- st_transform(canfax_raw, st_crs(sa))
elevator_reproj <- st_transform(elevator_raw, st_crs(sa))
popctr2011_reproj <- st_transform(popctr2011_raw, st_crs(sa))
popctr2021_reproj <- st_transform(popctr2021_raw, st_crs(sa))

roads_reproj <- st_transform(roads_raw, st_crs(sa))
roads_r <- rasterize(
  roads_reproj, me_r,
  field = 1,
  background = NA,
  touches = TRUE
)

#--- Write to disk ----
write_rast_tagged(crop_suit_reproj, "./processeddatasets", crop_suit_2019_filepath, crop_suit_2019_filename, "proj-naaea")
write_rast_tagged(plowprint_reproj, "./processeddatasets", plowprint_filepath, plowprint_filename, "proj-naaea")
write_sf_tagged(census_reproj, "./processeddatasets", census_subdivision_filepath, census_subdivision_filename, "proj-naaea")
write_sf_tagged(ecozones_reproj, "./processeddatasets", ecozones_filepath, ecozones_filename, "proj-naaea")
write_sf_tagged(ecoregions_reproj, "./processeddatasets", ecoregions_filepath, ecoregions_filename, "proj-naaea")
write_sf_tagged(ecodistricts_reproj, "./processeddatasets", ecodistricts_filepath, ecodistricts_filename, "proj-naaea")
write_sf_tagged(pa_reproj, "./processeddatasets", pa_filepath, pa_filename, "proj-naaea", driver = "GPKG")
write_sf_tagged(ncc_all, "./processeddatasets", ncc_easement_filepath, ncc_easement_filename, "proj-naaea", driver = "GPKG")
#write_sf_tagged(roads_reproj, "./processeddatasets", road_filepath, road_filename, "proj-naaea", layer = "roads_merge") #skip and write raster instead
write_rast_tagged(roads_r, "./processeddatasets", road_filepath, paste0(road_filename, ".tif"), "proj-naaea")
write_rast_tagged(lulc2020_prairies_reproj, file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2020_filename, "proj-naaea")
write_rast_tagged(lulc2015_prairies_reproj, file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2015_filename, "proj-naaea")
write_rast_tagged(lulc2010_prairies_reproj, file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2010_filename, "proj-naaea")
write_rast_tagged(lulc2000_prairies_reproj, file.path("./processeddatasets", aafc_filepath), "reprojected", aafc2000_filename, "proj-naaea")
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
write_rast_tagged(dem_reproj, "./processeddatasets", dem_filepath, dem_filename, "proj-naaea")
write_rast_tagged(psssm_mean_reproj, "./processeddatasets", psssm_filepath, psssm_mean_filename, "proj-naaea")
write_rast_tagged(psssm_sd_reproj, "./processeddatasets", psssm_filepath, psssm_sd_filename, "proj-naaea")
write_rast_tagged(tsm_planting_reproj, "./processeddatasets", tsm_filepath, tsm_planting_filename, "proj-naaea")
write_rast_tagged(tsm_growing_reproj, "./processeddatasets", tsm_filepath, tsm_growing_filename, "proj-naaea")
write_rast_tagged(gddgrow_5_hist_reproj, "./processeddatasets", climate_filepath, gddgrow_5_hist_filename, "proj-naaea")
write_rast_tagged(gddgrow_5_fut_reproj, "./processeddatasets", climate_filepath, gddgrow_5_fut_filename, "proj-naaea")
write_rast_tagged(hddheat_18_hist_reproj, "./processeddatasets", climate_filepath, hddheat_18_hist_filename, "proj-naaea")
write_rast_tagged(hddheat_18_fut_reproj, "./processeddatasets", climate_filepath, hddheat_18_fut_filename, "proj-naaea")
write_rast_tagged(prtot_hist_reproj, "./processeddatasets", climate_filepath, prtot_hist_filename, "proj-naaea")
write_rast_tagged(prtot_fut_reproj, "./processeddatasets", climate_filepath, prtot_fut_filename, "proj-naaea")
write_rast_tagged(tn_mean_hist_reproj, "./processeddatasets", climate_filepath, tn_mean_hist_filename, "proj-naaea")
write_rast_tagged(tn_mean_fut_reproj, "./processeddatasets", climate_filepath, tn_mean_fut_filename, "proj-naaea")
write_rast_tagged(tx_mean_hist_reproj, "./processeddatasets", climate_filepath, tx_mean_hist_filename, "proj-naaea")
write_rast_tagged(tx_mean_fut_reproj, "./processeddatasets", climate_filepath, tx_mean_fut_filename, "proj-naaea")
write_rast_tagged(nr_cdd_hist_reproj, "./processeddatasets", climate_filepath, nr_cdd_hist_filename, "proj-naaea")
write_rast_tagged(nr_cdd_fut_reproj, "./processeddatasets", climate_filepath, nr_cdd_fut_filename, "proj-naaea")
write_sf_tagged(canfax_reproj, "./processeddatasets", canfax_filepath, canfax_filename, "proj-naaea")
write_sf_tagged(elevator_reproj, "./processeddatasets", elevator_filepath, elevator_filename, "proj-naaea")
write_sf_tagged(popctr2011_reproj, "./processeddatasets", popctr_filepath, popctr2011_filename, "proj-naaea")
write_sf_tagged(popctr2021_reproj, "./processeddatasets", popctr_filepath, popctr2021_filename, "proj-naaea")
