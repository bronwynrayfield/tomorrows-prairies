
###################################################################################
##           PREPROCESSING: ALIGNING STUDY AREAS AND PILOT AREA                  ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 2: align all study areas and pilot area with the reference raster
#         created in step 1

# load packages
library(terra)
library(sf)
library(tidyverse)

# set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "./processeddatasets/pilot_area"

# study area
sa_dir <- "./processeddatasets/Studyarea"

####--- FILES ---####

# pilot area
pa_file <- file.path(pa_dir, "GNP_pilot_area_crop_CSD.shp")

# study area with which all other layers will be aligned
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")
  
# modelling study area
mod_sa_file <- file.path(sa_dir, "TP_study_area_prairies_ecozone.shp")

# modelling extent
mod_ext_file <- file.path(sa_dir, "TP_study_area_ecodistrict_buffer.shp")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

##############################################

# load reference raster
sa_rast <- rast(sa_file)

# print raster details 
origin(sa_rast)
res(sa_rast)
crs(sa_rast)

# visualize reference raster
plot(sa_rast)

# load study area shapefiles
mod_sa_vect <- vect(mod_sa_file)
mod_ext_vect <- vect(mod_ext_file)

# Reproject study areas to match reference raster CRS
mod_sa_vect_proj <- project(mod_sa_vect, crs(sa_rast))
mod_ext_vect_proj <- project(mod_ext_vect, crs(sa_rast))

# Rasterize study areas using reference raster as template
mod_sa_rast <- rasterize(mod_sa_vect_proj, sa_rast, field=1, touches=TRUE) # field=1 assigns constant value to rasterized cells
mod_ext_rast <- rasterize(mod_ext_vect_proj, sa_rast, field=1, touches=TRUE) # touches=TRUE ensures that any cell touched by the polygon is included

# print raster details 
origin(mod_sa_rast)
res(mod_sa_rast)
crs(mod_sa_rast)

# print raster details 
origin(mod_ext_rast)
res(mod_ext_rast)
crs(mod_ext_rast)

# visualize the new rasters
plot(mod_sa_rast)
plot(mod_ext_rast)

# Trim NA-only rows/columns
mod_sa_rast_trimmed <- trim(mod_sa_rast)
mod_ext_rast_trimmed <- trim(mod_ext_rast)

# plot again
plot(mod_sa_rast_trimmed)
plot(mod_ext_rast_trimmed)

# Check extents - they differ...
#ext(sa_rast)
#ext(sa_rast_trimmed)

#nrow(sa_rast); ncol(sa_rast)
#nrow(sa_rast_trimmed); ncol(sa_rast_trimmed)

# but origin the same!
#origin(sa_rast)
#res(sa_rast)
#crs(sa_rast)

# print raster details 
origin(mod_sa_rast_trimmed)
res(mod_sa_rast_trimmed)
crs(mod_sa_rast_trimmed)

# print raster details 
origin(mod_ext_rast_trimmed)
res(mod_ext_rast_trimmed)
crs(mod_ext_rast_trimmed)

# Save the raster as 8-bit unsigned integer
mod_sa_out_path <- file.path(sa_dir, "modelling_study_area_aligned_cnvc_Jan2026.tif")
writeRaster(mod_sa_rast_trimmed, mod_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

mod_ext_out_path <- file.path(sa_dir, "modelling_extent_aligned_cnvc_Jan2026.tif")
writeRaster(mod_ext_rast_trimmed, mod_ext_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Convert raster to polygons (excluding NA values)
mod_sa_polygons <- as.polygons(mod_sa_rast_trimmed, na.rm = TRUE, dissolve = TRUE)
mod_ext_polygons <- as.polygons(mod_ext_rast_trimmed, na.rm = TRUE, dissolve = TRUE)

# Export shapefile
mod_sa_out_shp <- file.path(sa_dir, "modelling_study_area_aligned_cnvc_Jan2026.shp")
writeVector(mod_sa_polygons, mod_sa_out_shp, overwrite = TRUE)

mod_ext_out_shp <- file.path(sa_dir, "modelling_extent_aligned_cnvc_Jan2026.shp")
writeVector(mod_ext_polygons, mod_ext_out_shp, overwrite = TRUE)
