
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS TO REFERENCE                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# STEP 3: Align all grassland layers to classification study area file from step 1
#         to prep for analysis

# load packages
library(terra)
library(sf)
library(tidyverse)
#library(rnaturalearth)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# nalcms directory
nalcms_root <- file.path(raw_grass_data_dir, "Landcover/CEC North American Land Change Monitoring System 2020/can_land_cover_2020v2_30m_tif/CAN_NALCMS_landcover_2020v2_30m/data")
nalcms_aligned_dir <- file.path(proc_grass_data_dir, "cec_nalcms/aligned")
nalcms_cropped_dir <- file.path(proc_grass_data_dir, "cec_nalcms/cropped")
nalcms_masked_dir <- file.path(proc_grass_data_dir, "cec_nalcms/masked")
nalcms_pilot_dir <- file.path(proc_grass_data_dir, "cec_nalcms/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# nalcms raster
nalcms_file <- file.path(nalcms_root, "CAN_NALCMS_landcover_2020v2_30m.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- StatCan Land Cover Register ----#

# Load rasters
r_nalcms <- rast(nalcms_file)

origin(r_nalcms) # -3.72529e-09 -1.00000e+01
r_nalcms

#class       : SpatRaster 
#size        : 152010, 188680, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2315310, 3345090, -250300, 4310000  (xmin, xmax, ymin, ymax)
#coord. ref. : WGS_1984_Lambert_Azimuthal_Equal_Area 
#source      : CAN_NALCMS_landcover_2020v2_30m.tif 
#color table : 1 
#categories  : Class_EN, Class_ES, Class_FR, Count 
#name        :     Class_EN 
#min value   :              
#max value   : Snow and Ice 

# print frequency table
freq_table <- freq(r_nalcms) 
print(freq_table) 

# Mark as categorical (factor)
r_nalcms <- as.factor(r_nalcms)

# crop before reprojecting to make manageable size 

# Define the bbox in lon/lat degrees (WGS84)
canada_ext_ll  <- ext(-140, -70, 30, 60)
canada_bbox_ll <- as.polygons(canada_ext_ll, crs = "EPSG:4326")

# Get CRS from the raster 
r_crs <- crs(r_nalcms)

# Project bbox polygon to the raster's CRS
canada_bbox_rcrs <- project(canada_bbox_ll, r_crs)

# Meter/CRS-based extent derived from projected polygon
canada_ext_rcrs <- ext(canada_bbox_rcrs)

# polygon for Canada
#canada_sf <- ne_countries(country = "Canada", returnclass = "sf")
#canada_v  <- vect(canada_sf)

# plot
#plot(canada_v, col = "grey90", border = "grey40", main = "Canada + bounding box")
#plot(canada_bbox, add = TRUE, border = "red", lwd = 2, lty = 2)

# crop and mask to bounding box
r_nalcms_crop <- crop(r_nalcms, canada_ext_rcrs)
r_nalcms_mask <- mask(r_nalcms_crop, canada_bbox_rcrs)

r_nalcms_mask
plot(r_nalcms_mask)

# Save the raster as 8-bit unsigned integer
nalcms_crop_out_path <- file.path(nalcms_cropped_dir, "nalcms_crop_bbox.tif")
writeRaster(r_nalcms_mask, nalcms_crop_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Set pixels with value 0 to NA - ocean, bay, US
r_nalcms_mask[r_nalcms_mask == 0] <- NA
plot(r_nalcms_mask)

# reproject before resampling
nalcms_reproj <- project(r_nalcms_mask, sa, method = "near") # better approach

freq_table2 <- freq(nalcms_reproj) 
print(freq_table2) 

origin(nalcms_reproj) # 
nalcms_reproj

# align to reference raster 
nalcms_reproj_algn <- resample(nalcms_reproj, sa, method = "near")
plot(nalcms_reproj_algn)

# Get values and counts
freq_table3 <- freq(nalcms_reproj_algn) 
print(freq_table3) # 

origin(nalcms_reproj_algn) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
nalcms_out_path <- file.path(nalcms_aligned_dir, "nalcms_reproj_algn.tif")
writeRaster(nalcms_reproj_algn, nalcms_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
nalcms_reproj_algn_trim <- trim(nalcms_reproj_algn)

# visualize the trimmed raster
plot(nalcms_reproj_algn_trim)

# print raster details 
origin(nalcms_reproj_algn_trim) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
nalcms_trim_out_path <- file.path(nalcms_aligned_dir, "nalcms_reproj_algn_trim.tif")
writeRaster(nalcms_reproj_algn_trim, nalcms_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (grass = 1)
nalcms_algn_binary <- classify(nalcms_reproj_algn, 
                            rcl = matrix(c(1, 0, # Temperate or sub-polar needleleaf forest
                                           2, 0, # Sub-polar taiga needleleaf forest
                                           5, 0, # Temperate or sub-polar broadleaf deciduous forest
                                           6, 0, # Temperate or sub-polar mixed forest
                                           8, 0, # Temperate or sub-polar shrubland
                                           10, 1, # Temperate or sub-polar grassland
                                           11, 0, # Sub-polar or polar shrubland-lichen-moss
                                           12, 0, # Sub-polar or polar grassland-lichen-moss
                                           13, 0, # Sub-polar or polar barren-lichen-moss
                                           14, 0, # Wetland
                                           15, 0, # Cropland
                                           16, 0, # Barren lands
                                           17, 0, # Urban
                                           18, 0, # Water
                                           19, 0), # Snow and ice
                                         ncol = 2, byrow = TRUE))
plot(nalcms_algn_binary)
nalcms_algn_binary

# Save the raster as 8-bit unsigned integer
nalcms_binary_out_path <- file.path(nalcms_aligned_dir, "nalcms_reproj_algn_binary_noshrub.tif")
writeRaster(nalcms_algn_binary, nalcms_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

#### ISOLATE CLASSES ####

# load raster
nalcms_file <- file.path(nalcms_masked_dir, "nalcms_reproj_algn_masked.tif")
nalcms <- rast(nalcms_file)

nalcms
plot(nalcms)

# isolate shrub
nalcms_algn_shrub <- classify(nalcms, 
                               rcl = matrix(c(1, 0, # Temperate or sub-polar needleleaf forest
                                              2, 0, # Sub-polar taiga needleleaf forest
                                              5, 0, # Temperate or sub-polar broadleaf deciduous forest
                                              6, 0, # Temperate or sub-polar mixed forest
                                              8, 1, # Temperate or sub-polar shrubland
                                              10, 0, # Temperate or sub-polar grassland
                                              11, 0, # Sub-polar or polar shrubland-lichen-moss
                                              12, 0, # Sub-polar or polar grassland-lichen-moss
                                              13, 0, # Sub-polar or polar barren-lichen-moss
                                              14, 0, # Wetland
                                              15, 0, # Cropland
                                              16, 0, # Barren lands
                                              17, 0, # Urban
                                              18, 0, # Water
                                              19, 0), # Snow and ice
                                            ncol = 2, byrow = TRUE))
plot(nalcms_algn_shrub)
nalcms_algn_shrub

# Save the raster as 8-bit unsigned integer
nalcms_shrub_out_path <- file.path(nalcms_aligned_dir, "nalcms_reproj_algn_shrub.tif")
writeRaster(nalcms_algn_shrub, nalcms_shrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

#### PILOT AREA RASTERS ####

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334

# crop to pilot area
crop_ext <- ext(pa)
nalcms_pa <- crop(nalcms_reproj_algn, crop_ext)
nalcms_pa
plot (nalcms_pa)


# mask to pilot area
nalcms_masked <- mask(nalcms_pa, pa)
nalcms_masked
plot(nalcms_masked)

# print raster details 
origin(nalcms_masked) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
nalcms_pa_out_path <- file.path(nalcms_pilot_dir, "nalcms_reproj_algn_crop_pa.tif")
writeRaster(nalcms_masked, nalcms_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
nalcms_bin_pa <- crop(nalcms_algn_binary, crop_ext)
nalcms_bin_pa
plot (nalcms_bin_pa)

# mask to pilot area
nalcms_bin_masked <- mask(nalcms_bin_pa, pa)
nalcms_bin_masked
plot(nalcms_bin_masked)

# Save the raster as 8-bit unsigned integer
nalcms_bin_pa_out_path <- file.path(nalcms_pilot_dir, "nalcms_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(nalcms_bin_masked, nalcms_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################

######### crop and mask to study area

# load raster
nalcms_file <- file.path(nalcms_aligned_dir, "nalcms_reproj_algn.tif")
nalcms <- rast(nalcms_file)

crop_ext <- ext(sa)

# crop raster to study area
nalcms_sa <- crop(nalcms, crop_ext)
nalcms_sa
plot (nalcms_sa)

# mask to study area
nalcms_masked <- mask(nalcms_sa, sa)
nalcms_masked
plot(nalcms_masked)

# Save the raster as 8-bit unsigned integer
nalcms_sa_out_path <- file.path(nalcms_mask_dir, "nalcms_reproj_algn_masked.tif")
writeRaster(nalcms_masked, nalcms_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load binary raster
nalcms_binary_file <- file.path(nalcms_aligned_dir, "nalcms_reproj_algn_binary_noshrub.tif")
nalcms_binary <- rast(nalcms_binary_file)

crop_ext <- ext(sa)

# crop raster to study area
nalcms_bin_sa <- crop(nalcms_binary, crop_ext)
nalcms_bin_sa
plot (nalcms_bin_sa)

# mask to study area
nalcms_bin_masked <- mask(nalcms_bin_sa, sa)
nalcms_bin_masked
plot(nalcms_bin_masked)

# Save the raster as 8-bit unsigned integer
nalcms_bin_sa_out_path <- file.path(nalcms_mask_dir, "nalcms_reproj_algn_binary_noshrub_masked.tif")
writeRaster(nalcms_bin_masked, nalcms_bin_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
