
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS TO REFERENCE                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 3: align all grassland layers to classification study area file from step 1
#         to prep for analysis

# load packages
library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)

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

# lcc directory
lcc_root <- file.path(raw_grass_data_dir, "./Land Cover/NRCan Land Cover of Canada 2020")
lcc_aligned_dir <- file.path(proc_grass_data_dir, "nrcan_lcc/aligned")
lcc_cropped_dir <- file.path(proc_grass_data_dir, "nrcan_lcc/cropped")
lcc_pilot_dir <- file.path(proc_grass_data_dir, "nrcan_lcc/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# lcc raster
lcc_file <- file.path(lcc_root, "landcover-2020-classification.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- NRCan Land Cover of Canada ----#

# Load rasters
r_lcc <- rast(lcc_file)

origin(r_lcc) # 0 0
plot(r_lcc)
r_lcc

#class       : SpatRaster 
#dimensions  : 160001, 190001, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2600010, 3100020, -885090, 3914940  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83(CSRS) / Canada Atlas Lambert (EPSG:3979) 
#source      : landcover-2020-classification.tif 
#color table : 1 
#name        : Canada2020 
#min value   :          1 
#max value   :         19 

# print frequency table
freq_table <- freq(r_lcc) 
print(freq_table) 

#  layer value      count
#1      1     1 2656885527
#2      1     2  267527073
#3      1     3          1
#4      1     5  222600062
#5      1     6  631127300
#6      1     8  585998755
#7      1    10  576963594
#8      1    11  600686833
#9      1    12  788852530
#10     1    13  731489329
#11     1    14  792439820
#12     1    15  500976019
#13     1    16  680173541
#14     1    17   71615576
#15     1    18 2463840984
#16     1    19  204148705

# Mark as categorical (factor)
r_lcc <- as.factor(r_lcc)
r_lcc
plot(r_lcc)

# crop before reprojecting to make manageable size 

# Define the bbox in lon/lat degrees (WGS84)
canada_ext_ll  <- ext(-140, -70, 30, 60)
canada_bbox_ll <- as.polygons(canada_ext_ll, crs = "EPSG:4326")

# Project bbox polygon to EPSG:3979 (meters)
canada_bbox_3979 <- project(canada_bbox_ll, "EPSG:3979")

# This is your meter-based extent (what you asked for)
canada_ext_3979 <- ext(canada_bbox_3979)

# polygon for Canada
#canada_sf <- ne_countries(country = "Canada", returnclass = "sf")
#canada_v  <- vect(canada_sf)

# plot
#plot(canada_v, col = "grey90", border = "grey40", main = "Canada + bounding box")
#plot(canada_bbox, add = TRUE, border = "red", lwd = 2, lty = 2)

# crop and mask to bounding box
r_lcc_crop <- crop(r_lcc, canada_ext_3979)
r_lcc_mask <- mask(r_lcc_crop, canada_bbox_3979)

r_lcc_mask
plot(r_lcc_mask)

# Save the raster as 8-bit unsigned integer
lcc_crop_out_path <- file.path(lcc_cropped_dir, "lcc_crop_bbox.tif")
writeRaster(r_lcc_mask, lcc_crop_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
#lcc_crop_trim <- trim(r_lcc_mask)

# visualize the trimmed raster
#plot(lcc_crop_trim)

# print raster details 
#origin(lcc_crop_trim) #
#lcc_crop_trim

# reproject before resampling
#lcc_reproj <- project(lcc_crop_trim, sa, method = "near") # better approach
lcc_reproj <- project(r_lcc_mask, sa, method = "near") # better approach

freq_table2 <- freq(lcc_reproj) 
print(freq_table2) 

origin(lcc_reproj) # 
lcc_reproj
plot(lcc_reproj)

# align to reference raster 
lcc_reproj_algn <- resample(lcc_reproj, sa, method = "near")
plot(lcc_reproj_algn)

# Get values and counts
freq_table3 <- freq(lcc_reproj_algn) 
print(freq_table3) 

origin(lcc_reproj_algn) # 6.696121 -2.550334
lcc_reproj_algn

# Save the raster as 8-bit unsigned integer
lcc_out_path <- file.path(lcc_aligned_dir, "lcc_reproj_algn.tif")
writeRaster(lcc_reproj_algn, lcc_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
lcc_reproj_algn_trim <- trim(lcc_reproj_algn)

# visualize the trimmed raster
plot(lcc_reproj_algn_trim)

# print raster details 
origin(lcc_reproj_algn_trim) # 6.696121 -2.550334
lcc_reproj_algn_trim

# Save the raster as 8-bit unsigned integer
lcc_trim_out_path <- file.path(lcc_aligned_dir, "lcc_reproj_algn_trim.tif")
writeRaster(lcc_reproj_algn_trim, lcc_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop and mask to study area
lcc_reproj_algn_crop <- crop(lcc_reproj_algn, ext(sa))
lcc_reproj_algn_mask <- mask(lcc_reproj_algn_crop, sa)

# Save the raster as 8-bit unsigned integer
lcc_sa_out_path <- file.path(lcc_cropped_dir, "lcc_reproj_algn_crop_sa.tif")
writeRaster(lcc_reproj_algn_mask, lcc_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (grass = 1)
lcc_algn_binary <- classify(lcc_reproj_algn, 
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
plot(lcc_algn_binary)
lcc_algn_binary

# Save the raster as 8-bit unsigned integer
lcc_binary_out_path <- file.path(lcc_aligned_dir, "lcc_reproj_algn_binary_noshrub.tif")
writeRaster(lcc_algn_binary, lcc_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334
res(pa) # 30 30

# crop to pilot area
crop_ext <- ext(pa)
lcc_pa <- crop(lcc_reproj_algn, crop_ext)
lcc_pa
plot (lcc_pa)


# mask to pilot area
lcc_masked <- mask(lcc_pa, pa)
lcc_masked
plot(lcc_masked)



# print raster details 
origin(lcc_masked) # 6.696121 -2.550334
res(lcc_masked) # 30 30

# Save the raster as 8-bit unsigned integer
lcc_pa_out_path <- file.path(lcc_cropped_dir, "lcc_reproj_algn_crop_pa.tif")
writeRaster(lcc_masked, lcc_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
lcc_bin_pa <- crop(lcc_algn_binary, crop_ext)
lcc_bin_pa
plot (lcc_bin_pa)

# mask to pilot area
lcc_bin_masked <- mask(lcc_bin_pa, pa)
lcc_bin_masked
plot(lcc_bin_masked)

# Save the raster as 8-bit unsigned integer
lcc_bin_pa_out_path <- file.path(lcc_pilot_dir, "lcc_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(lcc_bin_masked, lcc_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
