
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
#library(parallel)

# Use most cores, but leave 1 free (good practice)
#nthreads <- max(1, parallel::detectCores() - 1)

# Tell terra to use multiple threads
#terraOptions(threads = nthreads)

# (Optional) ensure on-disk temp files go somewhere fast
# terraOptions(tempdir = "D:/terra_tmp")  # example on Windows
#terraOptions(memfrac = 0.8)             # use more RAM if safe

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/class_inputs"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# AAFC LUTS directory
luts_tile_root <- file.path(raw_grass_data_dir, "aafc_luts_2020_30m_v5")
luts_merged_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/merged")
luts_cropped_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/cropped")
luts_aligned_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/aligned")
luts_pilot_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/pilot")


####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# AAFC LUTS 
luts_files <- file.path(luts_tile_root)

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
res(sa)
crs(sa)

####--- ALGIN RASTERS ---#### 

#----- AAFC LUTS ----#

# List all .tif files
luts_file_list <- list.files(luts_files, pattern = "\\.tif$", full.names = TRUE)

# Load them as a SpatRaster collection
r_luts <- lapply(luts_file_list, rast)

# Tiles do not share a common origin (they are not snapped to a common grid)
# Use the first raster as reference for grid alignment
ref_raster <- r_luts[[1]]

# Align all rasters to reference raster
# Use "near" for categorical data 
resample_method <- "near"

luts_aligned_luts <- lapply(r_luts, function(r) {
  project(r, ref_raster, method = "near", align = TRUE)
})

# Merge the aligned rasters
luts_merged <- Reduce(merge, luts_aligned_luts)

# print raster details 
origin(luts_merged) # -9.833142 -5.845974
res(luts_merged) # 30 30 
crs(luts_merged) # WGS 84 / UTM zone 7N (EPSG:32607)
ext(luts_merged) # SpatExtent : 607190.166858264, 6351740.16685826, 5421684.15402559, 9910794.15402559 (xmin, xmax, ymin, ymax)

# plot the merged raster
plot(luts_merged)

# Get values and counts
freq_table <- freq(luts_merged) 
print(freq_table) # 21: 13731441; 22: 6163427; 24: 21540208 etc.

# Save the raster as 8-bit unsigned integer
luts_merge_out_path <- file.path(luts_merged_dir, "luts_merged.tif")
writeRaster(luts_merged, luts_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load the merged raster
luts_merged_file <- file.path(luts_merged_dir, "luts_merged.tif")
luts_merged <- rast(luts_merged_file)
plot(luts_merged)
luts_merged

# Mark as categorical (factor)
luts_merged <- as.factor(luts_merged)

# reproject before resampling
luts_merged_reproj <- project(luts_merged, crs(sa)) # or should be project(ngi_reclass, sa) ??
luts_merged_reproj

# need to reclassify to avoid missing integers so new class doesn't get introduced when resampling
#luts_reclass <- classify(luts_merged, 
#                         rcl = matrix(c(21, 1,
#                                        22, 2,
#                                        24, 3,
#                                        25, 4,
#                                        28, 5,
#                                        29, 6,
#                                        31, 7,
#                                        41, 8,
#                                        42, 9,
#                                        43, 10,
#                                        44, 11), 
#                                      ncol = 2, byrow = TRUE))

# align to reference raster 
luts_merged_algn <- resample(luts_merged_reproj, sa, method = "near")
plot(luts_merged_algn)
luts_merged_algn
origin(luts_merged_algn) # 6.696121 -2.550334

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_18b847077116b_101252_80Q6udnsxVen2Mq.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : lyr.1 
#name        : lyr.1 
#min value   :    21 
#max value   :    92 

# Get values and counts
freq_table2 <- freq(luts_merged_algn) 
print(freq_table2) # 21: 8151488; 22: 3886913; 24: 9300748 etc.

# Save the raster as 8-bit unsigned integer
luts_algn_out_path <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn.tif")
writeRaster(luts_merged_algn, luts_algn_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
luts_merged_algn_trim <- trim(luts_merged_algn)
luts_merged_algn_trim

#class       : SpatRaster 
#size        : 86509, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 271047.4, 2866317  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_18b841dda2b92_101252_DvggOkjkTOQeAAI.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : lyr.1 
#name        : lyr.1 
#min value   :    21 
#max value   :    92

# visualize the trimmed raster
plot(luts_merged_algn_trim)

# print raster details 
origin(luts_merged_algn_trim) # 6.696121 -2.550334
res(luts_merged_algn_trim) # 30 30 
crs(luts_merged_algn_trim) # EPSG:9001
ext(luts_merged_algn_trim) # SpatExtent : -2321003.30387889, 1585686.69612111, 271047.4496659, 2866317.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
luts_trim_out_path <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn_trim.tif")
writeRaster(luts_merged_algn_trim, luts_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (values 2, 3, 5 become 1)
luts_algn_binary <- classify(luts_merged_algn, 
                            rcl = matrix(c(21, 0, # Settlement
                                           22, 0, # High Reflectance Settlement
                                           24, 0, # Settlement Forest
                                           25, 0, # Roads
                                           28, 0, # Vegetated Settlement
                                           29, 0, # Very High Reflectance Settlement
                                           31, 0, # Water
                                           41, 0, # Forest
                                           42, 0, # Forest Wetland
                                           43, 0, # Forest Regenerating after Harvest (<20 years)
                                           44, 0, # Forest Wetland Regenerating after Harvest (<20 years)
                                           47, 0, # Forest Regenerating after Harvest (20-29 years)
                                           48, 0, # Forest Wetland Regenerating after Harvest (20-29 years)
                                           49, 0, # Forest Regenerating after Fire (<20 years)
                                           51, 0, # Cropland
                                           52, 0, # Annual Cropland
                                           55, 0, # Land Converted to Cropland
                                           56, 0, # Land Converted to Annual Cropland
                                           63, 1, # Prairie Grassland
                                           64, 1, # Temperate Grassland
                                           65, 1, # Temperate Rainforest Grassland
                                           66, 1, # Boreal Grassland
                                           67, 1, # Montane Grassland
                                           68, 1, # Alpine Tundra Grassland
                                           69, 1, # Arctic Tundra Grassland
                                           73, 0, # Prairie Wetland
                                           74, 0, # Temperate Wetland
                                           75, 0, # Temperate Rainforest Wetland
                                           76, 0, # Boreal Wetland
                                           78, 0, # Alpine Tundra Wetland
                                           79, 0, # Arctic Tundra Wetland
                                           81, 0, # Newly-Detected Settlement (<10 years)
                                           82, 0, # Newly-Detected High Reflectance Settlement (<10 years)
                                           84, 0, # Newly-Detected Settlement Forest (<10 years)
                                           85, 0, # Newly-detected roads (<10 years)
                                           88, 0, # Newly-Detected Vegetated Settlement (<10 years)
                                           89, 0, # Newly-Detected Very High Reflectance Settlement (<10 years)
                                           91, 0, # Other Land
                                           92, 0), # Snow and Ice
                                         ncol = 2, byrow = TRUE))

# print raster details and plot
luts_algn_binary
plot(luts_algn_binary)

# Save the raster as 8-bit unsigned integer
luts_binary_out_path <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn_binary.tif")
writeRaster(luts_algn_binary, luts_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa)
pa

# crop to pilot area
crop_ext <- ext(pa)
luts_pa <- crop(luts_merged_algn, crop_ext)
luts_pa
plot (luts_pa)

# mask to pilot area
luts_masked <- mask(luts_pa, pa)
luts_masked
plot(luts_masked)

# print raster details 
origin(luts_masked)
res(luts_masked)
crs(luts_masked)

# Save the raster as 8-bit unsigned integer
luts_pa_out_path <- file.path(luts_pilot_dir, "aafc_luts_merged_reproj_algn_crop_pa.tif")
writeRaster(luts_masked, luts_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load binary raster
luts_binary_file <- file.path(luts_aligned_dir, "aafc_luts_merged_reproj_algn_binary.tif")
luts_algn_binary <- rast(luts_binary_file)

# crop binary raster to pilot area
luts_bin_pa <- crop(luts_algn_binary, crop_ext)
luts_bin_pa
plot (luts_bin_pa)

# mask to pilot area
luts_bin_masked <- mask(luts_bin_pa, pa)
luts_bin_masked
plot(luts_bin_masked)

# Save the raster as 8-bit unsigned integer
luts_bin_pa_out_path <- file.path(luts_pilot_dir, "luts_reproj_algn_crop_pa_binary.tif")
writeRaster(luts_bin_masked, luts_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

