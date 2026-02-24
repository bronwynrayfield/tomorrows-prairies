

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

# ECCC directory
eccc_root <- file.path(raw_grass_data_dir, "eccc_grass")
eccc_aligned_dir <- file.path(proc_grass_data_dir, "eccc_grass/aligned")
eccc_cropped_dir <- file.path(proc_grass_data_dir, "eccc_grass/cropped")
eccc_pilot_dir <- file.path(proc_grass_data_dir, "eccc_grass/pilot")

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# ECCC grass summary
eccc_file <- file.path(eccc_root, "grass_summary_a_class.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
res(sa) # 30 30
crs(sa) # EPSG:9001
ext(sa) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

####--- ALGIN RASTERS ---#### 

#----- ECCC ----#

# Load rasters
r_eccc <- rast(eccc_file)

res(r_eccc) # 30 30
origin(r_eccc) # 0 0
crs(r_eccc)
ext(r_eccc) # SpatExtent : -1452000, 24000, 6552000, 7524000 (xmin, xmax, ymin, ymax)
nrow(r_eccc); ncol(r_eccc) # 32400; 49200

plot(r_eccc)

# print frequency table
freq_table <- freq(r_eccc) 
print(freq_table) # 0: 737355633; 1: 1616626; 2: 1327389; 

# need to reclassify to avoid missing integers so new class doesn't get introduced when resampling
eccc_reclass <- classify(r_eccc, 
                        rcl = matrix(c(0, 0, # not currently grass or pasture
                                       1, 1, # currently grass, was forest
                                       2, 2, # currently grass, was shrub
                                       3, 3, # currently grass, was crop
                                       4, 4, # currently grass, was pasture
                                       5, 5, # grassland with no detected disturbance (back to 1984)
                                       6, 6, # currently grass, was built-up
                                       7, 7, # currently grass, was bare
                                       8, 8, # currently grass, was wetland
                                       9, 9, # currently grass, was water
                                       11, 10), # currently pasture (tame grasses)
                                     ncol = 2, byrow = TRUE))


unique_vals <- levels(as.factor(eccc_reclass))[[1]]$ID
unique_vals

# Mark as categorical (factor)
eccc_reclass <- as.factor(eccc_reclass)

# Add category labels
#levels(ngi_reclass) <- data.frame(
#  value = c(1, 2, 3),
#  label = c("Tame", "Native", "Mixed")
#)

res(eccc_reclass) # 30 30
origin(eccc_reclass) # 0 0 
crs(eccc_reclass) # EPSG:9001
ext(eccc_reclass) # SpatExtent : -1452000, 24000, 6552000, 7524000 (xmin, xmax, ymin, ymax)

# reproject before resampling
eccc_reclass_reproj <- project(eccc_reclass, crs(sa)) # or should be project(ngi_reclass, sa) ??

freq_table2 <- freq(eccc_reclass_reproj) 
print(freq_table2) # 0: 779640904; 1: 1714088; 2: 1393874

res(eccc_reclass_reproj) # 29.45835 29.45835
origin(eccc_reclass_reproj) # -3.814708 13.815769
crs(eccc_reclass_reproj) # EPSG:9001
ext(eccc_reclass_reproj) # SpatExtent : -1387963.62381047, 92171.3750936405, 919909.837478601, 2055293.71901432 (xmin, xmax, ymin, ymax)

# align to reference raster 
eccc_reclass_reproj_algn <- resample(eccc_reclass_reproj, sa, method = "near")
plot(eccc_reclass_reproj_algn)

# Get values and counts
freq_table3 <- freq(eccc_reclass_reproj_algn) 
print(freq_table3) # 0: 751735029; 1: 1653399; 2: 1343846

res(eccc_reclass_reproj_algn) # 30 30
origin(eccc_reclass_reproj_algn) # 6.696121 -2.550334
crs(eccc_reclass_reproj_algn) # EPSG:9001
ext(eccc_reclass_reproj_algn) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
eccc_out_path <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn.tif")
writeRaster(eccc_reclass_reproj_algn, eccc_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
eccc_reclass_reproj_algn_trim <- trim(eccc_reclass_reproj_algn)

# visualize the trimmed raster
plot(eccc_reclass_reproj_algn_trim)

# print raster details 
origin(eccc_reclass_reproj_algn_trim) # 6.696121 -2.550334
res(eccc_reclass_reproj_algn_trim) # 30 30 
crs(eccc_reclass_reproj_algn_trim) # EPSG:9001
ext(eccc_reclass_reproj_algn_trim) # SpatExtent : -1360253.30387889, 69336.6961211143, 1036197.4496659, 1963707.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
eccc_trim_out_path <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn_trim.tif")
writeRaster(eccc_reclass_reproj_algn_trim, eccc_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (values 1-10 become 1 for grass)
eccc_algn_binary <- classify(eccc_reclass_reproj_algn, 
                            rcl = matrix(c(0, 0,
                                           1, 1,
                                           2, 1,
                                           3, 1,
                                           4, 1,
                                           5, 1,
                                           6, 1,
                                           7, 1,
                                           8, 1,
                                           9, 1,
                                           10, 1), 
                                         ncol = 2, byrow = TRUE))

# Save the raster as 8-bit unsigned integer
eccc_binary_out_path <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn_binary.tif")
writeRaster(eccc_algn_binary, eccc_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
eccc_algn_nattam <- classify(eccc_algn, 
                             rcl = matrix(c(0, 0, # not currently grass or pasture
                                            1, 1, # currently grass, was forest
                                            2, 1, # currently grass, was shrub
                                            3, 1, # currently grass, was crop
                                            4, 1, # currently grass, was pasture
                                            5, 2, # grassland with no detected disturbance (back to 1984)
                                            6, 1, # currently grass, was built-up
                                            7, 1, # currently grass, was bare
                                            8, 1, # currently grass, was wetland
                                            9, 1, # currently grass, was water
                                            10, 1), # currently pasture (tame grasses)
                                          ncol = 2, byrow = TRUE))

# Save the raster as 8-bit unsigned integer
eccc_nattam_out_path <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn_natvtam.tif")
writeRaster(eccc_algn_nattam, eccc_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load eccc raster
eccc_algn_file <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn.tif")
eccc_algn <- rast(eccc_algn_file)
plot(eccc_algn)
eccc_algn

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa)
res(pa)
crs(pa)

# crop to pilot area
crop_ext <- ext(pa)
eccc_pa <- crop(eccc_algn, crop_ext)
eccc_pa
plot (eccc_pa)

# mask to pilot area
eccc_masked <- mask(eccc_pa, pa)
eccc_masked
plot(eccc_masked)

# print raster details 
origin(eccc_masked)
res(eccc_masked)
crs(eccc_masked)

# Save the raster as 8-bit unsigned integer
eccc_pa_out_path <- file.path(eccc_pilot_dir, "eccc_reproj_algn_crop_pa.tif")
writeRaster(eccc_masked, eccc_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load binary raster
eccc_binary_file <- file.path(eccc_aligned_dir, "eccc_reclass_reproj_algn_binary.tif")
eccc_algn_binary <- rast(eccc_binary_file)

# crop binary raster to pilot area
eccc_bin_pa <- crop(eccc_algn_binary, crop_ext)
eccc_bin_pa
plot (eccc_bin_pa)

# mask to pilot area
eccc_bin_masked <- mask(eccc_bin_pa, pa)
eccc_bin_masked
plot(eccc_bin_masked)

# Save the raster as 8-bit unsigned integer
eccc_bin_pa_out_path <- file.path(eccc_pilot_dir, "eccc_reproj_algn_crop_pa_binary.tif")
writeRaster(eccc_bin_masked, eccc_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
eccc_nt_pa <- crop(eccc_algn_nattam, crop_ext)
eccc_nt_pa
plot (eccc_nt_pa)

# mask to pilot area
eccc_nt_masked <- mask(eccc_nt_pa, pa)
eccc_nt_masked
plot(eccc_nt_masked)

# Save the raster as 8-bit unsigned integer
eccc_nt_pa_out_path <- file.path(eccc_pilot_dir, "eccc_reproj_algn_crop_pa_natvtam.tif")
writeRaster(eccc_nt_masked, eccc_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
