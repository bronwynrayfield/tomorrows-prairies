
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

# NGI tiles directory
ngi_tile_root <- file.path(raw_grass_data_dir, "cfga_ngi")
ngi_merged_dir <- file.path(proc_grass_data_dir, "cfga_ngi/merged")
ngi_cropped_dir <- file.path(proc_grass_data_dir, "cfga_ngi/cropped")
ngi_aligned_dir <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
ngi_pilot_dir <- file.path(proc_grass_data_dir, "cfga_ngi/pilot")

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# NGI rasters
ngi_ab_file <- file.path(ngi_tile_root, "AB-002.tif")
ngi_mb_file <- file.path(ngi_tile_root, "MB-003.tif")
ngi_sk_file <- file.path(ngi_tile_root, "SK-001.tif")

# NGI extent
ngi_extent_file <- "R:/CanadianGrasslands/WG_Classification/Classification_WG_intermediate_products/exported_shp_for_GEE/NGI_extent.shp"

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
res(sa)
crs(sa)
sa

####--- ALGIN RASTERS ---#### 

#----- NGI ----#

# Load rasters
r_ab <- rast(ngi_ab_file)
r_mb <- rast(ngi_mb_file)
r_sk <- rast(ngi_sk_file)

plot(r_ab)

res(r_ab) # 10 x 10
res(r_mb) # 10 x 10
res(r_sk) # 10 x 10
origin(r_ab) # 0 0
origin(r_mb) # 0 0
origin(r_sk) # 0 0 

ext(r_ab) # SpatExtent : 4622750, 5033770, 1616880, 2274690 (xmin, xmax, ymin, ymax)
ext(r_mb) # SpatExtent : 5507770, 5860000, 1443910, 1780640 (xmin, xmax, ymin, ymax)
ext(r_sk) # SpatExtent : 4890950, 5532640, 1480800, 2103970 (xmin, xmax, ymin, ymax)
nrow(r_ab); ncol(r_ab) # 65781; 41102
nrow(r_mb); ncol(r_mb) # 33673; 35223
nrow(r_sk); ncol(r_sk) # 62317; 64169

# Set pixels with value 0 to NA - NO! need to preserve zeroes
#r_ab[r_ab == 0] <- NA
#r_mb[r_mb == 0] <- NA
#r_sk[r_sk == 0] <- NA

# Merge 
merged_ngi <- mosaic(r_ab, r_mb, r_sk, fun = "max") # when rasters overlap, keep the largest value

res(merged_ngi) # 10 10
origin(merged_ngi) # 0 0 
crs(merged_ngi)
ext(merged_ngi) # SpatExtent : 4622750, 5860000, 1443910, 2274690 (xmin, xmax, ymin, ymax)

# plot the result
plot(merged_ngi, main = "Merged NGI Raster")

# print frequency table
freq_table <- freq(merged_ngi) 
print(freq_table) 

#  layer value      count
#1     1     0 5144984811
#2     1     2  894806781
#3     1     3  864399198
#4     1     5  214191378

# Save the raster as 8-bit unsigned integer
ngi_merge_out_path <- file.path(ngi_merged_dir, "cfga_ngi_mosaicked.tif")
writeRaster(merged_ngi, ngi_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## load NGI raster
ngi <- file.path(ngi_merged_dir, "cfga_ngi_mosaicked.tif")
ngi_merged <- rast(ngi)

# need to reclassify to 1, 2, 3 instead of 2, 3, 5 so new class (4) doesn't get introduced when resampling
#ngi_reclass <- classify(merged_ngi, 
#                            rcl = matrix(c(2, 1,
#                                           3, 2,
#                                           5, 3), 
#                                         ncol = 2, byrow = TRUE))

#freq_table3 <- freq(ngi_reclass) 
#print(freq_table3) 

#  layer value      count
#1     1     0 5144984811
#2     1     1  894806781
#3     1     2  864399198
#4     1     3  214191378

# Mark as categorical (factor)
ngi_merged <- as.factor(ngi_merged)
ngi_merged
plot(ngi_merged)

#class       : SpatRaster 
#size        : 83078, 123725, 1  (nrow, ncol, nlyr)
#resolution  : 10, 10  (x, y)
#extent      : 4622750, 5860000, 1443910, 2274690  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83 / Statistics Canada Lambert (EPSG:3347) 
#source      : cfga_ngi_mosaicked.tif 
#categories  : AB-002 
#name        : AB-002 
#min value   :      0 
#max value   :      5 

# Add category labels
#levels(ngi_reclass) <- data.frame(
#  value = c(0, 1, 2, 3),
#  label = c("Not grass", "Tame", "Native", "Mixed")
#)

#freq_table4 <- freq(ngi_reclass) 
#print(freq_table4) # 1: 894806781; 2: 864399198; 3: 214191378

#  layer  value      count
#1     1   <NA> 5144984811
#2     1   Tame  894806781
#3     1 Native  864399198
#4     1  Mixed  214191378

#res(ngi_reclass) # 10 10
#origin(ngi_reclass) # 0 0 
#crs(ngi_reclass) # EPSG:3347
#ext(ngi_reclass) # SpatExtent : 4622750, 5860000, 1443910, 2274690 (xmin, xmax, ymin, ymax)


# reproject before resampling
#ngi_reclass_merged_reproj <- project(ngi_reclass, crs(sa)) # or should be project(ngi_reclass, sa) ??
ngi_merged_reproj <- project(ngi_merged, crs(sa))

freq_table5 <- freq(ngi_merged_reproj) 
print(freq_table5) # 1: 954680046; 2: 917799068; 3: 227772560

#  layer value      count
#1     1     0 5500305860
#2     1     2  954680046
#3     1     3  917799068
#4     1     5  227772560

origin(ngi_merged_reproj) # 1.715237 -3.612827
ngi_merged_reproj

#class       : SpatRaster 
#size        : 105922, 130315, 1  (nrow, ncol, nlyr)
#resolution  : 9.754172, 9.754172  (x, y)
#extent      : -1307057, -35942.41, 908773.1, 1941954  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_21a6c5a8b200_137836_zxzxx9jyyx6UnnI.tif 
#categories  : AB-002 
#name        : AB-002 
#min value   :      0 
#max value   :      5 

# align to reference raster (NEAREST NEIGHBOR)
ngi_merged_reproj_algn <- resample(ngi_merged_reproj, sa, method = "near")
plot(ngi_merged_reproj_algn)

# Get values and counts
freq_table7 <- freq(ngi_merged_reproj_algn) 
print(freq_table7) # 1: 100921769; 1: 97023671; 3: 24084460

#  layer value     count
#1     1     0 581466091
#2     1     2 100921769
#3     1     3  97023671
#4     1     5  24084460

origin(ngi_merged_reproj_algn) # 6.696121 -2.550334
ngi_merged_reproj_algn

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_21a6c21c57717_137836_80Q6udnsxVen2Mq.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : AB-002 
#name        : AB-002 
#min value   :      0 
#max value   :      5 

# Save the raster as 8-bit unsigned integer
ngi_algn_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn.tif")
writeRaster(ngi_merged_reproj_algn, ngi_algn_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
ngi_merged_reproj_algn_trim <- trim(ngi_merged_reproj_algn)

# visualize the trimmed raster
plot(ngi_merged_reproj_algn_trim)

# print raster details 
origin(ngi_merged_reproj_algn_trim) # 6.696121 -2.550334
ngi_merged_reproj_algn_trim

#class       : SpatRaster 
#size        : 28401, 42370, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -1307063, -35963.3, 992127.4, 1844157  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source(s)   : memory
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : AB-002 
#name        : AB-002 
#min value   :      0 
#max value   :      5

# Save the raster as 8-bit unsigned integer
ngi_trim_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_trim.tif")
writeRaster(ngi_merged_reproj_algn_trim, ngi_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (values 2, 3, 5 become 1 = grass)
ngi_algn_binary <- classify(ngi_merged_reproj_algn, 
                               rcl = matrix(c(
                                 0, 0, # not grass
                                 2, 1, # tame
                                 3, 1, # native
                                 5, 1), # mixed
                               ncol = 2, byrow = TRUE))

plot(ngi_algn_binary)
ngi_algn_binary

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_21a6c522163a7_137836_DvggOkjkTOQeAAI.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#name        : AB-002 
#min value   :      0 
#max value   :      1

# Save the raster as 8-bit unsigned integer
ngi_binary_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_binary.tif")
writeRaster(ngi_algn_binary, ngi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load NGI raster
ngi <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn.tif")
ngi_algn <- rast(ngi)
plot(ngi_algn)

# use extent file to mask out zeroes that should be NA

# Load study area raster
ngi_extent <- vect(ngi_extent_file)
ngi_extent

plot(ngi_extent)

# project the vector to the raster CRS:
ngi_extent <- project(ngi_extent, crs(sa))

# Rasterize classification study area using CNVC raster as template
ngi_extent_rast <- rasterize(ngi_extent, sa, field=1, touches=TRUE) # field=1 assigns constant value to rasterized cells
# touches=TRUE ensures that any cell touched by the polygon is included

plot(ngi_extent_rast)

# crop to pilot area
ngi_crop_ext <- ext(ngi_extent_rast)
ngi_cropped <- crop(ngi_algn, ngi_crop_ext)
ngi_cropped
plot (ngi_cropped)

# mask to pilot area
ngi_masked <- mask(ngi_cropped, ngi_extent_rast)
ngi_masked
plot(ngi_masked)

# Save the raster as 8-bit unsigned integer
ngi_masked_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_good.tif")
writeRaster(ngi_masked, ngi_masked_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
ngi_algn_nattam <- classify(ngi_masked, 
                            rcl = matrix(c(
                              0, 0, # not grass
                              2, 1, # tame
                              3, 2, # native
                              5, 1), # mixed
                              ncol = 2, byrow = TRUE))
plot(ngi_algn_nattam)
ngi_algn_nattam

# Save the raster as 8-bit unsigned integer
ngi_nattam_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_natvtam_good.tif")
writeRaster(ngi_algn_nattam, ngi_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate each class
ngi_algn_native <- classify(ngi_masked, 
                            rcl = matrix(c(
                              0, 0, # not grass
                              2, 0, # tame
                              3, 1, # native
                              5, 0), # mixed
                              ncol = 2, byrow = TRUE))
plot(ngi_algn_native)
ngi_algn_native

# Save the raster as 8-bit unsigned integer
ngi_native_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_native.tif")
writeRaster(ngi_algn_native, ngi_native_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate each class
ngi_algn_tame <- classify(ngi_masked, 
                            rcl = matrix(c(
                              0, 0, # not grass
                              2, 1, # tame
                              3, 0, # native
                              5, 0), # mixed
                              ncol = 2, byrow = TRUE))
plot(ngi_algn_tame)
ngi_algn_tame

# Save the raster as 8-bit unsigned integer
ngi_tame_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_tame.tif")
writeRaster(ngi_algn_tame, ngi_tame_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate each class
ngi_algn_mixed <- classify(ngi_masked, 
                            rcl = matrix(c(
                              0, 0, # not grass
                              2, 0, # tame
                              3, 0, # native
                              5, 1), # mixed
                              ncol = 2, byrow = TRUE))
plot(ngi_algn_mixed)
ngi_algn_mixed

# Save the raster as 8-bit unsigned integer
ngi_mixed_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_mixed.tif")
writeRaster(ngi_algn_mixed, ngi_mixed_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
  
# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa)
res(pa)
crs(pa)

# crop to pilot area
crop_ext <- ext(pa)
ngi_pa <- crop(ngi_merged_reproj_algn, crop_ext)
ngi_pa
plot (ngi_pa)

# mask to pilot area
ngi_masked <- mask(ngi_pa, pa)
ngi_masked
plot(ngi_masked)

#class       : SpatRaster 
#size        : 1367, 1692, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -828263.3, -777503.3, 1111437, 1152447  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source(s)   : memory
#varname     : cfga_ngi_reclass_merged_reproj_algn_near 
#categories  : label 
#name        :  label 
#min value   :   Tame 
#max value   : Native 

# print raster details 
origin(ngi_masked) # 6.696121 -2.550334
res(ngi_masked)
crs(ngi_masked)

# Save the raster as 8-bit unsigned integer
ngi_pa_out_path <- file.path(ngi_pilot_dir, "ngi_merged_reproj_algn_crop_pa.tif")
writeRaster(ngi_masked, ngi_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
ngi_bin_pa <- crop(ngi_algn_binary, crop_ext)
ngi_bin_pa
plot (ngi_bin_pa)

# mask to pilot area
ngi_bin_masked <- mask(ngi_bin_pa, pa)
ngi_bin_masked
plot(ngi_bin_masked)

# Save the raster as 8-bit unsigned integer
ngi_bin_pa_out_path <- file.path(ngi_pilot_dir, "ngi_merged_reproj_algn_crop_pa_binary.tif")
writeRaster(ngi_bin_masked, ngi_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
ngi_nt_pa <- crop(ngi_algn_nattam, crop_ext)
ngi_nt_pa
plot (ngi_nt_pa)

# mask to pilot area
ngi_nt_masked <- mask(ngi_nt_pa, pa)
ngi_nt_masked
plot(ngi_nt_masked)

# Save the raster as 8-bit unsigned integer
ngi_nt_pa_out_path <- file.path(ngi_pilot_dir, "ngi_merged_reproj_algn_crop_pa_natvtam.tif")
writeRaster(ngi_nt_masked, ngi_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

