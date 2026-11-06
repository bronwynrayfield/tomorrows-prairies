
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
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# mgi tiles directory
mgi_tile_root <- file.path(raw_grass_data_dir, "Grassland/MB Manitoba Grassland Inventory")
mgi_cropped_dir <- file.path(proc_grass_data_dir, "mb_mgi/cropped")
mgi_aligned_dir <- file.path(proc_grass_data_dir, "mb_mgi/aligned")
mgi_pilot_dir <- file.path(proc_grass_data_dir, "mb_mgi/pilot")

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# mgi raster
mgi_file <- file.path(mgi_tile_root, "MGI.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
sa

####--- ALGIN RASTERS ---#### 

#----- MB Manitoba Grassland Inventory 2020 ----#

# Load rasters
r_mgi <- rast(mgi_file)

r_mgi
origin(r_mgi)
plot(r_mgi, main = "Raw MGI raster")

#class       : SpatRaster 
#dimensions  : 63100, 52547, 1  (nrow, ncol, nlyr)
#resolution  : 10, 10  (x, y)
#extent      : 263853.9, 789323.9, 5384950, 6015950  (xmin, xmax, ymin, ymax)
#coord. ref. : North_American_1983_CSRS_UTM_Zone_14N 
#source      : MGI.tif 
#categories  : CLASS 
#name        :   CLASS 
#min value   :      NA 
#max value   : Wetland 

# print frequency table
freq_table <- freq(r_mgi) 
print(freq_table)

#   layer    value      count
#1      1     <NA> 2096255987
#2      1   Forest  196719321
#3      1    Urban   32254912
#4      1    Water   83202237
#5      1 Cropland  472258047
#6      1   Native   76042255
#7      1    Shrub  134230939
#8      1     Tame   80524059
#9      1    Mixed   58914112
#10     1  Wetland   85313831

# reproject before resampling
mgi_reproj <- project(r_mgi, sa, method = "near") # better approach

# align to reference raster 
mgi_reproj_algn <- resample(mgi_reproj, sa, method = "near")

# inspect
plot(mgi_reproj_algn, main = "Aligned MGI raster")
origin(mgi_reproj_algn) 
mgi_reproj_algn

#class       : SpatRaster 
#dimensions  : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_244831913f9c_9288.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : CLASS 
#name        :   CLASS 
#min value   :      NA 
#max value   : Wetland 

# Get values and counts
freq_table2 <- freq(mgi_reproj_algn) 
print(freq_table2) 

#   layer    value     count
#1      1     <NA> 232944170
#2      1   Forest  21860412
#3      1    Urban   3585593
#4      1    Water   9249625
#5      1 Cropland  52499625
#6      1   Native   8456252
#7      1    Shrub  14920433
#8      1     Tame   8948515
#9      1    Mixed   6548895
#10     1  Wetland   9479581

# Save the raster as 8-bit unsigned integer
mgi_out_path <- file.path(mgi_aligned_dir, "mgi_reproj_algn.tif")
writeRaster(mgi_reproj_algn, mgi_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to fix
mgi_algn_fix <- classify(mgi_reproj_algn, 
                            rcl = matrix(c(
                              NA, 1, # NA
                              2, 2, # Forest
                              3, 3, # Urban
                              4, 4, # Water
                              5, 5, # Cropland
                              6, 6, # Native
                              7, 7, # Shrub
                              8, 8, # Tame
                              9, 9, # Mixed
                              10, 10), # Wetland
                              ncol = 2, byrow = TRUE))

mgi_algn_fix
plot(mgi_algn_fix)

# reclassify to fix
mgi_algn_fix2 <- classify(mgi_algn_fix, 
                         rcl = matrix(c(
                           0, NA, # NA
                           1, NA, # NA
                           2, 2, # Forest
                           3, 3, # Urban
                           4, 4, # Water
                           5, 5, # Cropland
                           6, 6, # Native
                           7, 7, # Shrub
                           8, 8, # Tame
                           9, 9, # Mixed
                           10, 10), # Wetland
                           ncol = 2, byrow = TRUE))

mgi_algn_fix2
plot(mgi_algn_fix2)

# Make categorical
mgi_algn_final <- as.factor(mgi_algn_fix2)
mgi_algn_final
plot(mgi_algn_final, main = "Corrected aligned MGI raster")

# Save the raster as 8-bit unsigned integer
mgi_out_path <- file.path(mgi_aligned_dir, "mgi_reproj_algn_fix.tif")
writeRaster(mgi_algn_final, mgi_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (0 = not grass, 1 = grass)
mgi_algn_binary <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 1, # Native
                              7, 0, # Shrub
                              8, 1, # Tame
                              9, 1, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

mgi_algn_binary
plot(mgi_algn_binary)

# Save the raster as 8-bit unsigned integer
mgi_binary_out_path <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_binary_noshrub.tif")
writeRaster(mgi_algn_binary, mgi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (0 = not grass, 1 = grass) - WITH SHRUB
mgi_algn_final_file <- file.path(mgi_aligned_dir, "mgi_reproj_algn_fix.tif")
mgi_algn_final <- rast(mgi_algn_final_file)
mgi_algn_final
plot(mgi_algn_final)

mgi_algn_binary_shrub <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 1, # Native
                              7, 1, # Shrub
                              8, 1, # Tame
                              9, 1, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

mgi_algn_binary_shrub
plot(mgi_algn_binary_shrub)

# Save the raster as 8-bit unsigned integer
mgi_binary_out_path <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_binary_withshrub.tif")
writeRaster(mgi_algn_binary_shrub, mgi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (0 = not grass, 1 = grass) - WITH SHRUB (SEP)
mgi_algn_binary_shrub_sep <- classify(mgi_algn_final, 
                                  rcl = matrix(c(
                                    2, 0, # Forest
                                    3, 0, # Urban
                                    4, 0, # Water
                                    5, 0, # Cropland
                                    6, 1, # Native
                                    7, 3, # Shrub
                                    8, 1, # Tame
                                    9, 1, # Mixed
                                    10, 0), # Wetland
                                    ncol = 2, byrow = TRUE))

mgi_algn_binary_shrub_sep
plot(mgi_algn_binary_shrub_sep)

# Save the raster as 8-bit unsigned integer
mgi_binary_out_path <- file.path(mgi_aligned_dir, "mb_mgi_reproj_algn_binary_withshrub_sep.tif")
writeRaster(mgi_algn_binary_shrub_sep, mgi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
mgi_algn_nattam <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 2, # Native
                              7, 0, # Shrub
                              8, 1, # Tame
                              9, 1, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

plot(mgi_algn_nattam, main = "Native vs tame MGI raster")
mgi_algn_nattam

# Save the raster as 8-bit unsigned integer
mgi_nattam_out_path <- file.path(mgi_aligned_dir, "mb_mgi_merged_reproj_algn_natvtam_noshrub.tif")
writeRaster(mgi_algn_nattam, mgi_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native) - WITH SHRUB
mgi_algn_nattam_shrub <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 2, # Native
                              7, 3, # Shrub
                              8, 1, # Tame
                              9, 1, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

plot(mgi_algn_nattam_shrub, main = "Native vs tame MGI raster with shrub")
mgi_algn_nattam_shrub

# Save the raster as 8-bit unsigned integer
mgi_nattam_out_path <- file.path(mgi_aligned_dir, "mb_mgi_merged_reproj_algn_natvtam_withshrub.tif")
writeRaster(mgi_algn_nattam_shrub, mgi_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# raster is outside of pilot area

############################################################################

# Isolate classes

# load mgi raster
mgi <- file.path(mgi_aligned_dir, "mgi_reproj_algn_fix.tif")
mgi_algn_final <- rast(mgi)

# Isolate native
mgi_algn_native <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 1, # Native
                              7, 0, # Shrub
                              8, 0, # Tame
                              9, 0, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

mgi_algn_native
plot(mgi_algn_native, main = "MGI native class")

# Save the raster as 8-bit unsigned integer
mgi_nat_out_path <- file.path(mgi_aligned_dir, "mgi_reproj_algn_native.tif")
writeRaster(mgi_algn_native, mgi_nat_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Isolate mixed
mgi_algn_mixed <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 0, # Native
                              7, 0, # Shrub
                              8, 0, # Tame
                              9, 1, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

mgi_algn_mixed
plot(mgi_algn_mixed, main = "MGI mixed class")

# Save the raster as 8-bit unsigned integer
mgi_mix_out_path <- file.path(mgi_aligned_dir, "mgi_reproj_algn_mixed.tif")
writeRaster(mgi_algn_mixed, mgi_mix_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Isolate tame
mgi_algn_tame <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 0, # Native
                              7, 0, # Shrub
                              8, 1, # Tame
                              9, 0, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

mgi_algn_tame
plot(mgi_algn_tame, main = "MGI tame class")

# Save the raster as 8-bit unsigned integer
mgi_tam_out_path <- file.path(mgi_aligned_dir, "mgi_reproj_algn_tame.tif")
writeRaster(mgi_algn_tame, mgi_tam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Isolate shrub
mgi_algn_shrub <- classify(mgi_algn_final, 
                            rcl = matrix(c(
                              2, 0, # Forest
                              3, 0, # Urban
                              4, 0, # Water
                              5, 0, # Cropland
                              6, 0, # Native
                              7, 1, # Shrub
                              8, 0, # Tame
                              9, 0, # Mixed
                              10, 0), # Wetland
                              ncol = 2, byrow = TRUE))

mgi_algn_shrub
plot(mgi_algn_shrub, main = "MGI shrub class")

# Save the raster as 8-bit unsigned integer
mgi_shrub_out_path <- file.path(mgi_aligned_dir, "mgi_reproj_algn_shrub.tif")
writeRaster(mgi_algn_shrub, mgi_shrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
