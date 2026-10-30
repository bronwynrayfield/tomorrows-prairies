
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

# pglc tiles directory
pglc_tile_root <- file.path(raw_grass_data_dir, "Grassland/AAFC_Prairie_Grassland")
#pglc_cropped_dir <- file.path(proc_grass_data_dir, "aafc_pglc/cropped")
pglc_aligned_dir <- file.path(proc_grass_data_dir, "aafc_pglc/aligned")
pglc_pilot_dir <- file.path(proc_grass_data_dir, "aafc_pglc/pilot")

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# pglc raster
pglc_file <- file.path(pglc_tile_root, "Grassland_classification/AAFC_Grassland_classification_2020_V1.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
sa

####--- ALGIN RASTERS ---#### 

#----- AAFC Prairie Grassland Land Cover 2020 ----#

# Load rasters
r_pglc <- rast(pglc_file)

r_pglc
origin(r_pglc) # 0 0 
plot(r_pglc, main = "Raw PGLC raster")

#class       : SpatRaster
#size        : 51505, 84760, 1  (nrow, ncol, nlyr)
#resolution  : 10, 10  (x, y)
#extent      : -1280140, -432540, 1016040, 1531090  (xmin, xmax, ymin, ymax)
#coord. ref. : Albers_Conic_Equal_Area
#source      : AAFC_Grassland_classification_2020_V1.tif
#color table : 1
#categories  : Class, Red, Green, Blue, Alpha
#name        : Class
#min value   :   100
#max value   :   122

# print frequency table
freq_table <- freq(r_pglc) 
print(freq_table)

#  layer value      count
#1     1   100 1613779554
#2     1   110  464732877
#3     1   122  454546095

# reproject before resampling
pglc_reproj <- project(r_pglc, sa, method = "near") # better approach

# align to reference raster 
pglc_reproj_algn <- resample(pglc_reproj, sa, method = "near")

# inspect
plot(pglc_reproj_algn, main = "Aligned PGLC raster")
origin(pglc_reproj_algn) 
pglc_reproj_algn

#class       : SpatRaster
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008)
#source      : spat_87281e703886_34600_zxzxx9jyyx6UnnI.tif
#color table : 1
#varname     : class_study_area_aligned_cnvc_Dec2025
#categories  : Class, Red, Green, Blue, Alpha
#name        : Class
#min value   :   100
#max value   :   122

# Get values and counts
freq_table2 <- freq(pglc_reproj_algn) 
print(freq_table2) 

#  layer value     count
#1     1   100 179307864
#2     1   110  51634502
#3     1   122  50508304

# Save the raster as 8-bit unsigned integer
pglc_out_path <- file.path(pglc_aligned_dir, "pglc_reproj_algn.tif")
writeRaster(pglc_reproj_algn, pglc_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (0 = not grass, 1 = grass)
pglc_algn_binary <- classify(pglc_reproj_algn, 
                            rcl = matrix(c(
                              100, 0, # Other land cover (aggregated non-target classes)
                              110, 1, # Native grassland
                              122, 1), # Seeded forage/high-disturbance grassland
                              ncol = 2, byrow = TRUE))

pglc_algn_binary
plot(pglc_algn_binary, main = "PGLC binary")

# Save the raster as 8-bit unsigned integer
pglc_binary_out_path <- file.path(pglc_aligned_dir, "pglc_reproj_algn_binary.tif")
writeRaster(pglc_algn_binary, pglc_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
pglc_algn_nattam <- classify(pglc_reproj_algn, 
                             rcl = matrix(c(
                               100, 0, # Other land cover (aggregated non-target classes)
                               110, 2, # Native grassland
                               122, 1), # Seeded forage/high-disturbance grassland
                               ncol = 2, byrow = TRUE))

plot(pglc_algn_nattam, main = "Native vs tame PGLC raster")
pglc_algn_nattam

# Save the raster as 8-bit unsigned integer
pglc_nattam_out_path <- file.path(pglc_aligned_dir, "pglc_reproj_algn_natvtam.tif")
writeRaster(pglc_algn_nattam, pglc_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

############################################################################

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334
res(pa) # 30 30
crs(pa)

# crop to pilot area
crop_ext <- ext(pa)
pglc_pa <- crop(pglc_reproj_algn, crop_ext)
pglc_pa
plot (pglc_pa)

# mask to pilot area
pglc_masked <- mask(pglc_pa, pa)
pglc_masked
plot(pglc_masked)

# print raster details 
origin(pglc_masked) # 6.696121 -2.550334
res(pglc_masked) # 30 30
crs(pglc_masked)

# Save the raster as 8-bit unsigned integer
pglc_pa_out_path <- file.path(pglc_pilot_dir, "pglc_reproj_algn_crop_pa.tif")
writeRaster(pglc_masked, pglc_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
pglc_bin_pa <- crop(pglc_algn_binary, crop_ext)
pglc_bin_pa
plot (pglc_bin_pa)

# mask to pilot area
pglc_bin_masked <- mask(pglc_bin_pa, pa)
pglc_bin_masked
plot(pglc_bin_masked)

# Save the raster as 8-bit unsigned integer
pglc_bin_pa_out_path <- file.path(pglc_pilot_dir, "pglc_reproj_algn_crop_pa_binary.tif")
writeRaster(pglc_bin_masked, pglc_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
pglc_nt_pa <- crop(pglc_algn_nattam, crop_ext)
pglc_nt_pa
plot (pglc_nt_pa)

# mask to pilot area
pglc_nt_masked <- mask(pglc_nt_pa, pa)
pglc_nt_masked
plot(pglc_nt_masked)

# Save the raster as 8-bit unsigned integer
pglc_nt_pa_out_path <- file.path(pglc_pilot_dir, "pglc_reproj_algn_crop_pa_natvtam.tif")
writeRaster(pglc_nt_masked, pglc_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

############################################################################

# Isolate classes

# load pglc raster
#pglc <- file.path(pglc_aligned_dir, "pglc_reproj_algn.tif")
#pglc_algn_final <- rast(pglc)

# Isolate native
pglc_algn_native <- classify(pglc_reproj_algn, 
                             rcl = matrix(c(
                               100, 0, # Other land cover (aggregated non-target classes)
                               110, 1, # Native grassland
                               122, 0), # Seeded forage/high-disturbance grassland
                               ncol = 2, byrow = TRUE))

pglc_algn_native
plot(pglc_algn_native, main = "PGLC native class")

# Save the raster as 8-bit unsigned integer
pglc_nat_out_path <- file.path(pglc_aligned_dir, "pglc_reproj_algn_native.tif")
writeRaster(pglc_algn_native, pglc_nat_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Isolate tame
pglc_algn_tame <- classify(pglc_reproj_algn, 
                             rcl = matrix(c(
                               100, 0, # Other land cover (aggregated non-target classes)
                               110, 0, # Native grassland
                               122, 1), # Seeded forage/high-disturbance grassland
                               ncol = 2, byrow = TRUE))

pglc_algn_tame
plot(pglc_algn_tame, main = "PGLC tame class")

# Save the raster as 8-bit unsigned integer
pglc_tam_out_path <- file.path(pglc_aligned_dir, "pglc_reproj_algn_tame.tif")
writeRaster(pglc_algn_tame, pglc_tam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
