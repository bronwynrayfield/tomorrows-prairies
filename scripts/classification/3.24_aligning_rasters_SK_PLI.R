
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

# PLI tiles directory
pli_tile_root <- file.path(raw_grass_data_dir, "Grassland/SK Prairie Landscape Inventory")
pli_merged_dir <- file.path(proc_grass_data_dir, "sk_pli/merged")
pli_cropped_dir <- file.path(proc_grass_data_dir, "sk_pli/cropped")
pli_aligned_dir <- file.path(proc_grass_data_dir, "sk_pli/aligned")
pli_pilot_dir <- file.path(proc_grass_data_dir, "sk_pli/pilot")

####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# pli rasters
pli_ap_file <- file.path(pli_tile_root, "Aspen Parkland/Aspen_Classified/Aspen_Classified.tif")
pli_cu_file <- file.path(pli_tile_root, "Cypress Upland/Cypress_Classified/Cypress_Classified.tif")
pli_mg_file <- file.path(pli_tile_root, "Mixed Grassland/PLI_Classification/PLI_Classification.tif")
pli_mm_file <- file.path(pli_tile_root, "Moist Mixed Grassland/MMG_Classified/MMG_Classified.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
sa

####--- ALGIN RASTERS ---#### 

#----- SK Prairie Landscape Inventory ----#

# Load rasters
r_ap <- rast(pli_ap_file)
r_cu <- rast(pli_cu_file)
r_mg <- rast(pli_mg_file)
r_mm <- rast(pli_mm_file)

r_ap
#class       : SpatRaster 
#size        : 67372, 84254, 1  (nrow, ncol, nlyr)
#resolution  : 7.989028, 7.989028  (x, y)
#extent      : 121403, 794510.5, 5422288, 5960525  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83 / UTM zone 13N (EPSG:26913) 
#source      : Aspen_Classified.tif 
#name        : Aspen_Classified 
#min value   :                1 
#max value   :                9

r_cu
#class       : SpatRaster 
#size        : 6495, 11211, 1  (nrow, ncol, nlyr)
#resolution  : 9.60993, 9.60993  (x, y)
#extent      : 137444.9, 245181.8, 5479745, 5542162  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83(CSRS) / UTM zone 13N (EPSG:2957) 
#source      : Cypress_Classified.tif 
#name        : Cypress_Classified 
#min value   :                  1 
#max value   :                  9 

r_mg
#class       : SpatRaster 
#size        : 87017, 75666, 1  (nrow, ncol, nlyr)
#resolution  : 9.999826, 9.999826  (x, y)
#extent      : 34979.03, 791625.9, 5350399, 6220554  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83(CSRS) / UTM zone 13N (EPSG:2957) 
#source      : PLI_Classification.tif 
#categories  : Name, Area_ha, Year, Count 
#name        :     Name 
#min value   : Cropland 
#max value   :    Trees

r_mm
#class       : SpatRaster 
#size        : 51175, 70478, 1  (nrow, ncol, nlyr)
#resolution  : 8.067117, 8.067117  (x, y)
#extent      : 154033.2, 722587.4, 5430163, 5842997  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83 / UTM zone 13N (EPSG:26913) 
#source      : MMG_Classified.tif 
#name        : MMG_Classified 
#min value   :              0 
#max value   :              9 

origin(r_ap) # 1.6935202 -0.7674738
origin(r_cu) # 3.658549 0.266043
origin(r_mg) # -0.360654  1.675645
origin(r_mm) # -0.3441006  0.8031349

plot(r_ap)
plot(r_cu)
plot(r_mg)
plot(r_mm)

# Set pixels with value 0 to NA 
r_mm[r_mm == 0] <- NA

plot(r_mm)
r_mm

# rasters have different origins and crs - need to align before merging
cu_algn <- project(r_cu, r_ap, method = "near", align = TRUE)
mg_algn <- project(r_mg, r_ap, method = "near", align = TRUE)
mm_algn <- project(r_mm, r_ap, method = "near", align = TRUE)

# check alignment
cu_algn
mg_algn
mm_algn

origin(r_ap) # 1.6935202 -0.7674738
origin(cu_algn) # 1.6935202 -0.7674738
origin(mg_algn) # 1.6935202 -0.7674738
origin(mm_algn) # 1.6935202 -0.7674738

# Merge 
#merged_pli <- mosaic(r_ap, r_cu, r_mg, r_mm, fun = "max") 
merged_pli <- merge(r_ap, cu_algn, mg_algn, mm_algn) 

origin(merged_pli)  # 1.6935202 -0.7674738
merged_pli

#class       : SpatRaster 
#size        : 108920, 95072, 1  (nrow, ncol, nlyr)
#resolution  : 7.989028, 7.989028  (x, y)
#extent      : 34977.66, 794510.5, 5350395, 6220560  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83 / UTM zone 13N (EPSG:26913) 
#source      : spat_3e688283376d_255624_aYjuhHgkOkYNCZW.tif 
#varname     : Aspen_Classified 
#name        : Aspen_Classified 
#min value   :                1 
#max value   :                9

# plot the result
plot(merged_pli, main = "Merged pli raster")

# print frequency table
freq_table <- freq(merged_pli) 
print(freq_table) 

#   layer value      count
#2      1     1 2054030976
#3      1     2  615138375
#4      1     3  142089106
#5      1     4  350088751
#6      1     5  276518375
#7      1     6   89878615
#8      1     7  156225101
#9      1     8   65178769
#10     1     9   17533578

# Save the raster as 8-bit unsigned integer
pli_merge_out_path <- file.path(pli_merged_dir, "sk_pli_merged.tif")
writeRaster(merged_pli, pli_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Mark as categorical (factor)
merged_pli <- as.factor(merged_pli)
plot(merged_pli)

# reproject before resampling
pli_merged_reproj <- project(merged_pli, sa, method = "near") 
pli_merged_reproj

freq_table2 <- freq(pli_merged_reproj) 
print(freq_table2) 

# align to reference raster 
pli_merged_reproj_algn <- resample(pli_merged_reproj, sa, method = "near")
plot(pli_merged_reproj_algn)
pli_merged_reproj_algn

# Get values and counts
freq_table3 <- freq(pli_merged_reproj_algn) 
print(freq_table3) 

# Save the raster as 8-bit unsigned integer
pli_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn.tif")
writeRaster(pli_merged_reproj_algn, pli_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
pli_merged_reproj_algn_trim <- trim(pli_merged_reproj_algn)

# visualize the trimmed raster
plot(pli_merged_reproj_algn_trim)

# Save the raster as 8-bit unsigned integer
pli_trim_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_trim.tif")
writeRaster(pli_merged_reproj_algn_trim, pli_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (0 = not grass, 1 = grass)
pli_algn_binary <- classify(pli_merged_reproj_algn, 
                            rcl = matrix(c(
                              1, 0, # Cropland
                              2, 1, # Native grassland
                              3, 1, # Mixed/altered grassland
                              4, 1, # Tame grassland
                              5, 0, # Water
                              6, 0, # Shrubs
                              7, 0, # Trees
                              8, 0, # Woody plants
                              9, 0), # Urban area
                              ncol = 2, byrow = TRUE))

pli_algn_binary

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_3e6887bb65f2c_255624_6Uxy02myoMrg5LO.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#name        : Aspen_Classified 
#min value   :                0 
#max value   :                1 

# Save the raster as 8-bit unsigned integer
pli_binary_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_binary_noshrub.tif")
writeRaster(pli_algn_binary, pli_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load pli raster
pli <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn.tif")
pli_algn <- rast(pli)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
pli_algn_nattam <- classify(pli_algn, 
                            rcl = matrix(c(
                              1, 0, # Cropland
                              2, 2, # Native grassland
                              3, 1, # Mixed/altered grassland
                              4, 1, # Tame grassland
                              5, 0, # Water
                              6, 0, # Shrubs
                              7, 0, # Trees
                              8, 0, # Woody plants
                              9, 0), # Urban area
                              ncol = 2, byrow = TRUE))
plot(pli_algn_nattam)
pli_algn_nattam

# Save the raster as 8-bit unsigned integer
pli_nattam_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_natvtam_noshrub.tif")
writeRaster(pli_algn_nattam, pli_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

#### ISOLATE CLASSES ####

# load pli raster
pli <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn.tif")
pli_algn <- rast(pli)

pli_algn
plot(pli_algn)

# isolate native
pli_algn_nat <- classify(pli_algn, 
                            rcl = matrix(c(
                              1, 0, # Cropland
                              2, 1, # Native grassland
                              3, 0, # Mixed/altered grassland
                              4, 0, # Tame grassland
                              5, 0, # Water
                              6, 0, # Shrubs
                              7, 0, # Trees
                              8, 0, # Woody plants
                              9, 0), # Urban area
                              ncol = 2, byrow = TRUE))
plot(pli_algn_nat)
pli_algn_nat

# Save the raster as 8-bit unsigned integer
pli_nat_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_native.tif")
writeRaster(pli_algn_nat, pli_nat_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# isolate tame (includes mixed)
pli_algn_tammix <- classify(pli_algn, 
                         rcl = matrix(c(
                           1, 0, # Cropland
                           2, 0, # Native grassland
                           3, 1, # Mixed/altered grassland
                           4, 1, # Tame grassland
                           5, 0, # Water
                           6, 0, # Shrubs
                           7, 0, # Trees
                           8, 0, # Woody plants
                           9, 0), # Urban area
                           ncol = 2, byrow = TRUE))
plot(pli_algn_tammix)
pli_algn_tammix

# Save the raster as 8-bit unsigned integer
pli_tammix_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_tameANDmixed.tif")
writeRaster(pli_algn_tammix, pli_tammix_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# isolate tame (no mixed)
pli_algn_tam <- classify(pli_algn, 
                         rcl = matrix(c(
                           1, 0, # Cropland
                           2, 0, # Native grassland
                           3, 0, # Mixed/altered grassland
                           4, 1, # Tame grassland
                           5, 0, # Water
                           6, 0, # Shrubs
                           7, 0, # Trees
                           8, 0, # Woody plants
                           9, 0), # Urban area
                           ncol = 2, byrow = TRUE))
plot(pli_algn_tam)
pli_algn_tam

# Save the raster as 8-bit unsigned integer
pli_tam_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_tame_nomixed.tif")
writeRaster(pli_algn_tam, pli_tam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# isolate mixed
pli_algn_mix <- classify(pli_algn, 
                         rcl = matrix(c(
                           1, 0, # Cropland
                           2, 0, # Native grassland
                           3, 1, # Mixed/altered grassland
                           4, 0, # Tame grassland
                           5, 0, # Water
                           6, 0, # Shrubs
                           7, 0, # Trees
                           8, 0, # Woody plants
                           9, 0), # Urban area
                           ncol = 2, byrow = TRUE))
plot(pli_algn_mix)
pli_algn_mix

# Save the raster as 8-bit unsigned integer
pli_mix_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_mixed.tif")
writeRaster(pli_algn_mix, pli_mix_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# isolate shrub
pli_algn_shrub <- classify(pli_algn, 
                         rcl = matrix(c(
                           1, 0, # Cropland
                           2, 0, # Native grassland
                           3, 0, # Mixed/altered grassland
                           4, 0, # Tame grassland
                           5, 0, # Water
                           6, 1, # Shrubs
                           7, 0, # Trees
                           8, 0, # Woody plants
                           9, 0), # Urban area
                           ncol = 2, byrow = TRUE))
plot(pli_algn_shrub)
pli_algn_shrub

# Save the raster as 8-bit unsigned integer
pli_shrub_out_path <- file.path(pli_aligned_dir, "sk_pli_merged_reproj_algn_shrub.tif")
writeRaster(pli_algn_shrub, pli_shrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

#### PILOT AREA RASTERS ####

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa)

# crop to pilot area
crop_ext <- ext(pa)
pli_pa <- crop(pli_merged_reproj_algn, crop_ext)
pli_pa
plot (pli_pa)

# mask to pilot area
pli_masked <- mask(pli_pa, pa)

# print raster details 
origin(pli_masked) # 
pli_masked
plot(pli_masked)

# Save the raster as 8-bit unsigned integer
pli_pa_out_path <- file.path(pli_pilot_dir, "pli_merged_reproj_algn_crop_pa.tif")
writeRaster(pli_masked, pli_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
pli_bin_pa <- crop(pli_algn_binary, crop_ext)
pli_bin_pa
plot (pli_bin_pa)

# mask to pilot area
pli_bin_masked <- mask(pli_bin_pa, pa)
pli_bin_masked
plot(pli_bin_masked)

# Save the raster as 8-bit unsigned integer
pli_bin_pa_out_path <- file.path(pli_pilot_dir, "pli_merged_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(pli_bin_masked, pli_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
pli_nt_pa <- crop(pli_algn_nattam, crop_ext)
pli_nt_pa
plot (pli_nt_pa)

# mask to pilot area
pli_nt_masked <- mask(pli_nt_pa, pa)
pli_nt_masked
plot(pli_nt_masked)

# Save the raster as 8-bit unsigned integer
pli_nt_pa_out_path <- file.path(pli_pilot_dir, "pli_merged_reproj_algn_crop_pa_natvtam_noshrub.tif")
writeRaster(pli_nt_masked, pli_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

############################################################################

# align each ecoregion raster to study area for use in accuracy assessment

# Load rasters
r_ap <- rast(pli_ap_file)
r_cu <- rast(pli_cu_file)
r_mg <- rast(pli_mg_file)
r_mm <- rast(pli_mm_file)

# Set pixels with value 0 to NA 
r_mm[r_mm == 0] <- NA

# Mark as categorical (factor)
r_ap <- as.factor(r_ap)
r_cu <- as.factor(r_cu)
r_mg <- as.factor(r_mg)
r_mm <- as.factor(r_mm)

plot(r_ap)
plot(r_cu)
plot(r_mg)
plot(r_mm)

# reproject before resampling
r_ap_reproj <- project(r_ap, sa, method = "near") 
r_cu_reproj <- project(r_cu, sa, method = "near")
r_mg_reproj <- project(r_mg, sa, method = "near")
r_mm_reproj <- project(r_mm, sa, method = "near")

r_ap_reproj
r_cu_reproj
r_mg_reproj
r_mm_reproj

# align to reference raster 
r_ap_reproj_algn <- resample(r_ap_reproj, sa, method = "near")
r_cu_reproj_algn <- resample(r_cu_reproj, sa, method = "near")
r_mg_reproj_algn <- resample(r_mg_reproj, sa, method = "near")
r_mm_reproj_algn <- resample(r_mm_reproj, sa, method = "near")

r_ap_reproj_algn
r_cu_reproj_algn
r_mg_reproj_algn
r_mm_reproj_algn

# Save the raster as 8-bit unsigned integer
pli_ap_out_path <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn.tif")
writeRaster(r_ap_reproj_algn, pli_ap_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_cu_out_path <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn.tif")
writeRaster(r_cu_reproj_algn, pli_cu_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_mg_out_path <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn.tif")
writeRaster(r_mg_reproj_algn, pli_mg_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_mm_out_path <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn.tif")
writeRaster(r_mm_reproj_algn, pli_mm_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (0 = not grass, 1 = grass)
pli_ap_algn_binary <- classify(r_ap_reproj_algn, 
                            rcl = matrix(c(
                              1, 0, # Cropland
                              2, 1, # Native grassland
                              3, 1, # Mixed/altered grassland
                              4, 1, # Tame grassland
                              5, 0, # Water
                              6, 0, # Shrubs
                              7, 0, # Trees
                              #8, 0, # Woody plants
                              9, 0), # Urban area
                              ncol = 2, byrow = TRUE))

pli_cu_algn_binary <- classify(r_cu_reproj_algn, 
                               rcl = matrix(c(
                                 1, 0, # Cropland
                                 2, 1, # Native grassland
                                 #3, 1, # Mixed/altered grassland
                                 4, 1, # Tame grassland
                                 5, 0, # Water
                                 6, 0, # Shrubs
                                 7, 0, # Trees
                                 #8, 0, # Woody plants
                                 9, 0), # Urban area
                                 ncol = 2, byrow = TRUE))

pli_mg_algn_binary <- classify(r_mg_reproj_algn, 
                               rcl = matrix(c(
                                 1, 0, # Cropland
                                 2, 1, # Native grassland
                                 3, 1, # Mixed/altered grassland
                                 4, 1, # Tame grassland
                                 5, 0, # Water
                                 6, 0, # Shrubs
                                 7, 0), # Trees
                                 #8, 0, # Woody plants
                                 #9, 0), # Urban area
                                 ncol = 2, byrow = TRUE))

pli_mm_algn_binary <- classify(r_mm_reproj_algn, 
                               rcl = matrix(c(
                                 1, 0, # Cropland
                                 2, 1, # Native grassland
                                 3, 1, # Mixed/altered grassland
                                 4, 1, # Tame grassland
                                 5, 0, # Water
                                 #6, 0, # Shrubs
                                 #7, 0, # Trees
                                 8, 0, # Woody plants
                                 9, 0), # Urban area
                                 ncol = 2, byrow = TRUE))

pli_ap_algn_binary
plot(pli_ap_algn_binary)
plot(pli_cu_algn_binary)
plot(pli_mg_algn_binary)
plot(pli_mm_algn_binary)

# Save the raster as 8-bit unsigned integer
pli_ap_binary_out_path <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_binary_noshrub.tif")
writeRaster(pli_ap_algn_binary, pli_ap_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_cu_binary_out_path <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_binary_noshrub.tif")
writeRaster(pli_cu_algn_binary, pli_cu_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_mg_binary_out_path <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_binary_noshrub.tif")
writeRaster(pli_mg_algn_binary, pli_mg_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_mm_binary_out_path <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_binary_noshrub.tif")
writeRaster(pli_mm_algn_binary, pli_mm_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
pli_ap_algn_nattam <- classify(r_ap_reproj_algn, 
                               rcl = matrix(c(
                                 1, 0, # Cropland
                                 2, 2, # Native grassland
                                 3, 1, # Mixed/altered grassland
                                 4, 1, # Tame grassland
                                 5, 0, # Water
                                 6, 0, # Shrubs
                                 7, 0, # Trees
                                 #8, 0, # Woody plants
                                 9, 0), # Urban area
                                 ncol = 2, byrow = TRUE))

pli_cu_algn_nattam <- classify(r_cu_reproj_algn, 
                               rcl = matrix(c(
                                 1, 0, # Cropland
                                 2, 2, # Native grassland
                                 #3, 1, # Mixed/altered grassland
                                 4, 1, # Tame grassland
                                 5, 0, # Water
                                 6, 0, # Shrubs
                                 7, 0, # Trees
                                 #8, 0, # Woody plants
                                 9, 0), # Urban area
                                 ncol = 2, byrow = TRUE))

pli_mg_algn_nattam <- classify(r_mg_reproj_algn, 
                               rcl = matrix(c(
                                 1, 0, # Cropland
                                 2, 2, # Native grassland
                                 3, 1, # Mixed/altered grassland
                                 4, 1, # Tame grassland
                                 5, 0, # Water
                                 6, 0, # Shrubs
                                 7, 0), # Trees
                                 #8, 0, # Woody plants
                                 #9, 0), # Urban area
                                 ncol = 2, byrow = TRUE))

pli_mm_algn_nattam <- classify(r_mm_reproj_algn, 
                               rcl = matrix(c(
                                 1, 0, # Cropland
                                 2, 2, # Native grassland
                                 3, 1, # Mixed/altered grassland
                                 4, 1, # Tame grassland
                                 5, 0, # Water
                                 #6, 0, # Shrubs
                                 #7, 0, # Trees
                                 8, 0, # Woody plants
                                 9, 0), # Urban area
                                 ncol = 2, byrow = TRUE))

pli_ap_algn_nattam
pli_cu_algn_nattam
pli_mg_algn_nattam
pli_mm_algn_nattam

# Save the raster as 8-bit unsigned integer
pli_ap_nattam_out_path <- file.path(pli_aligned_dir, "sk_pli_ap_reproj_algn_natvtam_noshrub.tif")
writeRaster(pli_ap_algn_nattam, pli_ap_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_cu_nattam_out_path <- file.path(pli_aligned_dir, "sk_pli_cu_reproj_algn_natvtam_noshrub.tif")
writeRaster(pli_cu_algn_nattam, pli_cu_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_mg_nattam_out_path <- file.path(pli_aligned_dir, "sk_pli_mg_reproj_algn_natvtam_noshrub.tif")
writeRaster(pli_mg_algn_nattam, pli_mg_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

pli_mm_nattam_out_path <- file.path(pli_aligned_dir, "sk_pli_mm_reproj_algn_natvtam_noshrub.tif")
writeRaster(pli_mm_algn_nattam, pli_mm_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
