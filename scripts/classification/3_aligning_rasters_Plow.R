
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

# plow directory
plow_root <- file.path(raw_grass_data_dir, "./Cropland/Plowprint WWF/WWF_Plowprint_2024/Spatial_data/WWF_PP_2024_report_full_extent")
plow_aligned_dir <- file.path(proc_grass_data_dir, "wwf_plowprint/aligned")
plow_cropped_dir <- file.path(proc_grass_data_dir, "wwf_plowprint/cropped")
plow_pilot_dir <- file.path(proc_grass_data_dir, "wwf_plowprint/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# plow raster
plow_file <- file.path(plow_root, "WWF_PP_2024_report_full_extent.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- WWF Plowprint ----#

# Load rasters
r_plow <- rast(plow_file)

origin(r_plow) # -6.991082 -5.858173
r_plow

#class       : SpatRaster 
#dimensions  : 126814, 153709, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2354947, 2256323, 311784.1, 4116204  (xmin, xmax, ymin, ymax)
#coord. ref. : NAD83 / Conus Albers (EPSG:5070) 
#source      : WWF_PP_2024_report_full_extent.tif 
#categories  : Display 
#name        :     Display 
#min value   :      Intact 
#max value   : Roads/Water

# print frequency table
freq_table <- freq(r_plow) 
print(freq_table) 

# Mark as categorical (factor)
r_plow <- as.factor(r_plow)

# crop before reprojecting to make manageable size 

# Define the bbox in lon/lat degrees (WGS84)
canada_ext_ll  <- ext(-140, -70, 30, 60)
canada_bbox_ll <- as.polygons(canada_ext_ll, crs = "EPSG:4326")

# Get CRS from the raster 
r_crs <- crs(r_plow)

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
r_plow_crop <- crop(r_plow, canada_ext_rcrs)
r_plow_mask <- mask(r_plow_crop, canada_bbox_rcrs)

r_plow_mask
plot(r_plow_mask)

# Save the raster as 8-bit unsigned integer
plow_crop_out_path <- file.path(plow_cropped_dir, "plow_crop_bbox.tif")
writeRaster(r_plow_mask, plow_crop_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reproject before resampling
plow_reproj <- project(r_plow, crs(sa)) # or should be project(ngi_reclass, sa) ??

freq_table2 <- freq(plow_reproj) 
print(freq_table2) 

origin(plow_reproj) # 1.8392815 0.4291859
plow_reproj

# align to reference raster 
plow_reproj_algn <- resample(plow_reproj, sa, method = "near")
plot(plow_reproj_algn)

# Get values and counts
freq_table3 <- freq(plow_reproj_algn) 
print(freq_table3) 

origin(plow_reproj_algn) # 6.696121 -2.550334
plow_reproj_algn

# Save the raster as 8-bit unsigned integer
plow_out_path <- file.path(plow_aligned_dir, "plow_reproj_algn.tif")
writeRaster(plow_reproj_algn, plow_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
plow_reproj_algn_trim <- trim(plow_reproj_algn)

# visualize the trimmed raster
plot(plow_reproj_algn_trim)

# print raster details 
origin(plow_reproj_algn_trim) # 6.696121 -2.550334
plow_reproj_algn_trim

# Save the raster as 8-bit unsigned integer
plow_trim_out_path <- file.path(plow_aligned_dir, "plow_reproj_algn_trim.tif")
writeRaster(plow_reproj_algn_trim, plow_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (grass = 1)
plow_algn_binary <- classify(plow_reproj_algn, 
                            rcl = matrix(c(0, 1, # Intact
                                           1, 1, # Intact
                                           2, 0, # New plowprint
                                           3, 0, # Plowprint
                                           40, 0, # Developed
                                           60, 0, # Developed
                                           100, 0), # Roads/water
                                         ncol = 2, byrow = TRUE))
plot(plow_algn_binary)
plow_algn_binary

# Save the raster as 8-bit unsigned integer
plow_binary_out_path <- file.path(plow_aligned_dir, "plow_reproj_algn_binary.tif")
writeRaster(plow_algn_binary, plow_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334

# crop to pilot area
crop_ext <- ext(pa)
plow_pa <- crop(plow_reproj_algn, crop_ext)
plow_pa
plot (plow_pa)

# mask to pilot area
plow_masked <- mask(plow_pa, pa)
plow_masked
plot(plow_masked)

# print raster details 
origin(plow_masked) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
plow_pa_out_path <- file.path(plow_pilot_dir, "plow_reproj_algn_crop_pa.tif")
writeRaster(plow_masked, plow_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
plow_bin_pa <- crop(plow_algn_binary, crop_ext)
plow_bin_pa
plot (plow_bin_pa)

# mask to pilot area
plow_bin_masked <- mask(plow_bin_pa, pa)
plow_bin_masked
plot(plow_bin_masked)

# Save the raster as 8-bit unsigned integer
plow_bin_pa_out_path <- file.path(plow_pilot_dir, "plow_reproj_algn_crop_pa_binary.tif")
writeRaster(plow_bin_masked, plow_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
