
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

# lcr directory
lcr_root <- file.path(raw_grass_data_dir, "Land Cover/StatCan Land Cover Register 2020/LCR_RCT_2020")
lcr_aligned_dir <- file.path(proc_grass_data_dir, "statcan_lcr/aligned")
lcr_cropped_dir <- file.path(proc_grass_data_dir, "statcan_lcr/cropped")
lcr_masked_dir <- file.path(proc_grass_data_dir, "statcan_lcr/masked")
lcr_pilot_dir <- file.path(proc_grass_data_dir, "statcan_lcr/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# lcr raster
lcr_file <- file.path(lcr_root, "LCR_RCT_2020.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- StatCan Land Cover Register ----#

# Load rasters
r_lcr <- rast(lcr_file)

origin(r_lcr) # -2.378897e-11  1.010392e-11
r_lcr

# print frequency table
freq_table <- freq(r_lcr) 
print(freq_table) 

# Mark as categorical (factor)
r_lcr <- as.factor(r_lcr)

# crop before reprojecting to make manageable size 

# Define the bbox in lon/lat degrees (WGS84)
canada_ext_ll  <- ext(-140, -70, 30, 60)
canada_bbox_ll <- as.polygons(canada_ext_ll, crs = "EPSG:4326")

# Get CRS from the raster 
r_crs <- crs(r_lcr)

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
r_lcr_crop <- crop(r_lcr, canada_ext_rcrs)
r_lcr_mask <- mask(r_lcr_crop, canada_bbox_rcrs)

r_lcr_mask
plot(r_lcr_mask)

# Save the raster as 8-bit unsigned integer
lcr_crop_out_path <- file.path(lcr_cropped_dir, "lcr_crop_bbox.tif")
writeRaster(r_lcr_mask, lcr_crop_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reproject before resampling
lcr_reproj <- project(r_lcr, sa, method = "near") # better approach

freq_table2 <- freq(lcr_reproj) 
print(freq_table2) 

origin(lcr_reproj) # 0.7826532 0.9658756
lcr_reproj

# align to reference raster 
lcr_reproj_algn <- resample(lcr_reproj, sa, method = "near")
plot(lcr_reproj_algn)

# Get values and counts
freq_table3 <- freq(lcr_reproj_algn) 
print(freq_table3) # 

origin(lcr_reproj_algn) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
lcr_out_path <- file.path(lcr_aligned_dir, "lcr_reproj_algn.tif")
writeRaster(lcr_reproj_algn, lcr_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
lcr_reproj_algn_trim <- trim(lcr_reproj_algn)

# visualize the trimmed raster
plot(lcr_reproj_algn_trim)

# print raster details 
origin(lcr_reproj_algn_trim) # 6.696121 -2.550334
res(lcr_reproj_algn_trim) # 30 30 
crs(lcr_reproj_algn_trim) # EPSG:9001
ext(lcr_reproj_algn_trim) # SpatExtent : -1565033.30387889, -10793.3038788857, 1062567.4496659, 2329887.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
lcr_trim_out_path <- file.path(lcr_aligned_dir, "lcr_reproj_algn_trim.tif")
writeRaster(lcr_reproj_algn_trim, lcr_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (grass = 1)
lcr_algn_binary <- classify(lcr_reproj_algn, 
                            rcl = matrix(c(1, 0, # Built-up and artificial surface
                                           2, 0, # Cropland
                                           3, 0, # Inland water body
                                           4, 0, # Treed (non-wetland)
                                           5, 0, # Treed wetland
                                           6, 0, # Treed area disturbance
                                           7, 1, # Grassland and shrubland
                                           8, 0, # Wetland (non-treed)
                                           9, 0, # Sparsely vegetated land
                                           10, 0, # Barren land
                                           11, 0), # Permanent snow and ice
                                         ncol = 2, byrow = TRUE))
plot(lcr_algn_binary)
lcr_algn_binary

# Save the raster as 8-bit unsigned integer
lcr_binary_out_path <- file.path(lcr_aligned_dir, "lcr_reproj_algn_binary.tif")
writeRaster(lcr_algn_binary, lcr_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334
res(pa) # 30 30
crs(pa)

# crop to pilot area
crop_ext <- ext(pa)
lcr_pa <- crop(lcr_reproj_algn, crop_ext)
lcr_pa
plot (lcr_pa)


# mask to pilot area
lcr_masked <- mask(lcr_pa, pa)
lcr_masked
plot(lcr_masked)



# print raster details 
origin(lcr_masked) # 6.696121 -2.550334
res(lcr_masked) # 30 30
crs(lcr_masked)

# Save the raster as 8-bit unsigned integer
lcr_pa_out_path <- file.path(lcr_cropped_dir, "lcr_reproj_algn_crop_pa.tif")
writeRaster(lcr_masked, lcr_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
lcr_bin_pa <- crop(lcr_algn_binary, crop_ext)
lcr_bin_pa
plot (lcr_bin_pa)

# mask to pilot area
lcr_bin_masked <- mask(lcr_bin_pa, pa)
lcr_bin_masked
plot(lcr_bin_masked)

# Save the raster as 8-bit unsigned integer
lcr_bin_pa_out_path <- file.path(lcr_pilot_dir, "lcr_reproj_algn_crop_pa_binary.tif")
writeRaster(lcr_bin_masked, lcr_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################

######### crop and mask to study area

# load raster
lcr_file <- file.path(lcr_aligned_dir, "lcr_reproj_algn.tif")
lcr <- rast(lcr_file)

crop_ext <- ext(sa)

# crop raster to study area
lcr_sa <- crop(lcr, crop_ext)
lcr_sa
plot (lcr_sa)

# mask to study area
lcr_masked <- mask(lcr_sa, sa)
lcr_masked
plot(lcr_masked)

# Save the raster as 8-bit unsigned integer
lcr_sa_out_path <- file.path(lcr_masked_dir, "lcr_reproj_algn_masked.tif")
writeRaster(lcr_masked, lcr_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load binary raster
lcr_binary_file <- file.path(lcr_aligned_dir, "lcr_reproj_algn_binary.tif")
lcr_binary <- rast(lcr_binary_file)

crop_ext <- ext(sa)

# crop raster to study area
lcr_bin_sa <- crop(lcr_binary, crop_ext)
lcr_bin_sa
plot (lcr_bin_sa)

# mask to study area
lcr_bin_masked <- mask(lcr_bin_sa, sa)
lcr_bin_masked
plot(lcr_bin_masked)

# Save the raster as 8-bit unsigned integer
lcr_bin_sa_out_path <- file.path(lcr_masked_dir, "lcr_reproj_algn_binary_masked.tif")
writeRaster(lcr_bin_masked, lcr_bin_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
