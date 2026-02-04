
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

# lcr directory
lcr_root <- file.path(raw_grass_data_dir, "./Land Cover/StatCan Land Cover Register 2020")
lcr_aligned_dir <- file.path(proc_grass_data_dir, "statcan_lcr/aligned")
lcr_cropped_dir <- file.path(proc_grass_data_dir, "statcan_lcr/cropped")

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

# reproject before resampling
lcr_reproj <- project(r_lcr, crs(sa)) # or should be project(ngi_reclass, sa) ??

freq_table2 <- freq(lcr_reproj) 
print(freq_table2) 

origin(lcr_reproj) # 0.7826532 0.9658756
lcr_reproj

# align to reference raster 
lcr_reproj_algn <- resample(lcr_reproj, sa, method = "near")
plot(lcr_reproj_algn)

# Get values and counts
freq_table3 <- freq(lcr_reproj_algn) 
print(freq_table3) # 0: 609557112; 1: 62790661; 2: 19747161; 3: 20456239

res(lcr_reproj_algn) # 30 30
origin(lcr_reproj_algn) # 6.696121 -2.550334
crs(lcr_reproj_algn) # EPSG:9001
ext(lcr_reproj_algn) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

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
#lcr_algn_binary <- classify(lcr_reproj_algn, 
                            rcl = matrix(c(0, 0, # Other land cover
                                           1, 1, # Cultivated grassland
                                           2, 1, # Natural/semi-natural grassland
                                           3, 0), # Open shrubland
                                         ncol = 2, byrow = TRUE))
plot(lcr_algn_binary)
lcr_algn_binary

# Save the raster as 8-bit unsigned integer
lcr_binary_out_path <- file.path(lcr_aligned_dir, "lcr_reclass_reproj_algn_binary.tif")
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
