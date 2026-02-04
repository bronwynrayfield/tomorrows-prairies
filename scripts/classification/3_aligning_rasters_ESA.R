
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

# esa directory
esa_root <- file.path(raw_grass_data_dir, "esa_worldcover")
esa_aligned_dir <- file.path(proc_grass_data_dir, "esa_worldcover/aligned")
esa_cropped_dir <- file.path(proc_grass_data_dir, "esa_worldcover/cropped")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# esa raster
esa_file <- file.path(esa_root, ".tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- ESA WorldCover ----#

# Load rasters
r_esa <- rast(esa_file)

origin(r_esa) # -2.378897e-11  1.010392e-11
r_esa


# mask global dataset to bounding box for Canada

canada_ext <- ext(
  -141,  # xmin (min longitude)
  -52,   # xmax (max longitude)
  40,    # ymin (min latitude)
  60     # ymax (max latitude)
)

# Make it a polygon (SpatVector) with CRS attached
canada_bbox <- as.polygons(canada_ext, crs = "EPSG:4326")

# polygon for Canada
canada_sf <- ne_countries(country = "Canada", returnclass = "sf")
canada_v  <- vect(canada_sf)

# plot
plot(canada_v, col = "grey90", border = "grey40", main = "Canada + bounding box")
plot(canada_bbox, add = TRUE, border = "red", lwd = 2, lty = 2)

# crop and mask to bounding box
r_esa_crop <- crop(r_esa, canada_ext)
r_esa_mask <- mask(r_esa_crop, canada_bbox)

r_esa_mask
plot(r_esa_mask)

# print frequency table
freq_table <- freq(r_esa_mask) 
print(freq_table) 

# Mark as categorical (factor)
r_esa_mask <- as.factor(r_esa_mask)

# reproject before resampling
esa_reproj <- project(r_esa_mask, crs(sa)) # or should be project(ngi_reclass, sa) ??

freq_table2 <- freq(esa_reproj) 
print(freq_table2) 

origin(esa_reproj) # 0.7826532 0.9658756
esa_reproj

# align to reference raster 
esa_reproj_algn <- resample(esa_reproj, sa, method = "near")
plot(esa_reproj_algn)

# Get values and counts
freq_table3 <- freq(esa_reproj_algn) 
print(freq_table3) # 0: 609557112; 1: 62790661; 2: 19747161; 3: 20456239

res(esa_reproj_algn) # 30 30
origin(esa_reproj_algn) # 6.696121 -2.550334
crs(esa_reproj_algn) # EPSG:9001
ext(esa_reproj_algn) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
esa_out_path <- file.path(esa_aligned_dir, "esa_reproj_algn.tif")
writeRaster(esa_reproj_algn, esa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
esa_reproj_algn_trim <- trim(esa_reproj_algn)

# visualize the trimmed raster
plot(esa_reproj_algn_trim)

# print raster details 
origin(esa_reproj_algn_trim) # 6.696121 -2.550334
res(esa_reproj_algn_trim) # 30 30 
crs(esa_reproj_algn_trim) # EPSG:9001
ext(esa_reproj_algn_trim) # SpatExtent : -1565033.30387889, -10793.3038788857, 1062567.4496659, 2329887.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
esa_trim_out_path <- file.path(esa_aligned_dir, "esa_reproj_algn_trim.tif")
writeRaster(esa_reproj_algn_trim, esa_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (grass = 1)
#esa_algn_binary <- classify(esa_reproj_algn, 
                            rcl = matrix(c(0, 0, # Other land cover
                                           1, 1, # Cultivated grassland
                                           2, 1, # Natural/semi-natural grassland
                                           3, 0), # Open shrubland
                                         ncol = 2, byrow = TRUE))
plot(esa_algn_binary)
esa_algn_binary

# Save the raster as 8-bit unsigned integer
esa_binary_out_path <- file.path(esa_aligned_dir, "esa_reclass_reproj_algn_binary.tif")
writeRaster(esa_algn_binary, esa_binary_out_path, datatype = "INT1U", 
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
esa_pa <- crop(esa_reproj_algn, crop_ext)
esa_pa
plot (esa_pa)


# mask to pilot area
esa_masked <- mask(esa_pa, pa)
esa_masked
plot(esa_masked)



# print raster details 
origin(esa_masked) # 6.696121 -2.550334
res(esa_masked) # 30 30
crs(esa_masked)

# Save the raster as 8-bit unsigned integer
esa_pa_out_path <- file.path(esa_cropped_dir, "esa_reproj_algn_crop_pa.tif")
writeRaster(esa_masked, esa_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)