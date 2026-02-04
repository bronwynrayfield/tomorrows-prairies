
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
raw_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/class_inputs"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# gpw directory
gpw_root <- file.path(raw_grass_data_dir, "wri_gpw_ann")
gpw_aligned_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/aligned")
gpw_cropped_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/cropped")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# gpw raster
gpw_file <- file.path(gpw_root, "gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- WRI Global Pasture Watch ----#

# Load rasters
r_gpw <- rast(gpw_file)

origin(r_gpw) # -2.378897e-11  1.010392e-11
r_gpw

#class       : SpatRaster 
#size        : 528004, 1436004, 1  (nrow, ncol, nlyr)
#resolution  : 0.00025, 0.00025  (x, y)
#extent      : -179.0005, 180.0005, -56.0005, 76.0005  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source      : gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2.tif 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2

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
r_gpw_crop <- crop(r_gpw, canada_ext)
r_gpw_mask <- mask(r_gpw_crop, canada_bbox)

r_gpw_mask
plot(r_gpw_mask)

#class       : SpatRaster 
#size        : 1e+05, 356000, 1  (nrow, ncol, nlyr)
#resolution  : 0.00025, 0.00025  (x, y)
#extent      : -141, -52, 40, 65  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source      : spat_6ef6c173924e2_454508_80Q6udnsxVen2Mq.tif 
#varname     : gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2 
#min value   :                                                             0 
#max value   :                                                             3 

# print frequency table
freq_table <- freq(r_gpw_mask) 
print(freq_table) 

#   layer value      count
#1     1     0 9707550923
#2     1     1  391054390
#3     1     2 6741738895
#4     1     3  155468746

# Mark as categorical (factor)
r_gpw_mask <- as.factor(r_gpw_mask)
plot(r_gpw_mask)
r_gpw_mask
cats(r_gpw_mask)
cats(gpw_reproj)

# reproject before resampling
gpw_reproj <- project(r_gpw_mask, crs(sa)) # or should be project(ngi_reclass, sa) ??

#class       : SpatRaster 
#size        : 204165, 414290, 1  (nrow, ncol, nlyr)
#resolution  : 16.64012, 16.64012  (x, y)
#extent      : -3482758, 3411077, 0.7167462, 3397331  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_6ef6c5aea1b5c_454508_hn1JHs296zpoh7J.tif 
#categories  : gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2 
#min value   :                                                             0 
#max value   :                                                             3 

# load reprojected file from scratch folder
#gpw_reproj_file <- file.path("R:/CanadianGrasslands/WG_Classification/github/tomorrows-prairies/scratch/spat_6ef6c5aea1b5c_454508_hn1JHs296zpoh7J.tif")
#gpw_reproj2 <- rast(gpw_reproj_file)


#file.exists(gpw_reproj_file)
#file.info(gpw_reproj_file)
#terra::describe(gpw_reproj_file)


#f <- gpw_reproj_file
#con <- file(f, "rb")
#hdr <- readBin(con, what="raw", n=16)
#close(con)

#hdr
#rawToChar(hdr[1:4])


freq_table2 <- freq(gpw_reproj) 
print(freq_table2) 

origin(gpw_reproj) # 2.6022461 0.7167462
gpw_reproj
plot(gpw_reproj)


# How many non-NA cells (approx; uses blocks)
global(!is.na(gpw_reproj), "sum")
global(is.na(gpw_reproj), "sum")

# Quick min/max computed from data (not just header)
minmax(gpw_reproj)


# Try reading a tiny window from the top-left corner
try(values(gpw_reproj, row=1, nrows=10, col=1, ncols=10))



# align to reference raster 
gpw_reproj_algn <- resample(gpw_reproj, sa, method = "near")
plot(gpw_reproj_algn)

# Get values and counts
freq_table3 <- freq(gpw_reproj_algn) 
print(freq_table3) # 0: 609557112; 1: 62790661; 2: 19747161; 3: 20456239

res(gpw_reproj_algn) # 30 30
origin(gpw_reproj_algn) # 6.696121 -2.550334
crs(gpw_reproj_algn) # EPSG:9001
ext(gpw_reproj_algn) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
gpw_out_path <- file.path(gpw_aligned_dir, "gpw_reproj_algn.tif")
writeRaster(gpw_reproj_algn, gpw_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
gpw_reproj_algn_trim <- trim(gpw_reproj_algn)

# visualize the trimmed raster
plot(gpw_reproj_algn_trim)

# print raster details 
origin(gpw_reproj_algn_trim) # 6.696121 -2.550334
res(gpw_reproj_algn_trim) # 30 30 
crs(gpw_reproj_algn_trim) # EPSG:9001
ext(gpw_reproj_algn_trim) # SpatExtent : -1565033.30387889, -10793.3038788857, 1062567.4496659, 2329887.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
gpw_trim_out_path <- file.path(gpw_aligned_dir, "gpw_reproj_algn_trim.tif")
writeRaster(gpw_reproj_algn_trim, gpw_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary (grass = 1)
gpw_algn_binary <- classify(gpw_reproj_algn, 
                             rcl = matrix(c(0, 0, # Other land cover
                                            1, 1, # Cultivated grassland
                                            2, 1, # Natural/semi-natural grassland
                                            3, 0), # Open shrubland
                                          ncol = 2, byrow = TRUE))
plot(gpw_algn_binary)
gpw_algn_binary

# Save the raster as 8-bit unsigned integer
gpw_binary_out_path <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_binary_noshrub.tif")
writeRaster(gpw_algn_binary, gpw_binary_out_path, datatype = "INT1U", 
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
gpw_pa <- crop(gpw_reproj_algn, crop_ext)
gpw_pa
plot (gpw_pa)
 

# mask to pilot area
gpw_masked <- mask(gpw_pa, pa)
gpw_masked
plot(gpw_masked)



# print raster details 
origin(gpw_masked) # 6.696121 -2.550334
res(gpw_masked) # 30 30
crs(gpw_masked)

# Save the raster as 8-bit unsigned integer
gpw_pa_out_path <- file.path(gpw_cropped_dir, "gpw_reproj_algn_crop_pa.tif")
writeRaster(gpw_masked, gpw_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

beep()
