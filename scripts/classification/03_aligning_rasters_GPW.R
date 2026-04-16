
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
raw_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/class_inputs"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# gpw directory
gpw_root <- file.path(raw_grass_data_dir, "wri_gpw_ann")
gpw_aligned_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/aligned")
gpw_cropped_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/cropped")
gpw_masked_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/masked")
gpw_pilot_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/pilot")

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
  -140,  # xmin (min longitude)
  -70,   # xmax (max longitude)
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
#size        : 80000, 240000, 1  (nrow, ncol, nlyr)
#resolution  : 0.00025, 0.00025  (x, y)
#extent      : -130, -70, 40, 60  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source      : spat_675c5c863b1_26460_80Q6udnsxVen2Mq.tif 
#varname     : gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2 
#min value   :                                                             0 
#max value   :                                                             3 

# print frequency table
freq_table <- freq(r_gpw_mask) 
print(freq_table) 

#  layer value      count
#1     1     0 8307959047
#2     1     1  391054021
#3     1     2 4647969785
#4     1     3  155468746

# Mark as categorical (factor)
r_gpw_mask2 <- as.factor(r_gpw_mask)

plot(r_gpw_mask)
r_gpw_mask

# reproject before resampling
gpw_reproj <- project(r_gpw_mask2, crs(sa)) # or should be project(ngi_reclass, sa) ??

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

# Try reading a tiny window from the top-left corner
#try(values(gpw_reproj, row=1, nrows=10, col=1, ncols=10))

freq_table2 <- freq(gpw_reproj) 
print(freq_table2) 

origin(gpw_reproj) # 2.789303 -4.316576
gpw_reproj
plot(gpw_reproj)

#class       : SpatRaster 
#size        : 147061, 261313, 1  (nrow, ncol, nlyr)
#resolution  : 18.13324, 18.13324  (x, y)
#extent      : -2674686, 2063765, -4.316576, 2666688  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_675c111b594b_26460_hn1JHs296zpoh7J.tif 
#categories  : gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2 
#min value   :                                                             0 
#max value   :                                                             3 

# align to reference raster 
gpw_reproj_algn <- resample(gpw_reproj, sa, method = "near")
plot(gpw_reproj_algn)

# Get values and counts
freq_table3 <- freq(gpw_reproj_algn) 
print(freq_table3) 

origin(gpw_reproj_algn) # 6.696121 -2.550334
gpw_reproj_algn

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_675c69f373fd_26460_lwU2fR7NDl3n2BP.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2 
#min value   :                                                             0 
#max value   :                                                             3 

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
gpw_reproj_algn_trim

#class       : SpatRaster 
#size        : 79870, 127419, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2236883, 1585687, 270327.4, 2666427  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_675c743f6362_26460_USRUenhvkoHhg0e.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2 
#min value   :                                                             0 
#max value   :                                                             3 

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

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : spat_675c63574560_26460_iAd7SUlxZLsaoD1.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#name        : gpw_grassland_rf.med.filt.bthr~40101_20241231_go_epsg.4326_v2 
#min value   :                                                             0 
#max value   :                                                             1

# Save the raster as 8-bit unsigned integer
gpw_binary_out_path <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_binary_noshrub.tif")
writeRaster(gpw_algn_binary, gpw_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load gpw raster
gpw_algn_file <- file.path(gpw_aligned_dir, "gpw_reproj_algn.tif")
gpw_algn <- rast(gpw_algn_file)
plot(gpw_algn)
gpw_algn

# reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native)
gpw_algn_nattam <- classify(gpw_algn, 
                            rcl = matrix(c(0, 0, # Other land cover
                                           1, 1, # Cultivated grassland
                                           2, 2, # Natural/semi-natural grassland
                                           3, 0), # Open shrubland
                                         ncol = 2, byrow = TRUE))
plot(gpw_algn_nattam)
gpw_algn_nattam

# Save the raster as 8-bit unsigned integer
gpw_nattam_out_path <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_natvtam_noshrub.tif")
writeRaster(gpw_algn_nattam, gpw_nattam_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate each class
gpw_algn_natseminat <- classify(gpw_algn, 
                            rcl = matrix(c(0, 0, # Other land cover
                                           1, 0, # Cultivated grassland
                                           2, 1, # Natural/semi-natural grassland
                                           3, 0), # Open shrubland
                                         ncol = 2, byrow = TRUE))
plot(gpw_algn_natseminat)
gpw_algn_natseminat

# Save the raster as 8-bit unsigned integer
gpw_natseminat_out_path <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_natseminat.tif")
writeRaster(gpw_algn_natseminat, gpw_natseminat_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate each class
gpw_algn_cult <- classify(gpw_algn, 
                                rcl = matrix(c(0, 0, # Other land cover
                                               1, 1, # Cultivated grassland
                                               2, 0, # Natural/semi-natural grassland
                                               3, 0), # Open shrubland
                                             ncol = 2, byrow = TRUE))
plot(gpw_algn_cult)
gpw_algn_cult

# Save the raster as 8-bit unsigned integer
gpw_cult_out_path <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_cultivated.tif")
writeRaster(gpw_algn_cult, gpw_cult_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate each class
gpw_algn_shrub <- classify(gpw_algn, 
                                rcl = matrix(c(0, 0, # Other land cover
                                               1, 0, # Cultivated grassland
                                               2, 0, # Natural/semi-natural grassland
                                               3, 1), # Open shrubland
                                             ncol = 2, byrow = TRUE))
plot(gpw_algn_shrub)
gpw_algn_shrub

# Save the raster as 8-bit unsigned integer
gpw_shrub_out_path <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_shrub.tif")
writeRaster(gpw_algn_shrub, gpw_shrub_out_path, datatype = "INT1U", 
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
gpw_pa_out_path <- file.path(gpw_pilot_dir, "gpw_reproj_algn_crop_pa.tif")
writeRaster(gpw_masked, gpw_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
gpw_bin_pa <- crop(gpw_algn_binary, crop_ext)
gpw_bin_pa
plot (gpw_bin_pa)

# mask to pilot area
gpw_bin_masked <- mask(gpw_bin_pa, pa)
gpw_bin_masked
plot(gpw_bin_masked)

# Save the raster as 8-bit unsigned integer
gpw_bin_pa_out_path <- file.path(gpw_pilot_dir, "gpw_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(gpw_bin_masked, gpw_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop native vs. tame raster to pilot area
gpw_nt_pa <- crop(gpw_algn_nattam, crop_ext)
gpw_nt_pa
plot (gpw_nt_pa)

# mask to pilot area
gpw_nt_masked <- mask(gpw_nt_pa, pa)
gpw_nt_masked
plot(gpw_nt_masked)

# Save the raster as 8-bit unsigned integer
gpw_nt_pa_out_path <- file.path(gpw_pilot_dir, "gpw_reproj_algn_crop_pa_natvtam_noshrub.tif")
writeRaster(gpw_nt_masked, gpw_nt_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################

######### crop and mask to study area

# load raster
gpw_file <- file.path(gpw_aligned_dir, "gpw_reproj_algn.tif")
gpw <- rast(gpw_file)

crop_ext <- ext(sa)

# crop raster to study area
gpw_sa <- crop(gpw, crop_ext)
gpw_sa
plot (gpw_sa)

# mask to study area
gpw_masked <- mask(gpw_sa, sa)
gpw_masked
plot(gpw_masked)

# Save the raster as 8-bit unsigned integer
gpw_sa_out_path <- file.path(gpw_masked_dir, "gpw_reproj_algn_masked.tif")
writeRaster(gpw_masked, gpw_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load binary raster
gpw_binary_file <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_binary_noshrub.tif")
gpw_binary <- rast(gpw_binary_file)

crop_ext <- ext(sa)

# crop raster to study area
gpw_bin_sa <- crop(gpw_binary, crop_ext)
gpw_bin_sa
plot (gpw_bin_sa)

# mask to study area
gpw_bin_masked <- mask(gpw_bin_sa, sa)
gpw_bin_masked
plot(gpw_bin_masked)

# Save the raster as 8-bit unsigned integer
gpw_bin_sa_out_path <- file.path(gpw_masked_dir, "gpw_reproj_algn_binary_noshrub_masked.tif")
writeRaster(gpw_bin_masked, gpw_bin_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load native vs tame raster
gpw_nattam_file <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_natvtam_noshrub.tif")
gpw_nattam <- rast(gpw_nattam_file)

crop_ext <- ext(sa)

# crop raster to study area
gpw_nt_sa <- crop(gpw_nattam, crop_ext)
gpw_nt_sa
plot (gpw_nt_sa)

# mask to study area
gpw_nt_masked <- mask(gpw_nt_sa, sa)
gpw_nt_masked
plot(gpw_nt_masked)

# Save the raster as 8-bit unsigned integer
gpw_nt_sa_out_path <- file.path(gpw_masked_dir, "gpw_reproj_algn_natvtam_noshrub_masked.tif")
writeRaster(gpw_nt_masked, gpw_nt_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load native raster
gpw_nat_file <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_natseminat.tif")
gpw_nat <- rast(gpw_nat_file)

crop_ext <- ext(sa)

# crop raster to study area
gpw_nat_sa <- crop(gpw_nat, crop_ext)
gpw_nat_sa
plot (gpw_nat_sa)

# mask to study area
gpw_nat_masked <- mask(gpw_nat_sa, sa)
gpw_nat_masked
plot(gpw_nat_masked)

# Save the raster as 8-bit unsigned integer
gpw_nat_sa_out_path <- file.path(gpw_masked_dir, "gpw_reproj_algn_natseminat_noshrub_masked.tif")
writeRaster(gpw_nat_masked, gpw_nat_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load native raster
gpw_cult_file <- file.path(gpw_aligned_dir, "gpw_reclass_reproj_algn_cultivated.tif")
gpw_cult <- rast(gpw_cult_file)

crop_ext <- ext(sa)

# crop raster to study area
gpw_cult_sa <- crop(gpw_cult, crop_ext)
gpw_cult_sa
plot (gpw_cult_sa)

# mask to study area
gpw_cult_masked <- mask(gpw_cult_sa, sa)
gpw_cult_masked
plot(gpw_cult_masked)

# Save the raster as 8-bit unsigned integer
gpw_cult_sa_out_path <- file.path(gpw_masked_dir, "gpw_reproj_algn_cultivated_noshrub_masked.tif")
writeRaster(gpw_cult_masked, gpw_cult_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
