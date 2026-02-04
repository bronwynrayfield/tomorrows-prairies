
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
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/github/tomorrows-prairies/scratch")

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

# ECCC directory
eccc_root <- file.path(raw_grass_data_dir, "eccc_grass")
eccc_cropped_dir <- file.path(proc_grass_data_dir, "eccc_grass/cropped")

# PHJV PUDL directory
pudl_root <- file.path(raw_grass_data_dir, "phjv_pudl")
pudl_cropped_dir <- file.path(proc_grass_data_dir, "phjv_pudl/cropped")

# Central Grasslands Roadmap directory
cgr_root <- file.path(raw_grass_data_dir, "cgr_gam_adv")
cgr_cropped_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/cropped")

# WRI annual grass directory
wri_root <- file.path(raw_grass_data_dir, "wri_gpw_ann")
wri_cropped_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/cropped")

# AAFC LUTS directory
luts_tile_root <- file.path(raw_grass_data_dir, "aafc_luts_2020_30m_v5")
luts_merged_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/merged")
luts_cropped_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/cropped")
luts_aligned_dir <- file.path(proc_grass_data_dir, "aafc_luts_2020_30m_v5/aligned")


####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# NGI rasters
ngi_ab_file <- file.path(ngi_tile_root, "AB-002.tif")
ngi_mb_file <- file.path(ngi_tile_root, "MB-003.tif")
ngi_sk_file <- file.path(ngi_tile_root, "SK-001.tif")

# AAFC LUTS 
luts_files <- file.path(luts_tile_root)

# ECCC grass summary
eccc_file <- file.path(eccc_root, "grass_summary_a_class.tif")

# WRI
wri_file <- file.path(wri_root, "gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2.tif")

# PUDL
pudl_file <- file.path(pudl_root, "PHJV_PUDLmask.tif")

# Central Grasslands Roadmap - Grassland Assessment Map (Advanced)
cgr_file <- file.path(cgr_root, "CGR_GAM_V2_separated.tif")

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

# Set pixels with value 0 to NA
r_ab[r_ab == 0] <- NA
r_mb[r_mb == 0] <- NA
r_sk[r_sk == 0] <- NA

# Merge 
merged_ngi <- merge(r_ab, r_mb, r_sk)

res(merged_ngi) # 10 10
origin(merged_ngi) # 0 0 
crs(merged_ngi)
ext(merged_ngi) # SpatExtent : 4622750, 5860000, 1443910, 2274690 (xmin, xmax, ymin, ymax)

# plot the result
plot(merged_ngi, main = "Merged NGI Raster")

# print frequency table
freq_table <- freq(merged_ngi) 
print(freq_table) # 2: 894806781; 3: 864399198; 5: 214191378

# Save the raster as 8-bit unsigned integer
ngi_merge_out_path <- file.path(ngi_merged_dir, "cfga_ngi_merged_no0s.tif")
writeRaster(merged_ngi, ngi_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# need to reclassify to 1, 2, 3 instead of 2, 3, 5 so new class (4) doesn't get introduced when resampling
ngi_reclass <- classify(merged_ngi, 
                            rcl = matrix(c(2, 1,
                                           3, 2,
                                           5, 3), 
                                         ncol = 2, byrow = TRUE))

freq_table3 <- freq(ngi_reclass) 
print(freq_table3) # 1: 894806781; 2: 864399198; 3: 214191378

# Mark as categorical (factor)
ngi_reclass <- as.factor(ngi_reclass)

# Add category labels
levels(ngi_reclass) <- data.frame(
  value = c(1, 2, 3),
  label = c("Tame", "Native", "Mixed")
)

freq_table4 <- freq(ngi_reclass) 
print(freq_table4) # 1: 894806781; 2: 864399198; 3: 214191378

res(ngi_reclass) # 10 10
origin(ngi_reclass) # 0 0 
crs(ngi_reclass) # EPSG:3347
ext(ngi_reclass) # SpatExtent : 4622750, 5860000, 1443910, 2274690 (xmin, xmax, ymin, ymax)


# reproject before resampling
ngi_reclass_merged_reproj <- project(ngi_reclass, crs(sa)) # or should be project(ngi_reclass, sa) ??

freq_table5 <- freq(ngi_reclass_merged_reproj) 
print(freq_table5) # 1: 954680046; 2: 917799068; 3: 227772560

res(ngi_reclass_merged_reproj) # 9.754172 9.754172
origin(ngi_reclass_merged_reproj) # 1.715237 -3.612827
crs(ngi_reclass_merged_reproj) # EPSG:9001
ext(ngi_reclass_merged_reproj) # SpatExtent : -1307057.29414187, -35942.4075207314, 908773.057216913, 1941954.4332725 (xmin, xmax, ymin, ymax)

# align to reference raster (OPTION 1: MODE)
#ngi_reclass_merged_reproj_algn_mode <- resample(ngi_reclass_merged_reproj, sa, method = "mode") # mode = good choice for categorical rasters
                                                                            # if output cells overlap with multiple input cells
#plot(ngi_reclass_merged_reproj_algn_mode)

# Get values and counts
#freq_table6 <- freq(ngi_reclass_merged_reproj_algn_mode) 
#print(freq_table6) # 1: 144602453; 2: 113848371; 3: 32231685

#res(ngi_reclass_merged_reproj_algn_mode) # 30 30
#origin(ngi_reclass_merged_reproj_algn_mode) # 6.696121 -2.550334
#crs(ngi_reclass_merged_reproj_algn_mode) # EPSG:9001
#ext(ngi_reclass_merged_reproj_algn_mode) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
#ngi_mode_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_reclass_merged_reproj_algn_mode.tif")
#writeRaster(ngi_reclass_merged_reproj_algn_mode, ngi_mode_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

# align to reference raster (OPTION 2: NEAREST NEIGHBOR)
ngi_reclass_merged_reproj_algn_near <- resample(ngi_reclass_merged_reproj, sa, method = "near")
plot(ngi_reclass_merged_reproj_algn_near)

# Get values and counts
freq_table7 <- freq(ngi_reclass_merged_reproj_algn_near) 
print(freq_table7) # 1: 100921769; 1: 97023671; 3: 24084460

res(ngi_reclass_merged_reproj_algn_near) # 30 30
origin(ngi_reclass_merged_reproj_algn_near) # 6.696121 -2.550334
crs(ngi_reclass_merged_reproj_algn_near) # EPSG:9001
ext(ngi_reclass_merged_reproj_algn_near) # SpatExtent : -2321003.30387889, 1585686.69612111, 270327.4496659, 2866497.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
ngi_near_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_reclass_merged_reproj_algn_near.tif")
writeRaster(ngi_reclass_merged_reproj_algn_near, ngi_near_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Trim NA-only rows/columns
ngi_reclass_merged_reproj_algn_near_trim <- trim(ngi_reclass_merged_reproj_algn_near)
#ngi_reclass_merged_reproj_algn_mode_trim <- trim(ngi_reclass_merged_reproj_algn_mode)

# visualize the trimmed raster
plot(ngi_reclass_merged_reproj_algn_near_trim)
#plot(ngi_reclass_merged_reproj_algn_mode_trim)

# print raster details 
origin(ngi_reclass_merged_reproj_algn_near_trim) # 6.696121 -2.550334
res(ngi_reclass_merged_reproj_algn_near_trim) # 30 30 
crs(ngi_reclass_merged_reproj_algn_near_trim) # EPSG:9001
ext(ngi_reclass_merged_reproj_algn_near_trim) # SpatExtent : -1242503.30387889, -41333.3038788857, 1062927.4496659, 1798587.4496659 (xmin, xmax, ymin, ymax)

#origin(ngi_reclass_merged_reproj_algn_mode_trim) # 6.696121 -2.550334
#res(ngi_reclass_merged_reproj_algn_mode_trim) # 30 30
#crs(ngi_reclass_merged_reproj_algn_mode_trim) # # EPSG:9001
#ext(ngi_reclass_merged_reproj_algn_mode_trim) # SpatExtent : -1242533.30387889, -41333.3038788857, 1062897.4496659, 1798587.4496659 (xmin, xmax, ymin, ymax)

# Save the raster as 8-bit unsigned integer
ngi_near_trim_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_reclass_merged_reproj_algn_near_trim.tif")
writeRaster(ngi_reclass_merged_reproj_algn_near_trim, ngi_near_trim_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Save the raster as 8-bit unsigned integer
#ngi_mode_trim_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_reclass_merged_reproj_algn_mode_trim.tif")
#writeRaster(ngi_reclass_merged_reproj_algn_mode_trim, ngi_mode_trim_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

# reclassify to binary (values 2, 3, 5 become 1)
ngi_algn_binary <- classify(ngi_reclass_merged_reproj_algn_near, 
                               rcl = matrix(c(2, 1,
                                 3, 1,
                                 5, 1), 
                               ncol = 2, byrow = TRUE))

# Save the raster as 8-bit unsigned integer
ngi_binary_out_path <- file.path(ngi_aligned_dir, "cfga_ngi_reclass_merged_reproj_algn_binary.tif")
writeRaster(ngi_algn_binary, ngi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

## NEED TO RUN
# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa)
res(pa)
crs(pa)

# crop to pilot area
crop_ext <- ext(pa)
luts_pa <- crop(luts_merged_algn, crop_ext)
luts_pa
plot (luts_pa)

# mask to pilot area
luts_masked <- mask(luts_pa, pa)
luts_masked
plot(luts_masked)

# print raster details 
origin(luts_masked)
res(luts_masked)
crs(luts_masked)

# Save the raster as 8-bit unsigned integer
luts_pa_out_path <- file.path(luts_cropped_dir, "aafc_luts_merged_reproj_algn_crop_pa.tif")
writeRaster(luts_masked, luts_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)







