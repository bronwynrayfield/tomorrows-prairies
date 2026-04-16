
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
library(googleCloudStorageR)

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

# dw directory
#dw_root <- file.path(raw_grass_data_dir, "Grassland/DynamicWorld/DW_2025_grass")
dw_merged_dir <- file.path(proc_grass_data_dir, "dw_2025_grass/merged")
dw_aligned_dir <- file.path(proc_grass_data_dir, "dw_2025_grass/aligned")
dw_cropped_dir <- file.path(proc_grass_data_dir, "dw_2025_grass/cropped")
dw_masked_dir <- file.path(proc_grass_data_dir, "dw_2025_grass/masked")
dw_pilot_dir <- file.path(proc_grass_data_dir, "dw_2025_grass/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

# ---- CONFIG ----
bucket   <- "tomorrowsprairies_processed"           
prefix   <- "DynamicWorld_2025_grass"       
pattern  <- "\\.tif$"
local_dir <- file.path(tempdir(), "gcs_download")
dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)

# ---- AUTH ----
# Option 1: Service account key file (recommended for automations)
Sys.setenv(GOOGLE_APPLICATION_CREDENTIALS = "R:/CanadianGrasslands/WG_Classification/Classification_WG_intermediate_products/tomorrows-prairies-cac51d6c98f0.json")

# Option 2: Interactive login
# googleCloudStorageR::gcs_auth()  # will open a browser for OAuth

# If using a service account:
gcs_auth(json_file = Sys.getenv("GOOGLE_APPLICATION_CREDENTIALS", unset = NA))

# ---- LIST OBJECTS IN THE FOLDER ----
# Set the default bucket
gcs_global_bucket(bucket)

# List all objects under the prefix (folder-like)
objs <- gcs_list_objects(prefix = prefix, detail = "full")
if (is.null(objs) || nrow(objs) == 0) {
  stop("No objects found under the specified prefix. Check 'bucket' and 'prefix'.")
}

# Filter to .tif (case-insensitive)
objs_tif <- subset(objs, grepl(pattern, name, ignore.case = TRUE))
if (nrow(objs_tif) == 0) {
  stop("No .tif files found under the prefix.")
}

message(sprintf("Found %d TIFF(s). Downloading to: %s", nrow(objs_tif), local_dir))

# ---- DOWNLOAD ----
local_files <- character(nrow(objs_tif))
for (i in seq_len(nrow(objs_tif))) {
  obj_path <- objs_tif$name[i]                  # full object path within the bucket
  # Create local subdirectories mirroring GCS structure
  local_path <- file.path(local_dir, basename(obj_path))
  gcs_get_object(object_name = obj_path, saveToDisk = local_path, overwrite = TRUE)
  local_files[i] <- local_path
}

# ---- READ WITH TERRA ----
# Build a VRT (best for many tiles)
r_vrt <- terra::vrt(local_files)
print(r_vrt)
plot(r_vrt)

# Or could stack/SpatRaster of multiple layers (if same extent/resolution)
# r_stack <- terra::rast(local_files)

# reproject before resampling
dw_reproj <- project(r_vrt, sa, method = "near") 

freq_table2 <- freq(dw_reproj) 
print(freq_table2) 

origin(dw_reproj) # 6.696121 -2.550334
dw_reproj
#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008) 
#source      : spat_33cec61e41183_212204_zxzxx9jyyx6UnnI.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#name        : spat_33cec12b76cb9_212204_o8up7yuDdJ6XHBX 
#min value   :                                         0 
#max value   :                                         2

# align to reference raster 
dw_reproj_algn <- resample(dw_reproj, sa, method = "near")
plot(dw_reproj_algn)

crop_ext <- ext(sa)

# crop raster to study area
dw_algn_sa <- crop(dw_reproj_algn, crop_ext)
dw_algn_sa
plot (dw_algn_sa)

# mask to study area
dw_algn_masked <- mask(dw_algn_sa, sa)
dw_algn_masked
plot(dw_algn_masked)

# reclassify to binary (grass = 1)
dw_algn_masked <- classify(dw_algn_masked, 
                           rcl = matrix(c(0, 0,  # not grassland
                                          2, 1), # grassland
                                        ncol = 2, byrow = TRUE))
plot(dw_algn_masked)
dw_algn_masked

# Get values and counts
freq_table3 <- freq(dw_algn_masked) 
print(freq_table3) # 

origin(dw_algn_masked) # 

# Save the raster as 8-bit unsigned integer
dw_out_path <- file.path(dw_masked_dir, "dw_reproj_algn_masked_binary.tif")
writeRaster(dw_algn_masked, dw_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)


# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334

# crop to pilot area
crop_ext <- ext(pa)
dw_pa <- crop(dw_algn_masked, crop_ext)
dw_pa
plot (dw_pa)


# mask to pilot area
dw_masked <- mask(dw_pa, pa)
dw_masked
plot(dw_masked)

# print raster details 
origin(dw_masked) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
dw_pa_out_path <- file.path(dw_pilot_dir, "dw_reproj_algn_crop_pa.tif")
writeRaster(dw_masked, dw_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
dw_bin_pa <- crop(dw_algn_binary, crop_ext)
dw_bin_pa
plot (dw_bin_pa)

# mask to pilot area
dw_bin_masked <- mask(dw_bin_pa, pa)
dw_bin_masked
plot(dw_bin_masked)

# Save the raster as 8-bit unsigned integer
dw_bin_pa_out_path <- file.path(dw_pilot_dir, "dw_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(dw_bin_masked, dw_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)