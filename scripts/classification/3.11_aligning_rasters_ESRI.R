
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

# esri directory
esri_merged_dir <- file.path(proc_grass_data_dir, "esri_rangelands/merged")
esri_aligned_dir <- file.path(proc_grass_data_dir, "esri_rangelands/aligned")
esri_masked_dir <- file.path(proc_grass_data_dir, "esri_rangelands/masked")
esri_pilot_dir <- file.path(proc_grass_data_dir, "esri_rangelands/pilot")

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
prefix   <- "ESRI_2024_rangeland"       
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
esri_vrt <- terra::vrt(local_files)
print(esri_vrt)
plot(esri_vrt)

#class       : SpatRaster 
#size        : 204060, 720377, 1  (nrow, ncol, nlyr)
#resolution  : 8.983153e-05, 8.983153e-05  (x, y)
#extent      : -139.0603, -74.34772, 41.66907, 60.00009  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source      : spat_33cec7dcf353a_212204_lwU2fR7NDl3n2BP.vrt 
#name        : spat_33cec7dcf353a_212204_lwU2fR7NDl3n2BP

# Or could stack/SpatRaster of multiple layers (if same extent/resolution)
# r_stack <- terra::rast(local_files)

# reproject before resampling
esri_reproj <- project(esri_vrt, sa, method = "near") 

freq_table2 <- freq(esri_reproj) 
print(freq_table2) 

origin(esri_reproj) # 
esri_reproj

# align to reference raster 
esri_reproj_algn <- resample(esri_reproj, sa, method = "near")
plot(esri_reproj_algn)

crop_ext <- ext(sa)

# crop raster to study area
esri_algn_sa <- crop(esri_reproj_algn, crop_ext)
esri_algn_sa
plot (esri_algn_sa)

# mask to study area
esri_algn_masked <- mask(esri_algn_sa, sa)
esri_algn_masked
plot(esri_algn_masked)

# Get values and counts
freq_table3 <- freq(esri_algn_masked) 
print(freq_table3) # 

origin(esri_algn_masked) # 

# Save the raster as 8-bit unsigned integer
esri_out_path <- file.path(esri_masked_dir, "esri_reproj_algn_masked_binary.tif")
writeRaster(esri_algn_masked, esri_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)


# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 

# crop to pilot area
crop_ext <- ext(pa)
esri_pa <- crop(esri_algn_masked, crop_ext)
esri_pa
plot (esri_pa)

# mask to pilot area
esri_masked <- mask(esri_pa, pa)
esri_masked
plot(esri_masked)

# print raster details 
origin(esri_masked) #

# Save the raster as 8-bit unsigned integer
esri_pa_out_path <- file.path(esri_pilot_dir, "esri_reproj_algn_crop_pa_binary.tif")
writeRaster(esri_masked, esri_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

