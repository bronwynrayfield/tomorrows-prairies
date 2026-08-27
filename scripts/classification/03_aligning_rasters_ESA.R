
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
library(reticulate)
library(geojsonio)
library(stars)
library(googleCloudStorageR)

# Pin reticulate to your conda Python (myenv)
Sys.setenv(RETICULATE_PYTHON = "E:/Users/alison.long/.conda/envs/myenv/python.exe")

py_config()  # confirm it shows path set in previous step

# Ensure pip and GEE API are installed into 'myenv'
reticulate::py_install("pip", pip = TRUE)
reticulate::py_install(c("earthengine-api", "numpy", "requests", "pillow"), pip = TRUE)

# Verify the module is visible
reticulate::py_module_available("ee")

# THEN load 
library(rgee)

# set python environment
#ee_install_set_pyenv("e:/Users/alison.long/.conda/envs/myenv/python.exe")

# authenticate and initialize rgee
ee_Authenticate()
ee_Initialize("alison.long", project = "tomorrows-prairies", gcs = TRUE)

SaK_file <- "R:/CanadianGrasslands/WG_Classification/Classification_WG_intermediate_products/tomorrows-prairies-cac51d6c98f0.json"

# Assign the SaK to a EE user.
ee_utils_sak_copy(
  sakfile =  SaK_file,
  users = "alison.long" # Unlike GD, we can use the same SaK for multiple users.
)

# Validate your SaK
ee_utils_sak_validate(bucket = "tp-rgee-backup")

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
esa_root <- file.path(raw_grass_data_dir, "Land Cover/ESA WorldCover 2021")
esa_merged_dir <- file.path(proc_grass_data_dir, "esa_worldcover/merged")
esa_aligned_dir <- file.path(proc_grass_data_dir, "esa_worldcover/aligned")
esa_cropped_dir <- file.path(proc_grass_data_dir, "esa_worldcover/cropped")
esa_masked_dir <- file.path(proc_grass_data_dir, "esa_worldcover/masked")
esa_pilot_dir <- file.path(proc_grass_data_dir, "esa_worldcover/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# ESA rasters
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W072_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W075_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W078_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W081_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W084_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W087_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W090_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W093_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W096_Map.tif")
esa_file_2 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W099_Map.tif")
esa_file_3 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W102_Map.tif")
esa_file_4 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W105_Map.tif")
esa_file_5 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W108_Map.tif")
esa_file_6 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W111_Map.tif")
esa_file_7 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W114_Map.tif")
esa_file_8 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W117_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W120_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W123_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W126_Map.tif")
esa_file_1 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N48W129_Map.tif")
esa_file_9 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W081_Map")
esa_file_9 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W084_Map")
esa_file_9 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W087_Map")
esa_file_9 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W090_Map")
esa_file_9 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W093_Map")
esa_file_9 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W096_Map")
esa_file_9 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W099_Map")
esa_file_10 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W102_Map")
esa_file_11 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W105_Map")
esa_file_12 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W108_Map")
esa_file_13 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W111_Map")
esa_file_14 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W114_Map")
esa_file_15 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W117_Map")
esa_file_16 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W120_Map")
esa_file_16 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W123_Map")
esa_file_16 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N51W126_Map")
esa_file_17 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W102_Map")
esa_file_18 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W105_Map")
esa_file_19 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W108_Map")
esa_file_20 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W111_Map")
esa_file_21 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W114_Map")
esa_file_22 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W117_Map")
esa_file_23 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W120_Map")
esa_file_24 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N54W123_Map")
esa_file_25 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N57W117_Map")
esa_file_26 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N57W120_Map")
esa_file_27 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N57W123_Map")
esa_file_28 <- file.path(esa_root, "ESA_WorldCover_10m_2021_V200_N57W126_Map")
esa_file_29 <- file.path(esa_root, "")
esa_file_30 <- file.path(esa_root, "")
esa_file_31 <- file.path(esa_root, "")

# Load the collection
wc_ic <- ee$ImageCollection("ESA/WorldCover/v200")


##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

# Ensure AOI CRS is geographic (WGS84), since WorldCover tiles are in EPSG:4326.
# (If sa has a projected CRS, reproject its geometry—not its values—just for the AOI.)
if (!is.lonlat(crs(sa))) {
  message("AOI raster is not lon/lat. Will transform geometry to EPSG:4326 for EE.")
}

# Get extent (SpatExtent) and convert it to a SpatVector polygon
ext_sa   <- ext(sa)
aoi_rect <- as.polygons(ext_sa, crs = crs(sa))  # <- returns a SpatVector polygon

# Convert to sf and reproject to WGS84 (EPSG:4326) for Earth Engine
aoi_rect_sf  <- st_as_sf(aoi_rect)
aoi_rect_wgs <- st_transform(aoi_rect_sf, 4326)

# Make it a single feature in EPSG:4326
aoi_rect_wgs_single <- aoi_rect_wgs |>
  st_make_valid() |>
  st_union() |>
  st_as_sf() |>
  st_set_crs(4326)

plot(aoi_rect_wgs_single)

# Use the geometry column (sfc), not the sf data frame:
aoi_sfc      <- st_geometry(aoi_rect_wgs_single)  # this is an sfc with 1 geometry
aoi_rect_ee  <- sf_as_ee(aoi_sfc)                 # <-- returns ee$Geometry

# The land-cover class band is "Map" in v200 tiles.
wc_img <- wc_ic$
  select("Map")$
  mosaic()$                # ensure a single image
  clip(aoi_rect_ee)

# Bring into R as terra::SpatRaster - GD/GCS more robust for larger AOIs 
wc_rast <- ee_as_rast(
  image     = wc_img,
  region    = aoi_rect_ee,
  scale     = 10,
  maxPixels = 1e13,
  fileDimensions = 8192,
  via       = "gcs",             # or "gcs" for Google Cloud Storage
  container = "tp-rgee-backup"        # Drive folder or GCS bucket
)

wc_rast

# moved files from temp folder to R:\CanadianGrasslands\existingdatasets\Land Cover\ESA WorldCover 2021\from_rgee

# Folder where the files landed 
tiles_dir <- "R:/CanadianGrasslands/existingdatasets/Land Cover/ESA WorldCover 2021/from_rgee"

# List & sort the tiles (ensures -0001, -0002, ..., -0065 order)
files <- list.files(
  tiles_dir, pattern = "^noid_image-\\d{4}\\.tif$", full.names = TRUE
)
files <- files[order(as.integer(sub(".*-(\\d+)\\.tif$", "\\1", files)))]
length(files)   # should be 65
head(files)

# Read all tiles as a multi-source SpatRaster
RastCollection <- sprc(files)  # stacks all tiles as separate sources (still lazy)

# Mosaic rasters
# If tiles overlap, 'first' keeps the first tile's pixel, 'max' picks the larger code.
esa_merged <- mosaic(RastCollection, fun = "first")

# Save the raster as 8-bit unsigned integer
esa_merge_out_path <- file.path(esa_merged_dir, "esa_worldcover_v200_10m_merged.tif")
writeRaster(esa_merged, esa_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

plot(esa_merged)
esa_merged

#class       : SpatRaster 
#size        : 302068, 819138, 1  (nrow, ncol, nlyr)
#resolution  : 8.983153e-05, 8.983153e-05  (x, y)
#extent      : -139.0607, -65.47631, 39.25826, 66.39349  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source      : spat_77ac656c9e9_30636_o8up7yuDdJ6XHBX.tif 
#varname     : noid_image-0001 
#name        : Map 
#min value   :   0 
#max value   : 100

# Set pixels with value 0 to NA - ocean, bay, US
esa_merged[esa_merged == 0] <- NA
plot(esa_merged)
esa_merged

#class       : SpatRaster 
#size        : 302068, 819138, 1  (nrow, ncol, nlyr)
#resolution  : 8.983153e-05, 8.983153e-05  (x, y)
#extent      : -139.0607, -65.47631, 39.25826, 66.39349  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source      : spat_77ac20941a8_30636_80Q6udnsxVen2Mq.tif 
#varname     : noid_image-0001 
#name        : Map 
#min value   :  10 
#max value   : 100 

# reproject before resampling
esa_reproj <- project(esa_merged, sa, method = "near") # better approach

freq_table2 <- freq(esa_reproj) 
print(freq_table2) 

origin(esa_reproj) # 
esa_reproj

# align to reference raster 
esa_reproj_algn <- resample(esa_reproj, sa, method = "near")
plot(esa_reproj_algn)

# Get values and counts
freq_table3 <- freq(esa_reproj_algn) 
print(freq_table3) # 

origin(esa_reproj_algn) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
esa_out_path <- file.path(esa_aligned_dir, "esa_reproj_algn.tif")
writeRaster(esa_reproj_algn, esa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

####--- ALIGN RASTERS ---#### 

#----- ESA WorldCover ----#

# Load rasters
r_esa_1 <- rast(esa_file_1)

origin(r_esa_1) # 
origin(r_esa_2) # 
r_esa_1
r_esa_2

# Merge 
merged_esa <- merge(r_esa_1, r_esa_2, r_esa_3)
merged_esa
plot(merged_esa)

# print frequency table
freq_table <- freq(merged_esa) 
print(freq_table) 

# Mark as categorical (factor)
merged_esa <- as.factor(merged_esa)

# Save the raster as 8-bit unsigned integer
esa_merge_out_path <- file.path(esa_merged_dir, "esa_worldcover_merged.tif")
writeRaster(merged_esa, esa_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop before reprojecting to make manageable size 

# Define the bbox in lon/lat degrees (WGS84)
canada_ext_ll  <- ext(-140, -70, 30, 60)
canada_bbox_ll <- as.polygons(canada_ext_ll, crs = "EPSG:4326")

# Get CRS from the raster 
r_crs <- crs(esa_merged)

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
r_esa_crop <- crop(esa_merged, canada_ext_rcrs)
r_esa_mask <- mask(r_esa_crop, canada_bbox_rcrs)

r_esa_mask
plot(r_esa_mask)

# Save the raster as 8-bit unsigned integer
esa_crop_out_path <- file.path(esa_cropped_dir, "esa_crop_bbox.tif")
writeRaster(r_esa_mask, esa_crop_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Set pixels with value 0 to NA - ocean, bay, US
r_esa_mask[r_esa_mask == 0] <- NA
plot(r_esa_mask)

# reproject before resampling
esa_reproj <- project(r_esa_mask, sa, method = "near") # better approach

freq_table2 <- freq(esa_reproj) 
print(freq_table2) 

origin(esa_reproj) # 
esa_reproj

# better way:

esa_out_path <- file.path(esa_aligned_dir, "esa_reproj.tif")

# Conservative and robust GDAL options
wopt <- list(
  datatype = "INT1U",            # or "INT1U"/"INT2U" 
  gdal     = c(
    "COMPRESS=LZW",
    "PREDICTOR=2",
    "TILED=YES",
    "BIGTIFF=YES",
    "BLOCKXSIZE=512",
    "BLOCKYSIZE=512"
  )
)

# Reproject to 'sa' with nearest neighbor
esa_reproj <- project(esa_merged, sa, method = "near",
                      filename = esa_out_path, overwrite = TRUE, wopt = wopt)
plot(esa_reproj)

# next step:

esa_algn_out_path <- file.path(esa_aligned_dir, "esa_merged_reproj_algn.tif")

# align to reference raster 
esa_reproj_algn <- resample(esa_reproj, sa, method = "near",
                            filename = esa_algn_out_path, overwrite = TRUE, wopt = wopt)
plot(esa_reproj_algn)
esa_reproj_algn
origin(esa_reproj_algn) # 6.696121 -2.550334

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic 
#source      : esa_merged_reproj_algn.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#name        : Map 
#min value   :   0 
#max value   : 100 



# reclassify to binary (grass = 1)
esa_algn_binary <- classify(esa_reproj_algn, 
                            rcl = matrix(c(0, NA,
                                           10, 0, # Tree cover
                                           20, 0, # Shrubland
                                           30, 1, # Grassland
                                           40, 0, # Cropland
                                           50, 0, # Built-up
                                           60, 0, # Bare/sparse vegetation
                                           70, 0, # Snow and ice
                                           80, 0, # Permanent water bodies
                                           90, 0, # Herbaceous wetland
                                           95, 0, # Mangroves
                                           100, 0), # Moss and lichen
                                         ncol = 2, byrow = TRUE))
plot(esa_algn_binary)
esa_algn_binary

# Save the raster as 8-bit unsigned integer
esa_binary_out_path <- file.path(esa_aligned_dir, "esa_reproj_algn_binary_noshrub.tif")
writeRaster(esa_algn_binary, esa_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to binary with shrub
esa_algn_binary <- classify(esa_reproj_algn, 
                            rcl = matrix(c(0, NA,
                                           10, 0, # Tree cover
                                           20, 1, # Shrubland
                                           30, 1, # Grassland
                                           40, 0, # Cropland
                                           50, 0, # Built-up
                                           60, 0, # Bare/sparse vegetation
                                           70, 0, # Snow and ice
                                           80, 0, # Permanent water bodies
                                           90, 0, # Herbaceous wetland
                                           95, 0, # Mangroves
                                           100, 0), # Moss and lichen
                                         ncol = 2, byrow = TRUE))
plot(esa_algn_binary)
esa_algn_binary

# Save the raster as 8-bit unsigned integer
esa_binary_out_path <- file.path(esa_aligned_dir, "esa_reproj_algn_binary_noshrub.tif")
writeRaster(esa_algn_binary, esa_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# isolate shrub class
esa_algn_shrub <- classify(esa_reproj_algn, 
                            rcl = matrix(c(0, NA,
                                           10, 0, # Tree cover
                                           20, 1, # Shrubland
                                           30, 0, # Grassland
                                           40, 0, # Cropland
                                           50, 0, # Built-up
                                           60, 0, # Bare/sparse vegetation
                                           70, 0, # Snow and ice
                                           80, 0, # Permanent water bodies
                                           90, 0, # Herbaceous wetland
                                           95, 0, # Mangroves
                                           100, 0), # Moss and lichen
                                         ncol = 2, byrow = TRUE))
plot(esa_algn_shrub)
esa_algn_shrub

# Save the raster as 8-bit unsigned integer
esa_shrub_out_path <- file.path(esa_aligned_dir, "esa_reproj_algn_shrub.tif")
writeRaster(esa_algn_shrub, esa_shrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334

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

# Save the raster as 8-bit unsigned integer
esa_pa_out_path <- file.path(esa_pilot_dir, "esa_reproj_algn_crop_pa.tif")
writeRaster(esa_masked, esa_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
esa_bin_pa <- crop(esa_algn_binary, crop_ext)
esa_bin_pa
plot (esa_bin_pa)

# mask to pilot area
esa_bin_masked <- mask(esa_bin_pa, pa)
esa_bin_masked
plot(esa_bin_masked)

# Save the raster as 8-bit unsigned integer
esa_bin_pa_out_path <- file.path(esa_pilot_dir, "esa_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(esa_bin_masked, esa_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

################

######### crop and mask to study area

# load raster
esa_file <- file.path(esa_aligned_dir, "esa_merged_reproj_algn.tif")
esa <- rast(esa_file)

crop_ext <- ext(sa)

# crop raster to study area
esa_sa <- crop(esa, crop_ext)
esa_sa
plot (esa_sa)

# mask to study area
esa_masked <- mask(esa_sa, sa)
esa_masked
plot(esa_masked)

# Save the raster as 8-bit unsigned integer
esa_sa_out_path <- file.path(esa_masked_dir, "esa_merged_reproj_algn_masked.tif")
writeRaster(esa_masked, esa_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

######### crop and mask to study area

# load binary raster
esa_binary_file <- file.path(esa_aligned_dir, "esa_reproj_algn_binary_noshrub.tif")
esa_binary <- rast(esa_binary_file)

crop_ext <- ext(sa)

# crop raster to study area
esa_bin_sa <- crop(esa_binary, crop_ext)
esa_bin_sa
plot (esa_bin_sa)

# mask to study area
esa_bin_masked <- mask(esa_bin_sa, sa)
esa_bin_masked
plot(esa_bin_masked)

# Save the raster as 8-bit unsigned integer
esa_bin_sa_out_path <- file.path(esa_masked_dir, "esa_reproj_algn_binary_noshrub_masked.tif")
writeRaster(esa_bin_masked, esa_bin_sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
