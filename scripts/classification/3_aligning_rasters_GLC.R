
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
library(rgee)

# authenticate and initialize rgee
#ee_Authenticate()
ee_Initialize("alison.long", project = "tomorrows-prairies", gcs = TRUE)

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

# glc directory
glc_root <- file.path(raw_grass_data_dir, "Land Cover/GLC_FCS30D/GLC_FCS30D_19852022maps_W95-W100")
glc_merged_dir <- file.path(proc_grass_data_dir, "glc_fcs30d/merged")
glc_aligned_dir <- file.path(proc_grass_data_dir, "glc_fcs30d/aligned")
glc_cropped_dir <- file.path(proc_grass_data_dir, "glc_fcs30d/cropped")
glc_pilot_dir <- file.path(proc_grass_data_dir, "glc_fcs30d/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# glc raster
glc_file_1 <- file.path(glc_root, "GLC_FCS30D_20002022_W95N50_Annual.tif")
glc_file_2 <- file.path(glc_root, "GLC_FCS30D_20002022_W95N55_Annual.tif")
glc_file_3 <- file.path(glc_root, "GLC_FCS30D_20002022_W95N60_Annual.tif")
glc_file_4 <- file.path(glc_root, "GLC_FCS30D_20002022_W100N50_Annual.tif")
glc_file_5 <- file.path(glc_root, "GLC_FCS30D_20002022_W100N55_Annual.tif")
glc_file_6 <- file.path(glc_root, "GLC_FCS30D_20002022_W100N60_Annual.tif")
glc_file_7 <- file.path(glc_root, "GLC_FCS30D_20002022_W105N50_Annual.tif")
glc_file_8 <- file.path(glc_root, "GLC_FCS30D_20002022_W105N55_Annual.tif")
glc_file_9 <- file.path(glc_root, "GLC_FCS30D_20002022_W105N60_Annual.tif")
glc_file_10 <- file.path(glc_root, "GLC_FCS30D_20002022_W110N50_Annual.tif")
glc_file_11 <- file.path(glc_root, "GLC_FCS30D_20002022_W110N55_Annual.tif")
glc_file_12 <- file.path(glc_root, "GLC_FCS30D_20002022_W110N60_Annual.tif")
glc_file_13 <- file.path(glc_root, "GLC_FCS30D_20002022_W115N50_Annual.tif")
glc_file_14 <- file.path(glc_root, "GLC_FCS30D_20002022_W115N55_Annual.tif")
glc_file_15 <- file.path(glc_root, "GLC_FCS30D_20002022_W115N60_Annual.tif")

# Load the collection
glc_ic <- ee$ImageCollection("projects/sat-io/open-datasets/GLC-FCS30D/annual")

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
aoi_rect_ee  <- sf_as_ee(aoi_sfc)                 # returns ee$Geometry

# Each band represents a year - band 23 represents data for 2022
glc_img <- glc_ic$
  select("b23")$
  mosaic()$                # ensure a single image
  clip(aoi_rect_ee)

# Bring into R as terra::SpatRaster - GD/GCS more robust for larger AOIs 
glc_rast <- ee_as_rast(
  image     = glc_img,
  region    = aoi_rect_ee,
  scale     = 10,
  maxPixels = 1e13,
  #fileDimensions = 8192,
  via       = "gcs",             # or "gcs" for Google Cloud Storage
  container = "tp-rgee-backup"        # Drive folder or GCS bucket
)

glc_rast

# moved files from temp folder to R:\CanadianGrasslands\existingdatasets\Land Cover\GLC_FCS30D\from_rgee

# Folder where the files landed (adjust if you've moved them)
tiles_dir <- "R:/CanadianGrasslands/existingdatasets/Land Cover/GLC_FCS30D/from_rgee"

# List & sort the tiles (ensures -0001, -0002, ..., -0065 order)
files <- list.files(
  tiles_dir, pattern = "^noid_image-\\d{4}\\.tif$", full.names = TRUE
)
files <- files[order(as.integer(sub(".*-(\\d+)\\.tif$", "\\1", files)))]
length(files)   # should be 65
head(files)

# Read all tiles as a multi-source SpatRaster
RastCollection <- sprc(files)  # stacks all tiles as separate sources 

# Conservative and robust GDAL options
wopt <- list(
  datatype = "INT1U",           
  gdal     = c(
    "COMPRESS=LZW",
    "PREDICTOR=2",
    "TILED=YES",
    "BIGTIFF=YES",
    "BLOCKXSIZE=512",
    "BLOCKYSIZE=512"
  )
)

# Mosaic rasters
# If tiles overlap, 'first' keeps the first tile's pixel, 'max' picks the larger value
glc_merge_out_path <- file.path(glc_merged_dir, "glc_merged.tif")
glc_merged <- mosaic(RastCollection, fun = "first",
                     filename = glc_merge_out_path, overwrite = TRUE, wopt = wopt)

plot(glc_merged, main = "glc_merged")
glc_merged

# Set pixels with value 0 to NA - ocean, bay, US
#glc_merged[glc_merged == 0] <- NA
#plot(glc_merged)
#glc_merged

# Reproject to 'sa' with nearest neighbor
glc_reproj_out_path <- file.path(glc_aligned_dir, "glc_merged_reproj.tif")
glc_reproj <- project(glc_merged, sa, method = "near",
                      filename = glc_reproj_out_path, overwrite = TRUE, wopt = wopt)

freq_table2 <- freq(glc_reproj) 
print(freq_table2) 

origin(glc_reproj) # 
glc_reproj

# align to reference raster 
glc_algn_out_path <- file.path(glc_aligned_dir, "glc_merged_reproj_algn.tif")
glc_algn <- resample(glc_reproj, sa, method = "near",
                            filename = glc_algn_out_path, overwrite = TRUE, wopt = wopt)
plot(glc_algn, main = "glc_algn")
cats(glc_algn)
glc_algn

# Mark as categorical (factor)
glc_algn <- as.factor(glc_algn)

# Get values and counts
freq_table3 <- freq(glc_reproj_algn) 
print(freq_table3) # 

origin(glc_reproj_algn) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
#glc_out_path <- file.path(glc_aligned_dir, "glc_reproj_algn.tif")
#writeRaster(glc_reproj_algn, glc_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

####--- ALIGN RASTERS ---#### 

#----- GLC_FCS30D ----#


# Read only band 23 (for 2022 data)
#b23_1 <- rast(glc_file_1, lyrs = 23)

# Load rasters
#r_glc <- rast(glc_file)

#origin(r_glc) #
#r_glc

# print frequency table
#freq_table <- freq(r_glc) 
#print(freq_table) 

# Mark as categorical (factor)
#r_glc <- as.factor(r_glc)

# crop before reprojecting to make manageable size 

# Define the bbox in lon/lat degrees (WGS84)
#canada_ext_ll  <- ext(-140, -70, 30, 60)
#canada_bbox_ll <- as.polygons(canada_ext_ll, crs = "EPSG:4326")

# Get CRS from the raster 
#r_crs <- crs(r_glc)

# Project bbox polygon to the raster's CRS
#canada_bbox_rcrs <- project(canada_bbox_ll, r_crs)

# Meter/CRS-based extent derived from projected polygon
#canada_ext_rcrs <- ext(canada_bbox_rcrs)

# polygon for Canada
#canada_sf <- ne_countries(country = "Canada", returnclass = "sf")
#canada_v  <- vect(canada_sf)

# plot
#plot(canada_v, col = "grey90", border = "grey40", main = "Canada + bounding box")
#plot(canada_bbox, add = TRUE, border = "red", lwd = 2, lty = 2)

# crop and mask to bounding box
#r_glc_crop <- crop(r_glc, canada_ext_rcrs)
#r_glc_mask <- mask(r_glc_crop, canada_bbox_rcrs)

#r_glc_mask
#plot(r_glc_mask)

# Save the raster as 8-bit unsigned integer
#glc_crop_out_path <- file.path(glc_cropped_dir, "glc_crop_bbox.tif")
#writeRaster(r_glc_mask, glc_crop_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

# Set pixels with value 0 to NA - ocean, bay, US
#r_glc_mask[r_glc_mask == 0] <- NA
#plot(r_glc_mask)

# reproject before resampling
#glc_reproj <- project(r_glc_mask, sa, method = "near") # better approach

#freq_table2 <- freq(glc_reproj) 
#print(freq_table2) 

#origin(glc_reproj) # 
#glc_reproj

# align to reference raster 
#glc_reproj_algn <- resample(glc_reproj, sa, method = "near")
#plot(glc_reproj_algn)

# Get values and counts
#freq_table3 <- freq(glc_reproj_algn) 
#print(freq_table3) # 

#origin(glc_reproj_algn) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
#glc_out_path <- file.path(glc_aligned_dir, "glc_reproj_algn.tif")
#writeRaster(glc_reproj_algn, glc_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

# Trim NA-only rows/columns
#glc_reproj_algn_trim <- trim(glc_reproj_algn)

# visualize the trimmed raster
#plot(glc_reproj_algn_trim)

# print raster details 
#origin(glc_reproj_algn_trim) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
#glc_trim_out_path <- file.path(glc_aligned_dir, "glc_reproj_algn_trim.tif")
#writeRaster(glc_reproj_algn_trim, glc_trim_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

# reclassify to binary (grass = 1)
glc_algn_binary <- classify(glc_algn, 
                               rcl = matrix(c(0, NA,
                                              10, 0, # Rainfed cropland
                                              11, 0, # Herbaceous cover cropland
                                              12, 0, # Tree or shrub cover (orchard) cropland
                                              20, 0, # Irrigated cropland
                                              51, 0, # Open evergreen broadleaved forest
                                              52, 0, # Closed evergreen broadleaved forest
                                              61, 0, # Open deciduous broadleaved forest
                                              62, 0, # Closed deciduous broadleaved forest
                                              71, 0, # Open evergreen needle-leaved forest
                                              72, 0, # Closed evergreen needle-leaved forest
                                              81, 0, # Open deciduous needle-leaved forest
                                              82, 0, # Closed deciduous needle-leaved forest
                                              91, 0, # Open mixed leaf forest
                                              92, 0, # Closed mixed leaf forest
                                              120, 0, # Shrubland
                                              121, 0, # Evergreen shrubland
                                              122, 0, # Deciduous shrubland
                                              130, 1, # Grassland
                                              140, 0, # Lichens and mosses
                                              150, 0, # Sparse vegetation
                                              152, 0, # Sparse shrubland
                                              153, 0, # Sparse herbaceous
                                              181, 0, # Swamp
                                              182, 0, # Marsh
                                              183, 0, # Flooded flat
                                              184, 0, # Saline
                                              185, 0, # Mangrove
                                              186, 0, # Salt marsh
                                              187, 0, # Tidal flat
                                              190, 0, # Impervious surfaces
                                              200, 0, # Bare areas
                                              201, 0, # Consolidated bare areas
                                              202, 0, # Unconsolidated bare areas
                                              210, 0, # Water body
                                              220, 0), # Permanent ice and snow
                                            ncol = 2, byrow = TRUE))
plot(glc_algn_binary)
glc_algn_binary

# Save the raster as 8-bit unsigned integer
glc_binary_out_path <- file.path(glc_aligned_dir, "glc_reproj_algn_binary_noshrub.tif")
writeRaster(glc_algn_binary, glc_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate classes
glc_algn_shrub <- classify(glc_algn, 
                            rcl = matrix(c(0, NA,
                                           10, 0, # Rainfed cropland
                                           11, 0, # Herbaceous cover cropland
                                           12, 0, # Tree or shrub cover (orchard) cropland
                                           20, 0, # Irrigated cropland
                                           51, 0, # Open evergreen broadleaved forest
                                           52, 0, # Closed evergreen broadleaved forest
                                           61, 0, # Open deciduous broadleaved forest
                                           62, 0, # Closed deciduous broadleaved forest
                                           71, 0, # Open evergreen needle-leaved forest
                                           72, 0, # Closed evergreen needle-leaved forest
                                           81, 0, # Open deciduous needle-leaved forest
                                           82, 0, # Closed deciduous needle-leaved forest
                                           91, 0, # Open mixed leaf forest
                                           92, 0, # Closed mixed leaf forest
                                           120, 1, # Shrubland
                                           121, 0, # Evergreen shrubland
                                           122, 0, # Deciduous shrubland
                                           130, 0, # Grassland
                                           140, 0, # Lichens and mosses
                                           150, 0, # Sparse vegetation
                                           152, 0, # Sparse shrubland
                                           153, 0, # Sparse herbaceous
                                           181, 0, # Swamp
                                           182, 0, # Marsh
                                           183, 0, # Flooded flat
                                           184, 0, # Saline
                                           185, 0, # Mangrove
                                           186, 0, # Salt marsh
                                           187, 0, # Tidal flat
                                           190, 0, # Impervious surfaces
                                           200, 0, # Bare areas
                                           201, 0, # Consolidated bare areas
                                           202, 0, # Unconsolidated bare areas
                                           210, 0, # Water body
                                           220, 0), # Permanent ice and snow
                                         ncol = 2, byrow = TRUE))
plot(glc_algn_shrub)
glc_algn_shrub

# Save the raster as 8-bit unsigned integer
glc_shrub_out_path <- file.path(glc_aligned_dir, "glc_reproj_algn_shrub.tif")
writeRaster(glc_algn_shrub, glc_shrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# reclassify to isolate classes - no deciduous shrub in study area!
#glc_algn_decidshrub <- classify(glc_algn, 
                           rcl = matrix(c(0, NA,
                                          10, 0, # Rainfed cropland
                                          11, 0, # Herbaceous cover cropland
                                          12, 0, # Tree or shrub cover (orchard) cropland
                                          20, 0, # Irrigated cropland
                                          51, 0, # Open evergreen broadleaved forest
                                          52, 0, # Closed evergreen broadleaved forest
                                          61, 0, # Open deciduous broadleaved forest
                                          62, 0, # Closed deciduous broadleaved forest
                                          71, 0, # Open evergreen needle-leaved forest
                                          72, 0, # Closed evergreen needle-leaved forest
                                          81, 0, # Open deciduous needle-leaved forest
                                          82, 0, # Closed deciduous needle-leaved forest
                                          91, 0, # Open mixed leaf forest
                                          92, 0, # Closed mixed leaf forest
                                          120, 0, # Shrubland
                                          121, 0, # Evergreen shrubland
                                          122, 1, # Deciduous shrubland
                                          130, 0, # Grassland
                                          140, 0, # Lichens and mosses
                                          150, 0, # Sparse vegetation
                                          152, 0, # Sparse shrubland
                                          153, 0, # Sparse herbaceous
                                          181, 0, # Swamp
                                          182, 0, # Marsh
                                          183, 0, # Flooded flat
                                          184, 0, # Saline
                                          185, 0, # Mangrove
                                          186, 0, # Salt marsh
                                          187, 0, # Tidal flat
                                          190, 0, # Impervious surfaces
                                          200, 0, # Bare areas
                                          201, 0, # Consolidated bare areas
                                          202, 0, # Unconsolidated bare areas
                                          210, 0, # Water body
                                          220, 0), # Permanent ice and snow
                                        ncol = 2, byrow = TRUE))
plot(glc_algn_decidshrub)
glc_algn_decidshrub

# Save the raster as 8-bit unsigned integer
glc_decidshrub_out_path <- file.path(glc_aligned_dir, "glc_reproj_algn_deciduousshrub.tif")
writeRaster(glc_algn_decidshrub, glc_decidshrub_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Load study area raster
pa <- rast(pa_file)

# print raster details 
origin(pa) # 6.696121 -2.550334

# crop to pilot area
crop_ext <- ext(pa)
glc_pa <- crop(glc_algn, crop_ext)
glc_pa
plot (glc_pa)


# mask to pilot area
glc_masked <- mask(glc_pa, pa)
glc_masked
plot(glc_masked)

# print raster details 
origin(glc_masked) # 6.696121 -2.550334

# Save the raster as 8-bit unsigned integer
glc_pa_out_path <- file.path(glc_pilot_dir, "glc_reproj_algn_crop_pa.tif")
writeRaster(glc_masked, glc_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# crop binary raster to pilot area
glc_bin_pa <- crop(glc_algn_binary, crop_ext)
glc_bin_pa
plot (glc_bin_pa)

# mask to pilot area
glc_bin_masked <- mask(glc_bin_pa, pa)
glc_bin_masked
plot(glc_bin_masked)

# Save the raster as 8-bit unsigned integer
glc_bin_pa_out_path <- file.path(glc_pilot_dir, "glc_reproj_algn_crop_pa_binary_noshrub.tif")
writeRaster(glc_bin_masked, glc_bin_pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
