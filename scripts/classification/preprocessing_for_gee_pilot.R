
###################################################################################
##         PREPROCESSING OF GRASSLAND PRODUCTS FOR GEE - PILOT AREA              ##
###################################################################################

##----------------------- Alison Long - October 2025 -------------------------

## Pixel alignment in GEE is approximate - abstracts away low-level raster alignment details.
## Use terra to align and preprocess rasters locally, 
## then upload them to GEE for large-scale analysis and visualization.

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "./processeddatasets/pilot_area"

# study area
sa_dir <- "./processeddatasets/Studyarea"

# grasslands
raw_grass_data_dir <- "./existingdatasets/class_inputs"
proc_grass_data_dir <- "./processeddatasets/grass"
out_grass_dir <- "./outputs/grass"

# NGI tiles directory
ngi_tile_root <- file.path(raw_grass_data_dir, "cfga_ngi")
ngi_merged_dir <- file.path(proc_grass_data_dir, "cfga_ngi/merged")
ngi_cropped_dir <- file.path(proc_grass_data_dir, "cfga_ngi/cropped")
ngi_aligned_dir <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")

# ECCC directory
eccc_root <- file.path(raw_grass_data_dir, "eccc_grass_summary")
eccc_cropped_dir <- file.path(proc_grass_data_dir, "eccc_grass_summary/cropped")

# PHJV PUDL directory
pudl_root <- file.path(raw_grass_data_dir, "phjv_pudl")
pudl_cropped_dir <- file.path(proc_grass_data_dir, "phjv_pudl/cropped")

# Central Grasslands Roadmap directory
cgr_root <- file.path(raw_grass_data_dir, "cgr_gam_adv")
cgr_cropped_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/cropped")

# WRI annual grass directory
wri_root <- file.path(raw_grass_data_dir, "wri_annual_grass")
wri_cropped_dir <- file.path(proc_grass_data_dir, "wri_annual_grass/cropped")


####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "GNP_pilot_area_crop_CSD.shp")

# Study Area
sa_file <- file.path(sa_dir, "classification_study_area.shp")

# CNVC raster
cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# NGI rasters
ngi_ab_file <- file.path(ngi_tile_root, "AB-002.tif")
ngi_mb_file <- file.path(ngi_tile_root, "MB-003.tif")
ngi_sk_file <- file.path(ngi_tile_root, "SK-001.tif")

# AAFC LUTS 

# ECCC grass summary
eccc_file <- file.path(eccc_root, "grass_summary_a_class.tif")

# WRI
wri_file <- file.path(wri_root, "gpw_grassland_rf.savgol.bthr_c_30m_20220101_20221231_go_epsg.4326_v1.tif")

# PUDL
pudl_file <- file.path(pudl_root, "PHJV_PUDLmask.tif")

# Central Grasslands Roadmap - Grassland Assessment Map (Advanced)
cgr_file <- file.path(cgr_root, "CGR_GAM_V2_separated.tif")

##############################################

# prep rasters for alignment

#----- merge NGI rasters ----

# Load rasters
r_ab <- rast(ngi_ab_file)
r_mb <- rast(ngi_mb_file)
r_sk <- rast(ngi_sk_file)

res(r_ab)
res(r_mb)
res(r_sk)
origin(r_ab)
origin(r_mb)
origin(r_sk)

ext(r_ab)
ext(r_mb)
ext(r_sk)
nrow(r_ab); ncol(r_ab)
nrow(r_mb); ncol(r_mb)
nrow(r_sk); ncol(r_sk)

# Set pixels with value 0 to NA
r_ab[r_ab == 0] <- NA
r_mb[r_mb == 0] <- NA
r_sk[r_sk == 0] <- NA

# Merge using mosaic with max function
merged_ngi <- mosaic(r_ab, r_mb, r_sk, fun = max)

res(merged_ngi)
origin(merged_ngi)
                        
# Optional: plot the result
plot(merged_ngi, main = "Merged NGI Raster")

# Save the raster as 8-bit unsigned integer
ngi_merge_out_path <- file.path(ngi_merged_dir, "cfga_ngi_merged_no0s.tif")
writeRaster(merged_ngi, ngi_merge_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

##############################################

# use CNVC file as snap raster

# Load CNVC raster
cnvc <- rast(cnvc_file)

# print raster details 
origin(cnvc)
res(cnvc)
crs(cnvc)

# visualize CNVC raster
plot(cnvc)

# Load pilot area shapefile
pa_vect <- vect(pa_file)

# Load classification study area shapefile
sa_vect <- vect(sa_file)

# Reproject pilot area to match CNVC raster CRS
pa_vect_proj <- project(pa_vect, crs(cnvc))

# Reproject classification study area to match CNVC raster CRS
sa_vect_proj <- project(sa_vect, crs(cnvc))

# Rasterize pilot area using CNVC raster as template
pa_rast <- rasterize(pa_vect_proj, cnvc, field=1, touches=TRUE) # field=1 assigns constant value to rasterized cells
# touches=TRUE ensures that any cell touched by the polygon is included

# Rasterize classification study area using CNVC raster as template
sa_rast <- rasterize(sa_vect_proj, cnvc, field=1, touches=TRUE) # field=1 assigns constant value to rasterized cells
# touches=TRUE ensures that any cell touched by the polygon is included

# print raster details 
origin(pa_rast)
res(pa_rast)
crs(pa_rast)

# print raster details 
origin(sa_rast)
res(sa_rast)
crs(sa_rast)

# Get frequency table of raster values
freq_table <- freq(pa_rast)

# Print the table - 1888287 pixels with value 1
print(freq_table)

# visualize the new raster
plot(pa_rast)

# Trim NA-only rows/columns
pa_rast_trimmed <- trim(pa_rast)
plot(pa_rast_trimmed)

ext(pa_rast)
ext(pa_rast_trimmed)

nrow(pa_rast); ncol(pa_rast)
nrow(pa_rast_trimmed); ncol(pa_rast_trimmed)

# Save the raster as 8-bit unsigned integer
pa_out_path <- file.path(pa_dir, "pilot_area_test_aligned_cnvc_trim_crop_CSD.tif")
writeRaster(pa_rast_trimmed, pa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Convert raster to polygons (excluding NA values)
pa_polygons <- as.polygons(pa_rast_trimmed, na.rm = TRUE, dissolve = TRUE)

# Export shapefile
pa_out_shp <- file.path(pa_dir, "pilot_area_test_aligned_cnvc_crop_CSD.shp")
writeVector(pa_polygons, pa_out_shp, overwrite = TRUE)


# rasterize study area shapefile

#------------------ CROP RASTERS TO PILOT AREA -------------------------

# designate file paths
merged_ngi_path <- file.path(ngi_merged_dir, "cfga_ngi_merged_no0s.tif")
pa_rast_path <- file.path(pa_dir, "pilot_area_test_aligned_cnvc.tif")

# read rasters
merged_ngi <- rast(merged_ngi_path)
pa_rast <- rast(pa_rast_path)

# reproject raster to desired crs
merged_ngi <- project(merged_ngi, crs(pa_rast))

# crop raster to pilot area 
cropped_ngi <- crop(merged_ngi, ext(pa_rast))

# align cropped raster to pilot area (resolution, origin, CRS)
crp_algn_ngi <- resample(cropped_ngi, pa_rast, method = "near") # nearest neighbor for categorical

# check alignment - returns TRUE if rasters are aligned in terms of CRS, resolution, origin, and extent
compareGeom(pa_rast, crp_algn_ngi, stopOnError = FALSE)

# print raster details 
origin(crp_algn_ngi)
res(crp_algn_ngi)
crs(crp_algn_ngi)

# visualize CNVC raster
plot(crp_algn_ngi)

# Save the raster as 8-bit unsigned integer
ngi_crp_algn_out_path <- file.path(ngi_cropped_dir, "cfga_ngi_crp_algn_cnvc.tif")
writeRaster(crp_algn_ngi, ngi_crp_algn_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

##### use SK raster only since mosaic is producing different resolution

# reproject NGI SK raster to match CRS of pilot area
ngi_sk_proj <- project(r_sk, crs(pa_rast))

# check output
origin(ngi_sk_proj)
res(ngi_sk_proj)
crs(ngi_sk_proj)

# snap origin of NGI raster to match pilot area while keeping its resolution
aligned_ngi_sk <- align(ngi_ski_proj, pa_rast, snap = "near")

# crop to extent of pilot area
cropped_ngi_sk <- crop(aligned_ngi_sk, ext(pa_rast))

# print raster details
origin(cropped_ngi_sk)
res(cropped_ngi_sk)
crs(cropped_ngi_sk)

# visualize
plot(cropped_ngi_sk)

# Save the raster
ngi_crp_algn_out_path <- file.path(ngi_cropped_dir, "cfga_ngi_crp_algn_native_res.tif")
writeRaster(cropped_ngi_sk, ngi_crp_algn_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)



#------------------- ALIGN RASTERS WITH CNVC ---------------------------



################################## OLD ############################################

####--- DIRECTORIES ---####
# study area
sa_dir <- "./processeddatasets/studyarea"

# pilot area
pa_dir <- "./processeddatasets/pilot_area"

# grasslands
raw_grass_data_dir <- "./existingdatasets"
proc_grass_data_dir <- "./processeddatasets/grass"
out_grass_dir <- "./outputs/grass"


####--- Files ---####
# Study Area
sa_file <- file.path(sa_dir, "classification_study_area.shp")

##############################################
# Rasterize study area
# Define projection, origin, resolution

# Load the shapefile
shp <- vect(sa_file) 

# Define raster parameters
res <- 10  # Resolution in map units 
origin <- c(0, 0)  # Origin coordinates
#crs <- "EPSG:10820"  # WGS 1984 Agriculture Canada Albers - did not work
crs <- "EPSG:5070"  # NAD83 / Conus Albers

#crs_wkt <- 'PROJCS["WGS 84 / Agriculture Canada Albers",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["latitude_of_center",40],PARAMETER["longitude_of_center",-96],PARAMETER["standard_parallel_1",44.75],PARAMETER["standard_parallel_2",55.75],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","10820"]]'

# Create an empty raster with desired extent, resolution, and CRS
r <- rast(ext(shp), resolution = res)
crs(r) <- crs
origin(r) <- origin

crs(r)

# Rasterize the shapefile
sa_raster <- rasterize(shp, r)  

# visualize the raster
plot(sa_raster)

# Save the rasteras 8-bit unsigned integer
sa_out_path <- file.path(sa_dir, "classification_study_area_aligned_10820.tif")
writeRaster(sa_raster, sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)


#--- Rasterize pilot area using sa_raster as template ---#

# Load pilot area shapefile
pilot_file <- file.path(pa_dir, "GNP_pilot_area_shifted.shp")
pilot_vect <- vect(pilot_file)

# Project pilot area to match sa_raster
pilot_vect_proj <- project(pilot_vect, sa_raster)
pilot_vect_proj <- project(pilot_vect, "EPSG:5070")
pilot_vect_proj <- project(pilot_vect, crs(sa_raster))

# Rasterize using sa_raster as template
# This will embed the smaller pilot area within the larger raster grid
pilot_raster <- rasterize(pilot_vect_proj, sa_raster, field = 1, background = NA, touches = TRUE)

summary(pilot_raster)

valid <- is.valid(pilot_vect_proj)
table(valid)


plot(sa_raster)
plot(pilot_vect_proj, add = TRUE, col = "red", lw = 6)


ext(sa_raster)
ext(pilot_vect_proj)

crs(pilot_vect)



# Load pilot shapefile using sf
pilot_sf <- st_read(file.path(pa_dir, "GNP_pilot_area_shifted.shp"))

# Reproject using sf to EPSG:5070
pilot_sf_proj <- st_transform(pilot_sf, crs = 5070)

# Convert back to terra
pilot_vect_proj <- vect(pilot_sf_proj)

# Check extent again
ext(pilot_vect_proj)



# visualize the raster
plot(pilot_raster)

# Save pilot area raster
pilot_out_path <- file.path(pa_dir, "GNP_pilot_area_shifted_aligned.tif")
writeRaster(pilot_raster, pilot_out_path, overwrite = TRUE)

message("Saved: ", pilot_out_path)

##############################################
# Merge tiles and then crop to study area
# Save merged and cropped rasters

# If RAM is limited, force temporary disk-backed raster operations
# terraOptions(memfrac = 0.5, tempdir = "D:/terra_temp")

#--- Directories for datasets that need to be merged ---
# AAFC LUTS 30m tiles directory
luts_30_tile_root <- file.path(raw_grass_data_dir, "luts_30m")
luts_30_merged_dir <- file.path(proc_grass_data_dir, "luts_30m/merged")
luts_30_cropped_dir <- file.path(proc_grass_data_dir, "luts_10m/cropped")

# AAFC LUTS 10m tiles directory
luts_10_tile_root <- file.path(raw_grass_data_dir, "luts_10m")
luts_10_merged_dir <- file.path(proc_grass_data_dir, "luts_10m/merged")
luts_10_cropped_dir <- file.path(proc_grass_data_dir, "luts_10m/cropped")

# AAFC ACI tiles directory
aci_tile_root <- file.path(raw_grass_data_dir, "aci")
aci_merged_dir <- file.path(proc_grass_data_dir, "aci/merged")
aci_cropped_dir <- file.path(proc_grass_data_dir, "aci/cropped")

# NGI tiles directory
ngi_tile_root <- file.path(raw_grass_data_dir, "ngi")
ngi_merged_dir <- file.path(proc_grass_data_dir, "ngi/merged")
ngi_cropped_dir <- file.path(proc_grass_data_dir, "ngi/cropped")

# GLC_FCS30D tiles directory
glc_tile_root <- file.path(raw_grass_data_dir, "glc")
glc_merged_dir <- file.path(proc_grass_data_dir, "glc/merged")
glc_cropped_dir <- file.path(proc_grass_data_dir, "glc/cropped")


#--- Merging then alignining and cropping to pilot area ---
# Create a tibble with file paths
tile_info <- tibble(
  path = all_tiles,
  file = basename(all_tiles)
)

# Check if there are at least two tiles
if (nrow(tile_info) < 2) {
  stop("Not enough tiles to mosaic.")
}

message("Processing all tiles together...")

# Read and mosaic tiles
rasters <- lapply(tile_info$path, rast)

# Use sa_aligned as reference for grid alignment
ref_raster <- rast(sa_aligned)

aligned <- lapply(rasters, function(r) {
  if (!compareGeom(r, ref_raster, stopOnError = FALSE)) {
    project(r, ref_raster, method = "near", align = TRUE)
  } else {
    r
  }
})

# mosaic together
mosaic_raster <- Reduce(merge, aligned)

# Crop and mask using ref_raster (sa_aligned)
cropped <- crop(mosaic_raster, ref_raster)
cropped_masked <- mask(cropped, ref_raster)

# Save outputs
writeRaster(mosaic_raster, file.path(luts_30_merged_dir, "LUTS_30m_merged.tif"), overwrite = TRUE)
writeRaster(cropped_masked, file.path(aafc_cropped_dir, "LUTS_30m_sa_aligned.tif"), overwrite = TRUE)

message("Saved: LUTS_30m_merged.tif and LUTS_30m_sa_aligned.tif")



#--- Load classification study area ---
#sa <- st_read(sa_file)

#--- Find all .tif files in LUTS 30m subdirectories ---
all_tiles <- list.files(luts_30_tile_root, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# Extract year from filenames like LU2000_u11_v5.tif
tile_info <- tibble(
  path = all_tiles,
  file = basename(all_tiles),
  year = str_extract(basename(all_tiles), "LU\\d{4}") %>% str_remove("LU")
) %>% 
  filter(!is.na(year))

#--- Loop over each year ---
for (yr in unique(tile_info$year)) {
  message("Processing year: ", yr)
  
  tile_paths <- tile_info %>%
    filter(year == yr) %>%
    pull(path)
  
  if (length(tile_paths) < 2) {
    warning("Skipping year ", yr, ": only ", length(tile_paths), " tile(s) found.")
    next
  }
  
  # Read and mosaic tiles
  rasters <- lapply(tile_paths, rast)
  
  # Tiles do not share a common origin (they are not snapped to a common grid)
  
  # Use r1 as reference for grid alignment
  ref_raster <- rasters[[1]]
  
  aligned <- lapply(rasters, function(r) {
    if (!compareGeom(r, ref_raster, stopOnError = FALSE)) {
      project(r, ref_raster, method = "near", align = TRUE)
    } else {
      r
    }
  })
  mosaic_raster <- Reduce(merge, aligned)
  
  # Reproject study area
  sa_vect <- vect(sa)
  sa_vect_proj <- project(sa_vect, mosaic_raster)
  
  
  # Crop and mask
  cropped <- crop(mosaic_raster, sa_vect_proj)
  cropped_masked <- mask(cropped, sa_vect_proj)
  
  # Save
  # Merged raster
  out_path <- file.path(luts_30_merged_dir, paste0("AAFC", yr, "_merged.tif"))
  writeRaster(mosaic_raster, out_path, overwrite = TRUE)
  
  # Cropped raster  
  out_path <- file.path(aafc_cropped_dir, paste0("AAFC", yr, "_sa.tif"))
  writeRaster(cropped_masked, out_path, overwrite = TRUE)
  
  message("Saved: ", out_path)
}


# ---- USER INPUT ----
input_dirs <- c("folder1/", "folder2/", "folder3/")  # Add as many as needed
output_folder <- "path/to/aligned_rasters/"
output_suffix <- "_aligned.tif"

# ---- LOAD REFERENCE RASTER ----

# reprojected in ArcGIS Pro
ref_arc <- rast("/Grassland/Prairie Habitat Joint Venture PUDL/PHJV_PUDLmask_CANAlbers.tif")

# ---- LOAD AND REPROJECT REFERENCE RASTER ----
#ref_raw <- rast("/Grassland/Prairie Habitat Joint Venture PUDL/PHJV/PHJV_PUDLmask.tif")

# Reproject reference raster to Canada Albers using PROJ.4 string from https://epsg.io/102001
# Albers Equal Area projection
# Standard parallels at 50°N and 70°N
# Central meridian at 96°W
# Latitude of origin at 40°N
# NAD83 datum
# Units in meters

#ref <- project(ref_raw, "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0
#                        +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")


# ---- FIND ALL RASTERS IN MULTIPLE FOLDERS ----
raster_files <- unlist(lapply(input_dirs, function(dir) {
  list.files(dir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
}))

# ---- ALIGN FUNCTION ----
align_raster <- function(r_path, ref, output_folder, suffix) {
  r <- rast(r_path)
  
  # Reproject and resample to match reference
  r_aligned <- project(r, ref, method = "near")
  r_aligned <- resample(r_aligned, ref, method = "near")
  
  # Create output filename
  out_name <- paste0(tools::file_path_sans_ext(basename(r_path)), suffix)
  out_path <- file.path(output_folder, out_name)
  
  # Save aligned raster
  writeRaster(r_aligned, out_path, overwrite = TRUE)
  return(out_path)
}

# ---- PROCESS ALL RASTERS ----
aligned_paths <- sapply(raster_files, align_raster, ref = ref, output_folder = output_folder, suffix = output_suffix)

cat("Aligned rasters saved to:\n")
print(aligned_paths)

# ---- CHECK ALIGNMENT ----
compareGeom(r1, r2_aligned, stopOnError = FALSE)
