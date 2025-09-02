
###################################################################################
##                  PREPROCESSING OF GRASSLAND PRODUCTS FOR GEE                  ##
###################################################################################

##----------------------- Alison Long - August 2025 -------------------------

## Pixel alignment in GEE is approximate - abstracts away low-level raster alignment details.
## Use terra to align and preprocess rasters locally, 
## then upload them to GEE for large-scale analysis and visualization.

# load packages
library(terra)
library(sf)
library(tidyverse)

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
crs <- "EPSG:10820"  # WGS 1984 Agriculture Canada Albers

# Create an empty raster with desired extent, resolution, and CRS
r <- rast(ext(shp), resolution = res)
crs(r) <- crs
origin(r) <- origin

# Rasterize the shapefile
r_rasterized <- rasterize(shp, r)  

# visualize the raster
plot(r_rasterized)

# Save the raster
sa_out_path <- file.path(sa_dir, "classification_study_area_aligned.tif")
writeRaster(r_rasterized, sa_out_path, overwrite = TRUE)


#--- Rasterize pilot area using r_rasterized as template ---#

# Load pilot area shapefile
pilot_file <- file.path(pa_dir, "GNP_pilot_area_shifted.shp")
pilot_vect <- vect(pilot_file)

# Project pilot area to match r_rasterized
pilot_vect_proj <- project(pilot_vect, r_rasterized)

# Rasterize using r_rasterized as template
# This will embed the smaller pilot area within the larger raster grid
pilot_raster <- rasterize(pilot_vect_proj, r_rasterized, field = 1, background = NA)

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
