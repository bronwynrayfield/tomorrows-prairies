
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

####--- PARAMETERS ---####
# resolution of rasters in m
res <- 10


####--- DIRECTORIES ---####
# study area
sa_dir <- "./processeddatasets/studyarea"

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
# Specify a field to rasterize
r_rasterized <- rasterize(shp, r, field = "Abbreviation")  

# Save the raster
writeRaster(r_rasterized, "classification_study_area_aligned.tif", overwrite = TRUE)


##############################################
# Merge tiles and then crop to study area
# Save merged and cropped rasters

# If RAM is limited, force temporary disk-backed raster operations
# terraOptions(memfrac = 0.5, tempdir = "D:/terra_temp")

#--- Directories ---
# AAFC LUTS tiles directory
tile_root <- file.path(raw_grass_data_dir, "luts")
luts_merged_dir <- file.path(proc_grass_data_dir, "luts/merged")
luts_cropped_dir <- file.path(proc_grass_data_dir, "luts/cropped")

# AAFC ACI tiles directory
tile_root <- file.path(raw_grass_data_dir, "aci")
aci_merged_dir <- file.path(proc_grass_data_dir, "aci/merged")
aci_cropped_dir <- file.path(proc_grass_data_dir, "aci/cropped")

# NGI tiles directory
tile_root <- file.path(raw_grass_data_dir, "ngi")
ngi_merged_dir <- file.path(proc_grass_data_dir, "ngi/merged")
ngi_cropped_dir <- file.path(proc_grass_data_dir, "ngi/cropped")

# GLC_FCS30D tiles directory
tile_root <- file.path(raw_grass_data_dir, "glc")
glc_merged_dir <- file.path(proc_grass_data_dir, "glc/merged")
glc_cropped_dir <- file.path(proc_grass_data_dir, "glc/cropped")

#--- Load prairie study area ---
sa <- st_read(sa_file)
# sa <- st_read(sa_40kmbuffer_noUS_file)

#--- Find all .tif files in subdirectories ---
all_tiles <- list.files(tile_root, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

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
  out_path <- file.path(aafc_merged_dir, paste0("AAFC", yr, "_merged.tif"))
  writeRaster(mosaic_raster, out_path, overwrite = TRUE)
  
  # Cropped raster  
  out_path <- file.path(aafc_cropped_dir, paste0("AAFC", yr, "_sa_no_buffer.tif"))
  #  out_path <- file.path(aafc_cropped_dir, paste0("AAFC", yr, "_40kmbuffer_noUS.tif"))
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
