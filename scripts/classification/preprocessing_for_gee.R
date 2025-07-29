
###################################################################################
##                  PREPROCESSING OF GRASSLAND PRODUCTS FOR GEE                  ##
###################################################################################

##----------------------- Alison Long - July 2025 -------------------------

## Pixel alignment in GEE is approximate - abstracts away low-level raster alignment details.
## Use terra to align and preprocess rasters locally, 
## then upload them to GEE for large-scale analysis and visualization.

# load packages
library(terra)

# set workign directory
setwd("R:/CanadianGrasslands/existingdatasets")

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
