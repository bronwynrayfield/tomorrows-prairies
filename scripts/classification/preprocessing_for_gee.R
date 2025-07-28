
###################################################################################
##                  PREPROCESSING OF GRASSLAND PRODUCTS FOR GEE                  ##
###################################################################################

##----------------------- Alison Long - July 2025 -------------------------

## Pixel alignment in GEE is approximate - abstracts away low-level raster alignment details.
## Use terra to align and preprocess rasters locally, 
## then upload them to GEE for large-scale analysis and visualization.

# load packages
library(terra)

# ---- USER INPUT ----
ref_path <- "path/to/reference.tif"
input_dirs <- c("folder1/", "folder2/", "folder3/")  # Add as many as needed
output_folder <- "path/to/aligned_rasters/"
output_suffix <- "_aligned.tif"

# ---- LOAD REFERENCE RASTER ----
#ref <- rast(ref_path)

# ---- LOAD AND REPROJECT REFERENCE RASTER ----
ref_raw <- rast("R:/CanadianGrasslands/existingdatasets/Grassland/Prairie Habitat Joint Venture PUDL/PHJV/PHJV_PUDLmask.tif")

# Define target CRS (e.g., Canada Albers Equal Area)
target_crs <- "EPSG:102001"  # Or use a custom PROJ string

# Reproject reference raster
ref <- project(ref_raw, target_crs, method = "near")


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
