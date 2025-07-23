##############################################
# Merge AAFC tiles and then crop to prairies
# Save merged and cropped rasters

# If RAM is limited, force temporary disk-backed raster operations
# terraOptions(memfrac = 0.5, tempdir = "D:/terra_temp")
#BR test# 

#--- Setup ---
source(file.path("scripts", "0-setup.R"))

#--- Directories ---
# AAFC tiles directory
tile_root <- file.path(raw_lulc_data_dir, "aafc")
aafc_merged_dir <- file.path(proc_lulc_data_dir, "aafc/merged")
aafc_cropped_dir <- file.path(proc_lulc_data_dir, "aafc/cropped")

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
