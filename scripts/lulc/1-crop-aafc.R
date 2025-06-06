##############################################
# Merge AAFC tiles and then crop to prairies
# Save merged and cropped rasters

# If RAM is limited, force temporary disk-backed raster operations
# terraOptions(memfrac = 0.5, tempdir = "D:/terra_temp")

#--- Setup ---
source(file.path("scripts", "0-setup.R"))

#--- Directories ---
# AAFC tiles directory
tile_root <- file.path(raw_lulc_data_dir, "aafc")
aafc_merged_dir <- file.path(proc_lulc_data_dir, "aafc/merged")
aafc_cropped_dir <- file.path(proc_lulc_data_dir, "aafc/cropped")

#--- Load prairie study area based on Vegetation Zones of Canada ---
veg_zones <- st_read(veg_zones_file)

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

  # Reproject veg zones boundary if needed
  prairie_match <- if (!st_crs(veg_zones) == crs(mosaic_raster)) {
    st_transform(prairie, crs(mosaic_raster))
  } else {
    veg_zones
  }
  
  prairie_vect <- vect(prairie_match)
  
  # Crop and mask
  cropped <- crop(mosaic_raster, prairie_vect)
  cropped_masked <- mask(cropped, prairie_vect)
  
  # Save
  # Merged raster
  out_path <- file.path(aafc_merged_dir, paste0("LU", yr, ".tif"))
  writeRaster(mosaic_raster, out_path, overwrite = TRUE)

  # Cropped raster  
  out_path <- file.path(aafc_cropped_dir, paste0("LU", yr, "_prairie.tif"))
  writeRaster(cropped_masked, out_path, overwrite = TRUE)

  message("Saved: ", out_path)
}
