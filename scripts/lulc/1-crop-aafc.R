library(terra)
library(sf)
library(dplyr)
library(stringr)

# If RAM is limited, force temporary disk-backed raster operations
# terraOptions(memfrac = 0.5, tempdir = "D:/terra_temp")

#--- Parameters ---
tile_root <- "D:/Projects/E104 - NU Grasslands/aafc"                            # Top-level directory
out_dir <- "D:/Projects/E104 - NU Grasslands/aafc_clipped"                      # Where to save mosaicked outputs
ecozone_file <- "D:/Projects/E104 - NU Grasslands/boundaries/prairies/vegetation_zones_of_canada_2020_prairies.shp"          # Ecozone boundary file

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#--- Load ecozone boundary ---
prairie <- st_read(ecozone_file)

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
for (yr in c("2010","2020")){#unique(tile_info$year)) {
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
  #save.image(file = paste0(yr, "_aafc_workspace_before_merge.RData"))
  mosaic_raster <- Reduce(merge, aligned)

  # Save
  out_path <- file.path(out_dir, paste0("LU", yr, ".tif"))
  writeRaster(mosaic_raster, out_path, overwrite = TRUE)
  
  # Reproject ecozone boundary if needed
  prairie_match <- if (!st_crs(prairie) == crs(mosaic_raster)) {
    st_transform(prairie, crs(mosaic_raster))
  } else {
    prairie
  }
  
  prairie_vect <- vect(prairie_match)
  
  # Crop and mask
  cropped <- crop(mosaic_raster, prairie_vect)
  clipped <- mask(cropped, prairie_vect)
  
  # Save
  out_path <- file.path(out_dir, paste0("LU", yr, "_prairie_clipped.tif"))
  writeRaster(clipped, out_path, overwrite = TRUE)
  
  message("✓ Saved: ", out_path)
}
