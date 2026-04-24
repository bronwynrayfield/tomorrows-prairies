# =============================================================================
# 1-merge-tiled-data.R
# =============================================================================
#
# Purpose:
#   Merges tiled raster datasets into single mosaics for downstream processing.
#   Currently handles AAFC LUTS only.
#   Slow to run — only needs to be re-run if raw tiles change.
#
#   If RAM is limited, force temporary disk-backed raster operations:
#   terraOptions(memfrac = 0.5, tempdir = "D:/terra_temp")
#
# Inputs:
#   - AAFC LUTS tiles (existingdatasets/landcover/AAFC LUTS 30m/)
#       Tiles named LU{year}_u{n}_v{n}.tif (e.g. LU2015_u11_v5.tif)
#       Years: 2000, 2010, 2015, 2020
#
# Outputs:
#   - AAFC LUTS merged rasters saved to:
#       processeddatasets/landcover/AAFC LUTS 30m/merged/AAFC{year}_merged.tif
#
# Notes:
#   - AAFC tiles do not share a common origin and are not snapped to a common
#     grid; tiles are aligned to the first tile before mosaicking using
#     project() with method = "near"
#   - Years with fewer than 2 tiles are skipped with a warning
#   - write_rast_tagged() not used here since outputs feed into
#     2-project-data.R which reads from the merged/ subfolder directly
#
# Dependencies:
#   - 0-setup.R (parameters, file paths, functions)
# =============================================================================

source(file.path("scripts/landcover", "0-setup.R"))


#--- AAFC LUTS ----
# AAFC tiles directory
tile_root <- file.path("./existingdatasets", aafc_filepath)
aafc_merged_dir <- file.path("./processeddatasets", aafc_filepath, "merged")

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

  # Save
  # Merged raster
  out_path <- file.path(aafc_merged_dir, paste0("AAFC", yr, "_merged.tif"))
  writeRaster(mosaic_raster, out_path, overwrite = TRUE)
  message("Saved: ", out_path)
}