
###################################################################################
##                   ANALYSIS: SPATIAL AGREEMENT SUMMARIES                       ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/class_inputs"
proc_grass_data_dir <- "./processeddatasets/grass"
out_grass_dir <- "./outputs/grass"

#-------- PIXELWISE AGREEMENT ------------

# List your raster files (modify pattern/path as needed)
files <- list.files("path/to/rasters", pattern = "\\.tif$", full.names = TRUE)

# Load as a SpatRaster stack
rstack <- rast(files)

# Count how many rasters have *non-NA* data for each pixel
count_nonNA <- app(rstack, fun = function(x) sum(!is.na(x)))

# Save result 
#writeRaster(count_nonNA, "pixel_data_count.tif", overwrite = TRUE)

#-------- AREA PER ECODISTRICT -----------------

# 1) Raster with data/NA
r <- rast("path/to/your_raster.tif")

# 2) Ecodistrict polygons
ed <- vect("path/to/ecodistricts.shp")

# Make sure CRS matches (reproject polygons to raster CRS if needed)
if (!same.crs(r, ed)) {
  ed <- project(ed, crs(r))
}

# (Optional but recommended) Create a stable ID field if none exists
if (is.null(ed$ecod_id)) {
  ed$ecod_id <- seq_len(nrow(ed))
}

# Create a raster of cell areas in hectares
area_ha <- cellSize(r, unit = "ha")   # SpatRaster of per-cell area in hectares

# Create a 1/0 mask of "has data"
has_data <- ifel(!is.na(r), 1, 0)

# Multiply mask by area to get hectares-with-data per cell
ha_with_data <- has_data * area_ha

# Zonal sum by ecodistrict; exact=TRUE weights partial cell overlaps by coverage
zsum <- extract(ha_with_data, ed, fun = sum, na.rm = TRUE, exact = TRUE)

# Join results back to polygons or keep as a table
out <- merge(as.data.frame(ed)[, c("ecod_id")], zsum, by.x = "ecod_id", by.y = "ID", all.x = TRUE)

# Clean up column names
names(out)[names(out) == names(out)[ncol(out)]] <- "hectares_with_data"

# Save as CSV (optional)
write.csv(out, "ecodistrict_hectares_with_data.csv", row.names = FALSE)

# If you also want the polygon layer with the attribute joined:
ed_out <- ed
ed_out$hectares_with_data <- out$hectares_with_data[match(ed_out$ecod_id, out$ecod_id)]
writeVector(ed_out, "ecodistricts_with_hectares.gpkg", filetype = "GPKG", overwrite = TRUE)

