

###################################################################################
##           PREPROCESSING: GENERATING STUDY AREA REFERENCE RASTER               ##
###################################################################################

##----------------------- Alison Long - December 2025 -------------------------

# step 1: create classification study area aligned to NatureServe CNVC raster
#         that all other study and pilot areas will nest within

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands")

####--- DIRECTORIES ---####

# study area
sa_dir <- "./processeddatasets/Studyarea"

# grasslands
raw_grass_data_dir <- "./existingdatasets/class_inputs"
proc_grass_data_dir <- "./processeddatasets/grass"
out_grass_dir <- "./outputs/grass"

####--- Files ---####

# Study Area
sa_file <- file.path(sa_dir, "classification_study_area.shp")

# CNVC raster
cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

####--- Use CNVC as snap raster ---####

# Load CNVC raster
cnvc <- rast(cnvc_file)

# print raster details 
origin(cnvc)
res(cnvc)
crs(cnvc)

# visualize CNVC raster
plot(cnvc)

# Load classification study area shapefile
sa_vect <- vect(sa_file)

# Reproject classification study area to match CNVC raster CRS
sa_vect_proj <- project(sa_vect, crs(cnvc))

# Rasterize classification study area using CNVC raster as template
sa_rast <- rasterize(sa_vect_proj, cnvc, field=1, touches=TRUE) # field=1 assigns constant value to rasterized cells
                                                                 # touches=TRUE ensures that any cell touched by the polygon is included

# print raster details 
origin(sa_rast)
res(sa_rast)
crs(sa_rast)

# Get frequency table of raster values
freq_table <- freq(sa_rast)

# Print the table - 4,434,655,083 pixels with value 1
print(freq_table)

# visualize the new raster
plot(sa_rast)

# Trim NA-only rows/columns
sa_rast_trimmed <- trim(sa_rast)
plot(sa_rast_trimmed)

# Check extents - they differ...
ext(sa_rast)
ext(sa_rast_trimmed)

nrow(sa_rast); ncol(sa_rast)
nrow(sa_rast_trimmed); ncol(sa_rast_trimmed)

# but origin the same!
origin(sa_rast)
res(sa_rast)
crs(sa_rast)

# Save the raster as 8-bit unsigned integer
sa_out_path <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")
writeRaster(sa_rast_trimmed, sa_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# Convert raster to polygons (excluding NA values)
sa_polygons <- as.polygons(sa_rast_trimmed, na.rm = TRUE, dissolve = TRUE)

# Export shapefile
sa_out_shp <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.shp")
writeVector(sa_polygons, sa_out_shp, overwrite = TRUE)

