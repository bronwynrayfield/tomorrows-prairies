
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS TO REFERENCE                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 3: align all grassland layers to classification study area file from step 1
#         to prep for analysis

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_data_dir  <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada"
proc_data_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts"

# Ecodistrict directory
ecodist_root <- file.path(raw_data_dir, "Ecodistricts")

####--- Files ---####

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# ecodistrict shapefile
ecodist_file <- file.path(ecodist_root, "ecodistricts_clean.shp")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
sa

####--- ALGIN RASTERS ---#### 

ecodist_vect <- vect(ecodist_file)
plot(ecodist_vect)
ecodist_vect

# Reproject ecodistrict shapefile to match reference raster CRS
ecodist_vect_proj <- project(ecodist_vect, crs(sa))

# Rasterize classification study area using CNVC raster as template
ecodist_rast <- rasterize(ecodist_vect_proj, sa, 
                          field = "ECODISTRIC", 
                          touches = TRUE) # ensures that any cell touched by the polygon is included

# print raster details 
origin(ecodist_rast)
ecodist_rast
plot(ecodist_rast)

# Make categorical
ecodist_rast <- as.factor(ecodist_rast)
plot(ecodist_rast)

# Save the raster as 8-bit unsigned integer
ecodist_out_path <- file.path(proc_data_dir, "ecodistricts_clean_rast.tif")
writeRaster(ecodist_rast, ecodist_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)
