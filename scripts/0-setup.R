library(terra)
library(sf)
library(dplyr)
library(stringr)

#--- Parameters ---


####--- Directories ---####
# Boundaries
veg_zones_dir <- "./data/processed/boundaries"

# LULC
raw_lulc_data_dir <- "./data/raw/lulc"
proc_lulc_data_dir <- "./data/processed/lulc"



####--- Files ---####
# Boundaries
veg_zones_file <- file.path(veg_zones_dir, "vegetation_zones_of_canada_2020_prairies.shp")



