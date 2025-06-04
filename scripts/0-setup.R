library(terra)
library(sf)
library(dplyr)
library(stringr)

#--- Parameters ---


####--- Directories ---####
# Boundaries
studyarea_dir <- "/data/processed/boundaries"

# LULC
raw_lulc_data_dir <- "/data/raw/lulc"
proc_lulc_data_dir <- "/data/processed/lulc"



####--- Files ---####
# Boundaries
studyarea_file <- file.path(studyarea_dir, "vegetation_zones_of_canada_2020_prairies.shp")



