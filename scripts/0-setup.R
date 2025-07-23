library(terra)
library(sf)
library(dplyr)
library(stringr)

#--- Parameters ---
# Resolution of rasters in m
res <- 30


####--- Directories ---####
# Study area
sa_dir <- "./data/processed/studyarea"

# Boundaries
boundaries_dir <- "./data/processed/boundaries"

# LULC
raw_lulc_data_dir <- "./data/raw/lulc"
proc_lulc_data_dir <- "./data/processed/lulc"
out_lulc_historic_dir <- "./outputs/lulc/historic"


####--- Files ---####
# Study Area
sa_file <- file.path(sa_dir, "study_area_nobuffer.shp")
sa_40kmbuffer_noUS_file <- file.path(sa_dir, "study_area_40kmbuffer_noUS.shp")

# Boundaries
veg_zones_file <- file.path(boundaries_dir, "vegetation_zones_of_canada_2020_sa_nobuffer.shp")

# LULC
aafc_crosswalk_file_2010 <- file.path(proc_lulc_data_dir, "aafc", "aafc-luts-crossawalk-2010.csv")
aafc_crosswalk_file_2015 <- file.path(proc_lulc_data_dir, "aafc", "aafc-luts-crossawalk-2015.csv")


