library(terra)
library(sf)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(maptiles)
library(ranger)
library(caret)
library(pROC)

#--- Parameters ---
# Resolution of rasters in m
res <- 30

####--- File names ---####
# Study Area
sa_file <- file.path("./processeddatasets/studyarea", "study_area_nobuffer.shp")

# File names
# Boundaries
pilot_s_filename <- "pilot_area_test_aligned_cnvc.shp"
pilot_r_filename <- "pilot_area_test_aligned_cnvc.tif"
vegzones_filename <- "vegetation_zones_of_canada_2020_sa_nobuffer.shp"
ecozones_filename <- "ecozone-prairies.shp"
ecoregions_filename <- "ecoregions-prairies.shp"
census_subdivision_filename <- "lcsd000a25a_e.shp"

# LULC
aafc2020_filename <- "AAFC2020_40kmbuffer_noUS.tif"
aafc2010_filename <- "AAFC2010_40kmbuffer_noUS.tif"
aafc2000_filename <- "AAFC2000_40kmbuffer_noUS.tif"
road_filename <- "NRN_SK_14_0_GPKG_en.gpkg"
#road_filename <- "roadless_data.mpk"
grass_loss_class_filename <- "recentchange_2022_loss_c05_grass_b01_class.tif"
grass_loss_year_filename <- "recentchange_2022_loss_c05_grass_b03_year.tif"
ipast_loss_class_filename <- "recentchange_2022_loss_c10_ipast_b01_class.tif"
ipast_loss_year_filename <- "recentchange_2022_loss_c10_ipast_b03_year.tif"
upast_loss_class_filename <- "recentchange_2022_loss_c04_upast_b01_class.tif"
upast_loss_year_filename <- "recentchange_2022_loss_c04_upast_b03_year.tif"

# Agriculture
crop_suit_2019_filename <- "gp_tillage"

# Protected Areas
pa_filename <- "ProtectedConservedArea_2024"

####--- File paths ---####
# partial path in either the existing or processed data directories

# Boundaries
pilot_filepath <- "./processeddatasets/pilot_area"
vegzones_filepath <- "boundaries/vegetation-zones-of-canada"
ecozones_filepath <- "boundaries/national-ecological-framework-for-canada"
census_subdivision_filepath <- "boundaries/census-subdivisions"

# LULC
aafc_filepath <- "landcover/aafc/cropped"
road_filepath <- "landcover/National Road Network/NRN_RRN_SK_GPKG"
#road_filepath <- "landcover/Roads - Poley et al/Data"
grass_loss_filepath <- "landcover/MostRecentChangeByClass/2022/loss/c05_grass"
ipast_loss_filepath <- "landcover/MostRecentChangeByClass/2022/loss/c10_ipast"
upast_loss_filepath <- "landcover/MostRecentChangeByClass/2022/loss/c04_upast"

# Agriculture
crop_suit_2019_filepath <- "agriculture/OlimbRobinson_GrassToGrain_EcologicalIndicators2019_datashare"

# Protected Areas
pa_filepath <- "management/ProtectedConservedArea_2024"
