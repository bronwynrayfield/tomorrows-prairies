library(terra)
library(sf)
library(dplyr)
#library(stringr)
library(tidyr)
library(ggplot2)
#library(maptiles)
library(ranger)
#library(caret)
#library(pROC)


#
# check available space on different mount points
system("df -h")

# set terra options for memory, temp directory, and progress bar
terraOptions(memfrac = 0.8,
             tempdir = "D:/terra-temp",
             progress = 1)

#--- Parameters ---
# Resolution of rasters in m
res <- 30

#--- File names ----
# Study Area
sa_filename <- "modelling_study_area_aligned_cnvc_Jan2026.shp"
sa_raster_filename <- "modelling_study_area_aligned_cnvc_Jan2026.tif"
me_filename <- "modelling_extent_aligned_cnvc_Jan2026.shp"
me_raster_filename <- "modelling_extent_aligned_cnvc_Jan2026.tif"

# File names
# Boundaries
pilot_filename <- "pilot_area_aligned_cnvc_Jan2026.shp"
vegzones_filename <- "vegetation_zones_of_canada_2020_sa_nobuffer.shp"
ecozones_filename <- "ecozones.shp"
ecoregions_filename <- "ecoregions_clean.shp"
ecodistricts_filename <- "ecodistricts_clean.shp"
census_subdivision_filename <- "lcsd000a25a_e.shp"

# LULC
aafc2020_filename <- "AAFC2020_merged.tif"
aafc2015_filename <- "AAFC2015_merged.tif"
aafc2010_filename <- "AAFC2010_merged.tif"
aafc2000_filename <- "AAFC2000_merged.tif"
#road_filename <- "NRN_SK_14_0_GPKG_en.gpkg"
road_filename <- "roads.gdb"
grass_loss_class_filename <- "recentchange_2022_loss_c05_grass_b01_class.tif"
grass_loss_year_filename <- "recentchange_2022_loss_c05_grass_b03_year.tif"
ipast_loss_class_filename <- "recentchange_2022_loss_c10_ipast_b01_class.tif"
ipast_loss_year_filename <- "recentchange_2022_loss_c10_ipast_b03_year.tif"
upast_loss_class_filename <- "recentchange_2022_loss_c04_upast_b01_class.tif"
upast_loss_year_filename <- "recentchange_2022_loss_c04_upast_b03_year.tif"
crop_loss_class_filename <- "recentchange_2022_loss_c03_crop_b01_class.tif"
crop_loss_year_filename <- "recentchange_2022_loss_c03_crop_b03_year.tif"
crop_gain_class_filename <- "recentchange_2022_gain_c03_crop_b01_class.tif"
crop_gain_year_filename <- "recentchange_2022_gain_c03_crop_b03_year.tif"
builtup_gain_class_filename <- "recentchange_2022_gain_c06_builtup_b01_class.tif"
builtup_gain_year_filename <- "recentchange_2022_gain_c06_builtup_b03_year.tif"

# Agriculture
crop_suit_2019_filename <- "gp_tillage"

# Protected Areas
pa_filename <- "ProtectedConservedArea_2024"

# Soil
oc_filename <- "oc_iso.10694.1995.mg.cm3_m_30m_b0cm..30cm_20050101_20101231_g_epsg.4326_v20250204_clip.tif"
bd_filename <- "bd.core_iso.11272.2017.g.cm3_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250204_clip.tif"
clay_filename <- "clay.tot_iso.11277.2020.wpct_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250523_clip.tif"
sand_filename <- "sand.tot_iso.11277.2020.wpct_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250523_clip.tif"
silt_filename <- "silt.tot_iso.11277.2020.wpct_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250523_clip.tif"
ph_filename <- "ph.h2o_iso.10390.2021.index_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250204_clip.tif"

# Topography
dem_filename <- "MRDEM Digital Terrain Model (COG).tif" 
slope_filename <- "slope.tif"
aspect_filename <- "aspect.tif"

#--- File paths ----
# partial path in either the existing or processed data directories

sa_filepath <- "./processeddatasets/studyarea"

# Boundaries
pilot_filepath <- "pilot_area"
vegzones_filepath <- "boundaries/vegetation-zones-of-canada"
ecozones_filepath <- "boundaries/national-ecological-framework-for-canada/Ecozones"
ecoregions_filepath <- "boundaries/national-ecological-framework-for-canada/Ecoregions"
ecodistricts_filepath <- "boundaries/national-ecological-framework-for-canada/Ecodistricts"
census_subdivision_filepath <- "boundaries/census-subdivisions"

# LULC
aafc_filepath <- "landcover/AAFC LUTS 30m"
#road_filepath <- "landcover/National Road Network/NRN_RRN_SK_GPKG"
road_filepath <- "landcover/Roads - Poley et al 2021"
grass_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c05_grass"
ipast_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c10_ipast"
upast_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c04_upast"
crop_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c03_crop"
crop_gain_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/gain/c03_crop"
builtup_gain_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/gain/c06_builtup"

# Agriculture
crop_suit_2019_filepath <- "cropland/OlimbRobinson_GrassToGrain_EcologicalIndicators2019_datashare"

# Protected Areas
pa_filepath <- "management/Canadian Protected and Conserved Areas Database/ProtectedConservedArea_2024"

# Soil
soil_filepath <- "soil"

# Topography
dem_filepath <- "topography/NRCan_MRDEM_30m"

#--- Scale ----
oc_scale <- 10
bd_scale <- 100
ph_scale <- 10

#--- No Data ----
oc_nodata <- 32767
bd_nodata <- 32767
clay_nodata <- 255
sand_nodata <- 255
silt_nodata <- 255
ph_nodata <- 255