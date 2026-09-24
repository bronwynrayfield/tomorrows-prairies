# set common libraries, parameters, functions, files, and directories

library(terra)
library(tidyterra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)
library(pROC)
library(usdm)
library(patchwork)
library(geodata)
library(foreign)
library(ncf)

# check available space on different mount points and memory
system("df -h")

# don't use scientific noation unless necessary
options(scipen = 999)

# set terra options for memory, temp directory, and progress bar
terraOptions(memfrac = 0.8,
             tempdir = "D:/terra-temp",
             progress = 1)

# --- Functions ----
source(file.path("scripts/landcover/functions", "functions-general.R"))
source(file.path("scripts/landcover/functions", "reclass-aafc-luts.R"))

#--- Parameters ----
# Resolution of rasters in m
res_m <- 30

#--- File names ----
# Study Area
sa_filename <- "modelling_study_area_aligned_cnvc_Jan2026.shp"
sa_raster_filename <- "modelling_study_area_aligned_cnvc_Jan2026.tif"
me_filename <- "modelling_extent_aligned_cnvc_Jan2026.shp"
me_raster_filename <- "modelling_extent_aligned_cnvc_Jan2026.tif"

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
road_filename <- "roads"
rail_filename <- "rails"
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
aci_filename <- "aci_merged.tif"
plowprint_filename <- "plow_reproj_algn_masked.tif"

# Protected Areas
pa_filename <- "ProtectedConservedArea_2024"
sk_cce_filename <- "CrownCEs_-4202437152155019093"

# Land tenure
ncc_easement_filename <- "NatureUnited_Ownership"

# Soil
oc_filename <- "oc_iso.10694.1995.mg.cm3_m_30m_b0cm..30cm_20050101_20101231_g_epsg.4326_v20250204_clip.tif"
bd_filename <- "bd.core_iso.11272.2017.g.cm3_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250204_clip.tif"
clay_filename <- "clay.tot_iso.11277.2020.wpct_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250523_clip.tif"
sand_filename <- "sand.tot_iso.11277.2020.wpct_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250523_clip.tif"
silt_filename <- "silt.tot_iso.11277.2020.wpct_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250523_clip.tif"
ph_filename <- "ph.h2o_iso.10390.2021.index_m_30m_b0cm..30cm_20200101_20221231_g_epsg.4326_v20250204_clip.tif"
psssm_mean_filename <- "psssm_mean_2011_2020.tif"
psssm_sd_filename <- "psssm_sd_2011_2020.tif"
tsm_planting_filename <- "tsm_planting_2020.tif"
tsm_growing_filename <- "tsm_growing_2020.tif"
dss_drainage_filename <- "dss_drainage.tif"
dss_watertbl_filename <- "dss_watertbl.tif"
dss_kp0_filename <- "dss_kp0.tif"
dss_kp33_filename <- "dss_kp33.tif"
dss_kp1500_filename <- "dss_kp1500.tif"
dss_ldepth_filename <- "dss_ldepth.tif"
dss_ldepth_max_filename <- "dss_ldepth_max.tif"

# Topography
w <- 27

dem_filename <- "MRDEM Digital Terrain Model (COG).tif" 
dss_slope_filename <- "dss_slope.tif"
tri_filename <- paste("tri", paste0(w, "x", w), sep = "-")
tpi_filename <- paste("tpi", paste0(w, "x", w), sep = "-")
srr_filename <- paste("srr", paste0(w, "x", w), sep = "-")
northness_filename <- paste("northness", paste0(w, "x", w), sep = "-")
eastness_filename <- paste("eastness", paste0(w, "x", w), sep = "-")

# Climate
hist_period <- "1991-2020"
fut_period <- "2021-2050"
ssp <- "ssp245"

gddgrow_5_hist_filename <- paste("gddgrow_5", ssp, hist_period, "p50.tif", sep="_") 
gddgrow_5_fut_filename <- paste("gddgrow_5", ssp, fut_period, "p50.tif", sep="_") 
hddheat_18_hist_filename <- paste("hddheat_18", ssp, hist_period, "p50.tif", sep="_") 
hddheat_18_fut_filename <- paste("hddheat_18", ssp, fut_period, "p50.tif", sep="_") 
prtot_hist_filename <- paste("prtot", ssp, hist_period, "p50.tif", sep="_") 
prtot_fut_filename <- paste("prtot", ssp, fut_period, "p50.tif", sep="_") 
tn_mean_hist_filename <- paste("tn_mean", ssp, hist_period, "p50.tif", sep="_") 
tn_mean_fut_filename <- paste("tn_mean", ssp, fut_period, "p50.tif", sep="_") 
tx_mean_hist_filename <- paste("tx_mean", ssp, hist_period, "p50.tif", sep="_") 
tx_mean_fut_filename <- paste("tx_mean", ssp, fut_period, "p50.tif", sep="_") 
nr_cdd_hist_filename <- paste("nr_cdd", ssp, hist_period, "p50.tif", sep="_") 
nr_cdd_fut_filename <- paste("nr_cdd", ssp, fut_period, "p50.tif", sep="_") 

# Socio-economic
canfax2009_filename <- "canfax_results_tbl_2005-2009.shp"
canfax2024_filename <- "canfax_results_tbl_2020-2024.shp"
elevator_filename <- "Prairies_ProcessElevators_Aug2026.shp"
popctr2011_filename <- "gpc_000a11a_e"
popctr2021_filename <- "lpc_000a21a_e"

#--- Directories ----
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
rail_filepath <- "landcover/National Railway Network"
grass_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c05_grass"
ipast_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c10_ipast"
upast_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c04_upast"
crop_loss_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/loss/c03_crop"
crop_gain_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/gain/c03_crop"
builtup_gain_filepath <- "landcover/ECCC_Prairies_Change/Prairies_Change_V1/MostRecentChangeByClass/2022/gain/c06_builtup"

# Agriculture
crop_suit_2019_filepath <- "cropland/OlimbRobinson_GrassToGrain_EcologicalIndicators2019_datashare"
aci_filepath <- "cropland/Annual Crop Inventory"
plowprint_filepath <- "cropland/Plowprint WWF"

# Protected Areas
pa_filepath <- "management/Canadian Protected and Conserved Areas Database/ProtectedConservedArea_2024"
sk_cce_filepath <- "management/SK Crown Conservation Easements"

# Land Tenure
ncc_easement_filepath <- "management/Tenure data (NCC)"

# Soil
soil_filepath <- "soil"

# Soil Moisture
psssm_filepath <- "soil/psssm"
tsm_filepath <- "soil/tsm"
dss_filepath <- "soil/detailed-soil-survey"

# Topography
dem_filepath <- "topography/NRCan_MRDEM_30m"

# Climate
climate_filepath <- "climate/CanDCS-M6/geotiff"

# Socio-economic
canfax_filepath <- "socio-economic/canfax"
elevator_filepath <- "socio-economic/process-elevators"
popctr_filepath <- "socio-economic/StatsCan Population Centres"

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

#--- LULC class codes ----
# Class scheme (identical in both years, so the maps cross-tabulate):
#   1 = native grassland   2 = tame grassland   3 = annual cropland
#   4 = settlement         5 = other            NA = outside modelling extent
cls_native <- 1
cls_tame   <- 2
cls_annual <- 3
cls_settlement <- 4
cls_other  <- 5

lulc_levels <- data.frame(
  value = c(cls_native, cls_tame, cls_annual, cls_settlement, cls_other),
  class = c("native grassland", "tame grassland", "annual cropland",
            "settlement", "other")
)