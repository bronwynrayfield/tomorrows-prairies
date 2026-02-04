
###################################################################################
##               PREPROCESSING: CROPPING RASTERS TO PILOT AREA                   ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------

# step 4: crop aligned grassland layers to pilot area file from step 2
#        to prep for pilot analysis

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

# NGI tiles directory
ngi_tile_root <- file.path(raw_grass_data_dir, "cfga_ngi")
ngi_merged_dir <- file.path(proc_grass_data_dir, "cfga_ngi/merged")
ngi_cropped_dir <- file.path(proc_grass_data_dir, "cfga_ngi/cropped")
ngi_aligned_dir <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")

# ECCC directory
eccc_root <- file.path(raw_grass_data_dir, "eccc_grass")
eccc_cropped_dir <- file.path(proc_grass_data_dir, "eccc_grass/cropped")

# PHJV PUDL directory
pudl_root <- file.path(raw_grass_data_dir, "phjv_pudl")
pudl_cropped_dir <- file.path(proc_grass_data_dir, "phjv_pudl/cropped")

# Central Grasslands Roadmap directory
cgr_root <- file.path(raw_grass_data_dir, "cgr_gam_adv")
cgr_cropped_dir <- file.path(proc_grass_data_dir, "cgr_gam_adv/cropped")

# WRI annual grass directory
wri_root <- file.path(raw_grass_data_dir, "wri_gpw_ann")
wri_cropped_dir <- file.path(proc_grass_data_dir, "wri_gpw_ann/cropped")


####--- Files ---####
# Pilot Area
pa_file <- file.path(pa_dir, "GNP_pilot_area_crop_CSD.shp")

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# NGI rasters
ngi_ab_file <- file.path(ngi_tile_root, "AB-002.tif")
ngi_mb_file <- file.path(ngi_tile_root, "MB-003.tif")
ngi_sk_file <- file.path(ngi_tile_root, "SK-001.tif")

# AAFC LUTS 

# ECCC grass summary
eccc_file <- file.path(eccc_root, "grass_summary_a_class.tif")

# WRI
wri_file <- file.path(wri_root, "gpw_grassland_rf.med.filt.bthr_c_30m_20240101_20241231_go_epsg.4326_v2.tif")

# PUDL
pudl_file <- file.path(pudl_root, "PHJV_PUDLmask.tif")

# Central Grasslands Roadmap - Grassland Assessment Map (Advanced)
cgr_file <- file.path(cgr_root, "CGR_GAM_V2_separated.tif")