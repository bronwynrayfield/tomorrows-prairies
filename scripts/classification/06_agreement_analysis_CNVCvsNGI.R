
###################################################################################
##             ANALYSIS: SPATIAL AGREEMENT SUMMARIES (CNVC vs NGI)               ##
###################################################################################

##----------------------- Alison Long - February 2026 -------------------------

# step 5: run descriptive data summaries to assess similarities and differences
#        among existing grassland datasets. 

# Native/mixed/tame – how do the two maps compare? CNVC map grouped into categories that align 
# with NGI. Pixel-level. Or hex level if CRS differences would inflate difference.
# Where there are mismatches – were some groups harder to map than others? Or at ecoregion level.
# NGI looks like it is mapping ditches as native prairie 
# Division level – do combine native and ruderal into natural grass vs tame. Could start analysis there if it more closely matches what I’m already doing.
# Ruderal = NGI mixed (expectation – want to see if this is the case)
# Plots from Aaron – want to know if those assignments match the group-level CNVC map. Could be mismatch because some plots are from 1980s but maybe hex would help. Their plots assign down to group level – different from what I would do with groundtruth data points (native vs tame etc)

## NGI native = CNVC native grass classes
## NGI mixed = CNVC ruderal
## NGI tame = CNVC ag pasture/hay

# load packages
library(terra)
library(ggplot2)
library(pheatmap)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)
library(ape)
library(sf)
library(viridis)
library(ggspatial)  # for scalebar / north arrow 
library(stringr)
library(patchwork)
library(lme4)
library(broom.mixed)
library(tools)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"

# grass raster directories
ngi_aligned_dir <- file.path(proc_grass_data_dir, "cfga_ngi/aligned")
cnvc_aligned_dir <- file.path(proc_grass_data_dir, "natureserve_cnvc/aligned")

# output 
out_dir <- "R:/CanadianGrasslands/WG_Classification/results/spatial_agreement_analysis"

####--- FILES ---####

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# ecodistricts
ecodist_file   <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean.shp"  

# grass raster files
ngi_file <- file.path(ngi_aligned_dir, "cfga_ngi_merged_reproj_algn_good.tif")
cnvc_file <- file.path(cnvc_aligned_dir, "cnvc_algn_good.tif")

####--- READ INTO TERRA ---####

# ecodistricts
#ecodist <- vect(ecodist_file)

# study area
#sa <- rast(sa_file)

# grass rasters
ngi <- rast(ngi_file)
cnvc <- rast(cnvc_file)


