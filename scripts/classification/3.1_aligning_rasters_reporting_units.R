
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

################### DIRECTORIES ###################

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
#raw_data_dir  <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada"
proc_data_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada"
proc_data_dir <- "R:/CanadianGrasslands/processeddatasets/boundaries"

# file directories
ecodist_root <- file.path(proc_data_dir, "Ecodistricts")
ecoreg_root  <- file.path(proc_data_dir, "Ecoregions")
ecoprov_root <- file.path(proc_data_dir, "Ecoprovinces")
ccs_root     <- file.path(proc_data_dir, "census-consolidated-subdivisions")

################### FILES ###################

# Study Area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# ecodistrict shapefile
ecodist_file <- file.path(ecodist_root, "ecodistricts_clean_aligned_Oct2026.shp")

# ecoregion shapefile
ecoreg_file <- file.path(ecoreg_root, "ecoregions_clean_classSA_Oct2026.shp")

# ecoprovince shapefile
ecoprov_file <- file.path(ecoprov_root, "ecoprovince_clean_classSA_Oct2026.shp")

# CCS shapefile
ccs_file <- file.path(ccs_root, "CCS_2021_aligned_class_sa.shp")

#------ use study area file as snap raster ------#

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa)
sa

################### ALIGN SHAPEFILES ###################

#------ ECODISTRICTS ------#
ecodist_vect <- vect(ecodist_file)
plot(ecodist_vect)
ecodist_vect

# Reproject ecodistrict shapefile to match reference raster CRS
ecodist_vect_proj <- project(ecodist_vect, sa)
ecodist_vect_proj
plot(ecodist_vect_proj)

# Export aligned shapefile
writeVector(ecodist_vect_proj, file.path(ecodist_root, "ecodistricts_clean_reproj.shp"))

#------ ECOREGIONS ------#
ecoreg_vect <- vect(ecoreg_file)
plot(ecoreg_vect)
ecoreg_vect

# Reproject ecoregrict shapefile to match reference raster CRS
ecoreg_vect_proj <- project(ecoreg_vect, sa)
ecoreg_vect_proj
plot(ecoreg_vect_proj)

# Export aligned shapefile
writeVector(ecoreg_vect_proj, file.path(ecoreg_root, "ecoregions_clean_aligned_classSA.shp"))

#------ ECOPROVINCES ------#
ecoprov_vect <- vect(ecoprov_file)
plot(ecoprov_vect)
ecoprov_vect

# Reproject ecoprovince shapefile to match reference raster CRS
ecoprov_vect_proj <- project(ecoprov_vect, sa)
ecoprov_vect_proj
plot(ecoprov_vect_proj)

# Export aligned shapefile
writeVector(ecoprov_vect_proj, file.path(ecoprov_root, "ecoprovinces_clean_aligned_classSA.shp"))

#------ CCS ------#
ccs_vect <- vect(ccs_file)
plot(ccs_vect)
ccs_vect

# Reproject CCS shapefile to match reference raster CRS
ccs_vect_proj <- project(ccs_vect, sa)
ccs_vect_proj
plot(ccs_vect_proj)

# Export aligned shapefile
writeVector(ccs_vect_proj, file.path(ccs_root, "ccs_clean_reproj.shp"))

################### ALIGN RASTERS ###################

#------ Rasterize ECODISTRICT file using CNVC raster as template ------#
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
ecodist_out_path <- file.path(proc_data_dir, "Ecodistricts/ecodistricts_clean_rast.tif")
writeRaster(ecodist_rast, ecodist_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

#------ Rasterize ECOREGION file using CNVC raster as template ------#
ecoreg_rast <- rasterize(ecoreg_vect_proj, sa, 
                          field = "ECOREGION_", 
                          touches = TRUE) # ensures that any cell touched by the polygon is included

# print raster details 
origin(ecoreg_rast)
ecoreg_rast
plot(ecoreg_rast)

# Make categorical
ecoreg_rast <- as.factor(ecoreg_rast)
plot(ecoreg_rast)

# Save the raster as 8-bit unsigned integer
ecoreg_out_path <- file.path(proc_data_dir, "Ecoregions/ecoregions_clean_aligned_classSA_rast.tif")
writeRaster(ecoreg_rast, ecoreg_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

#------ Rasterize ECOPROVINCE file using CNVC raster as template ------#

# Convert ecoprovince codes to integers (decimals won't be preserved in raster)
ecoprov_vect_proj$ECOPROV_INT <- round(ecoprov_vect_proj$ECOPROVINC * 10)
ecoprov_vect_proj

# Now rasterize
ecoprov_rast <- rasterize(ecoprov_vect_proj, sa, 
                         field = "ECOPROV_INT", 
                         touches = TRUE) # ensures that any cell touched by the polygon is included

# print raster details 
origin(ecoprov_rast)
ecoprov_rast
plot(ecoprov_rast)

# Make categorical
ecoprov_rast <- as.factor(ecoprov_rast)
plot(ecoprov_rast)

# Save the raster as 8-bit unsigned integer
ecoprov_out_path <- file.path(proc_data_dir, "Ecoprovinces/ecoprovinces_clean_aligned_classSA_rast.tif")
writeRaster(ecoprov_rast, ecoprov_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            #NAflag = 255,
            overwrite = TRUE)

#------ Rasterize CCS file using CNVC raster as template ------#
ccs_rast <- rasterize(ccs_vect_proj, sa, 
                          field = "CCSUID", 
                          touches = TRUE) # ensures that any cell touched by the polygon is included

# print raster details 
origin(ccs_rast)
ccs_rast
plot(ccs_rast)

# Make categorical
#ccs_rast <- as.factor(ccs_rast)
#plot(ccs_rast)

# Save the raster as 8-bit unsigned integer
ccs_out_path <- file.path(ccs_root, "ccs_aligned_class_sa_rast.tif")
writeRaster(ccs_rast, ccs_out_path, datatype = "INT2U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 65535,
            overwrite = TRUE)

############################################################
############################################################

###### Check for mismatches in ecological unit nesting ######

# This is a metadata validation, not a spatial validation of boundary alignment

# Read shapefile
#nef <- st_read("R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/archive/ecodistricts_clean_Jun2026.shp")
#nef <- st_read("R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_class_sa.shp")
nef <- st_read("R:/CanadianGrasslands/processeddatasets/boundaries/national-ecological-framework-for-canada/Ecodistricts/ecodistricts_clean_aligned_Oct2026.shp")

# Drop geometry for checks
df <- st_drop_geometry(nef)

str(df)

# Rename fields
df <- df %>%
  rename(
    ecodist = ECODISTRIC,
    ecoreg  = ECOREGION_,
    ecozone = ECOZONE_ID,
    ecoprov = ECOPROVINC
  )

# Address decimals in ecoprov
#ecoprov_vect_proj$ECOPROV_INT <- round(ecoprov_vect_proj$ECOPROVINC * 10)
df <- df %>% mutate(
    ecoprov = round(ecoprov, 1)
  )

# Use categorical version for logic
df <- df %>%
  mutate(ecoprov_cat = sprintf("%.1f", ecoprov))

# -------------------------------
# 1. Build expected lookup tables
# -------------------------------

# Expected mapping: ecodist -> ecoreg
lookup_ecodist <- df %>%
  count(ecodist, ecoreg) %>%
  group_by(ecodist) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup() %>%
  select(ecodist, expected_ecoreg = ecoreg)

# Expected mapping: ecoreg -> ecoprov
lookup_ecoreg <- df %>%
  count(ecoreg, ecoprov_cat) %>%
  group_by(ecoreg) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup() %>%
  select(ecoreg, expected_ecoprov = ecoprov_cat)

# Expected mapping: ecoprov -> ecozone
lookup_ecoprov <- df %>%
  count(ecoprov_cat, ecozone) %>%
  group_by(ecoprov_cat) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup() %>%
  select(ecoprov_cat, expected_ecozone = ecozone)

# -----------------------------------------
# 2. Join expected hierarchy back to dataset
# -----------------------------------------

df_check <- df %>%
  left_join(lookup_ecodist, by = "ecodist") %>%
  left_join(lookup_ecoreg, by = "ecoreg") %>%
  left_join(lookup_ecoprov, by = "ecoprov_cat")

# -----------------------------------------
# 3. Flag violations
# -----------------------------------------

df_check <- df_check %>%
  left_join(
    lookup_ecoreg %>%
      rename(expected_ecoreg_to_ecoprov = expected_ecoprov),
    by = c("expected_ecoreg" = "ecoreg")
  ) %>%
  left_join(
    lookup_ecoprov %>%
      rename(expected_ecoprov_to_ecozone = expected_ecozone),
    by = c("expected_ecoprov" = "ecoprov_cat")
  ) 

df_check <- df_check %>%
  mutate(
    # direct mismatches
    ecodist_mismatch = ecoreg != expected_ecoreg,
    ecoreg_mismatch  = ecoprov_cat != expected_ecoprov,
    ecoprov_mismatch = ecozone != expected_ecozone,
    
    # propagated hierarchy violations (cleaner + safer)
    ecodist_to_ecoprov_violation =
      ecoprov_cat != expected_ecoreg_to_ecoprov,
    
    ecodist_to_ecozone_violation =
      ecozone != expected_ecoprov_to_ecozone,
    
    any_violation =
      ecodist_mismatch |
      ecoreg_mismatch |
      ecoprov_mismatch |
      ecodist_to_ecoprov_violation |
      ecodist_to_ecozone_violation
  )

# -----------------------------------------
# 4. Extract problematic features
# -----------------------------------------

violations <- nef %>%
  slice(which(df_check$any_violation))

# Save results
st_write(violations, file.path(proc_data_dir, "nef_hierarchy_violations_Oct2026.shp"), delete_layer = TRUE)

# -----------------------------------------
# 5. Summary output
# -----------------------------------------

summary <- df_check %>%
  summarise(
    total = n(),
    violations = sum(any_violation),
    pct = 100 * violations / total
  )

print(summary)

# Optional: breakdown by type
colSums(df_check %>%
          select(ecodist_mismatch,
                 ecoreg_mismatch,
                 ecoprov_mismatch,
                 ecodist_to_ecoprov_violation,
                 ecodist_to_ecozone_violation))
