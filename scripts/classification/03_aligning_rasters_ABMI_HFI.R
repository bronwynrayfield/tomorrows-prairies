
###################################################################################
##               PREPROCESSING: ALIGNING RASTERS TO REFERENCE                    ##
###################################################################################

##----------------------- Alison Long - January 2026 -------------------------##

# step 3: align all grassland layers to classification study area file from step 1
#         to prep for analysis

# load packages
library(terra)
library(sf)
library(tidyverse)

# Set the temporary directory to a drive with sufficient space
terraOptions(tempdir = "R:/CanadianGrasslands/WG_classification/scratch")

####--- DIRECTORIES ---####

# pilot area
pa_dir <- "R:/CanadianGrasslands/processeddatasets/pilotarea"

# study area
sa_dir <- "R:/CanadianGrasslands/processeddatasets/studyarea"

# grasslands
raw_grass_data_dir <- "R:/CanadianGrasslands/existingdatasets"
proc_grass_data_dir <- "R:/CanadianGrasslands/WG_Classification/processed"
out_grass_dir <- "./outputs/grass"

# hfi directory
hfi_root <- file.path(raw_grass_data_dir, "Landcover/ABMI_Human_Footprint_Inventory")
hfi_aligned_dir <- file.path(proc_grass_data_dir, "abmi_hfi/aligned")
hfi_cropped_dir <- file.path(proc_grass_data_dir, "abmi_hfi/cropped")
hfi_rast_dir <- file.path(proc_grass_data_dir, "abmi_hfi/rasterized")
hfi_pilot_dir <- file.path(proc_grass_data_dir, "abmi_hfi/pilot")

####--- Files ---####
# Pilot area
pa_file <- file.path(pa_dir, "pilot_area_aligned_cnvc_Jan2026.tif")

# Study area
sa_file <- file.path(sa_dir, "class_study_area_aligned_cnvc_Dec2025.tif")

# CNVC raster
#cnvc_file <- file.path(raw_grass_data_dir, "ns_cnvc/NatureServe_CNVC_COMPRESSED_20251023.tif")

# hfi shapefile
hfi_file <- file.path(hfi_root, "ABMI_HFI_2022_cultivation.shp")

##### use study area file as snap raster #####

# Load study area raster
sa <- rast(sa_file)

# print raster details 
origin(sa) # 6.696121 -2.550334
sa

####--- ALGIN RASTERS ---#### 

#----- ABMI Human Footprint Inventory ----#

# Load vector file
v_hfi <- vect(hfi_file)
v_hfi

#class       : SpatVector 
#geometry    : polygons 
#dimensions  : 396434, 10  (geometries, attributes)
#extent      : 181243.9, 863395.3, 5426029, 6651749  (xmin, xmax, ymin, ymax)
#source      : ABMI_HFI_2022_cultivation.shp
#coord. ref. : NAD83 / Alberta 10-TM (Forest) (EPSG:3400) 
#names       : SOURCE          HFI_ID    FEATURE_TY  YEAR BNDRY_SOUR Modifier_F Modifier_Y YEAR_SOURC Shape_Leng
#type        :  <chr>           <chr>         <chr> <int>      <chr>      <chr>      <int>      <chr>      <num>
#values      : ABMI07 {ABCC15EA-5FAF~ ROUGH_PASTURE     0     ABMI07         NA          0         NA      533.6
#              ABMI12 {AACE1F39-1A91~  TAME_PASTURE     0     ABMI12         NA          0         NA       1042
#              ABMI07 {6CF796DD-A758~  TAME_PASTURE     0     ABMI07         NA          0         NA      480.7
#        Shape_Area
#        <num>
#          1.647e+04
#        4.631e+04
#        9219

# Print unique values in the FEATURE_TY column
unique_values <- unique(v_hfi$FEATURE_TY)
print(unique_values)

#[1] "ROUGH_PASTURE"         "TAME_PASTURE"          "CROP"                  "CULTIVATION_ABANDONED"
#[5] "FRUIT-VEGETABLES"

# Keep only rows where FEATURE_TY is NOT "FRUIT-VEGETABLES" or "CROP" or "CULTIVATION_ABANDONED"
v_hfi_filt <- v_hfi[!v_hfi$FEATURE_TY %in% c("FRUIT-VEGETABLES", "CROP", "CULTIVATION_ABANDONED"), ]

# Quick check
unique(v_hfi_filt$FEATURE_TY)
v_hfi_filt # 194641, 10  (geometries, attributes)

# Filter to include only records with year 2015 and later
v_hfi_filt_gte2015 <- v_hfi_filt[v_hfi_filt$YEAR >= 2015, ]
v_hfi_filt_gte2015 # 18952, 10  (geometries, attributes)

##---- Rasterize vector data ----##

# ---- 1) Ensure both layers use the same CRS as 'sa' ----
if (!same.crs(v_hfi, sa)) {
  v_hfi <- project(v_hfi, sa)  # reproject polygons to match 'sa'
}

# ---- 2) Build a categorical code for FEATURE_TY ----
# terra rasters can't store character values; we map each category to an integer
cats  <- sort(unique(v_hfi$FEATURE_TY))
codes <- seq_along(cats)
# attach integer code to the vector attributes
v_hfi$cult_type <- match(v_hfi$FEATURE_TY, cats)

# ---- 3) Create an empty raster aligned to 'sa' ----
r_tmpl <- rast(sa)  # inherits extent, resolution, origin, and CRS from 'sa'

# ---- 4) Rasterize polygons using the integer code ----
# For polygons, 'touches=TRUE' helps fill edge pixels; you can set it to FALSE if undesired.
r_code <- rasterize(v_hfi, y = r_tmpl,
                    field     = "cult_type",
                    touches   = TRUE,
                    background = NA)

plot(r_code)
r_code

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008) 
#source      : spat_33cec758535c2_212204_l9EFoqzFKoIF0vD.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#name        : cult_type 
#min value   :         1 
#max value   :         5 

# Mask to 'sa'’s valid footprint to enforce NAs:
#r_code <- mask(r_code, sa)

# ---- 5) Make the raster categorical (factor) and attach the lookup table (RAT) ----
r_cat <- as.factor(r_code)
# Set the level table so that cell values map back to the two-letter codes
levels(r_cat) <- data.frame(ID = codes, FEATURE_TY = cats)

# ---- 6) Write to disk ----
# Save the raster as 16-bit unsigned integer
hfi_out_path <- file.path(hfi_aligned_dir, "hfi_rast_reproj_algn.tif")
writeRaster(r_cat, hfi_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# ---- 7) Quick check ----
levels(r_cat)       # view the RAT
plot(r_cat)         # quick look (legend will show FEATURE_TY codes)
r_cat

# print frequency table
freq_table <- freq(r_cat) 
print(freq_table)


## reclassify to binary (1 = grass)

# 1) Confirm the levels table and the column name that holds the labels
levels(r_cat)[[1]] |> head()

# 2) Define which labels map to 1 (grass) and all others to 0
grass_codes <- c("ROUGH_PASTURE", "TAME_PASTURE")  

# Binary raster: 1 for grass classes, 0 otherwise (except NA)
hfi_algn_binary <- ifel(
  is.na(r_cat),                     # if the input cell is NA
  NA,                               # keep NA
  ifel(r_cat %in% grass_codes, 1, 0)  # otherwise 1 for grass, 0 for non-grass
)

# 5) Quick checks
hfi_algn_binary
freq(hfi_algn_binary)
plot(hfi_algn_binary) 

# Save the raster as 8-bit unsigned integer
hfi_binary_out_path <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary.tif")
writeRaster(hfi_algn_binary, hfi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# load raster
hfi <- file.path(hfi_aligned_dir, "hfi_rast_reproj_algn.tif")
hfi_algn <- rast(hfi)
plot(hfi_algn)
hfi_algn

#### reclassify to native vs. tame (0 = not grass, 1 = tame, 2 = native) - WILL BE SAME AS BINARY

# 1) Build lookup table
#lut_nt <- data.frame(
#  FEATURE_TY = labels_present,
#  recode    = ifelse(labels_present == "HG", 2,          # HG -> 2
#                     ifelse(labels_present %in% c("CP", "CPR"), 1,  # CP, CPR -> 1
#                            0))                                  # all others -> 0
#)

# 2) Substitute using FEATURE_TY
#hfi_algn_nattam <- terra::subs(
#  r_cat,
#  from       = lut_nt,
#  by         = "FEATURE_TY",
#  which      = "recode",
#  subsWithNA = TRUE   # unmatched labels become NA (but all are covered here)
#)

## another way

# Define mappings
#tame_codes <- c("CP", "CPR")   # -> 1
#native_codes <- c("HG")          # -> 2

#hfi_algn_nattam <- ifel(
#  is.na(hfi_algn), 
#  NA,
#  ifel(hfi_algn %in% native_codes, 2,    # HG -> 2
#       ifel(hfi_algn %in% tame_codes, 1,  # CP, CPR -> 1
#            0)                           # all others -> 0
#  )
#)

#plot(hfi_algn_nattam)
#hfi_algn_nattam

# Save the raster as 8-bit unsigned integer
#hfi_nattam_out_path <- file.path(hfi_aligned_dir, "hfi_reproj_algn_natvtam.tif")
#writeRaster(hfi_algn_nattam, hfi_nattam_out_path, datatype = "INT1U", 
#            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
#            NAflag = 255,
#            overwrite = TRUE)

## not present in pilot area!!

################################################

##---- Rasterize vector data ----##

# ---- 1) Ensure both layers use the same CRS as 'sa' ----
if (!same.crs(v_hfi_filt_gte2015, sa)) {
  v_hfi_filt_gte2015 <- project(v_hfi_filt_gte2015, sa)  # reproject polygons to match 'sa'
}

# ---- 2) Build a categorical code for FEATURE_TY ----
# terra rasters can't store character values; we map each category to an integer
cats  <- sort(unique(v_hfi_filt_gte2015$FEATURE_TY))
codes <- seq_along(cats)
# attach integer code to the vector attributes
v_hfi_filt_gte2015$cult_type <- match(v_hfi_filt_gte2015$FEATURE_TY, cats)

# ---- 3) Create an empty raster aligned to 'sa' ----
r_tmpl <- rast(sa)  # inherits extent, resolution, origin, and CRS from 'sa'

# ---- 4) Rasterize polygons using the integer code ----
# For polygons, 'touches=TRUE' helps fill edge pixels; you can set it to FALSE if undesired.
r_code <- rasterize(v_hfi_filt_gte2015, y = r_tmpl,
                    field     = "cult_type",
                    touches   = TRUE,
                    background = NA)

plot(r_code)
r_code

# Mask to 'sa'’s valid footprint to enforce NAs:
#r_code <- mask(r_code, sa)

# ---- 5) Make the raster categorical (factor) and attach the lookup table (RAT) ----
r_cat <- as.factor(r_code)
# Set the level table so that cell values map back to the two-letter codes
levels(r_cat) <- data.frame(ID = codes, FEATURE_TY = cats)
r_cat

#class       : SpatRaster 
#size        : 86539, 130223, 1  (nrow, ncol, nlyr)
#resolution  : 30, 30  (x, y)
#extent      : -2321003, 1585687, 270327.4, 2866497  (xmin, xmax, ymin, ymax)
#coord. ref. : North_America_Albers_Equal_Area_Conic (ESRI:102008) 
#source      : spat_62f546b7b217c_405332_o8up7yuDdJ6XHBX.tif 
#varname     : class_study_area_aligned_cnvc_Dec2025 
#categories  : FEATURE_TY 
#name        :    FEATURE_TY 
#min value   : ROUGH_PASTURE 
#max value   :  TAME_PASTURE 

# ---- 6) Write to disk ----
# Save the raster as 16-bit unsigned integer
hfi_out_path <- file.path(hfi_aligned_dir, "hfi_rast_reproj_algn_gte2015.tif")
writeRaster(r_cat, hfi_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)

# ---- 7) Quick check ----
levels(r_cat)       # view the RAT
plot(r_cat)         # quick look (legend will show FEATURE_TY codes)
r_cat

# print frequency table
freq_table <- freq(r_cat) 
print(freq_table)


## reclassify to binary (1 = grass)

# 1) Confirm the levels table and the column name that holds the labels
levels(r_cat)[[1]] |> head()

# 2) Define which labels map to 1 (grass) and all others to 0
grass_codes <- c("ROUGH_PASTURE", "TAME_PASTURE")  

# Binary raster: 1 for grass classes, 0 otherwise (except NA)
hfi_algn_binary_gte2015 <- ifel(
  is.na(r_cat),                     # if the input cell is NA
  NA,                               # keep NA
  ifel(r_cat %in% grass_codes, 1, 0)  # otherwise 1 for grass, 0 for non-grass
)

# 5) Quick checks
hfi_algn_binary_gte2015
freq(hfi_algn_binary_gte2015)
plot(hfi_algn_binary_gte2015) 

# Save the raster as 8-bit unsigned integer
hfi_binary_out_path <- file.path(hfi_aligned_dir, "hfi_reproj_algn_binary_gte2015.tif")
writeRaster(hfi_algn_binary_gte2015, hfi_binary_out_path, datatype = "INT1U", 
            gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
            NAflag = 255,
            overwrite = TRUE)
